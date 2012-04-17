#include "utils.h"

#include <cmath>
#include <cassert>
#include "AlignmentReader.h"

bool startPosComp(const Alignment& a, const Alignment& b)
{
	return a.getStart() < b.getStart();
}

void readAllAlignments(std::vector<Alignment>& alignments, const std::string& infile)
{
	AlignmentReader reader(AlignmentReader::input_tabdelimited, infile);
	Alignment a;

	while(reader.next(a))
		alignments.push_back(a);
}

void writeGammaCode(bit_file_c& out, long value)
{
	assert(out.good());

	if(value == 0) {
		out.PutBit(0);
		return;
	}

	unsigned length = (unsigned)floor(log2(value+1));		
	// The preceeding 1s
	for(int i = 0; i < length; i++)
		out.PutBit(1);

	// The delimiting 0
	out.PutBit(0);

	// The integer itself using log2(n) bits
	value = value - pow(2, length) + 1;

	out.PutBitsInt(&value, length, sizeof(long));
}

long readGammaCode(bit_file_c& in)
{
	long value = 0;

	// Reading the length of gamma code
	unsigned count = 0;

	int i;
	while((i = in.GetBit()) == 1) {
		count++;
	}

	if(i == 0 && count == 0)
		return 0;

	in.GetBitsInt(&value, count, sizeof(long));

	return value + pow(2,count) - 1;
}

// From readaligner
void revstr(std::string &t)
{
    char c;
    std::size_t n = t.size();
    for (std::size_t i = 0; i < n / 2; ++i) {
        c = t[i];
        t[i] = t[n - i - 1];
        t[n - i - 1] = c;
    }
}

// From readaligner
void complement(std::string &t)
{
    for (std::string::iterator it = t.begin(); it != t.end(); ++it)
    {
        switch (*it)
        {
        case('T'):
            *it = 'A';
            break;
        case('G'):
            *it = 'C';
            break;
        case('C'):
            *it = 'G';
            break;
        case('A'):
            *it = 'T';
            break;
        case('N'):
            *it = 'N';
            break;
        default:
            cerr << "Error: normalized sequence contains an invalid symbol: " << *it << std::endl;
            exit(1);
        }
    }
}


// Creates codes for chromosomes from given genomefile
// Note! Order matters, the file used to decompress must be the same genome file.
std::map<std::string, char> code_chromosomes(std::string genomefile) {

	ifstream in(genomefile.c_str());

	if(!in.is_open()) {
		std::cerr << "Failure to create the chromosome codes, file could not be opened." << std::endl;
		abort();
	}

	std::string row;

	std::vector<std::string> chromosome_names;

	while(getline(in, row)) {
		if(row.at(0) == '>')
			chromosome_names.push_back(row.substr(1, row.find_first_of(' ')-1));
	}

	std::map<string, char> codes;

	for(char i = 0; i < chromosome_names.size(); i++)
		codes[chromosome_names.at(i)] = i;

	return codes;
}

int getEditCode(char c)
{
	switch(c) {

		case 'A':
			return mismatch_A;
		case 'C':
			return mismatch_C;
		case 'G':
			return mismatch_G;
		case 'T':
			return mismatch_T;
		case 'N':
			return mismatch_N;
		case 'a':
			return insertion_A;
		case 'c':
			return insertion_C;
		case 'g':
			return insertion_G;
		case 't':
			return insertion_T;
		case 'D':
			return deletion;
		default:
			std::cerr << "ERROR AT GETEDITCODE() " << c << "\n";
			return -1;
	}
}


void writeEditOp(bit_file_c& out, long edPos, int edCode)
{
	writeGammaCode(out, edPos);
	out.PutBitsInt(&edCode, 4, sizeof(int));
}

std::pair<long, int> readEditOp(bit_file_c& in)
{
	long pos = readGammaCode(in);
	int code = 0;
	in.GetBitsInt(&code, 4, sizeof(int));
	return std::make_pair(pos, code);
}

long modifyString(int edCode, std::string& str, int index)
{
	if(index >= str.length())
		std::cerr << "Edit position " << index << " >= String length " << str.length() << '\n';
	switch(edCode) {
		case mismatch_A:
			str[index] = 'A';	
			return 0;
		case mismatch_C:
			str[index] = 'C';
			return 0;
		case mismatch_G:
			str[index] = 'A';	
			return 0;
		case mismatch_T:
			str[index] = 'C';
			return 0;
		case mismatch_N:
			str[index] = 'N';
			return 0;
		case insertion_A:
			str.insert(index, 1, 'A');
			return 1;
		case insertion_C:
			str.insert(index, 1, 'C');
			return 1;
		case insertion_G:
			str.insert(index, 1, 'G');
			return 1;
		case insertion_T:
			str.insert(index, 1, 'T');
			return 1;
		case deletion:
			str.erase(index, 1);
			return -1;
		default:
			std::cerr << edCode << " <- Unknown edit code\n";
			return 0;
	}
}
