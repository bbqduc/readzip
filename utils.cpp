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
		if(out.PutBit(value) == EOF) {
			std::cerr << "Error writing gamma code." << std::endl;
			out.Close();
			exit(-1);
		}
		return;
	}

	unsigned length = (unsigned)floor(log2(value+1));		
	// The preceeding 1s
	int temp = 1;
	for(int i = 0; i < length; i++) {
		if(out.PutBit(temp) == EOF) {
			std::cerr << "Error writing gamma code." << std::endl;
			out.Close();
			exit(-1);
		}
	}

	// The delimiting 0
	temp = 0;
	if(out.PutBit(temp) == EOF) {
		std::cerr << "Error writing gamme code." << std::endl;
		out.Close();
		exit(-1);
	}

	// The integer itself using log2(n) bits
	value = value - pow(2, length) + 1;

	if(out.PutBitsInt(&value, length, sizeof(long)) == EOF)
	{
		std::cerr << "Error writing gamma code." << std::endl;
		out.Close();
		exit(-2);
	}

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

	if(i == EOF) {
		std::cerr << "Error reading gamma code." << std::endl;
		in.Close();
	}

	if(in.GetBitsInt(&value, count, sizeof(long)) == EOF) {
		std::cerr << "Error reading gamma code." << std::endl;
		in.Close();
	}

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

