#include "MethodA.h"
#include "Alignment.h"
#include "AlignmentReader.h"
#include "bitfile.h"

#include <map>
#include <cmath>

using namespace std;


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

map<string, char> code_chromosomes(string genomefile) {

	ifstream in(genomefile.c_str());

	if(!in.is_open()) {
		cerr << "Failure to create the chromosome codes, file could not be opened." << endl;
		abort();
	}

	string row;

	vector<string> chromosome_names;

	while(getline(in, row)) {
		if(row.at(0) == '>')
			chromosome_names.push_back(row.substr(1, row.find_first_of(' ')-1));
	}

	map<string, char> codes;

	for(char i = 0; i < chromosome_names.size(); i++)
		codes[chromosome_names.at(i)] = i;

	return codes;
}

// Simplistic approach that only code mismatches for now
enum edit_codes {mismatch_A, mismatch_C, mismatch_G, mismatch_T};

// Compresses given alignment file.
// Returns true on success and false if there were any problems.
bool MethodA::compress_A(std::string inputfile, string outputfile, string genomefile) {

	AlignmentReader* reader = new AlignmentReader(AlignmentReader::input_tabdelimited, inputfile);

	Alignment a;

	bit_file_c out;

    /* open bit file for writing */
    try
    {
        out.Open(outputfile.c_str(), BF_WRITE);
    }
    catch (...)
    {
        return false;
    }

	map<string, char> chromosome_codes = code_chromosomes(genomefile);

	while(reader->next(a)) {

		// Output code for chromosome
		char chrom_code = chromosome_codes[a.getChromosome()];

		if(out.PutChar(chrom_code) == EOF) {
			cerr << "Error: writing chromosome code!" << endl;
			return false;
		}


		// Output code for strand
		int value;

		if(a.getStrand() == 'R')
			value = 0;
		else
			value = 1;

		if(out.PutBit(value) == EOF) {
			cerr << "Error: writing strand!" << endl;
			return false;
		}

		// Output code for start position (gamma code)
		// ( log2(x) = log(x)/log(2) )
		long start = a.getStart(); 
		unsigned start_length = (unsigned)floor(log(start+1)/log(2));		

		value = 1;

		// The preceeding 1s
		for(int i = 0; i < start_length; i++) {

			if(out.PutBit(value) == EOF) {
				cerr << "Error: writing start!" << endl;
				out.Close();
				return false;
			}
		}

		// The delimiting 0
		value = 0;
		if(out.PutBit(value) == EOF) {
			cerr << "Error: writing start!" << endl;
			out.Close();
			return false;
		}

		// The integer itself using log2(n) bits
		start = start - pow(2, start_length) + 1;

		if(out.PutBitsInt(&start, start_length, sizeof(long)) == EOF)
		{
			cerr << "Error: writing start!" << endl;
			out.Close();
			return false;
		}
		

		// Output code for length (gamma code)
		int length = a.getLength(); 
		unsigned length_length = (unsigned)floor(log(length+1)/log(2));

		value = 1;

		// The preceeding 1s
		for(int i = 0; i < length_length; i++) {

			if(out.PutBit(value) == EOF) {
				cerr << "Error: writing length!" << endl;
				out.Close();
				return false;
			}
		}

		// The delimiting 0
		value = 0;
		if(out.PutBit(value) == EOF) {
			cerr << "Error: writing length!" << endl;
			out.Close();
			return false;
		}
		length = length - pow(2, length_length) + 1;

		// The integer itself using log2(n) bits
		if(out.PutBitsInt(&length, length_length, sizeof(int)) == EOF)
		{
			cerr << "Error: writing length." << endl;
			out.Close();
			return false;
		}

		// Output codes for the edits (position with gamma code and the edits with fixed length)

		out.ByteAlign();
	}

	delete reader;
	out.Close();

	return true;

}

// Decompresses given compressed file using given genome file.
// Returns true on success and false if there were any problems.
bool MethodA::decompress_A(std::string inputfile, std::string outputfile, std::string genomefile){

	ofstream out(outputfile.c_str());
	bit_file_c in;
	ifstream in_genome(genomefile.c_str());

	try
	{
		in.Open(inputfile.c_str(), BF_READ);
	}
	catch (...)
	{
		cerr << "Failure to open the bitfile." << endl;
		return false;
	}

	if(!out.is_open()){
		cerr << "Failure to open the outputfile." << endl;
		return false;
	}

	// Reconstruct the codes for chromosomes
	map<string, char> chromosome_codes = code_chromosomes(genomefile);


	// Read chromosome content
	string info = "";
	string row;

	map<string, string> chromosomes;

	getline(in_genome, row);

	string id = row.substr(1, row.find_first_of(' ')-1);

	while(getline(in_genome, row)) {

		if(row.at(0) == '>') {
			chromosomes[id] = info;
			id = row.substr(1, row.find_first_of(' ')-1);
			info = "";
		}
		else {
			info.append(row);
		}


	}

	chromosomes[id] = info;

	while(!in.eof()) {
		// Get the values	
		char chromosome_code;

		// If this returns EOF, the file has been completely read, it's not an error.
		if((chromosome_code = in.GetChar()) == EOF) {
			in.Close();
			return true;
		}

		string chromosome;

		map<string, char>::iterator it;

		for(it = chromosome_codes.begin(); it != chromosome_codes.end(); it++) {
			if((it->second) == chromosome_code) {
				chromosome = it->first;
				break;
			}
		}

		char strand;

		int i;

		switch(i = in.GetBit()) {

			case 0:
				strand = 'R';
				break;
			case 1:
				strand = 'F';
				break;
			case EOF:
				cerr << "Failure to decompress strand." << endl;
				in.Close();
				return false;

		}

		long start = 0;

		// Reading the length of gamma code
		unsigned count = 0;

		while((i = in.GetBit()) == 1) {
			count++;
		}

		if(i == EOF) {
			cerr << "Failure to decompress start position (EOF)." << endl;
			in.Close();
			return false;
		}

		if(in.GetBitsInt(&start, count, sizeof(long)) == EOF) {
			cerr << "Failure to decompress start position (not enough bits)." << endl;
			in.Close();
			return false;
		}

		start = start + pow(2,count) - 1;


		int length = 0;

		// Reading the length of gamma code
		count = 0;

		while((i = in.GetBit()) == 1) {
			count++;
		}

		if(i == EOF) {
			cerr << "Failure to decompress length (EOF)." << endl;
			in.Close();
			return false;
		}

		if(in.GetBitsInt(&length, count, sizeof(int)) == EOF) {
			cerr << "Failure to decompress length (not enough bits)." << endl;
			in.Close();
			return false;
		}

		length = length + pow(2,count) - 1;





		vector<pair<int, char> > edits;


		// Reconstruct the data
		string data = chromosomes[chromosome].substr(start-1, length);

		for(int i = 0; i < edits.size(); i++) {

			data[edits.at(i).first-1] = edits.at(i).second;		

		}

		if(strand == 'R') {
			revstr(data);
			complement(data);
		}

		out << data << endl;

		in.ByteAlign();
	}
}
