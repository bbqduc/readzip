#include "MethodA.h"
#include "Alignment.h"
#include "AlignmentReader.h"
#include "bitfile.h"
#include "utils.h"

#include <map>
#include <cmath>

using namespace std;

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

		writeGammaCode(out, a.getStart());
		writeGammaCode(out, a.getLength());

		// Output codes for the edits (position with gamma code and the edits with fixed length)
		// Output of readaligner codes mismatches with ACGT and insertions with acgt
		vector<pair<int,char> > edits = a.getEdits();
		int size = edits.size();

		// Number of edits first that know how many to read back
		writeGammaCode(out, size);

		for(int i = 0; i < size; i++) {

			writeGammaCode(out, edits.at(i).first);

			int edit_value = 0;

			switch(edits.at(i).second) {

				case 'A':
					edit_value = mismatch_A;
					break;
				case 'C':
					edit_value = mismatch_C;
					break;
				case 'G':
					edit_value = mismatch_G;
					break;
				case 'T':
					edit_value = mismatch_T;
					break;
				case 'a':
					edit_value = insertion_A;
					break;
				case 'c':
					edit_value = insertion_C;
					break;
				case 'g':
					edit_value = insertion_G;
					break;
				case 't':
					edit_value = insertion_T;
					break;
				case 'D':
					edit_value = deletion;
					break;
			}

			if(out.PutBitsInt(&edit_value, 4, sizeof(edit_value)) == EOF) {
				cerr << "Error: writing edits!" << endl;
				return false;
			}

			edits.clear();

		}

		// Sets the buffer to nearest byte (decompression doesn't always work correctly if this isn't done)
		out.ByteAlign();

	}

	delete reader;
	out.Close();
	chromosome_codes.clear();

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

	while(true) {
		// Get the values	
		char chromosome_code;

		// If this returns EOF, the file has been completely read, it's not an error.
		if((chromosome_code = in.GetChar()) == EOF) {
			out.close();
			in.Close();
			chromosome_codes.clear();
			chromosomes.clear();
			return true;
		}

		string chromosome = "";

		map<string, char>::iterator it;

		for(it = chromosome_codes.begin(); it != chromosome_codes.end(); it++) {
			if((it->second) == chromosome_code) {
				chromosome = it->first;
				break;
			}
		}

		if(chromosome == "") {
			cerr << "Failure to decompress chromosome." << endl;
			in.Close();
			return false;
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

		long start = readGammaCode(in);
		long length = readGammaCode(in);

		vector<pair<int, char> > edits;

		// Read how many edits there are
		int edits_size = readGammaCode(in);

		// Read the edits
		for(int i = 0; i < edits_size; i++) {

			int pos = readGammaCode(in);
			int edit_number = 0;
			char edit;

			if(in.GetBitsInt(&edit_number, 4, sizeof(edit_number)) == EOF) {
				cerr << "Failure to decompress edits." << endl;
				return false;
			}

			// Returning the edits to form supported by Alignment object isn't necessary for reconstructing the sequence, 
			// but it will come in handy if user wants the alignment returned too (future development?)
			switch(edit_number) {

				case mismatch_A:
					edit = 'A';
					break;
				case mismatch_C:
					edit = 'C';
					break;
				case mismatch_G:
					edit = 'G';
					break;
				case mismatch_T:
					edit = 'T';
					break;
				case insertion_A:
					edit = 'a';
					break;
				case insertion_C:
					edit = 'c';
					break;
				case insertion_G:
					edit = 'g';
					break;
				case insertion_T:
					edit = 't';
					break;
				case deletion:
					edit = 'D';
					break;
			}

			edits.push_back(make_pair(pos, edit));
		}

		// Reconstruct the data
		string data = chromosomes[chromosome].substr(start-1, length);

		// Indels can mess up the indexes. Offset keeps track of them.
		int offset = 0;

		for(int i = 0; i < edits.size(); i++) {

			switch(edits.at(i).second) {

				case 'A':
					data[edits.at(i).first+offset] = 'A';	
					break;
				case 'C':
					data[edits.at(i).first+offset] = 'C';
					break;
				case 'G':
					data[edits.at(i).first+offset] = 'G';
					break;
				case 'T':
					data[edits.at(i).first+offset] = 'T';
					break;
				case 'a':
					data = data.substr(0,edits.at(i).first+offset) + "A" + data.substr(edits.at(i).first+offset, data.length());
					offset++;
					break;
				case 'c':
					data = data.substr(0,edits.at(i).first+offset) + "C" + data.substr(edits.at(i).first+offset, data.length());
					offset++;
					break;
				case 'g':
					data = data.substr(0,edits.at(i).first+offset) + "G" + data.substr(edits.at(i).first+offset, data.length());
					offset++;
					break;
				case 't':
					data = data.substr(0,edits.at(i).first+offset) + "T" + data.substr(edits.at(i).first+offset, data.length());
					offset++;
					break;
				case 'D':
					data = data.substr(0,edits.at(i).first+offset) + data.substr(edits.at(i).first+1+offset, data.length());
					offset--;
					break;
			}
		}

		if(strand == 'R') {
			revstr(data);
			complement(data);
		}

		out << data << endl;

		in.ByteAlign();

		edits.clear();
	}
}
