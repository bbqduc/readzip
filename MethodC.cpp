#include "MethodC.h"
#include "Alignment.h"
#include "AlignmentReader.h"
#include "bitfile.h"
#include "utils.h"

#include <map>
#include <cmath>

using namespace std;

// Compresses given alignment files.
// Returns true on success and false if there were any problems.
bool MethodC::compress_C(string first_inputfile, string second_inputfile, string outputfile, string genomefile) {

	AlignmentReader* first_reader = new AlignmentReader(AlignmentReader::input_tabdelimited, first_inputfile);
	AlignmentReader* second_reader = new AlignmentReader(AlignmentReader::input_tabdelimited, second_inputfile);

	Alignment a_1;
	Alignment a_2;

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

	map<string, int> chromosome_codes = code_chromosomes(genomefile);

	// Find out how many bits needed for fixed length
	int bits = ceil(log2(chromosome_codes.size()));

	while(first_reader->next(a_1)) {

		if(!(second_reader->next(a_2))) {

			cerr << "Second inputfile ended before the first, error in syncronizing the alignments." << endl;
			return false;

		}

		// Sanity check that pairs have been aligned correctly
		if((a_1.getChromosome() != a_2.getChromosome()) && (a_2.getStart() > a_1.getStart())) {

			cerr << "Mates have different chromosome, please check the alignment for read " << a_1.getName() << "." << endl;
			return false;

		}

		// Output code for chromosome
		int chrom_code = chromosome_codes[a_1.getChromosome()];

		if(out.PutBitsInt(&chrom_code, bits, sizeof(chrom_code)) == EOF) {
			cerr << "Error: writing chromosome code!" << endl;
			return false;
		}

		for(int mate = 1; mate <= 2; mate++) {

			// Output code for strand
			int value;

			if(mate == 1) {
				if(a_1.getStrand() == 'R')
					value = 0;
				else
					value = 1;
			}
			else {
				if(a_2.getStrand() == 'R')
					value = 0;
				else
					value = 1;
			}
			if(out.PutBit(value) == EOF) {
				cerr << "Error: writing strand!" << endl;
				return false;
			}

			vector<pair<int,char> > edits;


			if(mate == 1) {
				writeGammaCode(out, a_1.getStart());
				writeGammaCode(out, a_1.getLength());
			}
			else {
				writeGammaCode(out, a_2.getStart() - a_1.getStart());
				writeGammaCode(out, a_2.getLength());
			}

			// Output codes for the edits (position with gamma code and the edits with fixed length)
			// Output of readaligner codes mismatches with ACGT and insertions with acgt
			if(mate == 1)
				edits = a_1.getEdits();
			else
				edits = a_2.getEdits();

			int size = edits.size();

			int previous = 0;

			// Number of edits first that know how many to read back
			writeGammaCode(out, size);

			for(int i = 0; i < size; i++) {

				writeGammaCode(out, edits.at(i).first-previous);

				previous = edits.at(i).first;

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

			}

			edits.clear();

		}

	}

	// Check that there's nothing left in second inputfile
	if(second_reader->next(a_2)) {

		cerr << "First input file ended before the second, error in syncronizing the alignments." << endl;
		return false;

	}

	delete first_reader;
	delete second_reader;
	out.Close();

	return true;

}

// Decompresses given compressed file using given genome file.
// Returns true on success and false if there were any problems.
bool MethodC::decompress_C(std::string inputfile, std::string first_outputfile, std::string second_outputfile, std::string genomefile){

	ofstream out_1(first_outputfile.c_str());
	ofstream out_2(second_outputfile.c_str());
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

	if(!out_1.is_open() | !out_2.is_open()){
		cerr << "Failure to open the outputfile(s)." << endl;
		return false;
	}

	// Reconstruct the codes for chromosomes
	map<string, int> chromosome_codes = code_chromosomes(genomefile);

	// Find out how many bits needed for fixed length
	int bits = ceil(log2(chromosome_codes.size()));

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
		int chromosome_code = 0;

		// If this returns EOF, the file has been completely read, it's not an error.
		if((in.GetBitsInt(&chromosome_code, bits, sizeof(chromosome_code))) == EOF) {
			out_1.close();
			out_2.close();
			in.Close();
			chromosome_codes.clear();
			chromosomes.clear();
			return true;
		}

		string chromosome = "";

		map<string, int>::iterator it;

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

		long start = 0;
		long length = 0;

		for(int mate = 1; mate <= 2; mate++) {

			switch(i = in.GetBit()) {

				case 0:
					strand = 'R';
					break;
				case 1:
					strand = 'F';
					break;
				case EOF:
					out_1.close();
					out_2.close();
					in.Close();
					chromosome_codes.clear();
					chromosomes.clear();
					return true;

			}

			// For first mate, start is coded
			if(mate == 1)
				start = readGammaCode(in);
			// For second, it's the difference compared to first mate
			else {
				start += readGammaCode(in);
			}				

			if(start == EOF) {

				out_1.close();
				out_2.close();
				in.Close();
				chromosome_codes.clear();
				chromosomes.clear();
				return true;
			}


			length = readGammaCode(in);

			if(length == EOF) {

				out_1.close();
				out_2.close();
				in.Close();
				chromosome_codes.clear();
				chromosomes.clear();
				return true;
			}


			vector<pair<int, char> > edits;

			// Read how many edits there are
			int edits_size = readGammaCode(in);

			if(edits_size == -1) {

				out_1.close();
				out_2.close();
				in.Close();
				chromosome_codes.clear();
				chromosomes.clear();
				return true;
			}

			int pos = 0;

			// Read the edits
			for(int i = 0; i < edits_size; i++) {

				pos += readGammaCode(in);

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

			string data = "";

			if(chromosome != "*" && start != 0 && length != 0)
				data = chromosomes[chromosome].substr(start-1, length);

			// There's possibility that trailing zeros cause "valid" looking alignment, check!
			if(chromosome != "*" && (length == 0)) {
				out_1.close();
				out_2.close();
				in.Close();
				chromosome_codes.clear();
				chromosomes.clear();
				return true;
			}

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
					{
						if((data.length() > 0) && (edits.at(i).first != 0) && (edits.at(i).first < data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "A" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == 0)) {
							data = "A" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "A";
							offset++;						
						}
						else {
							data = "A";
						}
						break;
					}

					case 'c':
					{
						if((data.length() > 0) && (edits.at(i).first != 0) && (edits.at(i).first < data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "C" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == 0)) {
							data = "C" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "C";
							offset++;						
						}
						else {
							data = "C";
						}
						break;
					}

					case 'g':
					{
						if((data.length() > 0) && (edits.at(i).first != 0) && (edits.at(i).first < data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "G" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == 0)) {
							data = "G" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "G";
							offset++;						
						}
						else {
							data = "G";
						}
						break;
					}
					case 't':
					{
						if((data.length() > 0) && (edits.at(i).first != 0) && (edits.at(i).first < data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "T" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == 0)) {
							data = "T" + data.substr(edits.at(i).first+offset, data.length());
							offset++;
						}
						else if((data.length() > 0) && (edits.at(i).first == data.length())) {
							data = data.substr(0,edits.at(i).first+offset) + "T";
							offset++;						
						}
						else {
							data = "T";
						}
						break;
					}
					case 'D':
					{
						data = data.substr(0,edits.at(i).first+offset) + data.substr(edits.at(i).first+offset+1, data.length());
						offset--;
						break;
					}
				}
			}

			if(strand == 'R') {
				revstr(data);
				complement(data);
			}

			if(mate == 1)
				out_1 << data << endl;
			else
				out_2 << data << endl;

		}
	}
	return true;
}
