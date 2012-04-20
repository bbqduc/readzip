#include "utils.h"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "AlignmentReader.h"

bool startPosComp(const Alignment& a, const Alignment& b)
{
	return a.getStart() < b.getStart();
}

void writeAlignment(bit_file_c& out, Alignment& a, long prevPos)
{
	long posField = a.getStart() - prevPos;
	if(posField < 0)
		posField = -posField;
	long lengthField = a.getLength();
	long edField = a.getEdits().size();

	writeGammaCode(out, posField);
	writeGammaCode(out, lengthField);
	out.PutBit(a.getStrand() == 'F' ? 0 : 1);
	writeGammaCode(out, edField);

	// Write edit ops

	long prevEdPos = 0;
	for(int j = 0; j < edField; ++j)
	{
		long edPos = a.getEdits()[j].first;
		edPos -= prevEdPos; // Assuming here that edit ops come in increasing order by position
		prevEdPos = a.getEdits()[j].first;
		int edCode = getEditCode(a.getEdits()[j].second);
		writeEditOp(out, edPos, edCode);
	}
}

bool startPosPairComp(const std::pair<Alignment, Alignment>& a, const std::pair<Alignment, Alignment>& b)
{
	return startPosComp(a.first, b.first); 
}

void readAllAlignments(std::vector<Alignment>& alignments, const std::string& infile)
{
	AlignmentReader reader(AlignmentReader::input_tabdelimited, infile);
	Alignment a;

	while(reader.next(a))
		alignments.push_back(a);
}

void readAllPairAlignments(std::vector<std::pair<Alignment, Alignment> >& alignments, const std::string& infile1, const std::string& infile2)
{
	AlignmentReader first_reader(AlignmentReader::input_tabdelimited, infile1);
	AlignmentReader second_reader(AlignmentReader::input_tabdelimited, infile2);

	Alignment a, b;
	while(first_reader.next(a) && second_reader.next(b))
		alignments.push_back(std::make_pair(a,b));
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
	for(unsigned i = 0; i < length; i++)
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
std::map<std::string, int> code_chromosomes(std::string genomefile) {

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

	std::map<string, int> codes;

	sort(chromosome_names.begin(), chromosome_names.end());

	for(unsigned i = 0; i < chromosome_names.size(); i++)
		codes[chromosome_names.at(i)] = i;

	// Dummy for the "no alignment" cases
	codes["*"] = chromosome_names.size();

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
		case 'n':
			return insertion_N;
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

long modifyString(int edCode, std::string& str, size_t index)
{
	if(index > str.length())
		std::cerr << "Edit position " << index << " > String length " << str.length() << '\n';
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
		case insertion_N:
			str.insert(index, 1, 'N');
			return 1;
		case deletion:
			str.erase(index, 1);
			return -1;
		default:
			std::cerr << edCode << " <- Unknown edit code\n";
			return 0;
	}
}

bool align_single(std::string inputfile, std::string index, std::string outputfile, read_mode_t read_mode, bool maintainOrder) {

	std::string temp_file = outputfile + ".tmp";

	std::string callstring = "awk -f rnreads.awk " + inputfile + " > " + inputfile + ".tmp";

	system(callstring.c_str());

	callstring = "./readaligner/readaligner -P0 -i3 -v ";

	if(read_mode == read_mode_fasta)
		callstring += "--fasta ";
	else if(read_mode == read_mode_fastq)
		callstring += "--fastq ";
	else {
		std::cerr << "Unrecognized read mode. Exiting." << std::endl;
		exit(1);
	}

	system((callstring + " -o " + temp_file + " " + index + " " + inputfile + ".tmp").c_str());

	callstring = "sort -t 'd' -n +1 -2 " + temp_file + " > " + temp_file + ".sorted";
	system(callstring.c_str());
	system(("rm " + temp_file).c_str());

	//Create "insertion alignments" for unmapped reads
	ifstream in_reads((inputfile + ".tmp").c_str());
	ofstream out(outputfile.c_str());

	if(!in_reads.is_open() | !out.is_open()) {
		std::cerr << "Failure to open files. Exiting." << std::endl;
		exit(1);
	}

	AlignmentReader* alignment_reader = new AlignmentReader(AlignmentReader::input_tabdelimited, temp_file+".sorted");

	string id;
	string pattern;
	string row;

	Alignment a;

	bool no_more_alignments = false;

	int number_of_missing = 0;
	int total_number = 0;

	if(!(alignment_reader->next(a)))
		no_more_alignments = true;

	while(getline(in_reads, id)) {

		total_number++;

		id = id.substr(1, id.length());

		getline(in_reads,pattern);

		// Read quality rows off
		if(read_mode == read_mode_fastq) {
			getline(in_reads,row);
			getline(in_reads,row);
		}

		// Missing alignment
		if(no_more_alignments | (a.getName() != id)) {

			number_of_missing++;

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);

			edit_vector.push_back(make_pair(0, tolower(c)));

			for(unsigned i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				edit_vector.push_back(make_pair(maintainOrder ? (int)i : 0,tolower(c)));
			}

			Alignment new_a(id, strand, length, chromosome, start, edit_vector);

			out << new_a.toString() << endl;

		}

		else {
			out << a.toString() << endl;

			if(!alignment_reader->next(a))
				no_more_alignments = true;
		}

	}

	system(("rm " + inputfile + ".tmp").c_str());
	system(("rm " + temp_file + ".sorted").c_str());

	in_reads.close();
	out.close();


	if(number_of_missing > (0.5 * total_number))
		std::cerr << "Warning: more than half the reads failed to align, compression isn't good." << std::endl;

	return true;

}

bool align_pair(std::string inputfile_1, std::string inputfile_2, std::string index, std::string outputfile_1, std::string outputfile_2, read_mode_t read_mode, bool maintainOrder) {

	std::string temp_file_1 = outputfile_1 + ".tmp";
	std::string temp_file_2 = outputfile_2 + ".tmp";

	std::string callstring = "awk -f rnreads.awk " + inputfile_1 + " > " + inputfile_1 + ".tmp";

	system(callstring.c_str());
	callstring = "awk -f rnreads.awk " + inputfile_2 + " > " + inputfile_2 + ".tmp";
	system(callstring.c_str());

	// Ask for max 10 alignments per read, out of those bigger chance to find matching pair
	callstring = "./readaligner/readaligner -P0 -i3 -r10 -v ";

	if(read_mode == read_mode_fasta)
		callstring += "--fasta ";
	else if(read_mode == read_mode_fastq)
		callstring += "--fastq ";
	else {
		std::cerr << "Unrecognized read mode. Exiting." << std::endl;
		exit(1);
	}

	system((callstring + " -o " + temp_file_1 + " " + index + " " + inputfile_1 + ".tmp").c_str());
	system((callstring + " -o " + temp_file_2 + " " + index + " " + inputfile_2 + ".tmp").c_str());

	//Create "insertion alignments" for unmapped reads

	ifstream in_reads_1((inputfile_1 + ".tmp").c_str());
	ofstream out_1(outputfile_1.c_str());

	ifstream in_reads_2((inputfile_2 + ".tmp").c_str());
	ofstream out_2(outputfile_2.c_str());


	if(!in_reads_1.is_open() | !out_1.is_open() | !in_reads_2.is_open() | !out_2.is_open()) {
		std::cerr << "Failure to open files. Exiting." << std::endl;
		exit(1);
	}

	callstring = "sort -t 'd' -n +1 -2 " + temp_file_1 + " > " + temp_file_1 + ".sorted";
	system(callstring.c_str());
	system(("rm " + temp_file_1).c_str());
	callstring = "sort -t 'd' -n +1 -2 " + temp_file_2 + " > " + temp_file_2 + ".sorted";
	system(callstring.c_str());
	system(("rm " + temp_file_2).c_str());

	AlignmentReader* alignment_reader_1 = new AlignmentReader(AlignmentReader::input_tabdelimited, temp_file_1+".sorted");
	AlignmentReader* alignment_reader_2 = new AlignmentReader(AlignmentReader::input_tabdelimited, temp_file_2+".sorted");

	string id_1;
	string pattern_1;
	string row_1;

	Alignment a_1;

	string id_2;
	string pattern_2;
	string row_2;

	Alignment a_2;


	bool no_more_alignments_1 = false;
	bool no_more_alignments_2 = false;

	int number_of_missing = 0;
	int total_number = 0;

	if(!alignment_reader_1->next(a_1))
		no_more_alignments_1 = true;
	if(!alignment_reader_2->next(a_2))
		no_more_alignments_2 = true;

	while(getline(in_reads_1, id_1) && getline(in_reads_2, id_2)) {

		total_number++;

		id_1 = id_1.substr(1, id_1.length());

		getline(in_reads_1,pattern_1);

		id_2 = id_2.substr(1, id_1.length());

		getline(in_reads_2,pattern_2);

		// Read quality rows off
		if(read_mode == read_mode_fastq) {
			getline(in_reads_1,row_1);
			getline(in_reads_1,row_1);

			getline(in_reads_2,row_2);
			getline(in_reads_2,row_2);
		}

		// Missing alignment
		if(no_more_alignments_1 | (a_1.getName() != id_1) | no_more_alignments_2 | (a_2.getName() != id_2)) {

			number_of_missing++;

			// Create insertion alignments for the two
			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector_1;

			char c = pattern_1.at(0);

			edit_vector_1.push_back(make_pair(0, tolower(c)));

			for(unsigned i = 1; i < pattern_1.length(); i++) {
				c = pattern_1.at(i);
				edit_vector_1.push_back(make_pair(maintainOrder ? 1 : 0,tolower(c)));
			}

			Alignment new_a_1 = Alignment(id_1, strand, length, chromosome, start, edit_vector_1);

			out_1 << new_a_1.toString() << endl;


			std::vector<std::pair<int, char> > edit_vector_2;

			c = pattern_2.at(0);

			edit_vector_2.push_back(make_pair(0, tolower(c)));

			for(unsigned i = 1; i < pattern_2.length(); i++) {
				c = pattern_2.at(i);
				edit_vector_2.push_back(make_pair(maintainOrder ? 1 : 0,tolower(c)));
			}

			Alignment new_a_2 = Alignment(id_2, strand, length, chromosome, start, edit_vector_2);

			out_2 << new_a_2.toString() << endl;

			// Scroll the AlignmentReaders to next read if necessary
			while(!no_more_alignments_1 && (a_1.getName() == id_1))
				if(!alignment_reader_1->next(a_1))
					no_more_alignments_1 = true;

			while(!no_more_alignments_2 && (a_2.getName() == id_2))
				if(!alignment_reader_2->next(a_2))
					no_more_alignments_2 = true;
		}

		else {

			bool found_match = false;

			vector<Alignment> first_alignments;
			vector<Alignment> second_alignments;

			while(!no_more_alignments_1 && (a_1.getName() == id_1)) {
				first_alignments.push_back(a_1);
				if(!alignment_reader_1->next(a_1))
					no_more_alignments_1 = true;
			}

			while(!no_more_alignments_2 && (a_2.getName() == id_2)) {
				second_alignments.push_back(a_2);
				if(!alignment_reader_2->next(a_2))
					no_more_alignments_2 = true;
			}

			for(unsigned i = 0; i < first_alignments.size(); i++) {

				Alignment candidate_1 = first_alignments.at(i);

				for(unsigned j = 0; j < second_alignments.size(); j++) {

					Alignment candidate_2 = second_alignments.at(j);

					if((candidate_1.getChromosome() == candidate_2.getChromosome()) && (candidate_2.getStart() > candidate_1.getStart())) {
						out_1 << candidate_1.toString() << endl;
						out_2 << candidate_2.toString() << endl;
						found_match = true;
					}
					if(found_match)
						break;

				}
				if(found_match)
					break;
			}

			if(!found_match) {
				number_of_missing++;

				// Create insertion alignments for the two
				char strand = 'F';
				string chromosome = "*";
				long start = 0;
				int length = 0;

				std::vector<std::pair<int, char> > edit_vector_1;

				char c = pattern_1.at(0);

				edit_vector_1.push_back(make_pair(0, tolower(c)));

				for(unsigned i = 1; i < pattern_1.length(); i++) {
					c = pattern_1.at(i);
					edit_vector_1.push_back(make_pair(maintainOrder ? 1 : 0,tolower(c)));
				}

				Alignment new_a_1 = Alignment(id_1, strand, length, chromosome, start, edit_vector_1);

				out_1 << new_a_1.toString() << endl;


				std::vector<std::pair<int, char> > edit_vector_2;

				c = pattern_2.at(0);

				edit_vector_2.push_back(make_pair(0, tolower(c)));

				for(unsigned i = 1; i < pattern_2.length(); i++) {
					c = pattern_2.at(i);
					edit_vector_2.push_back(make_pair(maintainOrder ? 1 : 0,tolower(c)));
				}

				Alignment new_a_2 = Alignment(id_2, strand, length, chromosome, start, edit_vector_2);

				out_2 << new_a_2.toString() << endl;
			}
		}
	}

	in_reads_1.close();
	in_reads_2.close();
	out_1.close();
	out_2.close();

	system(("rm " + temp_file_1 + ".sorted").c_str());
	system(("rm " + temp_file_2 + ".sorted").c_str());

	system(("rm " + inputfile_1 + ".tmp").c_str());
	system(("rm " + inputfile_2 + ".tmp").c_str());

	if(number_of_missing > (0.5 * total_number))
		std::cerr << "Warning: more than half the reads failed to align, compression isn't good." << std::endl;

	return true;
}

long getRead(bit_file_c& in, const std::string& reference, std::string& out, long prevPos, bool decreasePos)
{
	long posField = readGammaCode(in);
	if(!in.good())
		return -1;
	if(decreasePos)
		posField = prevPos - posField;
	else
		posField += prevPos;
	long lengthField = readGammaCode(in);
	out = reference.substr(posField, lengthField);
	if(in.GetBit())
	{
		complement(out);
		std::reverse(out.begin(), out.end());
	}
	if(posField >= reference.length())
		std::cerr << posField << " >= " << reference.length() << '\n';

	long edField = readGammaCode(in);

	long lastEditPos = 0;
	long offset = 0;


	for(long i = 0; i < edField;++i)
	{
		std::pair<long, int> edOp = readEditOp(in);
		lastEditPos += edOp.first;
		offset += modifyString(edOp.second, out, lastEditPos+offset);
	}

	return posField;
}
