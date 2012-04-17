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

bool align_single(std::string inputfile, std::string index, std::string outputfile, read_mode_t read_mode) {

	std::string temp_file = outputfile + ".tmp";

	std::string callstring = "./readaligner/readaligner -i3 ";

	if(read_mode == read_mode_fasta)
		callstring += "--fasta ";
	else if(read_mode == read_mode_fastq)
		callstring += "--fastq ";
	else {
		std::cerr << "Unrecognized read mode. Exiting." << std::endl;
		exit(1);
	}

	system((callstring + " -o " + temp_file + " ./readaligner/" + index + " " + inputfile).c_str());


	//Create "insertion alignments" for unmapped reads

	ifstream in_reads(inputfile.c_str());
	ofstream out(outputfile.c_str());

	if(!in_reads.is_open() | !out.is_open()) {
		std::cerr << "Failure to open files. Exiting." << std::endl;
		exit(1);
	}

	AlignmentReader* alignment_reader = new AlignmentReader(AlignmentReader::input_tabdelimited, temp_file);

	string id;
	string pattern;
	string row;

	std::vector<std::pair<int, char> > dummy;

	Alignment a;

	bool no_more_alignments = false;

	if(!(alignment_reader->next(a)))
		no_more_alignments = true;

	while(getline(in_reads, id)) {

		id = id.substr(1, id.length());

		getline(in_reads,pattern);

		// Read quality rows off
		if(read_mode == read_mode_fastq) {
			getline(in_reads,row);
			getline(in_reads,row);
		}

		// Missing alignment
		if(no_more_alignments | (a.getName() != id)) {

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);

			edit_vector.push_back(make_pair(0, tolower(c)));

			for(unsigned i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				edit_vector.push_back(make_pair((int)i,tolower(c)));
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

	in_reads.close();
	out.close();

	system(("rm " + temp_file).c_str());

	return true;

}

// TODO Check for alignments being valid mates (same strand, same chromosome, correct order) and fix if necessary, 
// currently invalidity causes compression to exit
bool align_pair(std::string inputfile_1, std::string inputfile_2, std::string index, std::string outputfile_1, std::string outputfile_2, read_mode_t read_mode) {

	std::string temp_file_1 = outputfile_1 + ".tmp";
	std::string temp_file_2 = outputfile_2 + ".tmp";

	std::string callstring = "./readaligner/readaligner -i3 ";

	if(read_mode == read_mode_fasta)
		callstring += "--fasta ";
	else if(read_mode == read_mode_fastq)
		callstring += "--fastq ";
	else {
		std::cerr << "Unrecognized read mode. Exiting." << std::endl;
		exit(1);
	}

	system((callstring + " -o " + temp_file_1 + " ./readaligner/" + index + " " + inputfile_1).c_str());
	system((callstring + " -o " + temp_file_2 + " ./readaligner/" + index + " " + inputfile_2).c_str());

	//Create "insertion alignments" for unmapped reads


	// Process first file
	ifstream in_reads(inputfile_1.c_str());
	ofstream out(outputfile_1.c_str());

	if(!in_reads.is_open() | !out.is_open()) {
		std::cerr << "Failure to open files. Exiting." << std::endl;
		exit(1);
	}


	AlignmentReader* alignment_reader = new AlignmentReader(AlignmentReader::input_tabdelimited, temp_file_1);

	string id;
	string pattern;
	string row;

	Alignment a;

	bool no_more_alignments = false;

	if(!alignment_reader->next(a))
		no_more_alignments = true;

	while(getline(in_reads, id)) {

		id = id.substr(1, id.length());

		getline(in_reads,pattern);

		// Read quality rows off
		if(read_mode == read_mode_fastq) {
			getline(in_reads,row);
			getline(in_reads,row);
		}

		// Missing alignment
		if(no_more_alignments | (a.getName() != id)) {

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);

			edit_vector.push_back(make_pair(0, tolower(c)));

			for(unsigned i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				edit_vector.push_back(make_pair(1,tolower(c)));
			}

			Alignment new_a = Alignment(id, strand, length, chromosome, start, edit_vector);

			out << new_a.toString() << endl;

		}

		else {
			out << a.toString() << endl;

			if(!alignment_reader->next(a))
				no_more_alignments = true;

		}

	}

	in_reads.close();
	out.close();


	// Process second file
	in_reads.open(inputfile_2.c_str());
	out.open(outputfile_2.c_str());

	if(!in_reads.is_open() | !out.is_open()) {
		std::cerr << "Failure to open files. Exiting." << std::endl;
		exit(1);
	}


	alignment_reader = new AlignmentReader(AlignmentReader::input_tabdelimited, temp_file_2);

	no_more_alignments = false;

	if(!alignment_reader->next(a))
		no_more_alignments = true;

	while(getline(in_reads, id)) {

		id = id.substr(1, id.length());

		getline(in_reads,pattern);

		// Read quality rows off
		if(read_mode == read_mode_fastq) {
			getline(in_reads,row);
			getline(in_reads,row);
		}

		// Missing alignment
		if(no_more_alignments | (a.getName() != id)) {

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);

			edit_vector.push_back(make_pair(0, tolower(c)));

			for(unsigned i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				edit_vector.push_back(make_pair(1,tolower(c)));
			}

			Alignment new_a = Alignment(id, strand, length, chromosome, start, edit_vector);

			out << new_a.toString() << endl;

		}

		else {
			out << a.toString() << endl;
			if(!alignment_reader->next(a))
				no_more_alignments = true;
		}

	}

	in_reads.close();
	out.close();

	system(("rm " + temp_file_1).c_str());
	system(("rm " + temp_file_2).c_str());

	return true;
}
