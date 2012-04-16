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
		return -1;
	}

	if(in.GetBitsInt(&value, count, sizeof(long)) == EOF) {
		return -1;
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

	for(int i = 0; i < chromosome_names.size(); i++)
		codes[chromosome_names.at(i)] = i;

	// Dummy for the "no alignment" cases
	codes["*"] = chromosome_names.size();

	return codes;
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
		if(no_more_alignments | a.getName() != id) {

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);
			tolower(c);

			edit_vector.push_back(make_pair(0, c));

			for(int i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				tolower(c);
				edit_vector.push_back(make_pair(1,c));
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

}

// TODO Check for alignments being valid mates (same strand, same chromosome, correct order) and fix if necessary, 
// currently invalidity causes compression to exit
bool align_pair(std::string inputfile_1, std::string inputfile_2, std::string index, std::string outputfile_1, std::string outputfile_2, read_mode_t read_mode) {

	std::string temp_file_1 = outputfile_1 + ".tmp";
	std::string temp_file_2 = outputfile_2 + ".tmp";

	std::string callstring = "./readaligner/readaligner ";

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
		if(no_more_alignments | a.getName() != id) {

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);
			tolower(c);

			edit_vector.push_back(make_pair(0, c));

			for(int i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				tolower(c);
				edit_vector.push_back(make_pair(1,c));
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
		if(no_more_alignments | a.getName() != id) {

			char strand = 'F';
			string chromosome = "*";
			long start = 0;
			int length = 0;

			std::vector<std::pair<int, char> > edit_vector;

			char c = pattern.at(0);
			tolower(c);

			edit_vector.push_back(make_pair(0, c));

			for(int i = 1; i < pattern.length(); i++) {
				c = pattern.at(i);
				tolower(c);
				edit_vector.push_back(make_pair(1,c));
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

}
