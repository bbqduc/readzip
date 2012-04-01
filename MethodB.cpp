#include <cassert>
#include "bitfile.h"
#include <cmath>
#include "Alignment.h"
#include "AlignmentReader.h"
#include <algorithm>
#include <vector>

void writeGammaCode(bit_file_c& out, long value)
{
	assert(out.good());

	unsigned length = (unsigned)floor(log(value+1)/log(2));		
	// The preceeding 1s
	int temp = 1;
	for(int i = 0; i < length; i++) {
		if(out.PutBit(temp) == EOF) {
			std::cerr << "Error: writing start!" << std::endl;
			out.Close();
			exit(-1);
		}
	}

	// The delimiting 0
	temp = 0;
	if(out.PutBit(temp) == EOF) {
		std::cerr << "Error: writing start!" << std::endl;
		out.Close();
		exit(-1);
	}

	// The integer itself using log2(n) bits
	value = value - pow(2, length) + 1;

	if(out.PutBitsInt(&value, length, sizeof(long)) == EOF)
	{
		std::cerr << "Error: writing start!" << std::endl;
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

	if(i == EOF) {
		std::cerr << "Failure to decompress value position (EOF)." << std::endl;
		in.Close();
	}

	if(in.GetBitsInt(&value, count, sizeof(long)) == EOF) {
		std::cerr << "Failure to decompress value position (not enough bits)." << std::endl;
		in.Close();
	}

	return value + pow(2,count) - 1;
}

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

bool compress_B(std::string infile, string outputfile, string genomefile) {
	std::vector<Alignment> alignments;

	readAllAlignments(alignments, infile);
	std::cout << "Found " << alignments.size() << " alignments.\n";
	std::sort(alignments.begin(), alignments.end(), startPosComp);

	bit_file_c out;
	/* open bit file for writing */
	try { out.Open(outputfile.c_str(), BF_WRITE); }
	catch (...) { return false; }

	long prevPos = 0;
	for(size_t i = 0; i < alignments.size(); ++i)
	{
		long posField = alignments[i].getStart() - prevPos;
		long lengthField = alignments[i].getLength();
		long edField = alignments[i].getEdits().size();
		prevPos = posField;

		writeGammaCode(out, posField);
		std::cout << posField << "\n";
		writeGammaCode(out, lengthField);
		std::cout << lengthField << "\n";
		writeGammaCode(out, edField);
		std::cout << edField << "\n";

		// TODO : edit ops

		out.ByteAlign();
	}

	out.Close();
	return true;
}

bool decompress_B(std::string inputfile, std::string outputfile, std::string genomefile)
{
	ofstream out(outputfile.c_str());
	bit_file_c in;
	ifstream in_genome(genomefile.c_str());

	std::string refSeq;

	while(true)
	{
		std::string temp;
		getline(in_genome, temp);
		if(!in_genome) break;
		if(temp[0] != '>')
			refSeq += temp;
	}

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

	while(in.good())
	{
		long posField = readGammaCode(in);
		if(!in.good())
			break;
		long lengthField = readGammaCode(in);
		long edField = readGammaCode(in);
		in.ByteAlign();

		// TODO
	}


	return true;
}
