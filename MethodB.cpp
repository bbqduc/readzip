#include <cassert>
#include "bitfile.h"
#include <cmath>
#include "Alignment.h"
#include "AlignmentReader.h"
#include <algorithm>
#include <vector>

#include "utils.h"

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
