#include <cassert>
#include "bitfile.h"
#include <cmath>
#include "Alignment.h"
#include "AlignmentReader.h"
#include <algorithm>
#include <vector>

#include "utils.h"
#include "MethodB.h"


bool MethodB::compress(std::string infile, string outputfile, string genomefile) 
{
	std::vector<Alignment> alignments;

	readAllAlignments(alignments, infile);
	std::cout << "Found " << alignments.size() << " alignments.\n";
	std::sort(alignments.begin(), alignments.end(), startPosComp); // If pre-sorted wouldn't need so much memory

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
		prevPos = alignments[i].getStart();

		writeGammaCode(out, posField);
		writeGammaCode(out, lengthField);
		out.PutBit(alignments[i].getStrand() == 'F' ? 0 : 1);
		writeGammaCode(out, edField);

		// Write edit ops

		long prevEdPos = 0;
		for(int j = 0; j < edField; ++j)
		{
			long edPos = alignments[i].getEdits()[j].first;
			edPos -= prevEdPos; // Assuming here that edit ops come in increasing order by position
			prevEdPos = alignments[i].getEdits()[j].first;
			int edCode = getEditCode(alignments[i].getEdits()[j].second);
			writeEditOp(out, edPos, edCode);
		}
	}

	out.Close();
	return true;
}


bool MethodB::decompress(std::string inputfile, std::string outputfile, std::string genomefile)
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

	long prevPos = 0;
	long readNumber = 1;
	while(true)
	{
		long posField = readGammaCode(in);
		if(!in.good())
			break;
		posField += prevPos;
		prevPos = posField;
		long lengthField = readGammaCode(in);
		std::string read = refSeq.substr(posField, lengthField);
		if(in.GetBit())
		{
			complement(read);
			std::reverse(read.begin(), read.end());
		}
		if(posField >= refSeq.length())
			std::cerr << posField << " >= " << refSeq.length() << '\n';

		long edField = readGammaCode(in);
		std::cout << posField << '\t' << lengthField << '\t' << edField;

		long lastEditPos = 0;
		long offset = 0;
		for(long i = 0; i < edField;++i)
		{
			std::pair<long, int> edOp = readEditOp(in);
			lastEditPos += edOp.first;
			offset += modifyString(edOp.second, read, lastEditPos+offset);
			std::cout << '\t' << edOp.first << ' ' << edOp.second << ' ' << offset;
		}

		std::cout << '\n';

		if(read.length() > 0)
		{
			out << ">Read_" << readNumber++ << '\n';
			out << read << '\n';
		}
	}

	return true;
}
