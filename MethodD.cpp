#include <cassert>
#include "bitfile.h"
#include <cmath>
#include "Alignment.h"
#include "AlignmentReader.h"
#include <algorithm>
#include <vector>

#include "utils.h"
#include "MethodD.h"

// @author Johannes Ylinen

bool MethodD::compress(std::string inputfile, std::string inputfile2, std::string outputfile, std::string genomefile)
{
	std::vector<std::pair<Alignment, Alignment> > alignments;

	readAllPairAlignments(alignments, inputfile, inputfile2);
	std::cout << "Found " << alignments.size() << " alignments.\n";
	std::sort(alignments.begin(), alignments.end(), startPosPairComp); // If pre-sorted wouldn't need so much memory

	bit_file_c out;
	/* open bit file for writing */
	try { out.Open(outputfile.c_str(), BF_WRITE); }
	catch (...) { return false; }

	long prevPos = 0;
	for(size_t i = 0; i < alignments.size(); ++i)
	{
		writeAlignment(out, alignments[i].first, prevPos);
		prevPos = alignments[i].first.getStart();
		if(alignments[i].second.getStart() < prevPos)
			out.PutBit(1);
		else
			out.PutBit(0);
		writeAlignment(out, alignments[i].second, prevPos);
	}

	out.Close();
	return true;
}

bool MethodD::decompress(std::string inputfile, std::string first_outputfile, std::string second_outputfile, std::string genomefile)
{
	ofstream out1(first_outputfile.c_str());
	ofstream out2(second_outputfile.c_str());
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

	if(!out1 || !out2){
		cerr << "Failure to open the outputfile." << endl;
		return false;
	}

	long prevPos = 0;
	long readNumber = 1;
	while(true)
	{
		std::string read;
		prevPos = getRead(in, refSeq, read, prevPos);
		if(prevPos < 0)
			break;
		if(read.length() > 0)
		{
			out1 << ">Read_" << readNumber << '\n';
			out1 << read << '\n';
		}

		if(in.GetBit())
			getRead(in, refSeq, read, prevPos, true);
		else
			getRead(in, refSeq, read, prevPos, false);

		if(read.length() > 0)
		{
			out2 << ">Read_" << readNumber++ << '\n';
			out2 << read << '\n';
		}

	}

	return true;
}
