#include <cassert>
#include "bitfile.h"
#include <cmath>
#include "Alignment.h"
#include "AlignmentReader.h"
#include <algorithm>
#include <vector>

#include "utils.h"
#include "MethodB.h"

// @author Johannes Ylinen

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
		writeAlignment(out, alignments[i], prevPos);
		prevPos = alignments[i].getStart();
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
		std::string read;
		prevPos = getRead(in, refSeq, read, prevPos);
		if(prevPos < 0)
			break;
		if(read.length() > 0)
		{
			out << ">Read_" << readNumber++ << '\n';
			out << read << '\n';
		}
	}

	return true;
}
