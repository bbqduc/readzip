#include "Alignment.h"
#include "AlignmentReader.h"

AlignmentReader::AlignmentReader(input_format_t mode_, string file)
	: mode(mode_) {

	in.open(file.c_str());

	if(!in.is_open()) {

		cerr << "AlignmentReader: Failed to open the file " << file << "." << endl;
		abort();
	}
}


bool AlignmentReader::next(Alignment &a) {

	string row;

	getline(in, row);

	if(in.eof())
		return false;

	else {

		if(mode == input_tabdelimited) {

			vector<string> alignment_parts = split(row.c_str(), '\t');

			string name = alignment_parts.at(0);

			string chromosome = alignment_parts.at(1);

			long start = atol(alignment_parts.at(2).c_str());

			int length = atol(alignment_parts.at(3).c_str()) - start +1;

			char strand = alignment_parts.at(5).at(0);

			vector<pair<int,char> > edits;

			if(alignment_parts.at(6) != "") {

				vector<string> alignment_edits = split(alignment_parts.at(6).c_str(), ' ');

				for(unsigned i = 0; i < alignment_edits.size(); i = i+2) {

					edits.push_back(make_pair(atoi(alignment_edits.at(i).c_str()), alignment_edits.at(i+1).at(0)));

				}

			}

			a = Alignment(name, strand, length, chromosome, start, edits);

		}
	}

	return true;
}

