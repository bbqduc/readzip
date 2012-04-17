#include "Alignment.h"

#include <vector>
#include <iostream>
#include <sstream>

Alignment::Alignment(std::string name_, char strand_, int length_, std::string chromosome_, long start_, std::vector<std::pair<int, char> > edits_)
: name(name_), strand(strand_), length(length_), chromosome(chromosome_), start(start_), edits(edits_)
{}

Alignment::Alignment()
: name("NULL"), strand('0'), length(0), chromosome("NULL"), start(0), edits(0)
{}


std::string Alignment::toString() {

	std::stringstream sstm;

	sstm << this->getName() << '\t' << this->getChromosome() << '\t' << this->getStart() << '\t' << this->getStart() + this->getLength() -1
		<< '\t' << 1 << '\t' << this->getStrand() << '\t';

	std::vector<std::pair<int, char> > edit = this->edits;

	if(edit.size() > 0) {

		for(int i = 0; i < edit.size()-1; i++)
			sstm << edit.at(i).first << " " << edit.at(i).second << " ";

		sstm << edit.at(edit.size()-1).first << " " << edit.at(edits.size()-1).second;

	}
	return sstm.str();
}
