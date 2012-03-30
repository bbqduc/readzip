#include "Alignment.h"

#include <vector>

Alignment::Alignment(char strand_, int length_, std::string chromosome_, long start_, std::vector<std::pair<int, char> > edits_)
: strand(strand_), length(length_), chromosome(chromosome_), start(start_), edits(edits_)
{}

Alignment::Alignment()
: strand('0'), length(0), chromosome("NULL"), start(0), edits(0)
{}
