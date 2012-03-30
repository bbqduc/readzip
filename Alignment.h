#ifndef _Alignment_H_
#define _Alignment_H_

#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>

class Alignment {

public:

	enum input_mode_t {input_tabdelimited};

	Alignment(char strand_, int length_, std::string chromosome_, long start_, std::vector<std::pair<int, char> > edits_);
	Alignment();


	inline char getStrand() {
		return strand;
	}

	inline int getLength() {
		return length;
	}

	inline std::string getChromosome() {
		return chromosome;
	}

	inline long getStart() {
		return start;
	}

	inline std::vector<std::pair<int, char> > getEdits() {
		return edits;
	}

	


protected:

	char strand;
	int length;
	std::string chromosome;
	long start;
	std::vector<std::pair<int, char> > edits;


};

#endif // _Alignment_H_
