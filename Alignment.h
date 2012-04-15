#ifndef _Alignment_H_
#define _Alignment_H_

#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>

class Alignment {

public:

	enum input_mode_t {input_tabdelimited};

	Alignment(std::string name_, char strand_, int length_, std::string chromosome_, long start_, std::vector<std::pair<int, char> > edits_);
	Alignment();

	inline std::string getName() const{
		return name;
	}

	inline char getStrand() const{
		return strand;
	}

	inline int getLength() const{
		return length;
	}

	inline std::string getChromosome() const{
		return chromosome;
	}

	inline long getStart() const{
		return start;
	}

	inline const std::vector<std::pair<int, char> >& getEdits() const{
		return edits;
	}

	


protected:

	std::string name;
	char strand;
	int length;
	std::string chromosome;
	long start;
	std::vector<std::pair<int, char> > edits;


};

#endif // _Alignment_H_
