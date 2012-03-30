#ifndef _AlignmentReader_H_
#define _AlignmentReader_H_

#include <sstream>
#include <iostream>
#include <fstream>
#include "Alignment.h"

using namespace std;

class AlignmentReader {

public:
	enum input_format_t {input_tabdelimited };

	AlignmentReader(input_format_t mode_, std::string file);

	/* Reads the next alignment. */
	bool next(Alignment &alignment);

private:

	ifstream in;
	input_format_t mode;

	static vector<string> split(const char *str, char c = ' ')
	{
	    vector<string> result;

	    while(1)
	    {
	    	const char *begin = str;

	    	while(*str != c && *str)
	    		str++;

	    	result.push_back(string(begin, str));

	    	if(0 == *str++)
	    		break;
	    }

	    return result;
	}



};

#endif // _AlignmentReader_H_
