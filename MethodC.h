/* 
 * Implements method C from the project description: Compressing/decompressing paired-end reads where order needs to be maintained.
 *
 * @author Anna Kuosmanen
 *
 */

#ifndef _MethodC_H_
#define _MethodC_H_
#include <cstdlib>
#include <string>
#include <cstring>
#include "AlignmentReader.h"

class MethodC {

public:

	static bool compress_C(std::string first_inputfile, std::string second_inputfile, std::string outputfile, std::string genomefile);

	static bool decompress_C(std::string inputfile, std::string first_outputfile, std::string second_outputfile, std::string genomefile);

};

#endif //_MethodC_H_
