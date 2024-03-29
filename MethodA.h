/* 
 * Implements method A from the project description: Compressing/decompressing single-end reads where order needs to be maintained.
 *
 * @author Anna Kuosmanen
 *
 */

#ifndef _MethodA_H_
#define _MethodA_H_
#include <cstdlib>
#include <string>
#include <cstring>
#include "AlignmentReader.h"

class MethodA {

public:

	static bool compress_A(std::string inputfile, std::string outputfile, std::string genomefile);

	static bool decompress_A(std::string inputfile, std::string outputfile, std::string genomefile);

};

#endif //_MethodA_H_
