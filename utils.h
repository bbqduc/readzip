#include "bitfile.h"
#include <string>
#include <vector>
#include "Alignment.h"

void writeGammaCode(bit_file_c& out, long value);
long readGammaCode(bit_file_c& in);
bool startPosComp(const Alignment& a, const Alignment& b);
void readAllAlignments(std::vector<Alignment>& alignments, const std::string& infile);
void complement(std::string &t);
void revstr(std::string &t);
