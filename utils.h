/*
 * Generic tools for compressing/decompressing sequence data.
 *
 * @author Anna Kuosmanen & Johannes Ylinen
 *
 */
#include "bitfile.h"
#include <string>
#include <vector>
#include <map>
#include "Alignment.h"

/* Writes gamma code using bitfile. */
void writeGammaCode(bit_file_c& out, long value);

/* Reads gamma code using bitfile. */
long readGammaCode(bit_file_c& in);

/* Comparison operator for Alignments based on start positions. */
bool startPosComp(const Alignment& a, const Alignment& b);

/* Reads all alignments from the file to the vector using AlignmentReader. */
void readAllAlignments(std::vector<Alignment>& alignments, const std::string& infile);

/* Complements the sequence (source: readaligner). */
void complement(std::string &t);

/* Reverses the sequences (source: readaligner).*/
void revstr(std::string &t);

/* Creates codes for the chromosomes in the given genome file. */
std::map<std::string, int> code_chromosomes(std::string genomefile);

enum read_mode_t {read_mode_undef, read_mode_fasta, read_mode_fastq};

/* Prepares the reads for compression by aligning them (Single reads) */
bool align_single(std::string inputfile, std::string genome_file, std::string outputfile, read_mode_t read_mode) ;

/* Prepares the reads for compression by aligning them (Paired reads) */
bool align_pair(std::string input1, std::string input2, std::string genome_file, std::string outputfile_1, std::string outputfile_2, read_mode_t read_mode);

