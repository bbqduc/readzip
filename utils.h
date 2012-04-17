/*
 * Generic tools for compressing/decompressing sequence data.
 *
 * @author Anna Kuosmanen & Johannes Ylinen
 *
 */
#pragma once
#include "bitfile.h"
#include <string>
#include <vector>
#include <map>
#include "Alignment.h"

// Fixed length code (with 4 bits) can be used to display these
enum edit_codes_t {mismatch_A, mismatch_C, mismatch_G, mismatch_T, mismatch_N, insertion_A, insertion_C, insertion_G, insertion_T, deletion};
enum read_mode_t {read_mode_undef, read_mode_fasta, read_mode_fastq};

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

std::pair<long, int> readEditOp(bit_file_c& in);

void writeEditOp(bit_file_c& out, long edPos, int edCode);

int getEditCode(char c);

long modifyString(int edCode, std::string& str, size_t index);

/* Prepares the reads for compression by aligning them (Single reads) */
bool align_single(std::string inputfile, std::string genome_file, std::string outputfile, read_mode_t read_mode) ;

/* Prepares the reads for compression by aligning them (Paired reads) */
bool align_pair(std::string input1, std::string input2, std::string genome_file, std::string outputfile_1, std::string outputfile_2, read_mode_t read_mode);
