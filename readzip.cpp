#include <getopt.h>
#include <map>
#include "AlignmentReader.h"
#include "MethodA.h"
#include "MethodB.h"
#include "MethodC.h"
#include "MethodD.h"
#include "utils.h"

void print_help() {

	std::cerr << " Readzip: Tool for compressing/decompressing sequence reads." << std::endl << std::endl
			<< " Compression/decompression options:" << std::endl
			<< " -o                    Compress." << std::endl
			<< " -x                    Decompress." << std::endl << std::endl
			<< " Compression mode options: " << std::endl
			<< " -a                    Compress single-end reads, maintain order." << std::endl
			<< " -b                    Compress single-end reads, do not maintain order." << std::endl
			<< " -c                    Compress paired-end reads, maintain order." << std::endl
			<< " -d                    Compress paired-end reads, do not maintain order." << std::endl << std::endl
			<< " Input formats:" << std::endl
			<< " -f                    Fasta format." << std::endl
			<< " -q                    Fastq format." << std::endl << std::endl;
}


enum packing_mode_t {packing_mode_undef, packing_mode_a, packing_mode_b, packing_mode_c, packing_mode_d };
enum pack_unpack_mode_t {mode_undef, zip_mode, unzip_mode };

int main(int argc, char **argv) 
{

	packing_mode_t mode = packing_mode_undef;
	pack_unpack_mode_t xc_mode = mode_undef;
	read_mode_t read_mode = read_mode_undef;

	// Parse command line parameters
	int option_index = 0;
	int c;
	while((c = getopt(argc, argv, "abcdxofqh")) != -1)

	{

	switch(c)
		{
		case 'a':

			if(mode != packing_mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}

			mode = packing_mode_a; break;
		case 'b':
			if(mode != packing_mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}

			mode = packing_mode_b; break;
		case 'c':
			if(mode != packing_mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}

			mode = packing_mode_c; break;
		case 'd':
			if(mode != packing_mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}

			mode = packing_mode_d; break;
		case 'x':
			if(xc_mode != mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}
			xc_mode = unzip_mode; break;
		case 'o':
			if(xc_mode != mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}
			xc_mode = zip_mode; break;
		case 'f':
			if(read_mode != read_mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}			
			read_mode = read_mode_fasta;
			break;
		case 'q':
			if(read_mode != read_mode_undef) {
				std::cerr << "readzip: Conflicting mode parameters." << std::endl;
				exit(1);

			}			
			read_mode = read_mode_fastq;
			break;
		case 'h':
			print_help();
			exit(0);
		}
	}

	if (argc - optind == 0)
	{
		cerr << "readzip: for help on usage, use option -h." << endl;
		return 1;
	}


	// Sanity checks
	if(read_mode == read_mode_undef) {
		cerr << "readzip: Please specify either fasta or fastq format." << endl;
		exit(1);
	}
	if(xc_mode == mode_undef) {
		cerr << "readzip: Please specify either compression or decompression mode." << endl;
		exit(1);
	}
	if(mode == packing_mode_undef) {
		cerr << "readzip: Please specify method a, b, c or d." << endl;
		exit(1);
	}

	string genome_file, output_file, input_file;

	// TODO Calling each mode (fill up as modes are implemented)
	switch(mode) {

		case packing_mode_a:
		{

			// Parse filenames
			if (argc - optind < 3)
			{
				cerr << "readzip: missing input files!" << endl;
				return 1;
			}

			genome_file = string(argv[optind++]); 

			input_file = string(argv[optind++]);

			output_file = string(argv[optind++]);


			if(xc_mode == zip_mode) {

				string alignment_file = input_file + ".tab";

				std::string index = genome_file;

				// Call to align
				if(!(align_single(input_file, index, alignment_file, read_mode))) {
					std::cerr << "Error! Failure in aligning the reads." << std::endl;
					exit(1);
				}

				if(MethodA::compress_A(alignment_file, output_file, genome_file)) {
					std::cerr << "Done compressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the compression!" << std::endl;

				system(("rm " + alignment_file).c_str());

			}
			else {

				if(MethodA::decompress_A(input_file, output_file, genome_file)) {
					std::cerr << "Done decompressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the decompression!" << std::endl;
			}

			break;
		}			
		case packing_mode_b:

			genome_file = string(argv[optind++]); 
			input_file = string(argv[optind++]);
			output_file = string(argv[optind++]);

			if(xc_mode == zip_mode) {
				string alignment_file = input_file + ".tab";
				std::string index = genome_file;

				// Call to align
				if(!(align_single(input_file, index, alignment_file, read_mode, false))) {
					std::cerr << "Error! Failure in aligning the reads." << std::endl;
					exit(1);
				}

				if(MethodB::compress(alignment_file, output_file, genome_file)) {
					std::cerr << "Done compressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the compression!" << std::endl;
				system(("rm " + alignment_file).c_str());
			}
			else {
				if(MethodB::decompress(input_file, output_file, genome_file)) {
					std::cerr << "Done decompressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the compression!" << std::endl;
			}
			break;

		case packing_mode_c:

		{

			// Parse filenames
			if (argc - optind < 4)
			{
				cerr << "readzip: missing input files!" << endl;
				return 1;
			}

			string genome_file = string(argv[optind++]); 

			if(xc_mode == zip_mode) {

				string input_file_1 = string(argv[optind++]);

				string input_file_2 = string(argv[optind++]);

				string output_file = string(argv[optind++]);

				string alignment_file_1 = input_file_1 + ".tab";
				string alignment_file_2 = input_file_2 + ".tab";

				std::string index = genome_file;

				// Call to align
				if(!(align_pair(input_file_1, input_file_2, index, alignment_file_1, alignment_file_2, read_mode))) {
					std::cerr << "Error! Failure in aligning the reads." << std::endl;
					exit(1);
				}

				if(MethodC::compress_C(alignment_file_1, alignment_file_2, output_file, genome_file)) {
					std::cerr << "Done compressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the compression!" << std::endl;

				system(("rm " + alignment_file_1).c_str());
				system(("rm " + alignment_file_2).c_str());

			}
			else {

				string input_file = string(argv[optind++]);

				string output_file_1 = string(argv[optind++]);
				string output_file_2 = string(argv[optind++]);

				if(MethodC::decompress_C(input_file, output_file_1, output_file_2, genome_file)) {
					std::cerr << "Done decompressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the decompression!" << std::endl;
			}

			break;
		}			

		case packing_mode_d:

			string genome_file = string(argv[optind++]); 

			if(xc_mode == zip_mode) {
				string input_file_1 = string(argv[optind++]);

				string input_file_2 = string(argv[optind++]);

				string output_file = string(argv[optind++]);

				string alignment_file_1 = input_file_1 + ".tab";
				string alignment_file_2 = input_file_2 + ".tab";

				std::string index = genome_file;

				// Call to align
				if(!(align_pair(input_file_1, input_file_2, index, alignment_file_1, alignment_file_2, read_mode, false))) {
					std::cerr << "Error! Failure in aligning the reads." << std::endl;
					exit(1);
				}

				if(MethodD::compress(alignment_file_1, alignment_file_2, output_file, genome_file)) {
					std::cerr << "Done compressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the compression!" << std::endl;
//				system(("rm " + alignment_file_1).c_str());
//				system(("rm " + alignment_file_2).c_str());

			}
			else {
				string input_file = string(argv[optind++]);

				string output_file_1 = string(argv[optind++]);
				string output_file_2 = string(argv[optind++]);

				if(MethodD::decompress(input_file, output_file_1, output_file_2, genome_file)) {
					std::cerr << "Done decompressing." << std::endl;
				}
				else
					std::cerr << "Error! Something went wrong with the decompression!" << std::endl;
			}
			break;

	}

}
