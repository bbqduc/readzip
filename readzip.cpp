#include "Alignment.h"
#include "AlignmentReader.h"
#include "MethodA.h"
#include "MethodB.h"
#include "MethodC.h"
#include "MethodD.h"

#include <getopt.h>

enum packing_mode_t {packing_mode_a, packing_mode_b, packing_mode_c, packing_mode_d };
enum pack_unpack_mode_t {zip_mode, unzip_mode };

int main(int argc, char **argv) 
{

	packing_mode_t mode;
	pack_unpack_mode_t xc_mode;

	// Parse command line parameters
	int option_index = 0;
	int c;
	while((c = getopt(argc, argv, "abcdxo")) != -1)

	{

	// TODO sanity checks that no two conflicting options
	switch(c)
		{
		case 'a':
			mode = packing_mode_a; break;
		case 'b':
			mode = packing_mode_b; break;
		case 'c':
			mode = packing_mode_c; break;
		case 'd':
			mode = packing_mode_d; break;
		case 'x':
			xc_mode = unzip_mode; break;
		case 'o':
			xc_mode = zip_mode; break;
		}
	}


	// Parse filenames
	if (argc - optind == 0)
	{
		cerr << "readzip: missing input files!" << endl;
		return 1;
	}

	string first_file = string(argv[optind++]);
	string second_file;

	if(mode == packing_mode_c | mode == packing_mode_d) {
		if (argc - optind == 0)
		{
			cerr << "readzip: missing input files!" << endl;
			return 1;
		}
		
		second_file = string(argv[optind++]);
	}


	// TODO Calling each mode (fill up as modes are implemented)


}
