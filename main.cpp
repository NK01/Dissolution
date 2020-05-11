#include "Phase.h"
#include <iostream>
#include <fstream>
#include <string.h>

/*  This code will calculate dissolution in the two phase 
    system. We will utilize moving boundary algorithm
    and Concentration dependent Diffusivity.

    C++ standard used = c++17
    It is cross-platform capable code
    compiler used = clang (llvm 10.0)

    Copyright = Nishant Kumar (2020)
*/
int main()
{
    // FileStream Output for writing to the file
	std::ofstream out;

    // used to store the name of files for input or output
	char filename[256];
	char deltaFilename[256];

    // the number of partitions of the space grid
	int N{ (int)1e3 };

    strcpy(filename, "input_profile.txt");
	strcpy(deltaFilename, "deltax_profile.txt");


    return 0;
}