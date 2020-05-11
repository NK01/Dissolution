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
    char diffusivityFilename[256];

    strcpy(filename, "input_profile.txt");
	strcpy(deltaFilename, "deltax_profile.txt");
    strcpy(diffusivityFilename, "diffusivity_profile.txt");

    Phase FCC;

    FCC.lengthOfPhase = 100e-6;

    // FileStream Input for reading from the file
	std::ifstream in;
    in.open(filename);

    in >> FCC.numberOfControlVolumes;

    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        in >> FCC.concentration[0][i];
    }

    in.close();

    in.open(deltaFilename);

    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        in >> FCC.deltax[0][i];
    }

    in.close();

    in.open(diffusivityFilename);

    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        in >> FCC.diffusivity[0][i];
    }

    in.close();

    out.open("Initial_Conc.txt");
    
    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        out << "C[" << i << "]: " << FCC.concentration[0][i] << "\n";
    }

    out.close();

    double totalTime = 60 * 60;
    double t = 0;
    double dt = 1e-2;

    while (t < totalTime)
    {
        FCC.Diffusion(dt);

        t = t + dt;
    }

    out.open("Final_Conc.txt");
    
    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        out << "C[" << i << "]: " << FCC.concentration[0][i] << "\n";
    }

    out.close();

    return 0;
}