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

    // The FCC phase
    Phase FCC;

    // Length of FCC phase
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

    // The Laves phase
    Phase Laves;

    // Length of Laves phase
    Laves.lengthOfPhase = 100e-6;

    double totalTime{ 60 * 60 };
    double t{ 0 };
    double dt{ 1e-2 };
    double v{ 0 };

    while (t < totalTime)
    {
        if (FCC.lengthOfPhase > 0 && Laves.lengthOfPhase > 0)
		{
			v = (-Laves.diffusivity[0][0] * Laves.backGradient[0] * Laves.lengthOfPhase);
			v = v - (-FCC.diffusivity[0][0] * FCC.backGradient[0] * FCC.lengthOfPhase);
			//v = v / (y2 - y1);

			//x1old = x1;
			//x1 = x1 + (v * dt - 0);

			//x2old = x2;
			//x2 = x2 + (0 - v * dt);

			//bcc.set_length(x1);
			//hcp.set_length(x1);
			//fcc.set_length(x2);

			FCC.Diffusion(dt);
            Laves.Diffusion(dt);

			//hcp.add_precip(x1old);
		}
		else
		{
			if (FCC.lengthOfPhase <= 0)
			{
				FCC.lengthOfPhase = 0;
				//Laves.lengthOfPhase = L - x1;

				//bcc.set_length(x1);
				//hcp.set_length(x1);
				//fcc.set_length(x2);

				// Calculate Internal Diffusion
				Laves.Diffusion(dt);
			}
			else
			{
				Laves.lengthOfPhase = 0;
				//x1 = L - x2;

				//bcc.set_length(x1);
				//hcp.set_length(x1);
				//fcc.set_length(x2);

				// Calculate Internal Diffusion
				FCC.Diffusion(dt);
			}

		}

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