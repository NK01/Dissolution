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

    // The FCC phase
    Phase FCC(6, 172);

    // Length of FCC phase
    FCC.lengthOfPhase = 96.2;

    // Reading concentration
    FCC.ReadConcentration("fcc_conc.csv");

    // Reading deltax
    FCC.ReadDeltax("fcc_deltax.csv");

    // The Laves phase
    Phase Laves(6, 10);

    Laves.lengthOfPhase = 100 - FCC.lengthOfPhase;

    // Reading Concentration
    Laves.ReadConcentration("laves_conc.csv");

    // Reading deltax
    Laves.ReadDeltax("laves_deltax.csv");

    double totalTime{ 1.2e-2 };
    double t{ 0 };
    double dt{ 1e-2 };
    double v{ 0 };
    double tempVelocity{ 0 };

    while (t < totalTime)
    {
        if (FCC.lengthOfPhase > 0 && Laves.lengthOfPhase > 0)
		{
            v = (-Laves.diffusivity[0][0] * Laves.backGradient[0] * Laves.lengthOfPhase);
			v = v - (-FCC.diffusivity[0][0] * FCC.backGradient[0] * FCC.lengthOfPhase);
            
            for (int i = 1; i < FCC.numberOfSolutes; i++)
            {
                tempVelocity = (-Laves.diffusivity[i][0] * Laves.backGradient[i] * Laves.lengthOfPhase);
			    tempVelocity = tempVelocity - (-FCC.diffusivity[i][0] * FCC.backGradient[i] * FCC.lengthOfPhase);

                if (abs(tempVelocity) < abs(v))
                {
                    v = tempVelocity;
                }
                
            }
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
        out << "C[" << i << "]: \t";
        for (int j = 0; j < FCC.numberOfSolutes; j++)
        {
            out << FCC.concentration[j][i] << "\t\t";
        }
        
        out << "\n";
    }

    out.close();

    return 0;
}