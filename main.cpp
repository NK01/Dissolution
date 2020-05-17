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

    // Total Length
    double L{100};

    // The FCC phase
    Phase FCC(6, 172);

    // Length of FCC phase
    FCC.lengthOfPhase = 96.2;

    // Reading concentration
    FCC.ReadConcentration("fcc_conc.csv");

    // Reading deltax
    FCC.ReadDeltax("fcc_deltax.csv");

    // Reading diffusivities
    FCC.ReadDiffusivities("fcc_diffusivities.csv");

    // Adding Equilibrium concentrations
    FCC.backEquilibConc[0] = 17.07;
    FCC.backEquilibConc[1] = 14.93;
    FCC.backEquilibConc[2] = 8.71;
    FCC.backEquilibConc[3] = 5.19;
    FCC.backEquilibConc[4] = 1.94;
    FCC.backEquilibConc[5] = 0.57023;

    // The Laves phase
    Phase Laves(6, 10);

    Laves.lengthOfPhase = L - FCC.lengthOfPhase;

    // Reading Concentration
    Laves.ReadConcentration("laves_conc.csv");

    // Reading deltax
    Laves.ReadDeltax("laves_deltax.csv");

    // Adding Equilibrium concentrations
    Laves.backEquilibConc[0] = 15.41;
    Laves.backEquilibConc[1] = 18.879;
    Laves.backEquilibConc[2] = 37.05;
    Laves.backEquilibConc[3] = 4.93;
    Laves.backEquilibConc[4] = 0.35;
    Laves.backEquilibConc[5] = 0.07;

    double totalTime{ 60 };
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
            v = v / (Laves.backEquilibConc[0] - FCC.backEquilibConc[0]);
            
            for (int i = 1; i < FCC.numberOfSolutes; i++)
            {
                tempVelocity = (-Laves.diffusivity[i][0] * Laves.backGradient[i] * Laves.lengthOfPhase);
			    tempVelocity = tempVelocity - (-FCC.diffusivity[i][0] * FCC.backGradient[i] * FCC.lengthOfPhase);
                tempVelocity = tempVelocity / (Laves.backEquilibConc[i] - FCC.backEquilibConc[i]);

                if (abs(tempVelocity) < abs(v))
                {
                    v = tempVelocity;
                }
                
            }

			FCC.lengthOfPhase = FCC.lengthOfPhase + (v * dt - 0);

			Laves.lengthOfPhase = Laves.lengthOfPhase + (0 - v * dt);

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
				Laves.lengthOfPhase = L - FCC.lengthOfPhase;

				//bcc.set_length(x1);
				//hcp.set_length(x1);
				//fcc.set_length(x2);

				// Calculate Internal Diffusion
				Laves.Diffusion(dt);
			}
			else
			{
				Laves.lengthOfPhase = 0;
				FCC.lengthOfPhase = L - Laves.lengthOfPhase;

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