#include "Phase.h"
#include <iostream>
#include <fstream>
#include <string>

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
    double L{100e-6};

    // The FCC phase
    Phase FCC(6, 172);

    // Length of FCC phase
    FCC.lengthOfPhase = 96.2e-6;

    // Reading concentration
    FCC.ReadConcentration("fcc_conc.csv");

    // Reading deltax
    FCC.ReadDeltax("fcc_deltax.csv");

    // Reading diffusivities
    FCC.ReadDiffusivities("fcc_diffusivities.csv");

    // Adding Equilibrium concentrations
    FCC.frontEquilibConc[0] = 17.07;
    FCC.frontEquilibConc[1] = 14.93;
    FCC.frontEquilibConc[2] = 8.71;
    FCC.frontEquilibConc[3] = 5.19;
    FCC.frontEquilibConc[4] = 1.94;
    FCC.frontEquilibConc[5] = 0.57023;

    // The Laves phase
    Phase Laves(6, 10);

    Laves.lengthOfPhase = L - FCC.lengthOfPhase;

    // Reading Concentration
    Laves.ReadConcentration("laves_conc.csv");

    // Reading deltax
    Laves.ReadDeltax("laves_deltax.csv");

    // Adding Equilibrium concentrations
    Laves.frontEquilibConc[0] = 15.41;
    Laves.frontEquilibConc[1] = 18.879;
    Laves.frontEquilibConc[2] = 37.05;
    Laves.frontEquilibConc[3] = 4.93;
    Laves.frontEquilibConc[4] = 0.35;
    Laves.frontEquilibConc[5] = 0.07;

    double totalTime{ 60 * 60};
    double t{ 0 };
    double dt{ 1e-2 };
    double v{ 0 };
    double tempVelocity{ 0 };

    while (t < totalTime)
    {
        if (FCC.lengthOfPhase > 0 && Laves.lengthOfPhase > 0)
		{
            v = (-Laves.diffusivity[0][Laves.numberOfControlVolumes - 1] * Laves.frontGradient[0] * Laves.lengthOfPhase);
			v = v - (-FCC.diffusivity[0][FCC.numberOfControlVolumes - 1] * FCC.frontGradient[0] * FCC.lengthOfPhase);
            v = v / (Laves.frontEquilibConc[0] - FCC.frontEquilibConc[0]);
            
            for (int i = 1; i < FCC.numberOfSolutes; i++)
            {
                tempVelocity = (-Laves.diffusivity[i][Laves.numberOfControlVolumes - 1] * Laves.frontGradient[i] * Laves.lengthOfPhase);
			    tempVelocity = tempVelocity - (-FCC.diffusivity[i][FCC.numberOfControlVolumes - 1] * FCC.frontGradient[i] * FCC.lengthOfPhase);
                tempVelocity = tempVelocity / (Laves.frontEquilibConc[i] - FCC.frontEquilibConc[i]);

                if (abs(tempVelocity) < abs(v))
                {
                    v = tempVelocity;
                }
                
            }

            FCC.SetLength(v * dt - 0);

			Laves.SetLength(0 - v*dt);

			FCC.Diffusion(dt);
            Laves.Diffusion(dt);

            for (int i = 0; i < FCC.numberOfSolutes; i++)
            {
                FCC.concentration[i][FCC.numberOfControlVolumes - 1] = FCC.frontEquilibConc[i];
            }

            for (int i = 0; i < Laves.numberOfSolutes; i++)
            {
                Laves.concentration[i][Laves.numberOfControlVolumes - 1] = Laves.frontEquilibConc[i];
            }
		}
		else
		{
			if (FCC.lengthOfPhase <= 0)
			{
				FCC.lengthOfPhase = 0;
				Laves.lengthOfPhase = L - FCC.lengthOfPhase;

				// Calculate Internal Diffusion
				Laves.Diffusion(dt);
			}
			else
			{
				Laves.lengthOfPhase = 0;
				FCC.lengthOfPhase = L - Laves.lengthOfPhase;

				// Calculate Internal Diffusion
				FCC.Diffusion(dt);
			}

		}

        t = t + dt;
    }

    // Writing final concentration profile of FCC
    out.open("FCC_Conc.txt");
    
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

    // Writing final concentration profile of Laves
    out.open("Laves_Conc.txt");
    
    for (int i = 0; i < Laves.numberOfControlVolumes; i++)
    {
        out << "C[" << i << "]: \t";
        for (int j = 0; j < Laves.numberOfSolutes; j++)
        {
            out << Laves.concentration[j][i] << "\t\t";
        }
        
        out << "\n";
    }

    out.close();

    out.open("Interface_length.txt");

    out << "Final length of FCC: " << FCC.lengthOfPhase;

    out.close();

    return 0;
}