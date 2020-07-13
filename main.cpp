#include "Phase.h"
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <chrono>

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
    auto start = std::chrono::steady_clock::now();

    // FileStream Output for writing to the file
	std::ofstream output;

    std::vector<double> processConditions = std::vector <double>(4);

    std::ifstream processConditionInput;
    processConditionInput.open("ProcessCondition.csv");

    if (processConditionInput.is_open())
    {
        int row{0};
        int column{0};
        std::string value{"undefined"};
        while(processConditionInput)
        {
            std::string buff;
            std::getline(processConditionInput, buff);

            column = 0;
            int initialPos{0};
            int finalPos{0};

            if (row <= processConditions.size())
            {
                while (finalPos != std::string::npos)
                {
                    finalPos = buff.find_first_of(",", finalPos + 1);
                    value = buff.substr(initialPos, (finalPos - initialPos));
                    if (column == 1)
                    {
                        processConditions[row] = std::stod(value);
                    }
                    initialPos = finalPos + 1;
                    column++;
                }
            }

            row++;
        }
    }
    else
    {
        std::cout << "The file ProcessCondition.csv was not found! Using default values\n";
    }
    
    processConditionInput.close();

    // Total Length
    double L{ processConditions[2] };

    // The FCC phase
    Phase FCC;

    // Length of FCC phase
    FCC.lengthOfPhase = (100 - processConditions[3]) * L / 100;

    // Reading concentration
    FCC.ReadConcentration("fcc_conc.csv");

    // Reading deltax
    FCC.ReadDeltax("fcc_deltax.csv");

    // Reading diffusivities
    FCC.ReadDiffusivities("fcc_diffusivities.csv");

    // Adding Equilibrium concentrations 
    std::ifstream equilibConcInput;

    equilibConcInput.open("fcc_EquilibConc.csv");

    if (equilibConcInput.is_open())
    {
        int row{0};
        int column{0};
        std::string value{"undefined"};
        while(equilibConcInput)
        {
            std::string buff;
            std::getline(equilibConcInput, buff);

            column = 0;
            int initialPos{0};
            int finalPos{0};
            if (row != 0 && row <= 1)
            {
                while (finalPos != std::string::npos && column < FCC.numberOfSolutes)
                {
                    finalPos = buff.find_first_of(",", finalPos + 1);
                    value = buff.substr(initialPos, (finalPos - initialPos));
                    FCC.frontEquilibConc[column] = std::stod(value);
                    initialPos = finalPos + 1;
                    column++;
                }
            }

            row++;
        }
    }
    else
    {
        std::cout << "The file fcc_EquilibConc.csv was not found!\n";
    }
    
    equilibConcInput.close();

    // The Laves phase
    Phase Laves;

    Laves.lengthOfPhase = L - FCC.lengthOfPhase;

    // Reading Concentration
    Laves.ReadConcentration("laves_conc.csv");

    // Reading deltax
    Laves.ReadDeltax("laves_deltax.csv");

    // Reading diffusivities
    Laves.ReadDiffusivities("laves_diffusivities.csv");

    // Adding Equilibrium concentrations    
    equilibConcInput.open("laves_EquilibConc.csv");

    if (equilibConcInput.is_open())
    {
        int row{0};
        int column{0};
        std::string value{"undefined"};
        while(equilibConcInput)
        {
            std::string buff;
            std::getline(equilibConcInput, buff);

            column = 0;
            int initialPos{0};
            int finalPos{0};
            if (row != 0 && row <= 1)
            {
                while (finalPos != std::string::npos && column < Laves.numberOfSolutes)
                {
                    finalPos = buff.find_first_of(",", finalPos + 1);
                    value = buff.substr(initialPos, (finalPos - initialPos));
                    Laves.frontEquilibConc[column] = std::stod(value);
                    initialPos = finalPos + 1;
                    column++;
                }
            }

            row++;
        }
    }
    else
    {
        std::cout << "The file laves_EquilibConc.csv was not found!\n";
    }
    
    equilibConcInput.close();

    double totalTime{ 0 };
    double dt{ 1 };
    double t{ 0 };
    double v{ 0 };
    double tempVelocity{ 0 };

    output.open("Interface_length.txt");

    totalTime = processConditions[0];
    dt = processConditions[1];
    while (t < totalTime)
    {
        if (FCC.lengthOfPhase > 0 && Laves.lengthOfPhase > 0)
		{
            v = -(-Laves.frontGradient[0] * Laves.lengthOfPhase);
			v = v - (-FCC.frontGradient[0] * FCC.lengthOfPhase);
            v = v / (Laves.frontEquilibConc[0] - FCC.frontEquilibConc[0]);

            for (int i = 1; i < FCC.numberOfSolutes; i++)
            {
                tempVelocity = -(-Laves.frontGradient[i] * Laves.lengthOfPhase);
			    tempVelocity = tempVelocity - (-FCC.frontGradient[i] * FCC.lengthOfPhase);
                tempVelocity = tempVelocity / (Laves.frontEquilibConc[i] - FCC.frontEquilibConc[i]);

                if (abs(tempVelocity) < abs(v))
                {
                    v = tempVelocity;
                }
                
            }
            
            FCC.SetLength(v * dt - 0);

			Laves.SetLength(0 - v*dt);

			FCC.Diffusion(dt, 0);
            Laves.Diffusion(dt, 0);
		}
		else
		{
			if (FCC.lengthOfPhase <= 0)
			{
				FCC.lengthOfPhase = 0;
				Laves.lengthOfPhase = L - FCC.lengthOfPhase;

				// Calculate Internal Diffusion
				Laves.Diffusion(dt, -1);
			}
			else
			{
				Laves.lengthOfPhase = 0;
				FCC.lengthOfPhase = L - Laves.lengthOfPhase;

				// Calculate Internal Diffusion
				FCC.Diffusion(dt, -1);
			}

		}

        // progressbar display
        if (static_cast<int>(t/dt) % static_cast<int>(60 * 5 / dt) == 0)
        {
            double progress = t / totalTime;
            int barWidth = 70;
            std::cout << "[";

            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) 
            {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << static_cast<int>(progress * 100.0) << " %\r" << std::flush;

            output << FCC.lengthOfPhase << std::endl;
        }

        t = t + dt;
    }

    output << FCC.lengthOfPhase << std::endl;
    output.close();

    double totalLength{0.0};

    // Writing final concentration profile of FCC
    output.open("FCC_Conc.txt");
    
    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        totalLength += FCC.deltax[0][i];
        output << totalLength << "\t";
        for (int j = 0; j < FCC.numberOfSolutes; j++)
        {
            output << FCC.concentration[j][i] << "\t\t";
        }
        
        output << "\n";
    }

    output.close();

    // Writing final concentration profile of Laves
    totalLength = 0.0;
    output.open("Laves_Conc.txt");
    
    for (int i = 0; i < Laves.numberOfControlVolumes; i++)
    {
        totalLength += Laves.deltax[0][i];
        output << totalLength << "\t";
        for (int j = 0; j < Laves.numberOfSolutes; j++)
        {
            output << Laves.concentration[j][i] << "\t\t";
        }
        
        output << "\n";
    }

    output.close();

    std::cout << "\nSimulation completed\n";

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    return 0;
}