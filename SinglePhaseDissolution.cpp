#include "Phase.h"
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <chrono>

/*  This code will calculate dissolution of multi-component 
    system. We will utilize single phase diffusion mechanism
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
        std::cout << "The file ProcessCondition.csv was not found!\n"; 
        std::cout << "Using default values\n";
        for (int i = 0; i < 4; i++)
        {
            processConditions[i] = 0.0;
        }
        
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
        std::cout << "Using Default values (0.0) !\n";
        for (int i = 0; i < FCC.numberOfSolutes; i++)
        {
            FCC.frontEquilibConc[i] = 0.0;
            FCC.backEquilibConc[i] = 0.0;
        }
        
    }
    
    equilibConcInput.close();

    double totalTime{ 0 };
    double dt{ 1 };
    double t{ 0 };
    double v{ 0 };
    double tempVelocity{ 0 };

    // Writing initial concentration profile of FCC
    output.open("FCC_initial_Conc.txt");
    double currentLength{0.0};

    for (int i = 0; i < FCC.numberOfControlVolumes; i++)
    {
        currentLength += FCC.deltax[0][i];
        output << currentLength << "\t";
        for (int j = 0; j < FCC.numberOfSolutes; j++)
        {
            output << FCC.concentration[j][i] << "\t\t";
        }
        
        output << "\n";
    }

    output.close();

    output.open("Phase_Fraction.txt");

    double pf{0.0};
    double cvw{0.0};
    for (int i = 0; i < FCC.numberOfSolutes; i++)
    {
        for (int j = 0; j < FCC.numberOfControlVolumes; j++)
        {
            if (FCC.concentration[i][j] > 10.54)
            {
                cvw = FCC.concentration[i][j] - 10.54;
                cvw = cvw / 0.7269;
                pf += cvw*FCC.deltax[i][j]/100;
            }
        }
        pf = pf / FCC.lengthOfPhase;
    }
    output << pf << std::endl;

    totalTime = processConditions[0];
    dt = processConditions[1];
    while (t < totalTime)
    {
        if (FCC.lengthOfPhase <= 0)
		{
            std::cout << "The primary phase amount is zero!\n";
            break;
		}
		else
		{
			FCC.lengthOfPhase = L;
			// Calculate Internal Diffusion
			FCC.Diffusion(dt, -1);
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

            for (int i = 0; i < FCC.numberOfSolutes; i++)
            {
                pf = 0.0;
                for (int j = 0; j < FCC.numberOfControlVolumes; j++)
                {
                    if (FCC.concentration[i][j] > 10.7)
                    {
                        cvw = FCC.concentration[i][j] - 10.7;
                        cvw = cvw / 0.7269;
                        pf += cvw*FCC.deltax[i][j]/100;
                    }
                }
                pf = pf / FCC.lengthOfPhase;
            }
            output << pf << std::endl;
        }

        t = t + dt;
    }

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

    std::cout << "\nSimulation completed\n";

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    return 0;
}