#include "Phase.h"
#include <fstream>
#include <iostream>

// Default Constructor
Phase::Phase()
{
    numberOfSolutes = 1;
	numberOfControlVolumes = 3;
    lengthOfPhase = 100; // in microns
    xcumulative = 0;

    backGradient = std::vector <double>(numberOfSolutes);
    frontGradient = std::vector <double>(numberOfSolutes);
    backEquilibConc = std::vector <double>(numberOfSolutes);
    frontEquilibConc = std::vector <double>(numberOfSolutes);

    diffusivity = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    concentration = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    deltax = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));

    for (int i = 0; i < numberOfSolutes; i++)
    {
        for (int j = 0; j < numberOfControlVolumes; j++)
        {
            deltax[i][j] = 1.0;
        }
        
    }
}

// Constructor which takes number of control volumes as argument
Phase::Phase(int nElements)
{
    numberOfSolutes = nElements;
	numberOfControlVolumes = 3;
    lengthOfPhase = 100; // in microns
    xcumulative = 0;

    backGradient = std::vector <double>(numberOfSolutes);
    frontGradient = std::vector <double>(numberOfSolutes);
    backEquilibConc = std::vector <double>(numberOfSolutes);
    frontEquilibConc = std::vector <double>(numberOfSolutes);

    diffusivity = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    concentration = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    deltax = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    
    for (int i = 0; i < numberOfSolutes; i++)
    {
        for (int j = 0; j < numberOfControlVolumes; j++)
        {
            deltax[i][j] = 1.0;
        }
        
    }
}

// Constructor which takes number of control volumes and number of elements as argument
Phase::Phase(int nElements, int nCont)
{
    numberOfSolutes = nElements;
	numberOfControlVolumes = nCont;
    lengthOfPhase = 100; // in microns
    xcumulative = 0;

    backGradient = std::vector <double>(numberOfSolutes);
    frontGradient = std::vector <double>(numberOfSolutes);
    backEquilibConc = std::vector <double>(numberOfSolutes);
    frontEquilibConc = std::vector <double>(numberOfSolutes);

    diffusivity = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    concentration = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    deltax = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    
    for (int i = 0; i < numberOfSolutes; i++)
    {
        for (int j = 0; j < numberOfControlVolumes; j++)
        {
            deltax[i][j] = 1.0;
        }
    }
}

// Function to read concentration profile
void Phase::ReadConcentration(std::string concFilename)
{
    std::ifstream concInput;
    concInput.open(concFilename);

    if (concInput.is_open())
    {
        // number of rows in the excel file
        int rows{0};

        // number of column in the excel file
        int columns{0};

        // This loop will calculate the number of rows and columns in the excel file
        while(concInput)
        {
            std::string buff;
            std::getline(concInput, buff);
            int foundPos{0};
            while (foundPos != std::string::npos)
            {
                foundPos = buff.find_first_of(",", foundPos + 1);
                if (rows == 0)
                {
                    columns++;
                }
                
            }

            rows++;
        }
        rows = rows - 2;

        if (rows >=3 && columns >= 1)
        {
            numberOfSolutes = columns;
            numberOfControlVolumes = rows;

            backGradient = std::vector <double>(numberOfSolutes);
            frontGradient = std::vector <double>(numberOfSolutes);
            backEquilibConc = std::vector <double>(numberOfSolutes);
            frontEquilibConc = std::vector <double>(numberOfSolutes);

            diffusivity = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
            concentration = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
            deltax = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));

            for (int i = 0; i < numberOfSolutes; i++)
            {
                for (int j = 0; j < numberOfControlVolumes; j++)
                {
                    deltax[i][j] = 1.0;
                }
            }
        }
    }
    else
    {
        std::cout << "The file " << concFilename << " was not found!\n";
    }
    

    concInput.close();

    concInput.open(concFilename);

    if (concInput.is_open())
    {
        int row{0};
        int column{0};
        std::string value{"undefined"};
        while(concInput)
        {
            std::string buff;
            std::getline(concInput, buff);

            column = 0;
            int initialPos{0};
            int finalPos{0};
            if (row != 0 && row <= numberOfControlVolumes)
            {
                while (finalPos != std::string::npos && column < numberOfSolutes)
                {
                    finalPos = buff.find_first_of(",", finalPos + 1);
                    value = buff.substr(initialPos, (finalPos - initialPos));
                    concentration[column][row-1] = std::stod(value);
                    initialPos = finalPos + 1;
                    column++;
                }
            }

            row++;
        }
    }
    
    concInput.close();
}

// Function to read deltaX profile
void Phase::ReadDeltax(std::string deltaxFilename)
{
    std::ifstream deltaxInput;

    deltaxInput.open(deltaxFilename);

    if (deltaxInput.is_open())
    {
        int row{ 0 };
        int column{ 0 };
        std::string value{"undefined"};
        while(deltaxInput)
        {
            std::string buff;
            std::getline(deltaxInput, buff);

            column = 0;
            int initialPos{0};
            int finalPos{0};
            if (row != 0 && row <= numberOfControlVolumes)
            {
                while (finalPos != std::string::npos && column < numberOfSolutes)
                {
                    finalPos = buff.find_first_of(",", finalPos + 1);
                    value = buff.substr(initialPos, (finalPos - initialPos));
                    deltax[column][row-1] = std::stod(value);
                    initialPos = finalPos + 1;
                    column++;
                }
            }

            row++;
        }
    }
    else
    {
        std::cout << "The file " << deltaxFilename << " was not found!\n";
    }
    
    deltaxInput.close();
}

// Function to read diffusivity data
void Phase::ReadDiffusivities(std::string diffusivityFilename)
{
    std::ifstream diffusivityInput;

    diffusivityInput.open(diffusivityFilename);

    if (diffusivityInput.is_open())
    {
        int row{0};
        int column{0};
        std::string value{"undefined"};
        while(diffusivityInput)
        {
            std::string buff;
            std::getline(diffusivityInput, buff);

            column = 0;
            int initialPos{0};
            int finalPos{0};
            if (row != 0 && row <= numberOfControlVolumes)
            {
                while (finalPos != std::string::npos && column < numberOfSolutes)
                {
                    finalPos = buff.find_first_of(",", finalPos + 1);
                    value = buff.substr(initialPos, (finalPos - initialPos));
                    diffusivity[column][row-1] = std::stod(value);
                    initialPos = finalPos + 1;
                    column++;
                }
            }

            row++;
        }
    }
    else
    {
        std::cout << "The file " << diffusivityFilename << " was not found!\n";
    }
    
    diffusivityInput.close();
}

// Function to perform diffusion
void Phase::Diffusion(double dt, int index)
{
    double temp = dt;

	std::vector<double> a = std::vector <double>(numberOfControlVolumes);
	std::vector<double> b = std::vector <double>(numberOfControlVolumes);
	std::vector<double> c = std::vector <double>(numberOfControlVolumes);
	std::vector<double> d = std::vector <double>(numberOfControlVolumes);
	std::vector<double> cnew = std::vector <double>(numberOfControlVolumes);
	std::vector<double> dnew = std::vector <double>(numberOfControlVolumes);

    for (int j = 0; j < numberOfSolutes; j++)
    {
        for (int i = 0; i < numberOfControlVolumes; i++)
        {
            if (i == 0)
            {
                a[i] = 0;
                b[i] = 1 + 2 * diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i];
                c[i] = -2 * diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i];
                d[i] = concentration[j][i];
            }
            else if (i == numberOfControlVolumes - 1)
            {
                a[i] = -2 * diffusivity[j][i - 1] * temp / deltax[j][i] / deltax[j][i];
                b[i] = 1 + 2 * diffusivity[j][i - 1] * temp / deltax[j][i] / deltax[j][i];
                c[i] = 0;
                d[i] = concentration[j][i];
            }
            else
            {
                a[i] = -diffusivity[j][i - 1] * temp / deltax[j][i] / deltax[j][i];
                b[i] = 1 + (diffusivity[j][i - 1] + diffusivity[j][i + 1]) * temp / deltax[j][i] / deltax[j][i];
                c[i] = -diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i];
                d[i] = concentration[j][i];
            }
        }

        // Coefficient calculations
        cnew[0] = c[0] / b[0];
        dnew[0] = d[0] / b[0];

        for (int i = 1; i < numberOfControlVolumes; i++)
        {
            cnew[i] = c[i] / (b[i] - a[i] * cnew[i - 1]);
            dnew[i] = (d[i] - a[i] * dnew[i - 1]) / (b[i] - a[i] * cnew[i - 1]);
        }

        // Back Substitution
        concentration[j][numberOfControlVolumes - 1] = dnew[numberOfControlVolumes - 1];

        for (int i = numberOfControlVolumes - 2; i >= 0; i--)
        {
            concentration[j][i] = dnew[i] - cnew[i] * concentration[j][i + 1];
        }

        if (index != -1)
        {
            for (int i = 0; i < numberOfSolutes; i++)
            {
                concentration[i][numberOfControlVolumes - 1] = frontEquilibConc[i];
            }
        }

        frontGradient[j] = (3 * concentration[j][numberOfControlVolumes - 1] - 4 * concentration[j][numberOfControlVolumes - 2] + concentration[j][numberOfControlVolumes - 3]) / 2 / deltax[j][numberOfControlVolumes - 1] / lengthOfPhase;
        backGradient[j] = (-3 * concentration[j][0] + 4 * concentration[j][1] - concentration[j][2]) / 2 / deltax[j][0] / lengthOfPhase;

    }
}

// Adjusting the mesh according to new length
void Phase::SetLength(double dx)
{
    xcumulative += dx;
    lengthOfPhase += dx;

    if (xcumulative > 0)
    {
        while (xcumulative > deltax[0][numberOfControlVolumes - 1])
        {
            for (int i = 0; i < numberOfSolutes; i++)
            {
                concentration[i].push_back(backEquilibConc[i]);
                deltax[i].push_back(deltax[0][numberOfControlVolumes - 1]);
                diffusivity[i].push_back(diffusivity[i][numberOfControlVolumes - 1]);
            }
            
            xcumulative -= deltax[0][0] / 2;
            numberOfControlVolumes += 1;
        }
    }

    if (xcumulative < 0)
    {
        while (abs(xcumulative) > deltax[0][numberOfControlVolumes - 1])
        {
            if (numberOfControlVolumes > 3)
            {
                for (int i = 0; i < numberOfSolutes; i++)
                {
                    concentration[i].pop_back();
                    deltax[i].pop_back();
                    diffusivity[i].pop_back();
                }
                numberOfControlVolumes -= 1;
            }
            else
            {   
                if (lengthOfPhase > 0)
                {
                    for (int i = 0; i < numberOfSolutes; i++)
                    {
                        deltax[i][0] = lengthOfPhase / 3;
                        deltax[i][1] = lengthOfPhase / 3;
                        deltax[i][2] = lengthOfPhase / 3;
                    }
                }
            }
            xcumulative += deltax[0][0] / 2;
        }
    }
    
}