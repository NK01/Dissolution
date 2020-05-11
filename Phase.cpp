#include "Phase.h"

// Default Constructor
Phase::Phase()
{
    numberOfSolutes = 5;
	numberOfControlVolumes = 1000;
    lengthOfPhase = 100; // in microns

    backGradient = std::vector <double>(numberOfSolutes);
    frontGradient = std::vector <double>(numberOfSolutes);
    backEquilibConc = std::vector <double>(numberOfSolutes);
    frontEquilibConc = std::vector <double>(numberOfSolutes);

    diffusivity = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    concentration = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    deltax = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    
}

// Constructor which takes number of control volumes and number of elements as argument
Phase::Phase(int nElements, int nCont)
{
    numberOfSolutes = nElements;
	numberOfControlVolumes = nCont;
    lengthOfPhase = 100; // in microns

    backGradient = std::vector <double>(numberOfSolutes);
    frontGradient = std::vector <double>(numberOfSolutes);
    backEquilibConc = std::vector <double>(numberOfSolutes);
    frontEquilibConc = std::vector <double>(numberOfSolutes);

    diffusivity = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    concentration = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    deltax = std::vector<std::vector<double>>(numberOfSolutes, std::vector<double>(numberOfControlVolumes));
    
}

void Phase::Diffusion(double dt)
{
	double *a;
    double *b;
    double *c;
    double *d;
    double *cnew;
    double *dnew;

    double temp = dt / lengthOfPhase / lengthOfPhase;

	a = new double[numberOfControlVolumes];
	b = new double[numberOfControlVolumes];
	c = new double[numberOfControlVolumes];
	d = new double[numberOfControlVolumes];
	cnew = new double[numberOfControlVolumes];
	dnew = new double[numberOfControlVolumes];

    for (int j = 0; j < numberOfSolutes; j++)
    {
        for (int i = 0; i < numberOfControlVolumes; i++)
        {
            if (i == 0)
            {
                a[i] = 0;
                b[i] = 1 + 2 * diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i] / lengthOfPhase / lengthOfPhase;
                c[i] = -2 * diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i] / lengthOfPhase / lengthOfPhase;
                d[i] = concentration[j][i];
            }
            else if (i == numberOfSolutes - 1)
            {
                a[i] = -2 * diffusivity[j][i] * temp / deltax[j][i] / deltax[j][i] / lengthOfPhase / lengthOfPhase;
                b[i] = 1 + 2 * diffusivity[j][i] * temp / deltax[j][i] / deltax[j][i] / lengthOfPhase / lengthOfPhase;
                c[i] = 0;
                d[i] = concentration[j][i];
            }
            else
            {
                a[i] = -diffusivity[j][i] * temp / (deltax[j][i - 1] * (deltax[j][i - 1] + deltax[j][i]) / 2) / lengthOfPhase / lengthOfPhase;
                b[i] = 1 + (diffusivity[j][i + 1] + diffusivity[j][i]) * temp * (1 / deltax[j][i - 1] + 1 / deltax[j][i]) / ((deltax[j][i - 1] + deltax[j][i]) / 2) / lengthOfPhase / lengthOfPhase;
                c[i] = -diffusivity[j][i + 1] * temp / (deltax[j][i] * (deltax[j][i - 1] + deltax[j][i]) / 2) / lengthOfPhase / lengthOfPhase;
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

        for (int i = numberOfControlVolumes - 1; i >= 0; i--)
        {
            concentration[j][i] = dnew[i] - cnew[i] * concentration[j][i + 1];
        }

        backGradient[j] = (3 * concentration[j][numberOfControlVolumes - 1] - 4 * concentration[j][numberOfControlVolumes - 2] + concentration[j][numberOfControlVolumes - 3]) / 2 / deltax[j][numberOfControlVolumes - 1] / lengthOfPhase;
        frontGradient[j] = (-3 * concentration[j][0] + 4 * concentration[j][1] - concentration[j][2]) / 2 / deltax[j][0] / lengthOfPhase;

    }

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
	delete[] cnew;
	delete[] dnew;
}
