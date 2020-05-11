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

void Phase::Diffusion(double dt, double y, int index)
{
    //concentration[index] = y;

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
                b[i] = 1 + 2 * diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i];
                c[i] = -2 * diffusivity[j][i + 1] * temp / deltax[j][i] / deltax[j][i];
                d[i] = concentration[j][i];
            }
            else if (i == numberOfSolutes - 1)
            {
                a[i] = -2 * diffusivity[j][i] * temp / deltax[j][i] / deltax[j][i];
                b[i] = 1 + 2 * diffusivity[j][i] * temp / deltax[j][i] / deltax[j][i];
                c[i] = 0;
                d[i] = concentration[j][i];
            }
            else
            {
                a[i] = -D[i] * temp / (deltax[i - 1] * (deltax[i - 1] + deltax[i]) / 2);
                b[i] = 1 + (D[i + 1] + D[i]) * temp * (1 / deltax[i - 1] + 1 / deltax[i]) / ((deltax[i - 1] + deltax[i]) / 2);
                c[i] = -D[i + 1] * temp / (deltax[i] * (deltax[i - 1] + deltax[i]) / 2);
                d[i] = concentration[j][i];
            }
        }

        // Coefficient calculations
        cnew[0] = c[0] / b[0];
        dnew[0] = d[0] / b[0];

        for (int i = 1; i < N + 2; i++)
        {
            cnew[i] = c[i] / (b[i] - a[i] * cnew[i - 1]);
            dnew[i] = (d[i] - a[i] * dnew[i - 1]) / (b[i] - a[i] * cnew[i - 1]);
        }

        // Back Substitution
        phi[N + 1] = dnew[N + 1];

        for (int i = N; i >= 0; i--)
        {
            phi[i] = dnew[i] - cnew[i] * phi[i + 1];
        }

        back_gradient = (3 * phi[N + 1] - 4 * phi[N] + phi[N - 1]) / 2 / db;
        front_gradient = (-3 * phi[0] + 4 * phi[1] - phi[2]) / 2 / db;

    }

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
	delete[] cnew;
	delete[] dnew;
}
