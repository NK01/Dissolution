#ifndef PHASE_H_MAY_2020
#define PHASE_H_MAY_2020

#include <vector>

class Phase
{
    
public:

    // The number of Solute Elements in the phase
    int numberOfSolutes;

    // The number of control volumes phase is divided into
	int numberOfControlVolumes;

    // Total length of the phase
	double lengthOfPhase;

    // Concentration gradient at the back position of phase (i.e. at node 0)
	std::vector<double> backGradient;

    // Concentration gradient at the front position of phase (i.e. at final node)
	std::vector<double> frontGradient;

    // Equilibrium Concentration of the solute at the back position of phase (i.e. at node 0) 
    std::vector<double> backEquilibConc;

    // Equilibrium Concentration of the solute at the front position of phase (i.e. at final node)
	std::vector<double> frontEquilibConc;

    //Diffusivity of the solute inside the phase
	std::vector<std::vector<double>> diffusivity;

    // The Solute Profile inside the phase
	std::vector<std::vector<double>> concentration;

    // The length of the single control volume
	std::vector<std::vector<double>> deltax;

    // Default Constructor
	Phase();

    // Constructor which takes number of control volumes as argument
	Phase(int nElements);

    // Constructor which takes number of control volumes and number of elements as argument
	Phase(int nElements, int nCont);

    // Constructor which takes number of elements, number of control volumes and length of the phase as argument
	Phase(int nElements, int nCont, double length);

    // Default destructor
	//~Phase();

    // Function to perform diffusion for all the elements inside the phase
	void Diffusion(double dt);
};

#endif // !PHASE_H_MAY_2020
