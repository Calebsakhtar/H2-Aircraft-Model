#ifndef AIRCRAFT_MODEL_H
#define AIRCRAFT_MODEL_H

#include <math.h>

#include "../headers/MathTools.h"

namespace AircraftModel {

	// Implements equation (15-20) from:
	// https://www.sciencedirect.com/topics/engineering/friction-drag-coefficient
	// Which is the equation for a fully turbulent boundary layer over a flat plate
	// for air with compressibility accounted for.
	// 
	// Re is the reynolds number of the free-stream (needs checking)
	// M is the Mach number of the free-stream flow
	double compute_turb_frict_coeff(const double& Re, const double& M);

}

#endif