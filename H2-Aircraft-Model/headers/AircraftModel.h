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

	// ISA Returns the International Standard Atmosphere(ISA) values for ambient
	// temperature(op_T) in degrees kelvin, speed of sound(op_a) in m / s, pressure
	// (op_P)in kPa and density(op_rho) in Kg / m ^ 3 at a given input array height in
	// kilometers(ip_h).
	void ISA(const double& ip_h, double& op_T, double& op_a, double& op_P, double& op_rho);

}

#endif