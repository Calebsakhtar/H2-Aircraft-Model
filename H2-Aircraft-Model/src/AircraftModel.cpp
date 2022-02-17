
#include "../headers/AircraftModel.h"

namespace AircraftModel {

	double compute_turb_frict_coeff(const double& Re, const double& M) {
		// Implements equation (15-20) from:
		// https://www.sciencedirect.com/topics/engineering/friction-drag-coefficient
		// Which is the equation for a fully turbulent boundary layer over a flat plate
		// for air with compressibility accounted for
		// 
		// Re is the reynolds number of the free-stream (needs checking)
		// M is the Mach number of the free-stream flow

		double Cf_turb = 0.455;

		Cf_turb = Cf_turb / pow(log10(Re), 2.58);
		Cf_turb = Cf_turb / pow(1 + 0.144 * M * M, 0.65);

		return Cf_turb;
	}

}