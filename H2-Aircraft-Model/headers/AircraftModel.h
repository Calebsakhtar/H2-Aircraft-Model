#ifndef AIRCRAFT_MODEL_H
#define AIRCRAFT_MODEL_H

#include <math.h>

#include "../headers/MathTools.h"

namespace AircraftModel {

	// ISA Returns the International Standard Atmosphere(ISA) values for ambient
	// temperature(op_T) in degrees kelvin, speed of sound(op_a) in m / s, pressure
	// (op_P)in kPa, density(op_rho) in Kg / m ^ 3, and dynamic viscosity (op_visc) 
	// in Kg/m/s at a given input height in kilometers(ip_h).
	void ISA(const double& ip_h, double& op_T, double& op_a, double& op_P, double& op_rho,
		double& op_visc);
	
	// Implements equation (15-20) from:
	// https://www.sciencedirect.com/topics/engineering/friction-drag-coefficient
	// Which is the equation for a fully turbulent boundary layer over a flat plate
	// for air with compressibility accounted for.
	// 
	// Re is the reynolds number of the free-stream (needs checking)
	// M is the Mach number of the free-stream flow
	double compute_turb_frict_coeff(const double& Re, const double& M);

	double compute_delta_CD_fuselage(const double& current_d, const double& current_l,
		const double& next_d, const double& next_l, const double& M,
		const double& cruise_h);

	// Implements the Breguet Range Equation for prop aircraft, as seen in the following 
	// link: https://en.wikipedia.org/wiki/Range_(aeronautics)
	//
	// The inputs are propulsive efficiency (eta_prop), the Break Specific Fuel Consumption
	// (BSFC, in Kg/J), the Lift-Over-Drag ratio (L_D, assumed constant throughout the leg),
	// and the start and end total aircraft weights (w_start and w_end, both with the same 
	// units).
	//
	// The output is the range in m.
	double breguet_prop_range(const double& eta_prop, const double& BSFC, const double& L_D,
		const double& w_start, const double& w_end);

	// Implements the Breguet Range Equation for prop aircraft, as seen in the following 
	// link: https://en.wikipedia.org/wiki/Range_(aeronautics)
	//
	// The inputs are propulsive efficiency (eta_prop), the Break Specific Fuel Consumption
	// (BSFC, in Kg/J), the Lift-Over-Drag ratio (L_D, assumed constant throughout the leg),
	// and the range (range, in m).
	//
	// The output is the ratio of the weight of the aircraft at the start relative to the
	// end.
	double breguet_prop_wratio(const double& eta_prop, const double& BSFC, const double& L_D,
		const double& range);

	// Implements the engine performance table at ISA conditions given in the following link:
	// https://www.quora.com/At-cruise-speed-do-turboprops-run-at-their-maximal-rated-horse-power-If-not-how-much-less-is-that-given-horse-power-typically
	//
	// The input "h" is height in km, which must lie above 2.44 km and below 7.61 km. The 
	// returned output "PW127_BSFC" is the break-specific fuel consumption in kg/J.
	double compute_cruise_BSFC_PW127(const double& h);

	// Implements the take-off run estimation from Fig 5.4 in the Raymer book.
	//
	// The inputs are "MTOW" in tonnes, the wing area "S_wing" in m^2, the density ratio
	// relative to sea level "sl_rho_ratio", the take-off lift coefficient "CL_TO" and the
	// power of the engine in horsepower "BHP". The output "TO_run" is the take-off run in 
	// km.
	double compute_ground_run_raymer(const double& MTOW, const double& S_wing, const double&
		sl_rho_ratio, const double& CL_TO, const double& BHP);

	// Computes the propeller efficiency according to equation (13.15) from Raymer.
	//
	// The input "thrust" is the Thrust in kN, the "TAS" is the True Arispeed in m/s, and
	// "P" is the power in kW. If the thrust of one engine is given, the power given must
	// also correspond to one engine. The output is the propeller efficiency.
	double compute_eta_prop_raymer(const double& thrust, const double& TAS, const double& P);
}

#endif
