#ifndef AIRCRAFT_MODEL_H
#define AIRCRAFT_MODEL_H

#include <math.h>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "../headers/MathTools.h"

namespace AircraftModel {

	// ISA Returns the International Standard Atmosphere(ISA) values for ambient
	// temperature(op_T) in degrees kelvin, speed of sound(op_a) in m / s, pressure
	// (op_P)in kPa, density(op_rho) in Kg / m ^ 3, and dynamic viscosity (op_visc) 
	// in Kg/m/s at a given input height in kilometers(ip_h).
	//
	// Based on the following lecture slides: 
	// Engineering Tripos Part IIB, 4A7: Aircraft Aerodynamics and Design
	// Bill Dawes (with thanks to Dr Chez Hall), slides 20-22
	void ISA(const double& ip_h, double& op_T, double& op_a, double& op_P, double& op_rho,
		double& op_visc);
	
	// Calculates the IAS from the TAS and the altitude h
	// 
	// Input TAS must be in m/s
	// Input h must be in km
	// Output IAS is in m/s
	// 
	// This function adapts the following website:
	// https://aerotoolbox.com/airspeed-conversions/
	double TAS_to_IAS(const double& TAS, const double& h);

	// Calculates the TAS from the IAS and the altitude h
	// 
	// Input IAS must be in m/s
	// Input h must be in km
	// Output TAS is in m/s
	// 
	// This function adapts the following website:
	// https://aerotoolbox.com/airspeed-conversions/
	double IAS_to_TAS(const double& IAS, const double& h);
		
	// Implements equation (15-20) from:
	// https://www.sciencedirect.com/topics/engineering/friction-drag-coefficient
	// Which is the equation for a fully turbulent boundary layer over a flat plate
	// for air with compressibility accounted for.
	// 
	// Re is the reynolds number of the free-stream (needs checking)
	// M is the Mach number of the free-stream flow
	double compute_turb_frict_coeff(const double& Re, const double& M);

	// Computes a delta in the CD according to changes in the fuselage shape.
	//
	// This function adapts the method from http://wpage.unina.it/danilo.ciliberti/doc/Cusati.pdf
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
	// This information is also available in the ATR 72 FCOM: https://aviation-is.better-than.tv/atr72fcom.pdf
	// 
	// The tables assume a CG location of 25%.
	// 
	// The input "h" is height in km, which must lie above 2.44 km and below 7.61 km. The 
	// returned output "PW127_BSFC" is the break-specific fuel consumption in kg/J.
	double compute_cruise_BSFC_PW127(const double& h);

	// Implements the take-off run estimation from Fig 5.4 in the Raymer book.
	// Raymer, D. P., & American Institute of Aeronautics and Astronautics. (1989). 
	// Aircraft design: A conceptual approach. Washington, D.C: 
	// American Institute of Aeronautics and Astronautics
	//
	// The inputs are "MTOW" in tonnes, the wing area "S_wing" in m^2, the density ratio
	// relative to sea level "sl_rho_ratio", the take-off lift coefficient "CL_TO" and the
	// power of the engine in horsepower "BHP". The output "TO_run" is the take-off run in 
	// km.
	double compute_ground_run_raymer(const double& MTOW, const double& S_wing, const double&
		sl_rho_ratio, const double& CL_TO, const double& BHP);

	// Computes the propeller efficiency according to equation (13.15) from Raymer.
	// 
	// Raymer, D. P., & American Institute of Aeronautics and Astronautics. (1989). 
	// Aircraft design: A conceptual approach. Washington, D.C: 
	// American Institute of Aeronautics and Astronautics
	//
	// The input "thrust" is the Thrust in kN, the "TAS" is the True Arispeed in m/s, and
	// "P" is the power in kW. If the thrust of one engine is given, the power given must
	// also correspond to one engine. The output is the propeller efficiency.
	double compute_eta_prop_raymer(const double& thrust, const double& TAS, const double& P);

	// Computes the turboprop mass (in kg) for a given max. power requirement by using a
	// correlation from the data available in the following website under the 
	// specifications section:
	// https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_PW100
	// 
	// The input is the maximum power requirement of the engine in kW, and the output is
	// the likely engine mass in kg.
	double correl_turboprop_mass(const double& P_max);

	// Computes the turboprop Break-Specific Fuel Consumption at T/O for a given 
	// maximum power requirement by using a correlation from the data available in the 
	// following website under the specifications section:
	// https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_PW100
	// 
	// The input is the maximum power requirement of the engine in kW, and the output is
	// the likely BSFC in g/kWh.
	double correl_turboprop_TOBSFC(const double& P_max);

	// Performs a simple extrapolation to link the BSFC at TO with the BSFC at the chosen
	// cruise altitude h (which must be provided in km). Returns the BSFC at cruise with
	// the same units as the input units of the BSFC.
	double compute_new_engine_cruise_BSFC(const double& BSFC_TO, const double& h);

	// Given an input hydrogen fraction "H2_frac", a cruise altitude "h" in km, and a 
	// maximum power output for a single turboprop engine "P_max" in kW, return the hybrid
	// BSFC for a single engine in kg/J;
	// 
	// BSFC Equation: BSFC = Fuel Consumption/Power
	//
	// Please note that in this case "thermal efficiency" is the product of both cycle
	// and combustion efficiency.
	double calculate_hybrid_BSFC(const double& H2_frac, const double& h, const double& P_max);

	// Compute the total mass "op_calc_mass" in kg, the centre of gravity location "op_cg_loc"
	// in m, the payload mass "op_payload". It also states whether the volume and mass
	// constraints have been violated in "op_vio_vol" and "op_vio_mass" respectively. The nofuel
	// version of the outputs are those which assume a zero-fuel aircraft (they still include the
	// mass of hydrogen tanks). The mass of kerosene "op_M_JA1" in kg is also given.
	// 
	// The inputs are the mass of the engine "ip_M_engine" in kg, the ip TOTAL fuel mass 
	// "ip_M_fuel" in kg, and the H2 power fraction "ip_H2_frac" (power of hydrogen divided 
	// by power of kerosene).
	//
	// To compute the masses, first the volume of hydrogen needed is considered and the
	// remaining payload (to reach MTOW) is packed in the remaining space. This allows for the
	// aircraft to have a total mass that is LESS THAN MTOW, meaning that THIS FUNCTION NEEDS
	// TO BE ITERATED to achieve concordance to the payload fraction assumed in breguet and the
	// output of this program.
	//
	// This version loads a discrete number of passengers.
	bool compute_cg_loc_mass_pax(const double& ip_M_engine, const double& ip_M_fuel, const double& ip_H2_frac,
		double& op_cg_loc, double& op_calc_mass, double& op_cg_loc_nofuel, double& op_calc_mass_nofuel,
		double& op_payload, double& op_M_JA1, double& op_M_H2_net, int& op_num_pax, double& op_tank_l,
		bool& op_vio_mass, bool& op_vio_vol);

	// Compute the total mass "op_calc_mass" in kg, the centre of gravity location "op_cg_loc"
	// in m, the payload mass "op_payload". It also states whether the volume and mass
	// constraints have been violated in "op_vio_vol" and "op_vio_mass" respectively. The nofuel
	// version of the outputs are those which assume a zero-fuel aircraft (they still include the
	// mass of hydrogen tanks). The mass of kerosene "op_M_JA1" in kg is also given. The hydrogen
	// tank length "op_tank_l", mass of kerosene "op_M_JA1" total mass of hydrogen "op_M_H2_net",
	// and number of passengers on-board "op_num_pax" are also computed.
	// 
	// The inputs are the mass of the engine "ip_M_engine" in kg, the ip TOTAL fuel mass 
	// "ip_M_fuel" in kg, and the H2 power fraction "ip_H2_frac" (power of hydrogen divided 
	// by power of kerosene).
	//
	// To compute the masses, first the volume of hydrogen needed is considered and the
	// remaining payload (to reach MTOW) is packed in the remaining space. This allows for the
	// aircraft to have a total mass that is LESS THAN MTOW, meaning that THIS FUNCTION NEEDS
	// TO BE ITERATED to achieve concordance to the payload fraction assumed in breguet and the
	// output of this program.
	//
	// This version loads a continuous aircraft (passengers assumed to be evenly distributed)
	bool compute_cg_loc_mass(const double& ip_M_engine, const double& ip_M_fuel, const double& ip_H2_frac,
		double& op_cg_loc, double& op_calc_mass, double& op_cg_loc_nofuel, double& op_calc_mass_nofuel,
		double& op_payload, double& op_M_JA1, double& op_M_H2_net, int& op_num_pax, double& op_tank_l,
		bool& op_vio_mass, bool& op_vio_vol);

	// Modified version of compute_cg_loc_mass_pax that allows the weights to be distributed for an
	// aircraft off-design. The inputs are the "off-design" inputs.
	//
	// Compute the total mass "op_calc_mass" in kg, the centre of gravity location "op_cg_loc"
	// in m, the payload mass "op_payload". It also states whether the volume and mass
	// constraints have been violated in "op_vio_vol" and "op_vio_mass" respectively. The nofuel
	// version of the outputs are those which assume a zero-fuel aircraft (they still include the
	// mass of hydrogen tanks). The mass of kerosene "op_M_JA1" in kg is also given. The hydrogen
	// tank length "op_tank_l", mass of kerosene "op_M_JA1" total mass of hydrogen "op_M_H2_net",
	// and number of passengers on-board "op_num_pax" are also computed.
	// 
	// The inputs are the mass of the engine "ip_M_engine" in kg, the ip TOTAL fuel mass 
	// "ip_M_fuel" in kg, and the H2 power fraction "ip_H2_frac" (power of hydrogen divided 
	// by power of kerosene).
	//
	// To compute the masses, first the volume of hydrogen needed is considered and the
	// remaining payload (to reach MTOW) is packed in the remaining space. This allows for the
	// aircraft to have a total mass that is LESS THAN MTOW, meaning that THIS FUNCTION NEEDS
	// TO BE ITERATED to achieve concordance to the payload fraction assumed in breguet and the
	// output of this program.
	bool compute_cg_loc_mass_offdesign(const double& ip_M_engine, const double& ip_M_fuel, const double& ip_M_H2_design,
		const double& ip_H2_frac, const int& ip_num_seats, double& op_cg_loc, double& op_calc_mass,
		double& op_cg_loc_nofuel, double& op_calc_mass_nofuel, double& op_payload, double& op_M_JA1, double& op_M_H2_net,
		int& op_num_pax, bool& op_vio_mass, bool& op_vio_vol);
		
}

#endif
