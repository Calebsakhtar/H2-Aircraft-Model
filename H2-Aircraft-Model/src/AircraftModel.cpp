
#include "../headers/AircraftModel.h"

namespace AircraftModel {

	void ISA(const double& ip_h, double& op_T, double& op_a, double& op_P, double& op_rho,
		double& op_visc) {
		// ISA Returns the International Standard Atmosphere(ISA) values for ambient
		// temperature(op_T) in degrees kelvin, speed of sound(op_a) in m / s, pressure
		// (op_P)in kPa, density(op_rho) in Kg / m ^ 3, and dynamic viscosity (op_visc) 
		// in Kg/m/s at a given input height in kilometers(ip_h).
	
		// Initialize the relevant quantities at sea level. Viscosity data from:
		// http://www.aerodynamics4students.com/properties-of-the-atmosphere/sea-level-conditions.php
		const double T_sl = 288.15;
		const double a_sl = 340.3;
		const double P_sl = 101.325;
		const double rho_sl = 1.225;
		const double visc_sl = 1.789e-5;

		// Set the output quantities to equal to 0 (as initialization)
		op_T = 0;
		op_a = 0;
		op_P = 0;
		op_rho = 0;
		op_visc = 0;

		if (ip_h <= 0) {
			// Set the output quantities to equal to the sea level amount if below sea level
			op_T = T_sl;
			op_a = a_sl;
			op_P = P_sl;
			op_rho = rho_sl;
			op_visc = visc_sl;
		}
		else {

			// Initialize the temperature ratio
			double T_ratio = 1;

			if (ip_h < 11) {
				// Calculate ISA values below the tropopause
				op_T = 288.15 - 6.5 * ip_h;
				T_ratio = op_T / T_sl;

				op_a = a_sl * pow(T_ratio, 0.5);
				op_P = P_sl * pow(T_ratio, 5.256);
				op_rho = rho_sl * pow(T_ratio, 4.256);
			}
			else{
				// Calculate ISA values above the tropopause
				op_T = 216.65;
				T_ratio = op_T / T_sl;
				op_a = a_sl * pow(T_ratio, 0.5);
				op_P = P_sl * pow(T_ratio, 5.256);
				op_rho = rho_sl * pow(T_ratio, 4.256);

				double factor = exp(-0.1577 * (ip_h - 11));
				op_P = op_P * factor;
				op_rho = op_rho * factor;
			}

			// Sutherland's Law (from https://www.cfd-online.com/Wiki/Sutherland%27s_law)
			const double C1 = 1.458e-6;
			const double S = 110.4;

			op_visc = C1 * pow(op_T, 1.5) / (op_T + S);
		}
	}

	double TAS_to_IAS(const double& TAS, const double& h) {
		// Calculates the IAS from the TAS and the altitude h
		// 
		// Input TAS must be in m/s
		// Input h must be in km
		// Output IAS is in m/s
		// This function adapts the following website:
		// https://aerotoolbox.com/airspeed-conversions/

		// Initialize the outputs from ISA
		double T = 0;
		double a = 0;
		double P = 0;
		double rho = 0;
		double visc = 0;

		// Gamma Air
		double gamma = 1.4;
		const double rho_sl = 1.225;
		
		// Initialize the function output
		double IAS = 0;

		// Compute the atmospheric parameters
		ISA(h, T, a, P, rho, visc);

		// Compute the IAS
		const double M = TAS / a;
		const double P0_P = pow(1 + (gamma - 1) * M * M / 2, gamma / (gamma - 1));

		IAS = sqrt(2 * P * 1000 * (P0_P - 1) / rho_sl);

		return IAS;
	}

	double IAS_to_TAS(const double& IAS, const double& h) {
		// Calculates the TAS from the IAS and the altitude h
		// 
		// Input IAS must be in m/s
		// Input h must be in km
		// Output TAS is in m/s
		//
		// This function adapts the following website:
		// https://aerotoolbox.com/airspeed-conversions/


		// Initialize the outputs from ISA
		double T = 0;
		double a = 0;
		double P = 0;
		double rho = 0;
		double visc = 0;

		// Gamma Air
		double gamma = 1.4;
		const double rho_sl = 1.225;

		// Initialize the function output
		double TAS = 0;

		// Compute the atmospheric parameters
		ISA(h, T, a, P, rho, visc);

		// Compute the TAS
		TAS = 0.5 * rho_sl * IAS * IAS / (P * 1000) + 1;
		TAS = pow(TAS, (gamma - 1) / gamma);
		TAS = 2 * (TAS - 1) / (gamma - 1);
		TAS = a * sqrt(TAS);

		return TAS;
	}

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

	double compute_delta_CD_fuselage(const double& current_d, const double& current_l,
		const double& next_d, const double& next_l, const double& M, 
		const double& cruise_h) {

		// Calculate ISA values
		double T = 0;
		double a = 0;
		double P = 0;
		double rho = 0;
		double visc = 0;
		ISA(cruise_h, T, a, P, rho, visc);

		// Calculate TAS and the Reynolds number
		double TAS = a * M;
		const double Re = rho * TAS * current_l / visc;

		// Compute the fricion coefficient for a flat ptlate
		double Cf = compute_turb_frict_coeff(Re, M);

		// Compute the gradient of the Fuselage Skin Drag Coefficient
		double current_t = current_d / current_l;
		double dC_dt = Cf * (-pow(current_t, -2) - 0.05 * pow(current_t, -3) + 60);

		// Calculate change in the tube fineness ratio
		double next_t = next_d/next_l;
		double delta_t = next_t - current_t;

		// Calculate the delta in CD using a 1st order Taylor expansion
		return delta_t * dC_dt;
	}

	double breguet_prop_range(const double& eta_prop, const double& BSFC, const double& L_D,
		const double& w_start, const double& w_end) {
		// Implements the Breguet Range Equation for prop aircraft, as seen in the following 
		// link: https://en.wikipedia.org/wiki/Range_(aeronautics)
		//
		// The inputs are propulsive efficiency (eta_prop), the Break Specific Fuel Consumption
		// (BSFC, in Kg/J), the Lift-Over-Drag ratio (L_D, assumed constant throughout the leg),
		// and the start and end total aircraft weights (w_start and w_end, both with the same 
		// units).
		//
		// The output is the range in m.

		double range = eta_prop / (9.81 * BSFC);

		range = range * L_D * log(w_start / w_end); // log() is actually ln() in C++

		return range;
	}

	double breguet_prop_wratio(const double& eta_prop, const double& BSFC, const double& L_D,
		const double& range) {
		// Implements the Breguet Range Equation for prop aircraft, as seen in the following 
		// link: https://en.wikipedia.org/wiki/Range_(aeronautics)
		//
		// The inputs are propulsive efficiency (eta_prop), the Break Specific Fuel Consumption
		// (BSFC, in Kg/J), the Lift-Over-Drag ratio (L_D, assumed constant throughout the leg),
		// and the range (range, in m).
		//
		// The output is the ratio of the weight of the aircraft at the start relative to the
		// end.

		return exp(9.81 * BSFC * range / (L_D * eta_prop));
	}

	double compute_cruise_BSFC_PW127(const double& h) {
		// Implements the engine performance table at ISA conditions given in the following link:
		// https://www.quora.com/At-cruise-speed-do-turboprops-run-at-their-maximal-rated-horse-power-If-not-how-much-less-is-that-given-horse-power-typically
		//
		// The tables assume a CG location of 25%.
		// 
		// The input "h" is height in km, which must lie above 2.44 km and below 7.61 km. The 
		// returned output "PW127_BSFC" is the break-specific fuel consumption in kg/J.

		// Convert the height to flight level
		const double FL = 3.28084 * h * 10;

		// Check the flight level is reasonable
		assert(FL <= 250);
		assert(FL >= 80);

		// The maximum cruise power occurs when the engine is at 82% N2.
		const double P_cruise_max = 1589.8321; // kW (eq. to 2132 SHP, or 82% of 2750 SHP)

		// Data from the table
		const std::vector<double> FL_vec = 
			{ 80.0, 100., 120., 140., 160., 180., 200., 220., 240., 250. };
		const std::vector<double> Torque_percents =
			{ 94.5, 90.2, 86.1, 82.8, 78.9, 74.2, 68.9, 63.6, 58.3, 55.3 };
		const std::vector<double> fuel_rates =
			{ 464., 440., 418., 401., 381., 359., 335., 310., 286., 272. };

		// Interpolate between the table points to find the torque fraction and fuel consumption
		// per engine
		double T_frac = 0; // in percent
		MathTools::interpolate_fn(FL_vec, Torque_percents, FL, T_frac);
		double fuel_rate = 0; // in kg of fuel/hour/engine
		MathTools::interpolate_fn(FL_vec, fuel_rates, FL, fuel_rate);

		// Calculate and return the BSFC
		double PW127_BSFC = fuel_rate/(60 * 60 * T_frac * P_cruise_max * 10);
		return PW127_BSFC;
	}

	double compute_ground_run_raymer(const double& MTOW, const double& S_wing, const double&
		sl_rho_ratio, const double& CL_TO, const double& BHP) {
		// Implements the take-off run estimation from Fig 5.4 in the Raymer book.
		//
		// The inputs are "MTOW" in tonnes, the wing area "S_wing" in m^2, the density ratio
		// relative to sea level "sl_rho_ratio", the take-off lift coefficient "CL_TO" and the
		// power of the engine in horsepower "BHP". The output "TO_run" is the take-off run in 
		// km.

		// Compute the take-off parameter
		double TO_param = (MTOW / S_wing) / (sl_rho_ratio * CL_TO * BHP / MTOW);

		//// Check the take-off parameter is within bounds
		//assert(TO_param <= 600);
		//assert(TO_param >= 100);

		// Data from the figure
		const std::vector<double> TO_param_vec = { 100, 200, 300, 400, 500, 600 };
		const std::vector<double> run_vals = { 0.2189, 0.4275, 0.6759, 0.9641, 1.2523, 1.5604 };

		// Interpolate between the table points to find the torque fraction and fuel consumption
		// per engine
		double TO_run = 0; // in percent
		MathTools::interpolate_fn(TO_param_vec, run_vals, TO_param, TO_run);

		return TO_run;
	}

	double compute_eta_prop_raymer(const double& thrust, const double& TAS, const double& P) {
		// Computes the propeller efficiency according to equation (13.15) from Raymer.
		//
		// The input "thrust" is the Thrust in kN, the "TAS" is the True Arispeed in m/s, and
		// "P" is the power in kW. If the thrust of one engine is given, the power given must
		// also correspond to one engine. The output is the propeller efficiency.

		return thrust * TAS / P;
	}

	double correl_turboprop_mass(const double& P_max) {
		// Computes the turboprop mass (in kg) for a given max. power requirement by using a
		// correlation from the data available in the following website under the 
		// specifications section:
		// https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_PW100
		// 
		// The input is the maximum power requirement of the engine in kW, and the output is
		// the likely engine mass in kg.
	
		return 0.1391 * P_max + 200.75;
	}

	double correl_turboprop_TOBSFC(const double& P_max) {
		// Computes the turboprop Break-Specific Fuel Consumption at T/O for a given 
		// maximum power requirement by using a correlation from the data available in the 
		// following website under the specifications section:
		// https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_PW100
		// 
		// The input is the maximum power requirement of the engine in kW, and the output is
		// the likely BSFC in g/kWh.

		return 796.03 * pow(P_max, -0.136);
	}

	double compute_new_engine_cruise_BSFC(const double& BSFC_TO, const double& h) {
		// Performs a simple extrapolation to link the BSFC at TO with the BSFC at the chosen
		// cruise altitude h (which must be provided in km). Returns the BSFC at cruise with
		// the same units as the input units of the BSFC.

		const double PW127_BSFC_TO = 279; // g/kWh
		const double PW127_BSFC_cruise = compute_cruise_BSFC_PW127(h) * 3.6e+9; // g/kWh
		
		return BSFC_TO * PW127_BSFC_cruise / PW127_BSFC_TO;
	}

	double calculate_hybrid_BSFC(const double& H2_frac, const double& h, const double& P_max) {
		// Given an input hydrogen fraction "H2_frac", a cruise altitude "h" in km, and a 
		// maximum power output for a single turboprop engine "P_max" in kW, give the hybrid
		// BSFC for a single engine in kg/J;
		// 
		// Specific energy data:
		//    - JA1: https://en.wikipedia.org/wiki/Jet_fuel
		//    - H2: https://www.skai.co/hydrogen-details#:~:text=The%20specific%20energy%20of%20hydrogen,lithium%2Dion%20batteries%20(approximately%20
		//
		// BSFC Equation: BSFC = Fuel Consumption/Power
		//                     = ( H2_frac/c_H2 + 1/c_JA1 ) / ( H2_frac + 1 ) / eta_therm
		//
		// Please note that in this case "thermal efficiency" is the product of both thermal 
		// and combustion efficiency.
	
		const double c_JA1 = 43.15; // MJ/kg
		const double c_H2 = 142; // MJ/kg

		// First assume it's all kerosene
		const double BSFC_JA1_TO = correl_turboprop_TOBSFC(P_max); // g/kWh
		double BSFC_JA1_cruise = compute_new_engine_cruise_BSFC(BSFC_JA1_TO, h); // g/kWh

		// Convert BSFC from g/kWh to kg/MJ
		BSFC_JA1_cruise = BSFC_JA1_cruise / ( 3.6 * 1000 );

		// Using the BSFC equation, the thermal efficiency can be computed
		const double eta_therm = 1 / (BSFC_JA1_cruise * c_JA1);

		// Assuming the thermal efficiency is the same for both fuels, calculate the hybrid BSFC
		double BSFC_hybrid = (H2_frac/c_H2 + 1/c_JA1) / (H2_frac + 1) / eta_therm;

		// Convert the BSFC from kg/MJ to kg/J and return the result
		return BSFC_hybrid / (1e6);
	}

	bool compute_cg_loc_mass(const double& ip_M_engine, const double& ip_M_fuel,
		const double& ip_H2_frac, double& op_cg_loc, double& op_calc_mass, double& op_cg_loc_nofuel,
		double& op_calc_mass_nofuel, double& op_payload, double& op_M_JA1, double& op_num_pass,
		double& op_tank_l, bool& op_vio_mass, bool& op_vio_vol) {
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


		// Initialise all outputs
		op_cg_loc = 0;
		op_calc_mass = 0;
		op_payload = 0;
		op_vio_mass = false;
		op_vio_vol = false;

		// Initialize aircraft constants
		const double MTOW = 22000; // kg
		const double M_empty = 13500; // kg
		const double x_CG_empty = 12.202; // m
		const double x_CG_JA1 = 12.203; // m
		const double M_cargo_front_max = 928; // kg
		const double x_CG_cargo_front = 4.192; //m
		const double M_cargo_rear_max = 637; // kg
		const double x_CG_cargo_rear = 21.4555; // m
		const double pass_packing_density = 429.02; // kg/m
		const double x_per_row = 12.84 / 17. ; // m
		double M_pass_total = 0.; //kg

		// Initialize engine constants
		const double M_PW127 = 480; // kg
		const double x_CG_engine = 10.63; // m

		// Initialize fuel constants
		const double c_JA1 = 43.0; // MJ/kg
		const double c_H2 = 121.09; // MJ/kg
		const double rho_H2 = 67.3; // kg/m^3
		const double tank_eta = 0.63; // H2 kg req / system kg
		const double c = 1.121; // m of effective tank radius
		const double tank_vol_max = 47.11; //m^3 (in order to keep CG at same location as empty aircraft)

		// Calculate the base aircraft empty CG product and mass
		double CG_product = M_empty * x_CG_empty;
		double M_total = M_empty;

		// Account for the different engine
		const double delta_M_engines = 2 * ( ip_M_engine - M_PW127 ); // kg
		CG_product += delta_M_engines * x_CG_engine;
		M_total += delta_M_engines;

		// Calculate the amount of kerosene and hydrogen
		const double M_JA1 = ip_M_fuel / (ip_H2_frac * c_JA1 / c_H2 + 1); // kg
		const double M_H2 = M_JA1 * ip_H2_frac * c_JA1 / c_H2; // kg
		op_M_JA1 = M_JA1;

		// Account for the kerosene
		CG_product += M_JA1 * x_CG_JA1;
		M_total += M_JA1;

		// Calculate the volume and mass of the H2 system, and account for it
		const double Vol_H2sys = M_H2 / rho_H2; // m^3
		const double M_H2_system = M_H2 / tank_eta; // kg
		CG_product += M_H2_system * x_CG_empty;
		M_total += M_H2_system;

		// Check whether the volume constraint has been violated
		if (Vol_H2sys > tank_vol_max) {
			op_vio_vol = true;
		}
		else {
			op_vio_vol = false;
		}

		// Check whether the mass constraint has been violated
		double M_pay_remaining = MTOW - M_total;
		op_payload = M_pay_remaining;

		if (M_pay_remaining <= 0) {
			op_vio_mass = true;
		}
		else {
			op_vio_mass = false;
		}

		// First, fill out the rear cargo compartment
		if (M_pay_remaining >= M_cargo_rear_max) {
			M_pay_remaining += -M_cargo_rear_max;

			// Account for the rear cargo compartment
			CG_product += M_cargo_rear_max * x_CG_cargo_rear;
			M_total += M_cargo_rear_max;
		}
		else {
			// Put all the remaining payload in the rear cargo compartment 
			CG_product += M_pay_remaining * x_CG_cargo_rear;
			M_total += M_pay_remaining;
			
			op_payload = M_total - M_empty - delta_M_engines - M_H2_system - M_JA1;
			op_calc_mass = M_total;
			op_cg_loc = CG_product / M_total;
				
			return true;
		}

		// Secondly, fill out the front cargo compartment
		if (M_pay_remaining >= M_cargo_front_max) {
			M_pay_remaining += -M_cargo_front_max;

			// Account for the front cargo compartment
			CG_product += M_cargo_front_max * x_CG_cargo_front;
			M_total += M_cargo_front_max;
		}
		else {
			// Put all the remaining payload in the front cargo compartment 
			CG_product += M_pay_remaining * x_CG_cargo_front;
			M_total += M_pay_remaining;

			op_payload = M_total - M_empty - delta_M_engines - M_H2_system - M_JA1;
			op_calc_mass = M_total;
			op_cg_loc = CG_product / M_total;
			op_calc_mass_nofuel = M_total - M_JA1 - M_H2;
			op_cg_loc_nofuel = (CG_product - M_JA1 * x_CG_JA1 - M_H2 * x_CG_empty) 
				/ op_calc_mass_nofuel;

			return true;
		}

		// Fill out the cabin with the remaining payload
		const double l_tank = (Vol_H2sys - 4 * 3.14159265358979323846 * pow(c, 3)/3 ) /
			(3.14159265358979323846 * pow(c, 2)) + 2 * c;
		const double x_CG_payload1 = 5.86/2 + x_CG_empty/2 - l_tank / 4;
		const double x_CG_payload2 = 19.2426 / 2 + x_CG_empty / 2 + l_tank / 4;
		const double x_avail_payload1 = x_CG_empty - 5.86 - l_tank / 2; // 5.86 is the front start
																		// of the aircraft
		const double M_pass_max1 = x_avail_payload1 * pass_packing_density;
		const double x_avail_payload2 = 19.2426 - x_CG_empty - l_tank / 2; // 19.2426 is the end
																		   // coord of the plane
		const double M_pass_max2 = x_avail_payload2 * pass_packing_density;
		const double M_pass_max = M_pass_max1 + M_pass_max2;
			
		// Fill out the rear cabin first since it's larger
		if (M_pay_remaining >= M_pass_max2) {
			M_pay_remaining += -M_pass_max2;

			// Account for the front cargo compartment
			CG_product += M_pass_max2 * x_CG_payload2;
			M_total += M_pass_max2;
			M_pass_total += M_pass_max2;
		}
		else {
			// Put all the remaining payload in the front cargo compartment 
			CG_product += M_pay_remaining * x_CG_payload2;
			M_total += M_pay_remaining;
			M_pass_total += M_pay_remaining;
			
			op_payload = M_total - M_empty - delta_M_engines - M_H2_system - M_JA1;
			op_calc_mass = M_total;
			op_cg_loc = CG_product / M_total;
			op_calc_mass_nofuel = M_total - M_JA1 - M_H2;
			op_cg_loc_nofuel = (CG_product - M_JA1 * x_CG_JA1 - M_H2 * x_CG_empty)
				/ op_calc_mass_nofuel;

			return true;
		}

		// Then fill out the front cabin 
		if (M_pay_remaining >= M_pass_max1) {
			M_pay_remaining += -M_pass_max1;

			// Account for the front cargo compartment
			CG_product += M_pass_max1 * x_CG_payload1;
			M_total += M_pass_max1;
			M_pass_total += M_pass_max1;
		}
		else {
			// Put all the remaining payload in the front cargo compartment 
			CG_product += M_pay_remaining * x_CG_payload1;
			M_total += M_pay_remaining;
			M_pass_total += M_pay_remaining;
		}

		op_tank_l = l_tank;
		op_num_pass = floor(M_pass_total / 80.);
		op_payload = M_total - M_empty - delta_M_engines - M_H2_system - M_JA1;
		op_calc_mass = M_total;
		op_cg_loc = CG_product / M_total;
		op_calc_mass_nofuel = M_total - M_JA1 - M_H2;
		op_cg_loc_nofuel = (CG_product - M_JA1 * x_CG_JA1 - M_H2 * x_CG_empty)
			/ op_calc_mass_nofuel;

		return true;
	}

}
