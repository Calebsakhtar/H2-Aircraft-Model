
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
		// The input "h" is height in km. The output "PW127_BSFC" is the break-specific fuel
		// consumption in kg/J.

		// Convert the height to FL
		const double FL = 3.28084 * h * 10;

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

}
