
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

	void ISA(const double& ip_h, double& op_T, double& op_a, double& op_P, double& op_rho) {
		// ISA Returns the International Standard Atmosphere(ISA) values for ambient
		// temperature(op_T) in degrees kelvin, speed of sound(op_a) in m / s, pressure
		// (op_P)in kPa and density(op_rho) in Kg / m ^ 3 at a given input array height in
		// kilometers(ip_h).
	
		// Initialize the relevant quantities at sea level
		const double T_sl = 288.15;
		const double a_sl = 340.3;
		const double P_sl = 101.325;
		const double rho_sl = 1.225;

		// Set the output quantities to equal to 0 (as initialization)
		op_T = 0;
		op_a = 0;
		op_P = 0;
		op_rho = 0;

		if (ip_h <= 0) {
			// Set the output quantities to equal to the sea level amount if below sea level
			op_T = T_sl;
			op_a = a_sl;
			op_P = P_sl;
			op_rho = rho_sl;
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
		}
	}

	//double compute_delta_CD(const double& current_t, const double& current_length,
	//	const double& M, const double& cruise_alt) {
	//}
}