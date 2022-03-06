#ifndef AIRCRAFT_PLANEMAKERTOOLS_H
#define AIRCRAFT_PLANEMAKERTOOLS_H

#include <vector>
#include <string>
#include <fstream>

namespace PlaneMakerTools {

	// Writes a "quantity" to its corresponding "line" in the X-Plane acf file located in
	// "filepath".
	void write_quantity(const double& quantity, std::string line, const std::string& filepath);

	// Writes the relevant weight and balance data to the X-Plane acf file.
	// The inputs have the following details:
	// 
	//	- x_CG in m from the nose,
	//	- M_nofuel in kg
	void set_weight_data(const double& x_CG_nofuel, const double& M_nofuel,
		const std::string& acf_filepath);

	// Writes the relevant engine data to the X-Plane acf file.
	// The inputs have the following details:
	// 
	//	- All SFCs in kg/J, 
	//	- P_max in kW
	void set_engine_data(const double& ip_BSFC_full, const double& ip_BSFC_full_low, const double& P_max,
		const std::string& acf_filepath);

	// Writes the relevant fuel and tank data to the X-Plane acf file.
	// The inputs have the following details:
	// 
	//	 - M_fuel is in kg
	//	 - H2_M_prop ranges from 0 to 1
	//   - x_h2_tanks is in m, measured from the nose
	void set_fuel_data(double M_fuel, double x_h2_tanks, const double& H2_M_prop,
		const std::string& acf_filepath);

}

#endif
