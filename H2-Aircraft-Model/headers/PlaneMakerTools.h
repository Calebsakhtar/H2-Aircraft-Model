#ifndef AIRCRAFT_PLANEMAKERTOOLS_H
#define AIRCRAFT_PLANEMAKERTOOLS_H

#include <vector>
#include <string>
#include <fstream>

namespace PlaneMakerTools {

	// Find Line
	size_t find_line(const std::string& line, const std::string& acf_filepath,
		std::vector<std::string>& file_lines);

	void write_line(const size_t& idx, const std::string& line, std::vector<std::string> acf_lines,
		const std::string& acf_filepath);

	void write_quantity(const double& quantity, std::string line, const std::string& filepath);

	// Set weight data
	void set_weight_data(const double& x_CG_nofuel, const double& M_nofuel,
		const std::string& acf_filepath);

	// Set engine data
	void set_engine_data(const double& BSFC_Pmax_hybrid, const double& P_max,
		const std::string& acf_filepath);

	// Set fuel tank data

}

#endif
