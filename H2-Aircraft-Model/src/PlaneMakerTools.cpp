#include "../headers/PlaneMakerTools.h"

namespace PlaneMakerTools {

	size_t find_line(const std::string& line, const std::string& acf_filepath, 
		std::vector<std::string>& file_lines) {
	
		// Open the file containing the designs
		std::ifstream acf_file(acf_filepath);

		// Initialise a string to store the lines in
		std::string temp_line;

		// Make a vector to store the lines
		std::vector<std::string> acf_lines;

		// Initialize the indexes
		size_t i = 0;
		size_t idx = 0;

		// Read all lines
		while (std::getline(acf_file, temp_line)) {
			
			acf_lines.push_back(temp_line);

			// If we haven't found a match before, try to match
			if (idx == 0) {
				size_t temp_idx = temp_line.find(line);
				
				// If a match is found, store the current line index
				if (temp_idx != std::string::npos)
				{
					idx = i;
				}
			}

			// Increment the indez
			i++;
		}

		// Close the design file and clear the relevant vectors
		acf_file.close();

		file_lines = acf_lines;

		return idx;
	}

	void write_line(const size_t& idx, const std::string& line, std::vector<std::string> acf_lines,
		const std::string& acf_filepath) {
	
		// Replace the line in the file with the new line
		acf_lines[idx] = line;

		// Open the output file
		std::ofstream acf_file(acf_filepath);

		// Write the file
		for (size_t i = 0; i < acf_lines.size(); i++) {
		
			acf_file << acf_lines[i] << "\n";
		}

		// Close the file
		acf_file.close();
	}

	void write_quantity(const double& quantity, std::string line, const std::string& filepath) {

		// Initialise vector to hold file lines
		std::vector<std::string> file_lines;

		// Find the index of the desired line
		const size_t line_idx = find_line(line, filepath, file_lines);

		// Add the new quantity value to the BSFC line
		line += std::to_string(quantity);

		// Write the BSFC data
		write_line(line_idx, line, file_lines, filepath);
	}

	void set_weight_data(const double& x_CG_nofuel, const double& M_nofuel, 
		const std::string& acf_filepath) {
		// Input x_CG in m from the nose, input M_nofuel in kg
	
		const double x_origin = 10.505; // ft
		const double m_to_ft = 3.28084;
		const double kg_to_lb = 2.20462;
		
		// Convert length to feet and subtract origin
		const double x_CG = x_CG_nofuel * m_to_ft - x_origin;

		// Convert mass to lbs
		const double M = M_nofuel * kg_to_lb;

		// CG Line
		const std::string cg_line = "P acf/_cgZ ";
		write_quantity(x_CG, cg_line, acf_filepath);

		// Mass Line
		const std::string M_line = "P acf/_m_empty ";
		write_quantity(M, M_line, acf_filepath);
	}

	void set_engine_data(const double& ip_BSFC_full, const double& ip_BSFC_full_low, const double& P_max,
		const std::string& acf_filepath) {
		// Input BSFC in kg/J, input P_max in kW

		const double J_to_hp_hr = 3.725061361111e-7;
		const double kg_to_lb = 2.20462;
		const double high_to_low = 42/53;

		// Convert BSFC to lb/hp*hr
		const double BSFC_full = ip_BSFC_full * kg_to_lb / J_to_hp_hr;
		const double BSFC_half = BSFC_full * high_to_low;
		const double BSFC_full_low = ip_BSFC_full_low * kg_to_lb / J_to_hp_hr;
		const double BSFC_half_low = BSFC_full_low * high_to_low;

		// Convert power to hp
		const double P_max_hp = P_max * J_to_hp_hr * 1000;

		// BSFC_full High Altitude
		const std::string BSFC_full_line = "P acf/_SFC_full_hi_PRP ";
		write_quantity(BSFC_full, BSFC_full_line, acf_filepath);

		// BSFC_half High Altitude
		const std::string BSFC_half_line = "P acf/_SFC_half_hi_PRP ";
		write_quantity(BSFC_half, BSFC_half_line, acf_filepath);

		// BSFC_full Low Altitude
		const std::string BSFC_full_low_line = "P acf/_SFC_full_lo_PRP ";
		write_quantity(BSFC_full_low, BSFC_full_low_line, acf_filepath);

		// BSFC_half Low Altitude
		const std::string BSFC_half_low_line = "P acf/_SFC_half_lo_PRP ";
		write_quantity(BSFC_half_low, BSFC_half_low_line, acf_filepath);

		// Max Power
		const std::string Power_line = "P acf/_power_max_limit ";
		write_quantity(P_max_hp, Power_line, acf_filepath);
	}

	void set_fuel_data(double M_fuel, double x_h2_tanks, const double& H2_M_prop, 
		const std::string& acf_filepath) {
		// M_fuel is in kg
		// H2_M_prop ranges from 0 to 1
		// x_h2_tanks is in meters

		const double x_origin = 10.505; // ft
		const double m_to_ft = 3.28084;
		const double kg_to_lb = 2.20462;

		// Convert inputs to imperial units for X-Plane
		M_fuel *= kg_to_lb;
		x_h2_tanks = x_h2_tanks*m_to_ft - x_origin;

		// Calculate the proportion of fuel in each kerosene tank
		const double JA1_M_prop_single = 0.5 * (1 - H2_M_prop);

		// Tank 0 (Kerosene) Proportion
		const std::string tank0_line = "acf/_tank_rat/0 ";
		write_quantity(JA1_M_prop_single, tank0_line, acf_filepath);

		// Tank 1 (Kerosene) Proportion
		const std::string tank1_line = "acf/_tank_rat/1 ";
		write_quantity(JA1_M_prop_single, tank1_line, acf_filepath);

		// Tank 2 (LH2) Proportion
		const std::string tank2_line = "acf/_tank_rat/2 ";
		write_quantity(H2_M_prop, tank2_line, acf_filepath);

		// Tank 2 Lateral Locations
		const std::string tank2_lat1 = "P acf/_tank_xyz_full/2,0 ";
		write_quantity(0, tank2_lat1, acf_filepath);

		const std::string tank2_lat2 = "P acf/_tank_xyz_full/2,1 ";
		write_quantity(0, tank2_lat2, acf_filepath);

		// Empty Tank 2 Axial Locations
		const std::string tank2_ax = "P acf/_tank_xyz_full/2,2 ";
		write_quantity(x_h2_tanks, tank2_ax, acf_filepath);

		// Empty Tank 2 Lateral Locations
		const std::string tank2_lat1_empty = "P acf/_tank_xyz/2,0 ";
		write_quantity(0, tank2_lat1_empty, acf_filepath);

		const std::string tank2_lat2_empty = "P acf/_tank_xyz/2,1 ";
		write_quantity(0, tank2_lat2_empty, acf_filepath);

		// Empty Tank 2 Axial Locations
		const std::string tank2_ax_empty = "P acf/_tank_xyz/2,2 ";
		write_quantity(x_h2_tanks, tank2_ax_empty, acf_filepath);
	}


}