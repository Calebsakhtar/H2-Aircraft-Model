#include "../headers/MathTools.h"

namespace MathTools {

	bool interpolate_fn(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
		const double& input_x, double& result_y, const bool reverse) {
		// Given an input "input_x" that lies within the range of values of the sorted vector
		// "x_vals", return the corresponding output "result_y" by interpolating through the
		// output sample of points "y_vals".

		// Perform checks to the input data
		assert(x_vals.size() > 0);
		assert(x_vals.size() == y_vals.size()); // Check the inputs are of consistent size

		// If we are going from input to output (and not the reverse), perform these checks
		if (!reverse) {
			assert(std::is_sorted(x_vals.begin(), x_vals.end())); // Check x vector is sorted
			assert(x_vals[0] <= input_x); // Assert the input is at least equal to the smallest xval
			assert(x_vals.back() >= input_x); // Assert the input does not surpass the largest xval
		}

		// initialize the variable to hold the index 
		size_t j = 0;
	
		// Find the "slot" the input lies in
		for (size_t i = 0; i < x_vals.size() - 1; i++) {

			if (input_x <= x_vals[i + 1]) {
				// Perform Linear Interpolation within the found slot
				double x_frac = (input_x - x_vals[i]) / (x_vals[i + 1] - x_vals[i]);
				result_y = y_vals[i] + x_frac * (y_vals[i + 1] - y_vals[i]);

				// Return true since an answer has been found
				return true;
			}
		}

		// Return false since an answer has not been found
		return false;
	}

}