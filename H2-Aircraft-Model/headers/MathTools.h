#ifndef AIRCRAFT_MATHTOOLS_H
#define AIRCRAFT_MATHTOOLS_H

#include <vector>
#include <assert.h>
#include <algorithm>

namespace MathTools {

	// Given an input "input_x" that lies within the range of values of the sorted vector
	// "x_vals", return the corresponding output "result_y" by interpolating through the
	// output sample of points "y_vals".
	bool interpolate_fn(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
		const double& input_x, double& result_y, const bool reverse = false);

}

#endif
