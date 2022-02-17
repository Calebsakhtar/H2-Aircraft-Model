#ifndef AIRCRAFT_MATHTOOLS_H
#define AIRCRAFT_MATHTOOLS_H

#include <vector>
#include <assert.h>
#include <algorithm>

namespace MathTools {

	// Return the interpolation result
	bool interpolate_fn(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
		const double& input_x, double& result_y, const bool reverse = false);

}

#endif
