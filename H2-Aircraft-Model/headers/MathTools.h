#ifndef AIRCRAFT_MATHTOOLS_H
#define AIRCRAFT_MATHTOOLS_H

namespace MathTools {

	// Return the interpolation result
	bool interpolate_fn(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
		const double& input_x, double& result_y, const bool reverse = false);

}

#endif
