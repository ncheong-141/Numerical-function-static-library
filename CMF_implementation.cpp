/* Compile time variant implementation */

#include "SL_Comptime_Interface.h"

namespace CMF {

	/* isEven, isOdd*/
	bool isEven(const int& parameter) {

		// Calculate remainder of number divided by 2 (modulus) (if 0 then it is even)
		int modulus = parameter % 2;

		if (modulus == 0)
			return true;
		else
			return false;
	}
	bool isEven(const size_t& parameter) {

		// Calculate remainder of number divided by 2 (modulus) (if 0 then it is even)
		int modulus = parameter % 2;

		if (modulus == 0)
			return true;
		else
			return false;
	}
	bool isOdd(const int& parameter) {

		// Calculate remainder of number divided by 2 (modulus) (if 0 then it is even)
		int modulus = parameter % 2;

		if (modulus != 0)
			return true;
		else
			return false;
	}
	bool isOdd(const size_t& parameter) {

		// Calculate remainder of number divided by 2 (modulus) (if 0 then it is even)
		int modulus = parameter % 2;

		if (modulus != 0)
			return true;
		else
			return false;
	}

	/*Linear interpolation implementation */
	// Linear interp between two points which returns a reference array pt
	void l_interp(float(&pt0)[2], float(&pt1)[2], float(&out_x_y)[2]) {

		// Between point linear interp 
		// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);
		out_x_y[1] = ((pt0[1] * pt1[0]) - (pt0[1] * out_x_y[0]) + (pt1[1] * out_x_y[0]) - (pt1[1] * pt0[0])) / (pt1[0] - pt0[0]);
	}
	void l_interp(double(&pt0)[2], double(&pt1)[2], double(&out_x_y)[2]) {

		// Between point linear interp 
		// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);
		out_x_y[1] = ((pt0[1] * pt1[0]) - (pt0[1] * out_x_y[0]) + (pt1[1] * out_x_y[0]) - (pt1[1] * pt0[0])) / (pt1[0] - pt0[0]);

	}

	// Linear interpolation between two points which returns a value 
	float l_interp(float(&pt0)[2], float(&pt1)[2], float x) {

		// Between point linear interp 
		// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);
		return ((pt0[1] * pt1[0]) - (pt0[1] * x) + (pt1[1] * x) - (pt1[1] * pt0[0])) / (pt1[0] - pt0[0]);
	}
	double l_interp(double(&pt0)[2], double(&pt1)[2], double x) {

		// Between point linear interp 
		// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);
		return ((pt0[1] * pt1[0]) - (pt0[1] * x) + (pt1[1] * x) - (pt1[1] * pt0[0])) / (pt1[0] - pt0[0]);
	}
} // End of namespace CMF
