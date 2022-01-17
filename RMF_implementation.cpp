/* Implementation file for general func static library header*/

// Get interface files to link to console 
//#include "SL_Comptime_Interface.h"
#include "SL_Runtime_Interface.h"

namespace RMF {

	/*------------------ RMF Guards of the Static Library ---------------------------------------------*/
	// Note: - These are only used in the functions in this namespace. 
	//		 - For incorrect vector size of single, the mem_correcter checks for this (and corrects it). 

	/* ---- 1D array guards ----*/
	// - Prevent a zero size vector being used 
	// - Prevent two different vectors of a different length 
	template <typename T> void size_guard(std::vector<T>& dyn_arr_data, const char* fnName) {
		if (dyn_arr_data.size() == 0) {
			std::cout << "Error: Dynamic vector has zero size in function, " << fnName << ".\n";
			exit(-1);
		}
	}
	template <typename T> void size_guard(std::vector<T>& vectorName, std::vector<T>& vectorName2, const char* fnName) {

		// Guard to stop vectors/1D arrays of different lengths being multiplied together 
		if (vectorName.size() != vectorName2.size()) {
			std::cout << "Error: You cannot elementwise multiply two vectors/1D arrays of different lengths in function," << fnName << ".\n";
			exit(-1);
		}
	}

	/* ---- 2D array guards ---- */
	// - Prevent a zero size array being used 
	// - Since vector can have vectors of differnet lengths inside, check for same length. (Memory violation prevention)
	// - Prevent two different vectors of a different length 
	template <typename T> void size_guard(std::vector< std::vector<T>>& dyn_arr_data, const char* fnName) {

		size_t ref_rowsize = dyn_arr_data.size();
		size_t ref_colsize = dyn_arr_data[0].size();	// Calculate ref size only once (instead of in loop)

		if (ref_rowsize == 0) {
			std::cout << "Error: Dynamic vector has zero row dimension in function," << fnName << ".\n";
			exit(-1);
		}
		for (size_t i = 0; i < ref_rowsize; i++) {
			if (dyn_arr_data[i].size() != ref_colsize) {
				std::cout << "Error: Dynamic vector has different collumn dimensions in function," << fnName << ".\n";
				exit(-1);
			}
			else if (dyn_arr_data[i].size() == 0) {
				std::cout << "Error: Dynamic vector has zero collumn dimensions in function," << fnName << ".\n";
				exit(-1);
			}

		}
	}
	template <typename T> void size_guard(std::vector< std::vector<T>>& arrayName, std::vector< std::vector<T>>& arrayName2, const char* fnName) {
		// Row-wise checks 
		if (arrayName.size() != arrayName2.size()) {
			std::cout << "Error: You cannot elementwise multiply two arrays of different row lengths in function," << fnName << ".\n"; 			exit(-1);
			exit(-1);
		}

		// Establish reference lengths (can be any of the vectors since passed row-wise checks 
		size_t ref_rowsize = arrayName.size();

		// Col-wise checks (Ensures all collumns are the same size)
		for (size_t i = 0; i < ref_rowsize; i++) {
			if (arrayName[i].size() != arrayName2[i].size()) {
				std::cout << "Error: You cannot elementwise multiply two arrays of different collumn lengths in function," << fnName << " at row position: " << i << ".\n"; 			exit(-1);
				exit(-1);
			}
		}
	}

	// DYN_C2D variants. Note no need to check if cols are different dimensions.
	template <typename T> void size_guard(DYN_C2D<T>& dyn_arr_data, const char* fnName) {
		if (dyn_arr_data.size_cols() == 0 || dyn_arr_data.size_rows() == 0) {
			std::cout << "Error: Dynamic vector has zero row or collumn dimension in function," << fnName << ".\n";
			exit(-1);
		}
	}
	template <typename T> void size_guard(DYN_C2D<T>& arrayName, DYN_C2D<T>& arrayName2, const char* fnName) {
		// Row-wise checks 
		if ((arrayName.size_rows() != arrayName2.size_rows()) || (arrayName.size_cols() != arrayName2.size_cols())) {
			std::cout << "Error: Two arrays of different row or collumn lengths in function," << fnName << ".\n";
			exit(-1);
		}
	}


	/* ------------------ START OF LIBRARY FUNCTIONS: General math tools -------------------------------*/

	/* any() -> if any element is non-zero return true */
	template <typename T> bool any(std::vector<T>& dyn_arr_data) {
		for (size_t i = 0; i < dyn_arr_data.size(); i++) {
			if (dyn_arr_data[i] > 0)
				return true;
		}
		return false;
	}
	template <typename T> bool any(std::vector< std::vector<T>>& dyn_arr_data) {

		/* Internal guards -> Runtime checks*/
		// Since vector can have vectors of differnet lengths inside, check for same length.  Otherwise you will get a memory violation. 
		const char* fnName = "any()";
		size_guard<T>(dyn_arr_data, fnName);

		size_t ref_rowsize = dyn_arr_data.size();
		size_t ref_colsize = dyn_arr_data[0].size();	// Calculate ref size only once (instead of in loop)

		for (size_t i = 0; i < ref_rowsize; i++) {
			for (size_t j = 0; j < ref_colsize; j++) {
				if (dyn_arr_data[i][j] != 0)
					return false;
			}
		}

		return true;
	}

	/* all() -> if all elements are non-zero, return true */
	template <typename T> bool all(std::vector<T>& dyn_arr_data) {

		size_t ref_cols = dyn_arr_data.size();

		for (size_t i = 0; i < ref_cols; i++) {
			if (dyn_arr_data[i] != 0)
				return false;
		}
		return true;
	}
	template <typename T> bool all(std::vector< std::vector<T>>& dyn_arr_data) {


		/* Internal guards -> Runtime checks*/
		// Since vector can have vectors of differnet lengths inside, check for same length.  Otherwise you will get a memory violation. 
		const char* fnName = "all()";
		size_guard<T>(dyn_arr_data, fnName);

		size_t ref_rowsize = dyn_arr_data.size();
		size_t ref_colsize = dyn_arr_data[0].size();	// Calculate ref size only once (instead of in loop)

		for (size_t i = 0; i < ref_rowsize; i++) {
			for (size_t j = 0; j < ref_colsize; j++) {
				if (dyn_arr_data[i][j] != 0)
					return false;
			}
		}
		return true;
	}

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


	/* abs() -> Returns an absolute value*/
	template <typename T> T abs(T parameter) {

		if (parameter < 0)
			return (-1 * parameter);

		return parameter;
	}

	/* arange() fills array with incrementing values (+- 1) from specified start to end*/
	// Generates values between start and end according to discrete parameter (Note, cols must have an extra element for end point) 
	void				arange(std::vector<int>& dyn_arr_data, const int start, const int end) {

		// Pre-define collumn dimension for reducing amount of calculations 
		size_t ref_colsize = abs(end - start);

		/* Specific vector inputs
		- If the vector has no allocated memory/size -> Reserve, resize and populate vector
		- If vector is already allocated, check if dimensions match -> if not, re-define, if so, populate vector */

		// Optimize data structure for function
		mem_corrector(dyn_arr_data, ref_colsize);

		// Start function implementation where vector is populated
		int increment = 0;
		if (end > start) {
			for (size_t i = 0; i < ref_colsize; i++) {
				dyn_arr_data[i] = start + increment;
				increment++;
			}
		}
		else if (start > end) {
			for (size_t i = 0; i < ref_colsize; i++) {
				dyn_arr_data[i] = start + increment;
				increment--;
			}
		}

	}
	void				arange(std::vector<float>& dyn_arr_data, const float start, const float end, const size_t NP) {

		// Optimize data structure for function
		mem_corrector(dyn_arr_data, NP);

		float total_length = abs(end - start);			// Calculate total value between start and end
		float disc_inc = total_length / (NP - 1);       // Calculate the increments between each element
		float disc_val = 0;								// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;					// Write element
			disc_val += disc_inc;						// Add increment
		}
	}
	void				arange(std::vector<double>& dyn_arr_data, const double start, const double end, const size_t NP) {

		// Optimize data structure for function
		mem_corrector(dyn_arr_data, NP);

		double total_length = abs(end - start);			// Calculate total value between start and end
		double disc_inc = total_length / (NP - 1);      // Calculate the increments between each element
		double disc_val = 0;						    // Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;					// Write element
			disc_val += disc_inc;						// Add increment
		}
	}
	std::vector<int>	arange(const int start, const int end) {

		// Pre-define collumn dimension for reducing amount of calculations 
		size_t ref_colsize = abs(end - start);

		// Define return vector 
		std::vector<int> dyn_arr_data(ref_colsize);

		// Start function implementation where vector is populated
		int increment = 0;
		if (end > start) {
			for (size_t i = 0; i < ref_colsize; i++) {
				dyn_arr_data[i] = start + increment;
				increment++;
			}
		}
		else if (start > end) {
			for (size_t i = 0; i < ref_colsize; i++) {
				dyn_arr_data[i] = start + increment;
				increment--;
			}
		}

		return dyn_arr_data;
	}
	std::vector<float>	arange(const float start, const float end, const size_t NP) {

		std::vector<float> dyn_arr_data(NP);

		float total_length = abs(end - start);			// Calculate total value between start and end
		float disc_inc = total_length / (NP - 1);       // Calculate the increments between each element
		float disc_val = 0;								// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;					// Write element
			disc_val += disc_inc;						// Add increment
		}

		return dyn_arr_data;
	}
	std::vector<double> arange(const double start, const double end, const size_t NP) {

		std::vector<double> dyn_arr_data(NP);

		double total_length = abs(end - start);			// Calculate total value between start and end
		double disc_inc = total_length / (NP - 1);       // Calculate the increments between each element
		double disc_val = 0;								// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;					// Write element
			disc_val += disc_inc;						// Add increment
		}

		return dyn_arr_data;
	}

	/* tile() fills array with a specified value*/
	template <typename T> void tile(T tilevalue, std::vector<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "tile()");

		for (size_t i = 0; i < dyn_arr_data.size(); i++) {
			dyn_arr_data[i] = tilevalue;
		}

	}
	template <typename T> void tile(T tilevalue, std::vector< std::vector<T> >& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "tile()");

		for (size_t i = 0; i < dyn_arr_data.size(); i++) {
			for (size_t j = 0; j < dyn_arr_data[i].size(); j++) {
				dyn_arr_data[i][j] = tilevalue;
			}
		}
	}
	template <typename T> void tile(T tilevalue, DYN_C2D<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "tile()");

		// Loop over with iterator as much more efficient with DYN_C2D
		for (auto data_addr = dyn_arr_data.begin(); data_addr != dyn_arr_data.end(); data_addr++) {
			*data_addr = tilevalue;
		}
	}

	/* Elementwise multiplication of vectors */
	template <typename T> void elementwise(std::vector<T>& output, std::vector<T>& vectorName, std::vector<T>& vectorName2) {

		// Guard to stop vectors/1D arrays of different lengths being multiplied together 
		size_guard<T>(vectorName, vectorName2, "Elementwise");

		// Instantiate reserve space for vector for optimizing performance (reduces copying) 
		mem_corrector(output, vectorName.size());

		// Loop along elements of vector/ 1D arrays 
		for (int i = 0; i < vectorName.size(); i++) {
			output[i] = (vectorName[i] * vectorName2[i]);

			if (i == 0)
				std::cout << "Result: [" << output[i] << ", ";
			else if (i < vectorName.size() - 1)
				std::cout << output[i] << ", ";
			else
				std::cout << output[i] << "]" << std::endl;
		}
	}
	template <typename T> void elementwise(DYN_C2D<T>& output, DYN_C2D<T>& arrayName, DYN_C2D<T>& arrayName2) {

		// Guard to stop vectors/1D arrays of different lengths being multiplied together 
		size_guard<T>(arrayName, arrayName2, "Elementwise");

		// Instantiate reserve space for vector for optimizing performance (reduces copying) 
		output.mem_corrector(arrayName.size_rows(), arrayName.size_cols());

		// Loop along elements of array
		size_t contig_ind = 0;
		for (auto data_addr = output.begin(); data_addr != output.end(); data_addr++, contig_ind++) {
			*data_addr = arrayName(contig_ind) * arrayName2(contig_ind);
		}
	}

	/* max() find the maximum value of an array or between two values*/
	template <typename T> T max(std::vector<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "max()");

		// Set current max value as first element of vector for reference 
		T max_val = dyn_arr_data[0];

		// Loop over data, start at element 1. 
		for (auto data_addr = dyn_arr_data.begin() + 1; data_addr != dyn_arr_data.end(); data_addr++) {
			if (*data_addr > max_val) {
				max_val = *data_addr;
			}
		}

		return max_val;
	}
	template <typename T> T max(DYN_C2D<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "max()");

		// Set current max value as first element of vector for reference (using contig index)
		T max_val = dyn_arr_data(0);

		// Loop over data, start at element 1. 
		for (auto data_addr = dyn_arr_data.begin() + 1; data_addr != dyn_arr_data.end(); data_addr++) {
			if (*data_addr > max_val) {
				max_val = *data_addr;
			}
		}

		return max_val;
	}

	/* min() */
	template <typename T> T min(std::vector<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "min()");

		// Set current min value as first element of vector for reference 
		T min_val = dyn_arr_data[0];

		// Loop over data, start at element 1. 
		for (auto data_addr = dyn_arr_data.begin() + 1; data_addr != dyn_arr_data.end(); data_addr++) {
			if (*data_addr < min_val) {
				min_val = *data_addr;
			}
		}

		return min_val;
	}

	/* find() find an element in an array which equals a value and return an logical array */
	template <typename T> void find(T value, std::vector<T>& dyn_arr_data, std::vector<bool>& logic_arr) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "find()");

		// Optimize memory for locical vector 
		mem_corrector(logic_arr, dyn_arr_data.size());

		size_t index = 0;
		for (auto data_addr = dyn_arr_data.begin(); data_addr != dyn_arr_data.end(); data_addr++, index++) {
			if (*data_addr == value) {
				logic_arr[index] = true;
			}
		}
	}
	template <typename T> void find(T value, DYN_C2D<T>& dyn_arr_data, DYN_C2D<int>& logic_arr) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "find()");

		// Optimize memory for logical array in correspondance of dyn_arr_data 
		logic_arr.mem_corrector(dyn_arr_data.size_rows(), dyn_arr_data.size_cols());

		size_t contig_ind = 0;
		for (auto data_addr = logic_arr.begin(); data_addr != logic_arr.end(); data_addr++, contig_ind++) {
			if (dyn_arr_data(contig_ind) == value) {
				*data_addr = true;
			}
		}
	}

	/* sum () sums all element of an array (1 and 2D)*/
	template <typename T> T		sum(std::vector<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "sum()");

		T sum_val = 0;
		for (auto data_addr = dyn_arr_data.begin(); data_addr != dyn_arr_data.end(); data_addr++) {
			sum_val += *data_addr;
		}
		return sum_val;
	}
	template <typename T> void	sum(DYN_C2D<T>& dyn_arr_data, std::vector<T>& sum_result, bool sum_row_or_col) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "sum()");

		// Settign inputted for row-wise or col-wise summation; true = row-wise, false = colwise
		switch (sum_row_or_col) {
		case (false): {		// Col-wise

			// Optimize memory for output vector 
			mem_corrector(sum_result, dyn_arr_data.size_rows());

			for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
				T sum_val = 0;
				for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {
					sum_val += dyn_arr_data(i, j);
				}
				sum_result[i] = sum_val;
			}
			break;
		}
		case (true): {		// Row-wise

			// Optimize memory for output vector 
			mem_corrector(sum_result, dyn_arr_data.size_cols());

			for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {
				T sum_val = 0;
				for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
					sum_val += dyn_arr_data(i, j);
				}
				sum_result[j] = sum_val;
			}
			break;
		}
		}

	}

	/* mean() */
	template <typename T> double mean(std::vector<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "mean()");

		double mean_val = 0;

		// Loop over vector using an iterator 
		for (auto it = dyn_arr_data.begin(); it != dyn_arr_data.end(); it++) {
			mean_val += *it;
		}

		return mean_val = mean_val / (double)dyn_arr_data.size();
	}

	/* eye() -> generate an identity matrix */
	template <typename T> DYN_C2D<T> eye(size_t rows, size_t cols) {

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (rows != cols) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array in function, eye().\n";
			exit(-1);
		}

		DYN_C2D<T> dyn_arr_data(rows, cols);

		// Counter for diagonal element in relation to row index
		size_t diag_count = 0;

		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					dyn_arr_data(i, j) = 1;
				}
				else {
					dyn_arr_data(i, j) = 0;
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}

		return dyn_arr_data;
	}
	template <typename T> void eye(DYN_C2D<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "eye()");

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (dyn_arr_data.size_rows() != dyn_arr_data.size_cols()) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array in function, eye().\n";
			exit(-1);
		}

		// Counter for diagonal element in relation to row index
		size_t diag_count = 0;

		for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
			for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					dyn_arr_data(i, j) = 1;
				}
				else {
					dyn_arr_data(i, j) = 0;
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}
	}

	/* diag() -> generate a diagonal array from a vector/ or extract a vector from a diagonal array */
	template <typename T> void diag(DYN_C2D<T>& dyn_arr_data, T diag_val) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "diag()");

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (dyn_arr_data.size_rows() != dyn_arr_data.size_cols()) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array in function, diag().\n";
			exit(-1);
		}

		// Counter for diagonal element in relation to row index
		size_t diag_count = 0;

		for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
			for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					dyn_arr_data(i, j) = diag_val;
				}
				else {
					dyn_arr_data(i, j) = 0;
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}
	}
	template <typename T> void diag(DYN_C2D<T>& dyn_arr_data, std::vector<T>& vectorInput) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "diag()");

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (dyn_arr_data.size_rows() != dyn_arr_data.size_cols()) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array in function, diag().\n";
			exit(-1);
		}

		// Guard for input vector not being same length as square array (where rows = cols) 
		if (vectorInput.size() != dyn_arr_data.size_rows()) {
			std::cout << "Error: Input vector size not similar to array row/col size in function, diag().\n";
			exit(-1);
		}

		// Counter for diagonal element in relation to row index
		size_t diag_count = 0;

		for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
			for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					dyn_arr_data(i, j) = vectorInput[diag_count];
				}
				else {
					dyn_arr_data(i, j) = 0;
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}
	}
	template <typename T> void diag_extract(std::vector<T>& vectorOutput, DYN_C2D<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "diag()");

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (dyn_arr_data.size_rows() != dyn_arr_data.size_cols()) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array in function, diag().\n";
		}

		// Optimize memory for the output vector 
		mem_corrector(vectorOutput, dyn_arr_data.size_rows());

		// Counter for diagonal element in relation to row index
		size_t diag_count = 0;

		for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
			for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					vectorOutput[diag_count] = (dyn_arr_data(i, j));
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}

	}
	template <typename T> std::vector<T> diag_extract(DYN_C2D<T>& dyn_arr_data) {

		// Guard for zero vector 
		size_guard(dyn_arr_data, "diag()");

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (dyn_arr_data.size_rows() != dyn_arr_data.size_cols()) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array in function, diag().\n";
			exit(-1);
		}

		// Establish vector
		std::vector<T> vectorOutput(dyn_arr_data.size_rows());

		// Counter for diagonal element in relation to row index
		size_t diag_count = 0;

		for (size_t i = 0; i < dyn_arr_data.size_rows(); i++) {
			for (size_t j = 0; j < dyn_arr_data.size_cols(); j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					vectorOutput[diag_count] = (dyn_arr_data(i, j));
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}
		return vectorOutput;

	}

	/* flatten() -> flatten a 2D array (rows cols) onto a 1D array (cols) */
	template <typename T> void flatten(std::vector<T>& flattened_arr, DYN_C2D<T>& dyn_arr_data) {

		// Optomize/correct memory for flattened array vector 
		mem_corrector(flattened_arr, dyn_arr_data.size_rows() * dyn_arr_data.size_cols());
		size_t contig_ind = 0;

		// Loop over 2D data using iterator for speed and row and col info is required to be in contiguous format
		for (auto it = dyn_arr_data.begin(); it != dyn_arr_data.end(); it++, contig_ind++) {
			flattened_arr[contig_ind] = *it;
		}
	}
	template <typename T> std::vector<T> flatten(DYN_C2D<T>& dyn_arr_data) {

		// Establish vector
		std::vector<T> flattened_arr(dyn_arr_data.size_rows() * dyn_arr_data.size_cols());
		size_t contig_ind = 0;

		// Loop over 2D data using iterator for speed and row and col info is required to be in contiguous format
		for (auto it = dyn_arr_data.begin(); it != dyn_arr_data.end(); it++, contig_ind++) {
			flattened_arr[contig_ind] = *it;
		}

		return flattened_arr;
	}

	/* reshape() -> basically flips rows and collumn values; i.e. 2x3 array to a 3x2 array. */
	template <typename T> void reshape(DYN_C2D<T>& dyn_arr_data) {

		// Since data is contiguous, only need to adjust row/col values!
		size_t container = dyn_arr_data.size_rows();

		dyn_arr_data.setRefRowsize(dyn_arr_data.size_cols());
		dyn_arr_data.setRefcolsize(container);

	}

	/* l_interp -> linear interpolation between points. Note, integer not allowed as there is division in this function */
	// Linear interpolation of array data 
	void l_interp(DYN_C2D<float>& x_y_data, DYN_C2D<float>& x_y_newdata) {

		// Between point linear interp 
		// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);
	}
	void l_interp(DYN_C2D<double>& x_y_data, DYN_C2D<double>& x_y_newdata);

	/* flipdata -> flips data left to right, up to down etc - left to right -> settign 1, up to down setting 2 */
	template <typename T> void fliplr(std::vector<T>& dyn_arr_data) {

		const size_t cols = dyn_arr_data.size();

		if (isEven(cols) == true) {

			// Buffer to contain half the data (Instead of copying the full array proccessed to a new array and back again, 
			// copy half, overwrite the half on actual array and use rthe copy for the other half

			const size_t hd_cols = cols / 2;						// Pre-define cols/2 so only need to calculate once 
			std::vector<T> hd_buffer;		 						// Allocate buffer memory 
			hd_buffer.reserve(hd_cols);
			hd_buffer.resize(hd_cols);

			// hd_buffer[hd_cols];						

			for (int j = 0; j < hd_cols; j++) {						// Loop over half data

				hd_buffer[j] = dyn_arr_data[cols - j - 1];			// Copy data onto buffer (-1 for correct index number (e.g. cols = 10, last index = 9)
				dyn_arr_data[cols - j - 1] = dyn_arr_data[j];		// Overwrite end point in which buffer has just copied; 
				dyn_arr_data[j] = hd_buffer[j];						// Overwrite start point with buffer. 
			}
		}
		else {
			// For odd dimensions of array, mid point stays the same (e.g. array[9], index [4] does not change)

			const size_t hd_cols = (cols - 1) / 2;					// Pre-define cols/2 so only need to calculate once (-1 for element number on each side of midpt)
			std::vector<T> hd_buffer;		 						// Allocate buffer memory (size = elements on either side of mpt)
			hd_buffer.reserve(hd_cols);
			hd_buffer.resize(hd_cols);

			for (int j = 0; j < hd_cols; j++) {						// Loop over half data

				hd_buffer[j] = dyn_arr_data[cols - j - 1];				// Copy data onto buffer (-1 for correct index number (e.g. cols = 9, last index = 8)
				dyn_arr_data[cols - j - 1] = dyn_arr_data[j];				// Overwrite end point in which buffer has just copied; 
				dyn_arr_data[j] = hd_buffer[j];							// Overwrite start point with buffer. 
			}
		}
	}
	template <typename T> void fliplr(DYN_C2D<T>& dyn_arr_data);
	template <typename T> void flipud(DYN_C2D<T>& dyn_arr_data);

	/* concatenate() stitch together two arrays to produce one big array */
	// Row-wise concatenation 1D array (2* 1D array -> 2D array)
	template <typename T> void concatenate(std::vector<T>& dyn_arr_data, std::vector<T>& dyn_arr_data2, DYN_C2D<T>& concatenated_arr);

	// Column-wise concatenation (1D array[cols] +  1D array[cols_ -> 1D[cols + cols_] array)
	template <typename T> void concatenate(std::vector<T>& dyn_arr_data, std::vector<T>& dyn_arr_data2, std::vector<T>& concatenated_arr);

	// Row-wise concatentation 2D arrays (2D array[rows][cols] + 2D array[rows_][cols] -> 2D array[rows +rows_][cols]
	// Collumn-wise concatenation 2D arrays (2D array [rows][cols] + 2D array[rows][cols_] -> 2D array[rows][cols + cols_] 
	template <typename T> void concatenate(DYN_C2D<T>& dyn_arr_data, DYN_C2D<T>& dyn_arr_data2, DYN_C2D<T>& concatenated_arr, bool row_or_col_wise);


} // End of namespace RMF