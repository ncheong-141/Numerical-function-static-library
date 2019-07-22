#pragma once
/* Namespace for all compile time functions */

// External modules
#include <iostream>
#include <vector>
#include <chrono>

namespace instrumentation {

	class Timer {
	private:
		const char* stackframeName;	// Name of stackframe for timer

		// Assign start and end variables to memory with weird type
		std::chrono::time_point<std::chrono::steady_clock> start, end;
		std::chrono::duration<float> duration;		
	
	public:

		// Constructor (start timer here)
		Timer(const char* stackframe) : stackframeName(stackframe), duration(0) {
			start = std::chrono::high_resolution_clock::now();
		}

		// Destructor (end timer here) 
		~Timer() {
			end = std::chrono::high_resolution_clock::now();  	// end the timer
			duration = end - start;					// Calculate duration
			std::cout << "Time taken for " << stackframeName << ": " << duration.count() << " s\n";
		}
	};
}

/* Declare general function prototypes*/
// Template function definitions. Implementation and specification in .cpp file
/* NOTE: 		-- Non-type template parameters so impossible to forward declare/ put in implementation--*/
/* NOTE: - Passing in (&arr_data)[] is passing in a static array by reference. This avoids losing array dimension information
		   so the compiler can check, which, is safe since a flag will appear if user inputs wrong dimensions.
		   It passes a pointer to the whole array.	*/	
/* NOTE: - 3D arrays in C++ are an array of 2D arrays where the 2D array number (depth) is the first specified index => arr[depth][row][col]*/
/* NOTE: - Passing in return array by reference is good as it keeps memory on calling stackframe*/
/* NOTE: - template functions with non-type parameters is good cause it enforces compile time checks on array sizes inputted/ 
		   ensures they are the same dimensions so no mismatch and memory access violations. */
namespace CMF {

	/* Initializing stack array with parameter constructor (Not really but fills it in for you)*/
	// template<typename T, int cols>		void stack_array_params(T(&arr_data)[cols], )
	/* ----------------------- General math tools -----------------------------------------*/

	/* any() -> if any element is non-zero return true*/
	template<typename T, int cols> bool any(const T(&arr_data)[cols])
	{
		for (int i = 0; i < cols; i++) {
			if (arr_data[i] > 0)
				return true;
		}
		return false;
	}
	template<typename T, int rows, int cols> bool any(const T(&arr_data)[rows][cols]) 						// any() if any element is non-zero
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++) {
				if (arr_data[i][j] > 0)
					return true;
			}
		}

		return false;
	}
	template<typename T, int rows, int cols, int depth> bool any(const T(&arr_data)[depth][rows][cols]) {
		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < cols; k++) {
					if (arr_data[i][j][k] > 0)
						return true;
				}
			}
		}
		return false;
	}

	/* all() -> if all elements are non-zero return true*/
	template<typename T, int cols> bool all(const T(&arr_data)[cols])
	{
		for (int i = 0; i < cols; i++) {
			if (arr_data[i] != 0)
				return false;
		}
		return true;
	}
	template<typename T, int rows, int cols> bool all(const T(&arr_data)[rows][cols]) 						// any() if any element is non-zero
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++) {
				if (arr_data[i][j] != 0)
					return false;
			}
		}

		return true;
	}
	template<typename T, int rows, int cols, int depth> bool all(const T(&arr_data)[depth][rows][cols])
	{
		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < cols; k++) {
					if (arr_data[i][j][k] != 0)
						return false;
				}
			}
		}
		return true;
	}

	/*isEven, isOdd*/
	bool isEven(const int& parameter);
	bool isOdd(const int& parameter);

	/* abs() -> Returns an absolute value*/
	template<typename T> T abs(T parameter) {

		if (parameter < 0)
			return (-1 * parameter);

		return parameter;
	}
	
	/* arange() fills array with incrementing values (+- 1) from specified start to end*/
	template<typename T, int cols> void arange(T(&arr_data)[cols], const int start, const int end) {
		/* Guard for incorrect array */
		if (abs(end - start) != cols) {
			std::cout << "Error: Array dimensions do not match with start and end values\n"; exit(-1);
		}

		int increment = 0;
		if (end > start) {
			for (int i = 0; i < cols; i++) {
				arr_data[i] = start + increment;
				increment++;
			}
		}
		else if (start > end) {
			for (int i = 0; i < cols; i++) {
				arr_data[i] = start + increment;
				increment--;
			}
		}
	}
	// Generates values between start and end according to discrete parameter (Note, cols must have an extra element for end point) 
	template<int cols> void arange(float(&arr_data)[cols], const float start, const float end, const float NP) {

		// Guard for mismatch of array dimensions 
		if (cols != NP) {
			std::cout << "Error: Array dimensions mismatch in fn 'arange()'. Discretizing between start and end with current number of points (NP) yeilds array dimensions greater or lower\n";
			exit(-1);
		}

		float total_length = abs(end - start);		// Calculate total value between start and end
		float disc_inc = total_length / (NP - 1);       // Calculate the increments between each element
		float disc_val = 0;				// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (int j = 0; j < cols; j++) {

			arr_data[j] = disc_val;			// Write element
			disc_val += disc_inc;			// Add increment
		}
	}
	template<int cols> void arange(double(&arr_data)[cols], const double start, const double end, const double NP) {

		// Guard for mismatch of array dimensions 
		if (cols != NP) {
			std::cout << "Error: Array dimensions mismatch. Discretizing between start and NP yeilds array dimensions greater or lower\n";
			exit(-1);
		}

		double total_length = abs(end - start);		// Calculate total value between start and end
		double disc_inc = total_length / (NP - 1);      // Calculate the increments between each element
		double disc_val = 0;				// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (int j = 0; j < cols; j++) {

			arr_data[j] = disc_val;			// Write element
			disc_val += disc_inc;			// Add increment
		}
	}

	/* tile() fills array with a specified value*/
	template<typename T, int cols> void tile(T tilevalue, T(&arr_data)[cols]) {
		for (int i = 0; i < cols; i++) {
			arr_data[i] = tilevalue;
		}
	}
	template<typename T, int rows, int cols> void tile(T tilevalue, T(&arr_data)[rows][cols]) { 						// any() if any element is non-zero
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				arr_data[i][j] = tilevalue;
			}
		}
	}
	template<typename T, int rows, int cols, int depth> void tile(T tilevalue, T(&arr_data)[depth][rows][cols]) {
		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < cols; k++) {
					arr_data[i][j][k] = tilevalue;
				}
			}
		}
	}

	/* max() find the maximum value of an array or between two values*/
	template <typename T> T max(T& value, T& value2) {	//  Maximum value between two variables
		if (value2 > value) {
			return value2;
		}
		return value;
	}
	template <typename T, int cols> T max(T(&arr_data)[cols]) {

		T maxval = arr_data[0];			// Set current max value as first element

		for (int i = 1; i < cols; i++) {	// Start at 1 as first element is currently max value
			if (arr_data[i] > maxval) { maxval = arr_data[i]; }
		}
		return maxval;
	}
	template <typename T, int rows, int cols> T max(T(&arr_data)[rows][cols]) {

		T maxval = arr_data[0][0];	// Set current max value as first element

		for (int i = 0; i < rows; i++) {	// Start at 1 as first element is currently max value
			for (int j = 0; i < cols; j++) {
				if (arr_data[i][j] > maxval) { maxval = arr_data[i][j]; }
			}
		}
		return maxval;
	}
	template <typename T, int rows, int cols, int depth> T max(T(&arr_data)[depth][rows][cols]) {

		T maxval = arr_data[0][0][0];	// Set current max value as first element

		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < cols; k++) {
					if (arr_data[i][j][k] > maxval) { maxval = arr_data[i][j][k]; }
				}
			}
		}
		return maxval;
	}

	/* find() find an element in an array which equals a value and return an logical array */
	template <typename T, int cols>	void find(T& value, T(&arr_data)[cols], bool(&logic_arr)[cols]) {
		for (int i = 0; i < cols; i++) {
			if (arr_data[i] == value) { logic_arr[i] = true; }
			else { logic_arr[i] = false; }
		}
	}
	template <typename T, int rows, int cols> void find(T& value, T(&arr_data)[rows][cols], bool(&logic_arr)[rows][cols]) {

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (arr_data[i][j] == value) { logic_arr[i][j] = true; }
				else { logic_arr[i][j] = false; }
			}
		}

	}
	template <typename T, int rows, int cols, int depth> void find(T& value, T(&arr_data)[depth][rows][cols], bool(&logic_arr)[depth][rows][cols]) {
		for (int i = 0; i < depth; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < cols; k++) {
					if (arr_data[i][j][k] == value) { logic_arr[i][j][k] = true; }
					else { logic_arr[i][j][k] = false; }
				}
			}
		}

	}

	/* sum () sums all element of an array (1 and 2D)*/
	// 1D summation of all elements on collumns (return an integer) 
	template <typename T, int cols>	 T sum(T(&arr_data)[cols]) {

		T sum_val = 0; 					// Establish container for sum value
		for (int j = 0; j < cols; j++) {		// Loop over 1D array and sum components 
			sum_val += arr_data[j];
		}
		return sum_val;
	}
	// 2D summation of all elements along ->collumns<- (returns a 1D array)
	template <typename T, int rows, int cols> void sum(T(&arr_data)[rows][cols], T(&sum_result)[cols], int sum_row_or_col) {

		// Guard to stop invalid setting inputs
		if (sum_row_or_col != 0 && sum_row_or_col != 1) {
			std::cout << "Error: Invalid setting inputted; input 1 or 0 for collumn or row summation respectively.\n"; exit(-1);
		}
		// Switch statement for "row"(0) or "col"(1) input. 
		switch (sum_row_or_col) {
		case (1): { // Collumn summation
			for (int i = 0; i < rows; i++) {

				T sum_val = 0; 					// Establish container for sum value
				for (int j = 0; j < cols; j++) {		// Loop over 1D array and sum components 
					sum_val += arr_data[i][j];
				}
				sum_result[i] = sum_val;			// Write to the result
			}
			break;
		}
		case (0): { // Row summation 
			for (int j = 0; j < cols; j++) {

				T sum_val = 0; 					// Establish container for sum value
				for (int i = 0; i < rows; i++) {		// Loop over 1D array and sum components 
					sum_val += arr_data[i][j];
				}
				sum_result[j] = sum_val;			// Write to the result
			}
			break;
		}
		}
	}

	/* size_rows(), size_cols()- gets dimensions of array*/
	// 1D array size
	template<typename T, int cols> int size_cols(T(&arr_data)[cols]) {
		return cols;
	}
	// 2D array size 
	template<typename T, int rows, int cols> int size_cols(T(&arr_data)[rows][cols]) {
		return cols;
	}
	template<typename T, int rows, int cols> int size_rows(T(&arr_data)[rows][cols]) {
		return rows;
	}
	// 3D array size
	template <typename T, int rows, int cols, int depth> int size_cols(T(&arr_data)[depth][rows][cols]) {
		return cols;
	}
	template <typename T, int rows, int cols, int depth> int size_rows(T(&arr_data)[depth][rows][cols]) {
		return rows;
	}
	template <typename T, int rows, int cols, int depth> int size_depth(T(&arr_data)[depth][rows][cols]) {
		return depth;
	}

	/* mean() - the mean.. */
	template <typename T, int cols> double  mean(T(&arr_data)[cols]) {

		// Sum all elements
		double sum_val = 0; 				// Establish container for sum value
		for (int j = 0; j < cols; j++) {		// Loop over 1D array and sum components 
			sum_val += arr_data[j];
		}
		return sum_val / cols;
	}

	/* eye() -> generate an identity matrix */
	template <typename T, int rows, int cols> void eye(T(&output_eye_arr)[rows][cols]) {

		// Guard for incorrect array dimensions (rows must = cols for identity arrays) 
		if (rows != cols) {
			std::cout << "Error: Array dimensions mismatch. Must be a square array.\n";
		}
		// Counter for diagonal element in relation to row index
		int diag_count = 0;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {

				// Condition to set 1 for identity element
				if (j == diag_count) {
					output_eye_arr[i][j] = 1;
				}
				else {
					output_eye_arr[i][j] = 0;
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}
	}

	/* diag() -> generate a diagonal array from a vector/ or extract a vector from a diagonal array */
	template <typename T, int rows, int cols> void diag(T(&input_vect)[cols], T(&output_diag)[rows][cols]) {

		// Guard for incorrect array dimensions (rows must = cols for diagonal arrays) 
		if (rows != cols) {
			std::cout << "Error: Array dimensions mismatch. Diagonal array must be a square array.\n";
		}

		// Counter for diagonal element in relation to row index
		int diag_count = 0;

		// Note, for diag array both rows and cols must be equal 
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {

				// Condition to set diag element
				if (j == diag_count) {
					output_diag[i][j] = input_vect[i];
				}
				else {
					output_diag[i][j] = 0;
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}
	}
	template <typename T, int rows, int cols> void diag_extract(T(&input_diag)[rows][cols], T(&output_vect)[cols]) {

		// Guard for incorrect array dimensions (rows must = cols for diagonal arrays) 
		if (rows != cols) {
			std::cout << "Error: Array dimensions mismatch. Diagonal arrya must be a square array.\n";
		}

		// Counter for diagonal element in relation to row index
		int diag_count = 0;

		// Note, for diag array both rows and cols must be equal 
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {

				// Extract diagonal
				if (j == diag_count) {
					output_vect[i] = input_diag[i][j];
				}
			}

			// Increment diag counter for next row 
			diag_count++;
		}

	}

	/* flatten() -> flatten a 2D array (rows cols) onto a 1D array (cols) */
	template <typename T, int rows, int cols> void flatten(T(&input_flatten_arr)[rows][cols], T(&output_vect)[rows * cols]) {

		// Offset for establishing index position for each row on the vector 
		int vect_offset = 0;

		// Loop over array to be flattened 
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {

				output_vect[j + vect_offset] = input_flatten_arr[i][j];
			}

			// Add the collumn value to offset data of next implementation of row data.
			vect_offset += cols;
		}
	}

	/* reshape() -> basically flips rows and collumn values; i.e. 2x3 array to a 3x2 array. */
	// 1D arrays to 2D
	template <typename T, int cols> void reshape(T(&arr_data)[cols], T(&reshaped_arr)[cols][1]) {

		// Loop over array
		for (int j = 0; j < cols; j++) {
			reshaped_arr[j][0] = arr_data[j];
		}
	}
	template <typename T, int cols> void reshape(T(&arr_data)[cols], T(&reshaped_arr)[1][cols]) {

		// Loop over array
		for (int j = 0; j < cols; j++) {
			reshaped_arr[0][j] = arr_data[j];
		}
	}
	// 2D array reshape
	template <typename T, int rows, int cols> void reshape(T(&arr_data)[rows][cols], T(&reshaped_arr)[cols][rows]) {

		// Loop over array
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				reshaped_arr[j][i] = arr_data[i][j];
			}
		}
	}

	/* l_interp -> linear interpolation between points. Note, integer not allowed as there is division in this function */
	// Linear interpolation between two points
	void l_interp(float(&pt0)[2], float(&pt1)[2], float(&out_x_y)[2]);
	void l_interp(double(&pt0)[2], double(&pt1)[2], double(&out_x_y)[2]);

	// Linear interpolation which returns a value 
	float  l_interp(float(&pt0)[2], float(&pt1)[2], float x);
	double l_interp(double(&pt0)[2], double(&pt1)[2], double x);

	// Linear interpolation of array data 
	template <int cols, int cols_> void l_interp(float(&xdat)[cols], float(&ydat)[cols], float(&xnew)[cols_], float(&ynew)[cols_]) {

		// Floating pt precision is minimum 6 significant digits by the IEEE 754 standard (ref, learncpp.com)
		// Guard if discretized data is smaller/larger than the precision
		if (((abs(xdat[cols - 1] - xdat[0]) / cols) < 1e-6) || ((abs(xnew[cols_ - 1] - xnew[0]) / cols) < 1e-6))
			std::cout << "Warning: discretized data in function 'l_interp' is smaller than minimum finite precision (6 significant digits) of floating point numbers.";

		// Guard if larger
		if (((abs(xdat[cols - 1] - xdat[0]) / cols) > 1e6) || ((abs(xnew[cols_ - 1] - xnew[0]) / cols) > 1e6))
			std::cout << "Warning: discretized data in function 'l_interp' is larger than minimum finite precision (6 significant digits) of floating point numbers.";

		// Visit every pt 0 and 1 on xdat and find corresponding pts of xnew between the points.
		for (int j = 0; j < cols; j++) {

			//float x_0 = xdat[j]; 
			//float x_1 = xdat[j + 1];
			//float y_0 = ydat[j];
			//float y_1 = ydat[j + 1];

			// Bool variable to establish if point has been found for loop efficiency 
			// (dont keep looping if you don't expect anymore points as xdat is monotonically increasing)
			bool found = false;

			// Loop over xnew data to establish pts inbetween x0 and x1 and calculate 
			for (int f_ind = 0; f_ind < cols_; f_ind++) {
				// If xnew is inbetween x0 and x1 
				if ((xnew[f_ind] < xdat[j + 1]) && (xnew[f_ind] > xdat[j])) {

					// Linear interp and store 
					// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);

					ynew[f_ind] = ((ydat[j] * xdat[j + 1]) - (ydat[j] * xnew[f_ind]) + (ydat[j + 1] * xnew[f_ind]) - (ydat[j + 1] * xdat[j])) / (xdat[j + 1] - xdat[j]);

					// Set variable to 
					found = true;
				}
				// if xnew = x0
				else if (xnew[f_ind] == xdat[j]) {
					ynew[f_ind] = ydat[j];

					// Set variable to 
					found = true;
				}
				// if xnew = x1
				else if (xnew[f_ind] > 0.99f * xdat[j + 1] && xnew[f_ind] < 1.01f * xdat[j + 1]) {
					ynew[f_ind] = ydat[j + 1];

					// Set variable to 
					found = true;
				}
				else if (found == true) {

					// Break from loop to stop executing it if pts have been found
					break;
				}


			}

		}
	}
	template <int cols, int cols_> void l_interp(double(&xdat)[cols], double(&ydat)[cols], double(&xnew)[cols_], double(&ynew)[cols_]) {

		// Visit every pt 0 and 1 on xdat and find corresponding pts of xnew between the points.
		for (int j = 0; j < cols; j++) {

			//float x_0 = xdat[j]; 
			//float x_1 = xdat[j + 1];
			//float y_0 = ydat[j];
			//float y_1 = ydat[j + 1];

			// Bool variable to establish if point has been found for loop efficiency 
			// (dont keep looping if you don't expect anymore points as xdat is monotonically increasing)
			bool found = false;

			// Loop over xnew data to establish pts inbetween x0 and x1 and calculate 
			for (int f_ind = 0; f_ind < cols_; f_ind++) {
				// If xnew is inbetween x0 and x1 
				if ((xnew[f_ind] < xdat[j + 1]) && (xnew[f_ind] > xdat[j])) {

					// Linear interp and store 
					// y = ((y0 * x1) - (y0 * x) + (y1 * x) - (y1 * x0)) / (x1 - x0);

					ynew[f_ind] = ((ydat[j] * xdat[j + 1]) - (ydat[j] * xnew[f_ind]) + (ydat[j + 1] * xnew[f_ind]) - (ydat[j + 1] * xdat[j])) / (xdat[j + 1] - xdat[j]);

					// Set variable to 
					found = true;
				}
				// if xnew = x0
				else if (xnew[f_ind] == xdat[j]) {
					ynew[f_ind] = ydat[j];

					// Set variable to 
					found = true;
				}
				// if xnew = x1
				else if (xnew[f_ind] > 0.99 * xdat[j + 1] && xnew[f_ind] < 1.01 * xdat[j + 1]) {
					ynew[f_ind] = ydat[j + 1];

					// Set variable to 
					found = true;
				}
				else if (found == true) {

					// Break from loop to stop executing it if pts have been found
					break;
				}


			}

		}
	}

	/* flipdata -> flips data left to right, up to down etc - left to right -> settign 1, up to down setting 2 */
	template <typename T, int cols>	void fliplr(T(&arr_data)[cols]) {

		if (isEven(cols) == true)
		{
			// Buffer to contain half the data (Instead of copying the full array proccessed to a new array and back again, 
			// copy half, overwrite the half on actual array and use rthe copy for the other half

			const int hd_cols = cols / 2;					// Pre-define cols/2 so only need to calculate once 
			T		  hd_buffer[hd_cols];				// Allocate buffer memory 

			for (int j = 0; j < hd_cols; j++) {				// Loop over half data

				hd_buffer[j] = arr_data[cols - j - 1];	// Copy data onto buffer (-1 for correct index number (e.g. cols = 10, last index = 9)
				arr_data[cols - j - 1] = arr_data[j];			// Overwrite end point in which buffer has just copied; 
				arr_data[j] = hd_buffer[j];				// Overwrite start point with buffer. 
			}


		}
		else { // (Odd amount of collumns)

			// For odd dimensions of array, mid point stays the same (e.g. array[9], index [4] does not change)

			const int hd_cols = (cols - 1) / 2;				// Pre-define cols/2 so only need to calculate once (-1 for element number on each side of midpt)
			T		  hd_buffer[hd_cols];				// Allocate buffer memory (size = elements on either side of mpt)

			for (int j = 0; j < hd_cols; j++) {				// Loop over half data

				hd_buffer[j] = arr_data[cols - j - 1];			// Copy data onto buffer (-1 for correct index number (e.g. cols = 9, last index = 8)
				arr_data[cols - j - 1] = arr_data[j];			// Overwrite end point in which buffer has just copied; 
				arr_data[j] = hd_buffer[j];				// Overwrite start point with buffer. 
			}

		}
	}
	template <typename T, int rows, int cols> void fliplr(T(&arr_data)[rows][cols]) {

		if (isEven(cols) == true)
		{
			// Buffer to contain half the data (Instead of copying the full array proccessed to a new array and back again, 
			// copy half, overwrite the half on actual array and use rthe copy for the other half

			const int hd_cols = cols / 2;					// Pre-define cols/2 so only need to calculate once 
			T		  hd_buffer[hd_cols];				// Allocate buffer memory 

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < hd_cols; j++) {			// Loop over half data

					hd_buffer[j] = arr_data[i][cols - j - 1];	// Copy data onto buffer (-1 for correct index number (e.g. cols = 10, last index = 9)
					arr_data[i][cols - j - 1] = arr_data[i][j];	// Overwrite end point in which buffer has just copied; 
					arr_data[i][j] = hd_buffer[j];			// Overwrite start point with buffer. 
				}
			}
		}
		else { // (Odd amount of collumns)

			// For odd dimensions of array, mid point stays the same (e.g. array[9], index [4] does not change)

			const int hd_cols = (cols - 1) / 2;					// Pre-define cols/2 so only need to calculate once (-1 for element number on each side of midpt)
			T		  hd_buffer[hd_cols];					// Allocate buffer memory (size = elements on either side of mpt)

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < hd_cols; j++) {				// Loop over half data

					hd_buffer[j] = arr_data[i][cols - j - 1];		// Copy data onto buffer (-1 for correct index number (e.g. cols = 9, last index = 8)
					arr_data[i][cols - j - 1] = arr_data[i][j];		// Overwrite end point in which buffer has just copied; 
					arr_data[i][j] = hd_buffer[j];				// Overwrite start point with buffer. 
				}
			}
		}

	}
	template <typename T, int rows, int cols> void flipud(T(&arr_data)[rows][cols]) {


		if (isEven(rows) == true)
		{
			// Buffer to contain half the data (Instead of copying the full array proccessed to a new array and back again, 
			// copy half, overwrite the half on actual array and use rthe copy for the other half

			const int hd_rows = rows / 2;					// Pre-define cols/2 so only need to calculate once 
			T		  hd_buffer[hd_rows];				// Allocate buffer memory 

			for (int i = 0; i < hd_rows; i++) {
				for (int j = 0; j < cols; j++) {			// Loop over half data

					hd_buffer[i] = arr_data[rows - i - 1][j];	// Copy data onto buffer (-1 for correct index number (e.g. cols = 10, last index = 9)
					arr_data[rows - i - 1][j] = arr_data[i][j];	// Overwrite end point in which buffer has just copied; 
					arr_data[i][j] = hd_buffer[i];			// Overwrite start point with buffer. 
				}
			}
		}
		else { // (Odd amount of collumns)

			// For odd dimensions of array, mid point stays the same (e.g. array[9], index [4] does not change)

			const int hd_rows = (rows - 1) / 2;				// Pre-define cols/2 so only need to calculate once 
			T		  hd_buffer[hd_rows];				// Allocate buffer memory 

			for (int i = 0; i < hd_rows; i++) {
				for (int j = 0; j < cols; j++) {			// Loop over half data

					hd_buffer[i] = arr_data[rows - i - 1][j];	// Copy data onto buffer (-1 for correct index number (e.g. cols = 10, last index = 9)
					arr_data[rows - i - 1][j] = arr_data[i][j];	// Overwrite end point in which buffer has just copied; 
					arr_data[i][j] = hd_buffer[i];			// Overwrite start point with buffer. 
				}
			}
		}

	}
	
	/* TO DO:
	inverse, transpose, determinant, eigenvalue, reshape, macro for concatenating [a,b] [a;b], fliplr, flipud,
	size, length, flatten, eye, diag, mean, median, mode, standard dev, varience, polyfit, linear interp, c*/

	/* concatenate() stitch together two arrays to produce one big array */
	// Row-wise concatenation 1D array (2* 1D array -> 2D array)
	template <typename T, int cols>	void concatenate(T(&arr_data)[cols], T(&arr_data2)[cols], T(&concatted_arr)[2][cols]) {
		for (int j = 0; j < cols; j++) {
			concatted_arr[0][j] = arr_data[j];
			concatted_arr[1][j] = arr_data2[j];
		}
	}

	// Row-wise concatentation 2D arrays (2D array[rows][cols] + 2D array[rows_][cols] -> 2D array[rows +rows_][cols]
	template <typename T, int rows, int rows_, int cols>	void concatenate(T(&arr_data)[rows][cols], T(&arr_data2)[rows_][cols], T(&concatted_arr)[rows + rows_][cols]) {

		// Establish which array has the largest row dimension. 
		// First array has biggest row dimension (arr_data)
		if (rows > rows_) {

			// Set up code structure depending on larger dimension to facilitate use of only one loop
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {

					// Fill in array data for first array with loop that goes to "rows" (first array row dimension)
					concatted_arr[i][j] = arr_data[i][j];

					// Condition since row_ < rows, do not exceed array dimensions. 
					if (i < rows_)
						concatted_arr[i + rows][j] = arr_data2[i][j];
				}
			}
		}
		// Second array has biggest row dimension (arr_data2)
		else if (rows_ > rows) {

			// Set up code structure depending on larger dimension to facilitate use of only one loop
			for (int i = 0; i < rows_; i++) {
				for (int j = 0; j < cols; j++) {

					// Condition since rows < rows_, do not exceed array dimensions. 
					if (i < rows)
						concatted_arr[i][j] = arr_data[i][j];

					// Fill in array data for first array with loop that goes to "rows_" (second array row dimension)
					concatted_arr[i + rows][j] = arr_data2[i][j];
				}
			}
		}
		// If first and second array dimensions are the same 
		else if (rows == rows_) {

			// Doesnt matter which dimension you use for row loop since same size
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					concatted_arr[i][j] = arr_data[i][j];
					concatted_arr[i + rows][j] = arr_data2[i][j];
				}
			}
		}
	}

	// Column-wise concatenation (1D array[cols] +  1D array[cols_ -> 1D[cols + cols_] array)
	template <typename T, int cols, int cols_> void concatenate(T(&arr_data)[cols], T(&arr_data2)[cols_], T(&concatted_arr)[cols + cols_]) {

		// Establish which array has the largest collumn dimension. 
		// First array has biggest collumn dimension (arr_data)
		if (cols > cols_) {

			// Set up code structure depending on larger dimension to facilitate use of only one loop
			for (int j = 0; j < cols; j++) {

				// Fill in array data for first array with loop that goes to "cols" (first array row dimension)
				concatted_arr[j] = arr_data[j];

				// Condition since col_ < cols, do not exceed array dimensions. 
				if (j < cols_)
					concatted_arr[j + cols] = arr_data2[j];
			}
		}
		// Second array has biggest collumn dimension (arr_data2) -> cols_
		else if (cols_ > cols) {
			for (int j = 0; j < cols_; j++) {

				if (j < cols)
					concatted_arr[j] = arr_data[j];

				concatted_arr[j + cols] = arr_data2[j];
			}
		}
		// If first and second array have same length of collumns 
		else if (cols == cols_) {

			// Doesnt matter which col dimension used for loop
			for (int j = 0; j < cols; j++) {
				concatted_arr[j] = arr_data[j];
				concatted_arr[j + cols] = arr_data2[j];
			}
		}
	}

	// Collumn-wise concatenation 2D arrays (2D array [rows][cols] + 2D array[rows][cols_] -> 2D array[rows][cols + cols_] 
	template <typename T, int rows, int cols, int cols_> void concatenate(T(&arr_data)[rows][cols], T(&arr_data2)[rows][cols_], T(&concatted_arr)[rows][cols + cols_]) {

		// Establish which array has the largest collumn dimension. 
		// First array has biggest collumn dimension (arr_data)
		if (cols > cols_) {

			// Set up code structure depending on larger dimension to facilitate use of only one loop
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {

					// Fill in array data for first array with loop that goes to "rows" (first array row dimension)
					concatted_arr[i][j] = arr_data[i][j];

					// Condition since row_ < rows, do not exceed array dimensions. 
					if (j < cols_)
						concatted_arr[i][j + cols] = arr_data2[i][j];
				}
			}
		}
		// Second array has biggest cols_ dimension (arr_data2)
		else if (cols_ > cols) {

			// Set up code structure depending on larger dimension to facilitate use of only one loop
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols_; j++) {

					// Condition since cols < cols_, do not exceed array dimensions. 
					if (j < cols)
						concatted_arr[i][j] = arr_data[i][j];

					// Fill in array data for first array with loop that goes to "cols_" (second array row dimension)
					concatted_arr[i][j + cols] = arr_data2[i][j];
				}
			}
		}
		// If first and second array dimensions are the same 
		else if (cols == cols_) {

			// Doesnt matter which dimension you use for collumn loop since same size
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					concatted_arr[i][j] = arr_data[i][j];
					concatted_arr[i][j + cols] = arr_data2[i][j];
				}
			}
		}
	}

	/* slice() -> Extracts data from an array within given bounds */
	template<typename T, int start_ind, int end_ind, int cols> void slice(T(&arr_data)[cols], T(&out_arr)[end_ind - start_ind + 1]) {

		// Guard for exceeding array dimensions
		if ((start_ind < 0) || (end_ind >= cols)) {
			std::cout << "Error: slice() exceeding array dimensions.\n";
			exit(-1);
		}

		for (int j = start_ind; j <= end_ind; j++) {
			out_arr[j - start_ind] = arr_data[j];
		}
	}

	/* ---------------------------- Transformations ----------------------------------------*/
	/* Elementwise vector multiplication*/
	// Create new vector from function 
	template<typename T, int cols> void elementwise(T (&vectorName)[cols], T (&vectorName2)[cols], T (&output)[cols]) {

		for (int i = 0; i < cols; i++)
		{
			output[i] = vectorName[i] * vectorName2[i];

			//if (i == 0)
			//	std::cout << "Result: [" << output[i] << ", ";
			//else if (i < sizeV1 - 1)
			//	std::cout << output[i] << ", ";
			//else
			//	std::cout << output[i] << "]\n";
		}
	}
	// Keep existing vector where first inputted vector is overwritten
	template<typename T, int cols> void elementwise(T(&vectorName_main)[cols], T(&vectorName_multiplier)[cols]) {

		for (int i = 0; i < cols; i++)
		{
			vectorName_main[i] = vectorName_main[i] * vectorName_multiplier[i];

			//if (i == 0)
			//	std::cout << "Result: [" << output[i] << ", ";
			//else if (i < sizeV1 - 1)
			//	std::cout << output[i] << ", ";
			//else
			//	std::cout << output[i] << "]\n";
		}
	}
	
	// Macro for elementwise (.*)

	/* Cross product*/

	/* Determinant */

	/* Eigenvalue calculation */



	/* ---------------------------- Geometry tools -----------------------------------------*/
	/*Generate points of a line*/
	template<int cols> void genpt_straightline(float(&x_data)[cols], float (&y_data)[cols], float grad, float intercept) {

		// Loop over x_data and calculate y points 
		for (int j = 0; j < cols; j++) {
			y_data[j] = grad * x_data[j] + intercept;
		}
	}
	
}

