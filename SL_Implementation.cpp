/* Implementation file for general func static library header*/

// Get interface files to link to console 
#include "SL_Runtime_Interface.h"
#include "SL_Comptime_Interface.h"

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
		if (( arrayName.size_rows() != arrayName2.size_rows()) || (arrayName.size_cols() != arrayName2.size_cols())) {
			std::cout << "Error: You cannot elementwise multiply two arrays of different row or collumn lengths in function," << fnName << ".\n";
			exit(-1);
		}
	}


	/*------------------ RMF Vector Engineers of the Static Library ------------------------------------*/
	// Note: - These are only used in the functions in this namespace. 

	/* ---- 1D vectors memory optimizing functions ---- */

	// mem_reserver()
	// - If a vector is implemented with no size,  this function allocates correct memory
	// - It also checks the capacity for excess memory
	template <typename T> void mem_reserver(std::vector<T>& dyn_arr_data, size_t ref_colsize) {

		if (dyn_arr_data.capacity() != ref_colsize) {		// If incorrect memory/no memory allocated
			dyn_arr_data.reserve(ref_colsize);		// Allocate memory to avoid unneccessary copying
		}
	}
	
	// mem_corrector()  
	// - If a vector is inputted with incorrect size, this function corrects 
	//   the size and allocates correct memory
	// - It also checks the capacity for excess memory if the vector is the correct size
	template <typename T> void mem_corrector(std::vector<T>& dyn_arr_data, size_t ref_colsize) {

		if (dyn_arr_data.size() > ref_colsize) {		// If vector is bigger than needed size 
			dyn_arr_data.resize(ref_colsize);		// Re-size vector to required size
			dyn_arr_data.shrink_to_fit();			// Deallocate any un-used memory
		}
		else if (dyn_arr_data.size() < ref_colsize) {		// If vector is smaller than needed size
			dyn_arr_data.reserve(ref_colsize);		// Allocate correct memory (note, this does not add onto memory, end-start is the value of allocation)
									// Furthermore, with testing, if you only re-size it allocates extra memory!!
			dyn_arr_data.resize(ref_colsize);		// Re-size vector to required size
		}
		else {							// If vector is same size
			if (dyn_arr_data.capacity() != ref_colsize) {	// If incorrect memory allocated
				dyn_arr_data.shrink_to_fit();		// Deallocate any un-used memory
			}
		}
	}

	/* ---- 2D vectors memory optimizing functions ---- */
	/* Note, this is for std:;vector< std::vector<T> > variants of 2D arrays. 
			 For DYN_C2D arrays, the mem corrector is a member function */

	// mem_reserver 
	template <typename T> void mem_reserver(std::vector< std::vector<T> >& dyn_arr_data, size_t rows, size_t cols) {
		
		dyn_arr_data.reserve(rows);			// Reserve rows	
		for (size_t j = 0; j < cols; j++) {
			dyn_arr_data[j].reserve(cols);		// Reserve collumns 
		}
	}

	// mem_corrector 
	template <typename T> void mem_corrector(std::vector < std::vector<T> >& dyn_arr_data, size_t rows, size_t cols) {

		if (dyn_arr_data.size() > rows) {	// If array is bigger rowwise 
			dyn_arr_data.resize(rows);	// Re-size vector to required size 
			dyn_arr_data.shrink_to_fit();	// Deallocate any unused memory 
		}
		else if (dyn_arr_data.size() < rows) {	// If array rows is smaller 
			dyn_arr_data.reserve(rows);		
			dyn_arr_data.resize(rows);	// Re-size vector to required size 
		}
		else {					// If array is same size rowwise
			
			// Trim excess data 
			if (dyn_arr_data.capacity() != rows) {	
				dyn_arr_data.shrink_to_fit(); 
			}
		}

		// Check collumn dimensions 
		for (size_t j = 0; j < cols; j++) {
			if (dyn_arr_data[j].size() > cols) {
				dyn_arr_data[j].resize(rows);
				dyn_arr_data[j].shrink_to_fit();
			}
			else if (dyn_arr_data[j].size() < cols) {
				dyn_arr_data[j].reserve(cols);
				dyn_arr_data[j].resize(cols);
			}
			else {
				if (dyn_arr_data[j].capacity() != cols) {
					dyn_arr_data[j].reserve(cols);
				}
			}
		}
	}


	/* ------------------ START OF LIBRARY FUNCTIONS: General math tools -------------------------------*/

	/* any() -> if any element is non-zero return true */
	template <typename T> bool any(std::vector<T>& dyn_arr_data){
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
		
		for (size_t i = 0; i < ref_rowsize; i++){
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
	bool isOdd(const int& parameter) {

		// Calculate remainder of number divided by 2 (modulus) (if 0 then it is even)
		int modulus = parameter % 2;

		if (modulus != 0)
			return true;
		else
			return false;
	}

	/* abs() -> Returns an absolute value*/
	template <typename T> T abs(T parameter){

		if (parameter < 0)
			return (-1 * parameter);

		return parameter;
	}

	/* arange() fills array with incrementing values (+- 1) from specified start to end*/
	// Generates values between start and end according to discrete parameter (Note, cols must have an extra element for end point) 
	void arange(std::vector<int>& dyn_arr_data, const int start, const int end) {

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
	void arange(std::vector<float>& dyn_arr_data, const float start, const float end, const size_t NP) {

		// Optimize data structure for function
		mem_corrector(dyn_arr_data, NP);

		float total_length = abs(end - start);		// Calculate total value between start and end
		float disc_inc = total_length / (NP - 1);       // Calculate the increments between each element
		float disc_val = 0;				// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;		// Write element
			disc_val += disc_inc;			// Add increment
		}
	}
	void arange(std::vector<double>& dyn_arr_data, const double start, const double end, const size_t NP) {

		// Optimize data structure for function
		mem_corrector(dyn_arr_data, NP);
	
		double total_length = abs(end - start);		// Calculate total value between start and end
		double disc_inc = total_length / (NP - 1);      // Calculate the increments between each element
		double disc_val = 0;				// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;		// Write element
			disc_val += disc_inc;			// Add increment
		}
	}
	std::vector<int> arange(const int start, const int end) {

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

		float total_length = abs(end - start);		// Calculate total value between start and end
		float disc_inc = total_length / (NP - 1);       // Calculate the increments between each element
		float disc_val = 0;				// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;		// Write element
			disc_val += disc_inc;			// Add increment
		}

		return dyn_arr_data;
	}
	std::vector<double> arange(const double start, const double end, const size_t NP) {

		std::vector<double> dyn_arr_data(NP);

		double total_length = abs(end - start);		// Calculate total value between start and end
		double disc_inc = total_length / (NP - 1);      // Calculate the increments between each element
		double disc_val = 0;				// Initialize discrete value of elements 

		// Loop over collumns and write values 
		for (size_t j = 0; j < NP; j++) {

			dyn_arr_data[j] = disc_val;		// Write element
			disc_val += disc_inc;			// Add increment
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
		for (int i = 0; i < vectorName.size(); i++)	{
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
	template <typename T> T	sum(std::vector<T>& dyn_arr_data) {
		
		// Guard for zero vector 
		size_guard(dyn_arr_data, "sum()"); 

		T sum_val = 0; 
		for (auto data_addr = dyn_arr_data.begin(); data_addr != dyn_arr_data.end(); data_addr++) {
			sum_val += *data_addr;
		}
		return sum_val; 
	}
	template <typename T> void sum(DYN_C2D<T>& dyn_arr_data, std::vector<T>& sum_result, bool sum_row_or_col) {

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
						sum_val += dyn_arr_data(i,j);
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
	/* ------------- Below is to be implemented, however, see compile time implementation for these functions --------------------*/
	/* mean() */
	template <typename T> double mean(std::vector<T>& dyn_arr_data);

	/* eye() -> generate an identity matrix */
	template <typename T> void eye(DYN_C2D<T>& dyn_arr_output, size_t rows, size_t cols);

	/* diag() -> generate a diagonal array from a vector/ or extract a vector from a diagonal array */
	template <typename T> void diag(DYN_C2D<T>& dyn_arr_output, std::vector<T>& vectorInput);
	template <typename T> void diag_extract(std::vector<T>& vectorOutput, DYN_C2D<T>& dyn_arr_input);

	/* flatten() -> flatten a 2D array (rows cols) onto a 1D array (cols) */
	template <typename T> void flatten(std::vector<T>& flattened_arr, DYN_C2D<T>& dyn_arr_input);

	/* reshape() -> basically flips rows and collumn values; i.e. 2x3 array to a 3x2 array. */
	template <typename T> void reshape(DYN_C2D<T>& dyn_arr_data);

	/* l_interp -> linear interpolation between points. Note, integer not allowed as there is division in this function */
	// Linear interpolation of array data 
	void l_interp(DYN_C2D<float>& x_y_data, DYN_C2D<float>& x_y_newdata);
	void l_interp(DYN_C2D<double>& x_y_data, DYN_C2D<double>& x_y_newdata);

	/* flipdata -> flips data left to right, up to down etc - left to right -> settign 1, up to down setting 2 */
	template <typename T> void fliplr(std::vector<T>& dyn_arr_data);
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

	/* -------------------- Explicit template instansiations ------------------------*/
	// Forward declarations which are needed for the linker to link the template function prototype
	// in header with that in implementation! (Explicit template instantiation)						

	/* any() */
	template bool any<int>(std::vector<int>& dyn_arr_data);
	template bool any<int>(std::vector< std::vector<int>>& dyn_arr_data);
	template bool any<float>(std::vector<float>& dyn_arr_data);
	template bool any<float>(std::vector< std::vector<float>>& dyn_arr_data);
	template bool any<double>(std::vector<double>& dyn_arr_data);
	template bool any<double>(std::vector< std::vector<double>>& dyn_arr_data);

	/* all() */
	template bool all<int>(std::vector<int>& dyn_arr_data);
	template bool all<int>(std::vector< std::vector<int>>& dyn_arr_data);
	template bool all<float>(std::vector<float>& dyn_arr_data);
	template bool all<float>(std::vector< std::vector<float>>& dyn_arr_data);
	template bool all<double>(std::vector<double>& dyn_arr_data);
	template bool all<double>(std::vector< std::vector<double>>& dyn_arr_data);

	/* abs() */
	template int RMF::abs<int>(int);
	template size_t RMF::abs<size_t>(size_t);
	template float RMF::abs<float>(float);
	template double RMF::abs<double>(double);

	/* tile()*/
	template void RMF::tile<int>(int, std::vector<int>&);
	template void RMF::tile<float>(float, std::vector<float>&);
	template void RMF::tile<double>(double, std::vector<double>&);
	template void RMF::tile<int>(int, std::vector< std::vector<int> >&);
	template void RMF::tile<float>(float, std::vector< std::vector<float> >&);
	template void RMF::tile<double>(double,  std::vector< std::vector<double> >&);
	template void RMF::tile<int>(int, DYN_C2D<int>&);
	template void RMF::tile<float>(float, DYN_C2D<float>&);
	template void RMF::tile<double>(double, DYN_C2D<double>&);

	/* max() */
	template int RMF::max(std::vector<int>&);
	template float RMF::max(std::vector<float>&);
	template double RMF::max(std::vector<double>&);
	template int RMF::max(DYN_C2D<int>&);
	template float RMF::max(DYN_C2D<float>&);
	template double RMF::max(DYN_C2D<double>&);

	/* find() */
	template void RMF::find(int, std::vector<int>&, std::vector<bool>&);
	template void RMF::find(float, std::vector<float>&, std::vector<bool>&);
	template void RMF::find(double, std::vector<double>&, std::vector<bool>&);
	template void RMF::find(int, DYN_C2D<int>&, DYN_C2D<int>&);
	template void RMF::find(float, DYN_C2D<float>&, DYN_C2D<int>&);
	template void RMF::find(double, DYN_C2D<double>&, DYN_C2D<int>&);

	/* sum() */
	template int    RMF::sum(std::vector<int>&);
	template float  RMF::sum(std::vector<float>&);
	template double RMF::sum(std::vector<double>&);
	template void	RMF::sum(DYN_C2D<int>&, std::vector<int>&, bool);
	template void	RMF::sum(DYN_C2D<float>&, std::vector<float>&, bool);
	template void	RMF::sum(DYN_C2D<double>&, std::vector<double>&, bool);



	/* elementwise multiplication */
	template void RMF::elementwise(std::vector<int>&, std::vector<int>&, std::vector<int>&);
	template void RMF::elementwise(std::vector<float>&, std::vector<float>&, std::vector<float>&);
	template void RMF::elementwise(std::vector<double>&, std::vector<double>&, std::vector<double>&);
	template void elementwise<int>(DYN_C2D<int>&, DYN_C2D<int>&, DYN_C2D<int>&);
	template void elementwise<float>(DYN_C2D<float>&, DYN_C2D<float>&, DYN_C2D<float>&);
	template void elementwise<double>(DYN_C2D<double>&, DYN_C2D<double>&, DYN_C2D<double>&);


} // End of namespace RMF



/* Compile time function variant implementations
Note: most in header file since using non-type template functions */
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
	bool isOdd(const int& parameter) {

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
