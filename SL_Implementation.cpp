/* Implementation file for general func static library header*/

// Get interface files to link to console 
#include "SL_Runtime_Interface.h"
#include "SL_Comptime_Interface.h"

namespace RMF {

	/*------------------ RMF Guards of the Static Library ---------------------------------------------*/
	// Note: - These are only used in the functions in this namespace. 
	
	/* ---- 1D array guards ----*/
	// Prevent a zero size vector being used (as it is pointless or wrong) 
	template <typename T> void colsize_guard(std::vector<T>& dyn_arr_data, const char* fnName) {
		if (dyn_arr_data.size() == 0) {
			std::cout << "Error: Dynamic vector has zero size in function, " << fnName << ".\n";
			exit(-1);
		}
	}

	/* ---- 2D array guards ---- */
	// - Prevent a zero size vector being used (as it is pointless or wrong) 
	// - Since vector can have vectors of differnet lengths inside, check for same length. (Memory violation prevention) 
	template <typename T> void colsize_guard(std::vector< std::vector<T>>& dyn_arr_data, const char* fnName) {

		size_t ref_rowsize = dyn_arr_data.size();
		size_t ref_colsize = dyn_arr_data[0].size();	// Calculate ref size only once (instead of in loop)

		for (size_t i = 0; i < ref_rowsize; i++) {
			if (dyn_arr_data[i].size() != ref_colsize) {
				std::cout << "Error: Dynamic vector has different collumn dimensions in function," << fnName << ".\n";
				exit(-1);
			}
		}
	}
	template <typename T> void rowsize_guard(std::vector< std::vector<T>>& dyn_arr_data, const char* fnName) {

		if (dyn_arr_data.size() == 0) {
			std::cout << "Error: Dynamic vector has zero row dimension in function," << fnName << ".\n";
			exit(-1);
		}
	}


	/*------------------ RMF Vector Engineers of the Static Library ------------------------------------*/
	// (Optomizers) 

	/* ---- 1D vectors ---- */

	// mem_reserver()
	// - If a vector is implemented with no size,  this function allocates correct memory
	// - It also checks the capacity for excess memory
	template <typename T> void mem_reserver(std::vector<T>& dyn_arr_data, size_t ref_colsize) {

		if (dyn_arr_data.capacity() != ref_colsize) {		// If incorrect memory/no memory allocated
			dyn_arr_data.reserve(ref_colsize);				// Allocate memory to avoid unneccessary copying
		}
	}
	
	// mem_corrector()  
	// - If a vector is inputted with incorrect size, this function corrects 
	//   the size and allocates correct memory
	// - It also checks the capacity for excess memory if the vector is the correct size
	template <typename T> void mem_corrector(std::vector<T>& dyn_arr_data, size_t ref_colsize) {

		if (dyn_arr_data.size() > ref_colsize) {		// If vector is bigger than needed size 
			dyn_arr_data.resize(ref_colsize);			// Re-size vector to required size
			dyn_arr_data.shrink_to_fit();				// Deallocate any un-used memory
		}
		else if (dyn_arr_data.size() < ref_colsize) {	// If vector is smaller than needed size
			dyn_arr_data.reserve(ref_colsize);			// Allocate correct memory (note, this does not add onto memory, end-start is the value of allocation)
														// Furthermore, with testing, if you only re-size it allocates extra memory!!
			dyn_arr_data.resize(ref_colsize);			// Re-size vector to required size
		}
		else {												// If vector is same size
			if (dyn_arr_data.capacity() != ref_colsize) {	// If incorrect memory allocated
				dyn_arr_data.shrink_to_fit();				// Deallocate any un-used memory
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
		colsize_guard<T>(dyn_arr_data, fnName);

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
		colsize_guard<T>(dyn_arr_data, fnName);

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
	template <typename T> void tile(T tilevalue, const size_t cols, std::vector<T>& dyn_arr_data) {
		
		// Optimize (and correct if neccessary) the data structure for function
		mem_corrector(dyn_arr_data, cols);
		
		for (int i = 0; i < cols; i++) {
			dyn_arr_data[i] = tilevalue;
		}

	}
	template <typename T> void tile(T tilevalue, const size_t rows, const size_t cols, std::vector< std::vector<T> >& dyn_arr_data) {
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				dyn_arr_data[i][j] = tilevalue;
			}
		}
	}

	template <typename T> void elementwise(std::vector<T>& output, std::vector<T>& vectorName, std::vector<T>& vectorName2)
	{

		// Guard to stop vectors/1D arrays of different lengths being multiplied together 
		if (vectorName.size() != vectorName2.size())
		{
			std::cout << "Error: You cannot elementwise multiply two vectors/1D arrays of different lengths" << std::endl;
			exit(-1);
		}
		else
		{
			// Instantiate reserve space for vector for optimizing performance (reduces copying) 
			mem_corrector(output, vectorName.size());

			// Loop along elements of vector/ 1D arrays 

			for (int i = 0; i < vectorName.size(); i++)
			{
				output[i] = (vectorName[i] * vectorName2[i]);

				if (i == 0)
					std::cout << "Result: [" << output[i] << ", ";
				else if (i < vectorName.size() - 1)
					std::cout << output[i] << ", ";
				else
					std::cout << output[i] << "]" << std::endl;
			}
		}
	}	

	/* -------------------- Explicit instanciations ------------------------*/

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
	template void RMF::tile<int>(int, const size_t, std::vector<int>&);
	template void RMF::tile<float>(float, const size_t, std::vector<float>&);
	template void RMF::tile<double>(double, const size_t, std::vector<double>&);
	template void RMF::tile<int>(int, const size_t, const size_t, std::vector< std::vector<int> >&);
	template void RMF::tile<float>(float, const size_t, const size_t, std::vector< std::vector<float> >&);
	template void RMF::tile<double>(double, const size_t, const size_t, std::vector< std::vector<double> >&);

	/* elementwise multiplication */
	template void RMF::elementwise(std::vector<int>&, std::vector<int>&, std::vector<int>&);
	template void RMF::elementwise(std::vector<float>&, std::vector<float>&, std::vector<float>&);
	template void RMF::elementwise(std::vector<double>&, std::vector<double>&, std::vector<double>&);

} // End of namespace RUNTIME

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




	// Explicit template specification
	// Forward declarations which are needed for the linker to link the template function prototype
	// in header with that in implementation! (Explicit template instantiation)						

	// Required to split the implementation of the static library with the public interface (API)
	// Use of extern keyword to tell the linker to find the function declarations in another .obj file 
	// (created by the .cpp file)


	
}
