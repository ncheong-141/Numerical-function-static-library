#pragma once
/* Namespace for all runtime functions */

// External modules
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>


/* IDEA!
 Make guards active and not-active on debug/release mode using a macro like the error log system.
i.e., you run your program in debug mode, check if everything is all good, then deactive the guards on release mode
with the same code.
Note: this is a garbage idea if your data size is time dependent*/

/* Need a plotter*/

/* --------------------------------------- Macros ---------------------------------------------------*/
/* Log system macro for debugging*/

/* Active/In-active guards*/


/*  ---------------------------------Start of RMF API -----------------------------------------------*/
/*  Note:- When using std::vector loops must be indexed with size_t as vector.size(),.capacity() is size_t datatype
		 - When reading an element of a vector, [], is faster than at(i) except there is no bound checks
		   which shouldnt really be needed if you use size() of the vectors.
		   inline - this function will be defined in multiple translation units, don't worry about it. The linker needs to make sure all translation units use a single instance of the variable/function.

	Note:- inline - this function will be defined in multiple translation units, don't worry about it.
		   The linker needs to make sure all translation units use a single instance of the variable/function.
		  -Generally, declaring templates inline is pointless, as they have the linkage semantics of inline
		   already.
		  -However, explicit specialization and instantiation of templates require inline to be used.
	Note:- The standard has a specific feature to improve the efficiency of returning by value.
		   It's called "copy elision", and more specifically in this case the
		   "named return value optimization (NRVO). The vector will be moved not copied*/


namespace RMF {

	/*------------------ Dynamic array generators ------------------------------------------*/

	// 2D array
	template <typename T> void generate_dyn_2D(std::vector< std::vector<T> >& vectorName, const int rows, const int cols) {

		// Reserve memory/pre-allocate to avoid unnesseccary copying and deleting if the vector 
		// is not the correct length initially. 
		vectorName.reserve(rows);

		for (int i = 0; i < rows; i++) {

			// Construct vectors within element to further avoid copying (use emplace_back).
			// i.e. avoid constructing vector on function stack frame and copying it into the vector. 
			// emplace_back is constructing a vector with parameters (cols) where cols is size. 
			vectorName.emplace_back(cols);
		}
	}

	// 3D array. Index similar to static arrays array[depth][rows][cols]
	template <typename T> void generate_dyn_3D(std::vector< std::vector< std::vector<T> > >& vectorName, const int depth, const int rows, const int cols) {

		// Reserve memory/pre-allocate to avoid unnesseccary copying and deleting if the vector 
		// is not the correct length initially. 
		vectorName.reserve(depth);

		// 3D array is an array of 2D arrays (or planes)
		std::vector < std::vector<int> > plane;
		generate_dyn_2D(plane, rows, cols);

		for (int k = 0; k < depth; k++) {

			// Construct vectors within element to further avoid copying (use emplace_back).
			// i.e. avoid constructing vector on function stack frame and copying it into the vector. 
			// emplace_back is constructing the plane into each element of the vector (thus making it 3D) 
			vectorName.emplace_back(plane);
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

	/* ---- 2D vectors memory optimizing functions ---- */
	/* Note, this is for std:;vector< std::vector<T> > variants of 2D arrays.
			 For DYN_C2D arrays, the mem corrector is a member function */

			 // mem_reserver 
	template <typename T> void mem_reserver(std::vector< std::vector<T> >& dyn_arr_data, size_t rows, size_t cols) {

		dyn_arr_data.reserve(rows);				// Reserve rows	
		for (size_t j = 0; j < cols; j++) {
			dyn_arr_data[j].reserve(cols);		// Reserve collumns 
		}
	}

	// mem_corrector 
	template <typename T> void mem_corrector(std::vector < std::vector<T> >& dyn_arr_data, size_t rows, size_t cols) {

		if (dyn_arr_data.size() > rows) {	// If array is bigger rowwise 
			dyn_arr_data.resize(rows);		// Re-size vector to required size 
			dyn_arr_data.shrink_to_fit();	// Deallocate any unused memory 
		}
		else if (dyn_arr_data.size() < rows) {	// If array rows is smaller 
			dyn_arr_data.reserve(rows);
			dyn_arr_data.resize(rows);			// Re-size vector to required size 
		}
		else {									// If array is same size rowwise

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

	/* ---------------- Contiguous multi-dimensional arrays with std::vector---------------*/
	/* In a 2D array, the rows are not stored contiguously as they are they're own std::vector,
	   this leads to performance decreases so a class containing a single std::vector, where
	   the indexes are treated as a 2D array */

	template <typename T> class DYN_C2D {

	private:
		size_t ref_rowsize;					// Size of rows (upper limit)
		size_t ref_colsize;					// Size of cols (upper limit)
		std::vector<T> contiguous_data;		// The data for the array 

	public:

		/* Default object constructor*/
		DYN_C2D() : ref_rowsize(0), ref_colsize(0) {}

		// Constructor for objects (e.g. class) with parameter instansiaton 
		DYN_C2D(size_t cols) : ref_rowsize(0), ref_colsize(cols) {
		}

		/* 2D array Constructor */
		DYN_C2D(size_t rows, size_t cols) : ref_rowsize(rows), ref_colsize(cols) {

			// Allocate 1D contiguous memory for vector
			contiguous_data.reserve(rows * cols);
			contiguous_data.resize(rows * cols);
		}

		/* Destructor */
		~DYN_C2D() {
			//std::cout << "2D contig array destroyed\n";
			contiguous_data.~vector();
		}

		/* ------------------------- Operator overloads --------------------------------------------- */

		/* Accessing contiguous data through 2D indices(reading and writing) */
		T& operator()(size_t i, size_t j) {

			/* Bounds check */
			if (i >= ref_rowsize || j >= ref_colsize) {
				std::cout << "Error: Exceeding 2D array dimensions\n" << "Value of i: " << i << "\nValue of j: " << j << "\n";
				std::cin.get(); exit(-1);
			}

			return contiguous_data[(i * ref_colsize) + j];
		}

		// Const variant for working with const data types
		const T& operator()(size_t i, size_t j) const {

			/* Bounds check */
			if (i >= ref_rowsize || j >= ref_colsize) {
				std::cout << "Error: Exceeding 2D array dimensions\n"; exit(-1);
			}
			return contiguous_data[(i * ref_colsize) + j];
		}

		/* Accessing contiguous data through 1D indices (i.e. if you loop rows*cols) (reading and writing) */
		T& operator()(size_t ij) {

			/* Bounds check */
			if (ij > ref_rowsize* ref_colsize) {
				std::cout << "Error: Exceeding 2D array dimensions\n" << "Value of ij: " << ij << "\nValue of rows*cols: " << ref_rowsize * ref_colsize << "\n";
				std::cin.get(); exit(-1);
			}

			return contiguous_data[ij];
		}

		// Const variant for working with const data types
		const T& operator()(size_t ij) const {

			/* Bounds check */
			if (ij > ref_rowsize* ref_colsize) {
				std::cout << "Error: Exceeding 2D array dimensions\n" << "Value of ij: " << ij << "\nValue of rows*cols: " << ref_rowsize * ref_colsize << "\n";
				std::cin.get(); exit(-1);
			}

			return contiguous_data[ij];
		}

		/* ----------------------- Vector functionality interface ------------------------------- */
		// This allows the std::vector functions to be accessed from the class

		// Access first and end element
		T front() { return contiguous_data.front(); }
		T back() { return contiguous_data.back(); }

		/* data() -> Pointer to the underlying element storage. For non-empty containers, the returned
					 pointer compares equal to the address of the first element.*/
		T* data() { return contiguous_data.data(); }
		const T* data() const { return contiguous_data.data(); }

		/* begin(), end, rbegin, rend and const variants (bloody weird datatypes)*/

		// Defining these weird types as an easier and more understandable semantic 
		using _iterator = std::_Vector_iterator<std::_Vector_val<std::_Simple_types<T> > >;
		using _c_iterator = std::_Vector_const_iterator<std::_Vector_val<std::_Simple_types<T> > >;
		using _r_iterator = std::reverse_iterator<std::_Vector_iterator<std::_Vector_val<std::_Simple_types<T> > > >;
		using _cr_iterator = std::reverse_iterator<std::_Vector_const_iterator<std::_Vector_val<std::_Simple_types<T> > > >;

		// All returns an iterator to the start or end (with const variants). (Efficient way to loop)
		_iterator			begin() { return contiguous_data.begin(); }
		_c_iterator			cbegin() { return contiguous_data.cbegin(); }
		_iterator			end() { return contiguous_data.end(); }
		_c_iterator			cend() { return contiguous_data.cend(); }

		_r_iterator			rbegin() { return contiguous_data.rbegin(); }
		_cr_iterator		crbegin() { return contiguous_data.crbegin(); }
		_r_iterator			rend() { return contiguous_data.rend(); }
		_cr_iterator		crend() { return contiguous_data.crend(); }

		/* empty() -> checks whether container is empty */
		bool empty() { return contiguous_data.empty(); }

		/* size() -> returns the total size of the vector */
		size_t size() { return contiguous_data.size(); }
		size_t size_rows() const { return ref_rowsize; }
		size_t size_cols() const { return ref_colsize; }

		/* reserve -> reserve/pre-allocate storage for vector */
		void reserve(size_t rows, size_t cols) { contiguous_data.reserve(rows * cols); }
		void reserve(size_t reserve_val) { contiguous_data.reserve(reserve_val); }

		/* max_size() -> Returns the maximum number of elements the container is able to hold due
						 to system or library implementation limitations. */
		size_t max_size() { return contiguous_data.max_size(); }

		/* Capacity() -> returns number of elements that can be held in currently allocated storage */
		size_t capacity() { return contiguous_data.capacity(); }

		/* shrink_to_fit -> reduce memory usage by freeing unused memory*/
		void shrink_to_fit() { contiguous_data.shrink_to_fit(); }

		/* clear() -> Erases all elements from the container. After this call, size() returns zero.
					  Invalidates any references, pointers, or iterators referring to contained elements.
					  Any past-the-end iterators are also invalidated.
					  Note, capacity, i.e. memory, is NOT freed*/
		void clear() { contiguous_data.clear(); }

		/* insert(), emplace(), erase() not being added. */

		/* resize() -> */
		// Need to account for change of rows/cols 
		void resize(size_t rows, size_t cols) {

			// Simply readjust member variables 
			ref_colsize = cols;
			ref_rowsize = rows;

			// Re-size contiguous data and trim any excess memory 
			contiguous_data.resize(rows * cols);
			contiguous_data.shrink_to_fit();

		}

		/* emplace_back -> constructs an element in-place at the end */
		// Note, empalce_back, push_back and pop_back all work on ROWS, not single elements (except for element).
		// Element_emplace_back allows you to input objects with different parameters for each element

		// empalce_back a single row (all parameters the same in the row/or default constructor)) 
		void emplace_back(T object) {

			// Update reference rows and cols 			
			ref_rowsize += 1;		// Increment the row size

			// Note, no reserve as this can hugely detriment performance if this function
			// is called alot of times. 
			for (size_t j = 0; j < ref_colsize; j++) {
				contiguous_data.emplace_back(object);
			}
		}

		// Efficient variant of emplace_back for multiple rows growth.  (all parameters the same/or default constructor)) 
		void emplace_back(T object, size_t rows) {

			// Update reference rows and cols 			
			ref_rowsize += rows;		// Increment the row size

			// Allocate an extra collumn of memory
			/* Notes on allocating a collumn of memory:
				- When you are growing an array (adding rows), the memory is deallocated and re-allocatted
				  in a differnet place. This happens when you use reserve() too.
				- Therefore, manually allocating the memory (e.g. for each row as no. of cols) is bad and
				  performs much worse than not reserving memory.
				  This is because emplace_back in the std::vector pre-allocates extra memory, likely a
				  function of how many times emplace_back has been called back to back.
				  This means that manually reserving memory actually has ALOT more re-allocations and
				  thus the performance loss.
				- Furthermore, as the array grows the time it takes to re-allocate also grows (as more data
				  is required to be moved! Thus making the additional re-allocations quadratically more costly*/
				  /* However, for simply adding a small amount of rows, reserving is actually faster as the cost of
					 additional reallocations is negligable compared to the performance gain of allocating the correct
					 memory*/

			contiguous_data.reserve(contiguous_data.capacity() + (ref_colsize * rows));

			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < ref_colsize; j++) {
					contiguous_data.emplace_back(object);
				}
			}
		}

		/* NOTE NEED TO ENSURE THIS DOES NOT MAKE INCOMPLETE ROWS*/
		// emplace_back a single element, allows for different parameters in each element but requires index to check for new row. 
		void element_emplace_back(T object, size_t index) {

			// Update reference rows and cols 	
			if (index == 0)
				ref_rowsize += 1;		// Increment the row size

			if (index >= ref_colsize) {
				std::cout << "Error, exceeding DYN_C2D column dimensions.\n";
				exit(-1);
			}

			// Note, no reserve as this can hugely detriment performance if this function
			// is called alot of times. Reserve should be called outside this function.
			contiguous_data.emplace_back(object);
		}

		// Integral check of data (function should always be called after element emplace back to ensure rows are full/not zero too. 
		void DYN_C2D_complete_rows_integral_check() {

			// If rowsize is zero, first row was not completer 
			if (ref_rowsize == 0) {
				std::cout << "Error, zero row size in DYN_C2D data structure.\n";
				exit(-1);
			}
			if (this->size() != ref_rowsize * ref_colsize) {
				std::cout << "Error, row*col does not equate to DYN_C2D size. Row*col = " << ref_colsize * ref_rowsize << " and size of DYN_C2D is: " << this->size();
				std::cout << ".\nMust be missing elements in last row or constant amount in every row.\n";
				exit(-1);
			}


		}

		/* push_back -> adds an element to the end*/
		// This variant of push_back works for all objects too. 
		void push_back(T object) {
			// Update reference rows and cols 			
			ref_rowsize += 1;		// Increment the row size

			// Note, no reserve as this can hugely detriment performance if this function
			// is called alot of times. 
			for (size_t j = 0; j < ref_colsize; j++) {
				contiguous_data.push_back(object);
			}
		}

		// Efficient variant of push_back for multiple rows growth. 
		void push_back(T object, size_t rows) {

			// Update reference rows and cols 			
			ref_rowsize += rows;		// Increment the row size

			contiguous_data.reserve(contiguous_data.capacity() + (ref_colsize * rows));

			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < ref_colsize; j++) {
					contiguous_data.emplace_back(object);
				}
			}
		}

		/* pop_back-> removes last element*/
		// Note, this works in ROWS 
		void pop_back() {

			// Loop over one collumn
			for (size_t j = 0; j < ref_colsize; j++) {
				contiguous_data.pop_back();
			}

			// Update the rowsize
			ref_rowsize -= 1;

			// Trim the excess memory since it is not deallocated automatically. 
			contiguous_data.shrink_to_fit();
		}

		void pop_back(size_t rows) {

			// Loop over rows*cols ( same as loop i rows, loop j cols but faster). 
			for (size_t ij = 0; ij < rows * ref_colsize; ij++) {
				contiguous_data.pop_back();
			}

			// Update the rowsize
			ref_rowsize -= rows;

			// Trim the excess memory since it is not deallocated automatically. 
			contiguous_data.shrink_to_fit();
		}


		/* Additional vector functionality*/

		// Object initialization (as with an object (e.g. class) you cant predefine the array dims 
		// if the object takes parameters+)
		void param_object_instansiation(T object, size_t rows, size_t cols) {

			// Allocate memory for objects to avoid unneccesary copying
			contiguous_data.reserve(rows * cols);

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {

					// Construct object directly in element to avoid a copy/move command (e.g. use emplace_back) 
					contiguous_data.emplace_back(object);
				}
			}

			// Update referece dimension variables
			ref_rowsize = rows;
			ref_colsize = cols;
		}

		// mem_corrector().
		// This function corrects the memory/size of the array and adjusts the row, col parameters 
		// accordingly. 
		void mem_corrector(size_t rows, size_t cols) {

			// Maximum contiguous index inputed
			size_t max_contig_index = rows * cols;

			// Compare max contiguous index with array
			// resize member function for DYN_C2D already trims excess data and adjusts the 
			// member variables ref_colsize and ref_rowsize. 
			if (contiguous_data.size() > max_contig_index) {		// If array size is larger than intended
				this->resize(rows, cols);
			}
			else if (contiguous_data.size() < max_contig_index) {	// If array size is smaller.. 
				contiguous_data.reserve(max_contig_index);
				this->resize(rows, cols);
			}
			else {	// If array is same size, check capcity for excess memory 
				if (contiguous_data.capacity() != max_contig_index) {
					contiguous_data.shrink_to_fit();
				}
			}
		}

		void setRefcolsize(size_t set_val) {
			ref_colsize = set_val;
		}

		void setRefRowsize(size_t set_val) {
			ref_rowsize = set_val;
		}
	};


	/* ----------------------- General math tools -----------------------------------------*/

	/* any() -> if any element is non-zero return true*/
	template <typename T> bool any(std::vector<T>& dyn_arr_data);
	template <typename T> bool any(std::vector< std::vector<T>>& dyn_arr_data);

	/* all() -> if all elements are non-zero, return true */
	template <typename T> bool all(std::vector<T>& dyn_arr_data);
	template <typename T> bool all(std::vector< std::vector<T>>& dyn_arr_data);

	/*isEven, isOdd*/
	bool isEven(const int& parameter);
	bool isEven(const size_t& parameter);
	bool isOdd(const int& parameter);
	bool isOdd(const size_t& parameter);

	/* abs() -> Returns an absolute value*/
	template<typename T> T abs(T parameter);

	/* arange() fills array with incrementing values (+- 1) from specified start to end*/
	// Generates values between start and end according to discrete parameter (Note, cols must have an extra element for end point) 
	void				arange(std::vector<int>& dyn_arr_data, const int start, const int end);
	void				arange(std::vector<float>& dyn_arr_data, const float start, const float end, const size_t NP);
	void				arange(std::vector<double>& dyn_arr_data, const double start, const double end, const size_t NP);
	std::vector<int>	arange(const int start, const int end);
	std::vector<float>	arange(const float start, const float end, const size_t NP);
	std::vector<double> arange(const double start, const double end, const size_t NP);

	/* tile() fills array with a specified value*/
	template <typename T> void tile(T tilevalue, std::vector<T>& dyn_arr_data);
	template <typename T> void tile(T tilevalue, std::vector< std::vector<T> >& dyn_arr_data);
	template <typename T> void tile(T tilevalue, DYN_C2D<T>& dyn_arr_data);

	/* std-style elemen-wise vector multiplcation*/
	template <typename T> void elementwise(std::vector<T>& output, std::vector<T>& vectorName, std::vector<T>& vectorName2);
	template <typename T> void elementwise(std::vector<T>& output, std::vector<T>& vectorName, std::vector<T>& vectorName2);
	template <typename T> void elementwise(DYN_C2D<T>& output, DYN_C2D<T>& arrayName, DYN_C2D<T>& arrayName2);

	/* max() find the maximum value of an array or between two values*/
	template <typename T> T max(std::vector<T>& dyn_arr_data);
	template <typename T> T max(DYN_C2D<T>& dyn_arr_data);

	/* min()*/
	template <typename T> T min(std::vector<T>& dyn_arr_data);

	/* find() find an element in an array which equals a value and return an logical array */
	template <typename T> void find(T value, std::vector<T>& dyn_arr_data, std::vector<bool>& logic_arr);
	template <typename T> void find(T value, DYN_C2D<T>& dyn_arr_data, DYN_C2D<int>& logic_arr);

	/* sum () sums all element of an array (1 and 2D)*/
	template <typename T> T		sum(std::vector<T>& dyn_arr_data);
	template <typename T> void	sum(DYN_C2D<T>& dyn_arr_data, std::vector<T>& sum_result, bool sum_row_or_col);

	/* mean() */
	template <typename T> double mean(std::vector<T>& dyn_arr_data);

	/* eye() -> generate an identity matrix */
	template <typename T> DYN_C2D<T> eye(size_t rows, size_t cols);
	template <typename T> void		 eye(DYN_C2D<T>& dyn_arr_output);

	/* diag() -> generate a diagonal array from a vector/ or extract a vector from a diagonal array */
	template <typename T> void diag(DYN_C2D<T>& dyn_arr_data, T diag_val);
	template <typename T> void diag(DYN_C2D<T>& dyn_arr_data, std::vector<T>& vectorInput);
	template <typename T> void diag_extract(std::vector<T>& vectorOutput, DYN_C2D<T>& dyn_arr_data);
	template <typename T> std::vector<T> diag_extract(DYN_C2D<T>& dyn_arr_data);

	/* flatten() -> flatten a 2D array (rows cols) onto a 1D array (cols) */
	template <typename T> void flatten(std::vector<T>& flattened_arr, DYN_C2D<T>& dyn_arr_data);
	template <typename T> std::vector<T> flatten(DYN_C2D<T>& dyn_arr_data);

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

	/* slice() -> Extracts data from an array within given bounds */




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
	template void RMF::tile<double>(double, std::vector< std::vector<double> >&);

	template void RMF::tile<int>(int, DYN_C2D<int>&);
	template void RMF::tile<float>(float, DYN_C2D<float>&);
	template void RMF::tile<double>(double, DYN_C2D<double>&);

	/* elementwise multiplication */
	template void RMF::elementwise(std::vector<int>&, std::vector<int>&, std::vector<int>&);
	template void RMF::elementwise(std::vector<float>&, std::vector<float>&, std::vector<float>&);
	template void RMF::elementwise(std::vector<double>&, std::vector<double>&, std::vector<double>&);

	template void elementwise<int>(DYN_C2D<int>&, DYN_C2D<int>&, DYN_C2D<int>&);
	template void elementwise<float>(DYN_C2D<float>&, DYN_C2D<float>&, DYN_C2D<float>&);
	template void elementwise<double>(DYN_C2D<double>&, DYN_C2D<double>&, DYN_C2D<double>&);

	/* max() */
	template int	RMF::max(std::vector<int>&);
	template float	RMF::max(std::vector<float>&);
	template double RMF::max(std::vector<double>&);

	template int	RMF::max(DYN_C2D<int>&);
	template float	RMF::max(DYN_C2D<float>&);
	template double RMF::max(DYN_C2D<double>&);

	/* min() */
	template int	RMF::min(std::vector<int>&);
	template float	RMF::min(std::vector<float>&);
	template double RMF::min(std::vector<double>&);

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

	/* mean() */
	template double RMF::mean(std::vector<int>&);
	template double RMF::mean(std::vector<float>&);
	template double RMF::mean(std::vector<double>&);

	/* eye() */
	template DYN_C2D<int> RMF::eye(size_t, size_t);
	template DYN_C2D<float> RMF::eye(size_t, size_t);
	template DYN_C2D<double> RMF::eye(size_t, size_t);

	template void RMF::eye(DYN_C2D<int>&);
	template void RMF::eye(DYN_C2D<float>&);
	template void RMF::eye(DYN_C2D<double>&);

	/* diag() */
	template void RMF::diag(DYN_C2D<int>&, int);
	template void RMF::diag(DYN_C2D<float>&, float);
	template void RMF::diag(DYN_C2D<double>&, double);

	template void RMF::diag(DYN_C2D<int>&, std::vector<int>&);
	template void RMF::diag(DYN_C2D<float>&, std::vector<float>&);
	template void RMF::diag(DYN_C2D<double>&, std::vector<double>&);

	/* diag_extract() */
	template void RMF::diag_extract(std::vector<int>&, DYN_C2D<int>&);
	template void RMF::diag_extract(std::vector<float>&, DYN_C2D<float>&);
	template void RMF::diag_extract(std::vector<double>&, DYN_C2D<double>&);

	template std::vector<int> diag_extract(DYN_C2D<int>&);
	template std::vector<float> diag_extract(DYN_C2D<float>&);
	template std::vector<double> diag_extract(DYN_C2D<double>&);

	/* flatten*/
	template void RMF::flatten(std::vector<int>&, DYN_C2D<int>&);
	template void RMF::flatten(std::vector<float>&, DYN_C2D<float>&);
	template void RMF::flatten(std::vector<double>&, DYN_C2D<double>&);
	template std::vector<int> RMF::flatten(DYN_C2D<int>&);
	template std::vector<float> RMF::flatten(DYN_C2D<float>&);
	template std::vector<double> RMF::flatten(DYN_C2D<double>&);

	/* reshape() */
	template void RMF::reshape(DYN_C2D<int>&);
	template void RMF::reshape(DYN_C2D<float>&);
	template void RMF::reshape(DYN_C2D<double>&);

	/* fliplr() */
	template void RMF::fliplr(std::vector<int>&);
	template void RMF::fliplr(std::vector<float>&);
	template void RMF::fliplr(std::vector<double>&);

}


