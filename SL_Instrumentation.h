#pragma once
/* =================================== Instrumentation interface file =================================== */

// Contains data acquisition, timers and datatypes (except for DYN_C2D)

// Internal modules 
#include "SL_Runtime_interface.h"

// External modules 
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>


/* Error macro */
#define Error_msg(X) std::cout << X << ".\n"; exit(-1);

/* Add
New file for instrumentaation and datatypes
ADD Tuple
update max, min to return index values.
update find to return an index.
update data acquisition to accept more datatypes*/
namespace instrumentation {



	/*=================================== Data acquisition instrumentation ===================================*/
	// Max row establisher 
	size_t Data_sizerow_establish(const char* filename) {

		// Get data from resource files. 
		std::ifstream inputfile(filename);
		std::string line;

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		size_t count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();
			count_i++;

		}
		inputfile.close();

		return count_i;
	}

	// Num cols establisher ( assuming constant collumn size in every row)
	size_t Data_sizecol_establish(const char* filename, const char delimiter) {
		
		
		// Get data from resource files. 
		std::ifstream inputfile(filename);
		std::string line;

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Instantiate string stream with row data from text file. 
		std::getline(inputfile, line);

		// Count number of values based on constant delimiter position (since 2 digit numbers count as 2 values in strings)
		size_t count = 1; 
		for (size_t i = 1; i < line.size(); i++) {
			if (line[i-1] == delimiter) {
				count++; 
			}
		}

		// Close file 
		inputfile.close();

		// Return the size of the row
		return count; 
	}


	// Stack memory 
	/* One dimension - row-wise in data file */
	template <int cols>
	void Data_acquisition(const char* filename, int(&stack_array)[cols], bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;
		std::string element; 

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		bool end_line = false;
		unsigned int count = 0; 
		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);
			std::getline(ss, element);

			// Convert string into integer values
			stack_array[count] = std::stoi(element);
			count++;
		}
		inputfile.close();
	};

	/* One dimension - column wise in data file with constant delimiter */
	template <int cols>
	void Data_acquisition(const char* filename, int(&stack_array)[cols], const char delimiter) {

		// Get data from resource files. 
		std::ifstream inputfile(filename);
		std::string line;

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Instantiate string stream with row data from text file. 
		std::getline(inputfile, line);
		std::stringstream ss(line);
		std::string data_row[cols];

		// Loop over row and obtain collumn data into the string stream. 
		for (unsigned int col = 0; col < cols; col++) {

			if (col != cols - 1)
				std::getline(ss, data_row[col], delimiter);
			else
				std::getline(ss, data_row[col]);

			// Convert string into integer values
			stack_array[col] = std::stoi(data_row[col]);
		}

		// Close file 
		inputfile.close();
	};

	/* Colwise but no delimiter*/
	template <int cols>
	void Data_acquisition(const char* filename, int(&stack_array)[cols]) {

		// Get data from resource files. 
		std::ifstream inputfile(filename);
		std::string line;

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Instantiate string stream with row data from text file. 
		std::getline(inputfile, line);

		// Check if data is the same size
		if (line.size() != cols) {
			std::cout << "Error: data file and structure not the same size. \nData file collumn size: " << line.size() << "\nData structure collumn size: " << cols << "\n";
			exit(-1);
		} 

		// Loop over row and obtain collumn data into the string stream. 
		for (size_t col = 0; col < line.size(); col++) {

			// Convert char into integer values
			stack_array[col] = line[col] - 48; 
		}

		// Close file 
		inputfile.close();
	};

	/* 2 dimensions */ 
	// Multi delimiters : - int overload 
	template <int rows, int cols>
	void Data_acquisition(const char* filename, int(&stack_array)[rows][cols], const char multi_delim[cols - 1], bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;
		std::string data_row[cols];

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		unsigned int count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end... (Coult use a do while loop instead)
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			for (unsigned int col = 0; col < cols; col++) {

				if (col != col - 1)
					std::getline(ss, data_row[col], multi_delim[col]);
				else
					std::getline(ss, data_row[col]);

				// Convert string into integer values
				stack_array[count_i][col] = std::stoi(data_row[col]);
			}

			// Incremenet the row count
			count_i++;

			// Check that data size matches with input stack array (avoid accessing unallocated memory.)
			if (count_i > rows) {
				std::cout << "\nError: input stack_array does not match text file data dimensions. Rows in text file: " << count_i << ".\n" << "Rows in stack_array: " << rows << ".\n";
				exit(-1);
			}

		}
		inputfile.close();
	};

	// Multi-delimiter: - string overload
	template <int rows, int cols>
	void Data_acquisition(const char* filename, std::string(&stack_array)[rows][cols], const char multi_delim[cols - 1], bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		unsigned int count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end... (Coult use a do while loop instead)
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			for (unsigned int col = 0; col < cols; col++) {

				if (col != col - 1)
					std::getline(ss, stack_array[count_i][col], multi_delim[col]);
				else
					std::getline(ss, stack_array[count_i][col]);
			}

			// Incremenet the row count
			count_i++;

			// Check that data size matches with input stack array (avoid accessing unallocated memory.)
			if (count_i > rows) {
				std::cout << "\nError: input stack_array does not match text file data dimensions. Rows in text file: " << count_i << ".\n" << "Rows in stack_array: " << rows << ".\n";
				exit(-1);
			}

		}
		inputfile.close();
	};

	// 2D Constant delimiter 
	template<int rows, int cols>
	void Data_acquisition(const char* filename, int(&stack_array)[rows][cols], const char delimiter, bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;
		std::string data_row[cols];

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		unsigned int count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			for (unsigned int col = 0; col < cols; col++) {

				if (col != cols - 1)
					std::getline(ss, data_row[col], delimiter);
				else
					std::getline(ss, data_row[col]);

				// Convert string into integer values
				stack_array[count_i][col] = std::stoi(data_row[col]);
			}

			// Incremenet the row count
			count_i++;

			// Check that data size matches with input stack array (avoid accessing unallocated memory.)
			if (count_i > rows) {
				std::cout << "\nError: input stack_array does not match text file data dimensions. Rows in text file: " << count_i << ".\n" << "Rows in stack_array: " << rows << ".\n";
				exit(-1);
			}

		}
		inputfile.close();
	};

	template<int rows, int cols>
	void Data_acquisition(const char* filename, std::string(&stack_array)[rows][cols], const char delimiter, bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;
		std::string data_row[cols];

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		unsigned int count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			for (unsigned int col = 0; col < cols; col++) {

				if (col != cols - 1)
					std::getline(ss, data_row[col], delimiter);
				else
					std::getline(ss, data_row[col]);

				// Convert string into integer values
				stack_array[count_i][col] = data_row[col];
			}

			// Incremenet the row count
			count_i++;

			// Check that data size matches with input stack array (avoid accessing unallocated memory.)
			if (count_i > rows) {
				std::cout << "\nError: input stack_array does not match text file data dimensions. Rows in text file: " << count_i << ".\n" << "Rows in stack_array: " << rows << ".\n";
				exit(-1);
			}

		}
		inputfile.close();
	};

	// 2D no delimiter
	template<int rows, int cols>
	void Data_acquisition(const char* filename, std::string(&stack_array)[rows][cols]) {

		// Initialize variables data type for extracting data. 
		std::string line;
		std::string data_row;

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		unsigned int count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);
			std::getline(ss, data_row);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			for (unsigned int col = 0; col < cols; col++) {

				// Convert string into integer values
				stack_array[count_i][col] = data_row[col];
			}

			// Incremenet the row count
			count_i++;

			// Check that data size matches with input stack array (avoid accessing unallocated memory.)
			if (count_i > rows) {
				std::cout << "\nError: input stack_array does not match text file data dimensions. Rows in text file: " << count_i << ".\n" << "Rows in stack_array: " << rows << ".\n";
				exit(-1);
			}

		}
		inputfile.close();
	};

	// ----------------------------- Dynamic memory loader

	// Single string input
	void Data_acquisition(const char* filename, std::string& string, const char delimiter, bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			std::getline(ss, string, delimiter);
		}
		inputfile.close();
	};

	// Multiple strings in file
	void Data_acquisition(const char* filename, std::vector<std::string>& string_vect) {

		// Initialize variables data type for extracting data. 
		std::string line;

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		while (std::getline(inputfile, line)) {
			string_vect.push_back(line);
		}
		inputfile.close();
	};

	// 1D collumn wise for strings 
	void Data_acquisition(const char* filename, std::vector<std::string>& string_vect, const char delimiter) {

		// Get data from resource files. 
		std::ifstream inputfile(filename);
		std::string line;

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Instantiate string stream with row data from text file. 
		std::getline(inputfile, line);
		std::stringstream ss(line);
		
		// Count number of values based on constant delimiter position (since 2 digit numbers count as 2 values in strings)
		size_t count = 1;
		for (size_t i = 1; i < line.size(); i++) {
			if (line[i - 1] == delimiter) {
				count++;
			}
		}

		// Check to ensure string_vect is the same size as the input file. 
		if (string_vect.size() != count) {
			string_vect.resize(count);
			string_vect.shrink_to_fit(); 
		}

		// Loop over row and obtain collumn data into the string stream. 
		for (unsigned int col = 0; col < count; col++) {

			if (col != count - 1)
				std::getline(ss, string_vect[col], delimiter);
			else
				std::getline(ss, string_vect[col]);
		}

		// Close file 
		inputfile.close();
	};

	// Using DYN_C2D 
	void Data_acquisition(const char* filename, RMF::DYN_C2D<int>& dyn_arr, const char delimiter, bool del_first_char) {

		// Initialize variables data type for extracting data. 
		std::string line;
		std::vector<std::string> data_row(dyn_arr.size_cols());

		// Get data from resource files. 
		std::ifstream inputfile(filename);

		// Check if file is open
		if (!inputfile.is_open()) {
			std::cout << "\n\nNo file found...\n\n";
			exit(-1);
		}

		// Count for while loop index for row tracking
		unsigned int count_i = 0;
		bool end_line = false;

		while (std::getline(inputfile, line) && end_line == false) {

			//  If this line is the end...
			end_line = inputfile.eof();

			// If delete first character of line desired
			if (del_first_char == true) { line.erase(line.begin()); }

			// Instantiate string stream with row data from text file. 
			std::stringstream ss(line);

			// Loop over collumns of data (note, need to make a check for this that you are not going over limits) 
			for (unsigned int col = 0; col < dyn_arr.size_cols(); col++) {

				if (col != col - 1)
					std::getline(ss, data_row[col], delimiter);
				else
					std::getline(ss, data_row[col]);

				// Convert string into integer values
				dyn_arr(count_i, col) = std::stoi(data_row[col]);
			}

			// Incremenet the row count
			count_i++;

			// Check that data size matches with input stack array (avoid accessing unallocated memory.)
			if (count_i > dyn_arr.size_rows()) {
				std::cout << "\nError: input DYN_C2D array does not match text file data dimensions. Rows in text file: " << count_i << ".\n" << "Rows in DYN_C2D array: " << dyn_arr.size_rows() << ".\n";
				exit(-1);
			}

		}
		inputfile.close();
	};


	//Convert std::string vectors to int, float and doubles. 
	std::vector<std::vector<int>> convert_string_vector_int(std::vector<std::string>& string_vect, char delimiter) {

		// Set up output vector (with appropiate dtype)
		std::vector< std::vector<int> > out_vect; 
		out_vect.resize(string_vect.size()); 

		bool warning_decimal_number_to_int_floor = false; 

		// Loop down vector and convert each string individually 
		for (size_t i = 0; i < string_vect.size(); i++) {

			// String container for each seperate numeric value
			std::string container;

			for (size_t j = 0; j < string_vect[i].size(); j++) {

				// While looping down the string, extract each figure up until the delimiter
				if (string_vect[i][j] == ' ') {}									// Do nothing if it is a space
				else if (string_vect[i][j] == '.' || string_vect[i][j] =='f') {
					warning_decimal_number_to_int_floor = true; 
					container.push_back(string_vect[i][j]);
				}
				else if (string_vect[i][j] != delimiter) {							// Add on any character which is not the delimiter, space or 'f'

					// Add character which is part of the numeric value
					container.push_back(string_vect[i][j]);
				}
				else {																// Delimiter found
					
					// Add numeric value to the output vector in corresponding row, i.
					out_vect[i].push_back(std::stoi(container));

					// Clear the current container
					container.clear();
				}
			}

			// Accoutn for last index
			// Add numeric value to the output vector in corresponding row, i.
			out_vect[i].push_back(std::stoi(container));
			container.clear();
		}

		if (warning_decimal_number_to_int_floor == true) {
			std::cout << "Warning: Converted decimal value to int, numeric values have been floored.\n";
		}

		// Return the converted vector
		return out_vect;
	}

	std::vector<std::vector<float>> convert_string_vector_float(std::vector<std::string>& string_vect, char delimiter) {

		// Set up output vector (with appropiate dtype)
		std::vector< std::vector<float> > out_vect;
		out_vect.resize(string_vect.size());

		// Loop down vector and convert each string individually 
		for (size_t i = 0; i < string_vect.size(); i++) {

			// String container for each seperate numeric value
			std::string container;

			for (size_t j = 0; j < string_vect[i].size(); j++) {

				// While looping down the string, extract each figure up until the delimiter

				if (string_vect[i][j] == ' ' || string_vect[i][j] ==  'f') {}		// Do nothing if it is a space or float char 'f'
				else if (string_vect[i][j] != delimiter) {							// Add on any char to container which is no delimiter, space or 'f'

					// Add character which is part of the numeric value
					container.push_back(string_vect[i][j]);
				}
				else {																// Delimiter found

					// Add numeric value to the output vector in corresponding row, i.
					out_vect[i].push_back(std::stof(container));

					// Clear the current container
					container.clear();
				}
			}
			// Accoutn for last index
			// Add numeric value to the output vector in corresponding row, i.
			out_vect[i].push_back(std::stof(container));
			container.clear();
		}

		// Return the converted vector
		return out_vect;
	}

	std::vector<std::vector<double>> convert_string_vector_double(std::vector<std::string>& string_vect, char delimiter) {

		// Set up output vector (with appropiate dtype)
		std::vector< std::vector<double> > out_vect;
		out_vect.resize(string_vect.size());

		// Loop down vector and convert each string individually 
		for (size_t i = 0; i < string_vect.size(); i++) {

			// String container for each seperate numeric value
			std::string container;

			for (size_t j = 0; j < string_vect[i].size(); j++) {

				// While looping down the string, extract each figure up until the delimiter
				if (string_vect[i][j] == ' ' || string_vect[i][j] == 'f') {}		// Do nothing if it is a space or 'f' char
				else if (string_vect[i][j] != delimiter) {							// Add on any char which isnt delimiter, space or 'f'

					// Add character which is part of the numeric value
					container.push_back(string_vect[i][j]);
				}
				else {																// Delimiter found

					// Add numeric value to the output vector in corresponding row, i.
					out_vect[i].push_back(std::stod(container));

					// Clear the current container
					container.clear();
				}
			}

			// Accoutn for last index
			// Add numeric value to the output vector in corresponding row, i.
			out_vect[i].push_back(std::stod(container));
			container.clear();
		}

		// Return the converted vector
		return out_vect;
	}


	/*=================================== Data types ===================================*/
	/* Datatime class for operating with dates*/
	class datetime {
	private:
		std::string datetime_str;

	public:
		/* Member variables*/
		int year;
		int month;
		int day;
		int hour;
		int minute;
		int second;

		/* Constructor */
		datetime() : year(0), month(0), day(0), hour(0), minute(0), second(0) {};

		// Constructor with vars;  
		datetime(const char* datetime, const char* fmt) {

			// Convert datetimes to string for standard lib string operations. (const char* allows more flexible input)
			datetime_str = datetime;

			/* Different formats for datetime (loads to add!) */
			if (fmt == "YYYY:MM:DD HH:MM") {

				// Check if dates are in CORRECT format. 
				if (datetime[4] != ':' || datetime[7] != ':' || datetime[10] != ' ' || datetime[13] != ':') {
					std::cout << "Error: format of string in Datetime incorrect. Must be in format 'YYYY:MM:DD HH:MM'.\n";
					exit(-1);
				}

				// Extract variables
				year = std::stoi(datetime_str.substr(0, 4));
				month = std::stoi(datetime_str.substr(5, 2));
				day = std::stoi(datetime_str.substr(8, 2));
				hour = std::stoi(datetime_str.substr(11, 2));
				minute = std::stoi(datetime_str.substr(14, 2));

				// Unspecified variables = 0 
				second = 0;
			}
			else if (fmt == "YYYY-MM-DD HH:MM") {

				// Check if dates are in CORRECT format. 
				if (datetime[4] != '-' || datetime[7] != '-' || datetime[10] != ' ' || datetime[13] != '-') {
					std::cout << "Error: format of string in Datetime incorrect. Must be in format 'YYYY:MM:DD HH:MM'.\n";
					exit(-1);
				}

				// Extract variables
				year = std::stoi(datetime_str.substr(0, 4));
				month = std::stoi(datetime_str.substr(5, 2));
				day = std::stoi(datetime_str.substr(8, 2));
				hour = std::stoi(datetime_str.substr(11, 2));
				minute = std::stoi(datetime_str.substr(14, 2));

				// Unspecified variables = 0 
				second = 0;
			}
			else {
				std::cout << "Error: please input a correct format (fmt)\n";
				exit(-1);
			}
		}

		/* Member functions */

		// Print datetime 
		void print(const char* fmt) {

			if (fmt == "YYYY:MM:DD HH:MM") {
				if (hour < 10 && minute < 10)
					std::cout << year << ":" << month << ":" << day << " 0" << hour << ":0" << minute;
				else if (hour < 10 && minute >= 10)
					std::cout << year << ":" << month << ":" << day << " 0" << hour << ":" << minute;
				else
					std::cout << year << ":" << month << ":" << day << " " << hour << ":" << minute;
			}
			else if (fmt == "YYYY-MM-DD HH:MM") {


				if (hour < 10 && minute < 10)
					std::cout << year << "-" << month << "-" << day << " 0" << hour << ":0" << minute;
				else if (hour < 10 && minute >= 10)
					std::cout << year << "-" << month << "-" << day << " 0" << hour << ":" << minute;
				else
					std::cout << year << "-" << month << "-" << day << " " << hour << ":" << minute;
			}
			else {
				std::cout << "Error, please enter a format YYYY:MM:DD HH:MM or YYYY-MM-DD HH:MM\n";
				exit(-1);
			}
		}

		// Calculate difference in datetimes (itself and another datetime. 
		// Need to account for different number of days in each month, particularly february! Leaving it simple for now till i need this functionality. 
		double difference(datetime& different_date, const char* unit) {

			// Difference vars
			double diff_year = std::abs(year - different_date.year);
			double diff_month = std::abs(month - different_date.month);
			double diff_day = std::abs(day - different_date.day);
			double diff_hour = std::abs(hour - different_date.hour);
			double diff_minute = std::abs(minute - different_date.minute);
			double diff_second = std::abs(second - different_date.second);

			if (unit == "second") {
				return ((diff_second)+(diff_minute * 60) + (diff_hour * 3600) + (diff_day * 86400) + (diff_month * 2678400) + (diff_year * 32140800));
			}
			if (unit == "minute") {
				return ((diff_second / 60) + (diff_minute)+(diff_hour * 60) + (diff_day * 1440) + (diff_month * 44640) + (diff_year * 535680));
			}
		}

		// Filling an empty datetime instance with a date (allows for default constructor to be used and not be pointless) 
		void fromstring(const char* datetime, const char* fmt) {

			// Convert datetimes to string for standard lib string operations. (const char* allows more flexible input)
			datetime_str = datetime;

			/* Different formats for datetime (loads to add!) */
			if (fmt == "YYYY:MM:DD HH:MM") {

				// Check if dates are in CORRECT format. 
				if (datetime_str[4] != ':' || datetime_str[7] != ':' || datetime_str[10] != ' ' || datetime_str[13] != ':') {
					std::cout << "Error: format of string in Datetime incorrect. Must be in format 'YYYY:MM:DD HH:MM'.\n";
					exit(-1);
				}

				// Extract variables
				year = std::stoi(datetime_str.substr(0, 4));
				month = std::stoi(datetime_str.substr(5, 2));
				day = std::stoi(datetime_str.substr(8, 2));
				hour = std::stoi(datetime_str.substr(11, 2));
				minute = std::stoi(datetime_str.substr(14, 2));

				// Unspecified variables = 0 
				second = 0;
			}
			else if (fmt == "YYYY-MM-DD HH:MM") {

				// Check if dates are in CORRECT format. 
				if (datetime_str[4] != '-' || datetime_str[7] != '-' || datetime_str[10] != ' ' || datetime_str[13] != ':') {
					std::cout << "Error: format of string in Datetime incorrect. Must be in format 'YYYY:MM:DD HH:MM'.\n";
					exit(-1);
				}

				// Extract variables
				year = std::stoi(datetime_str.substr(0, 4));
				month = std::stoi(datetime_str.substr(5, 2));
				day = std::stoi(datetime_str.substr(8, 2));
				hour = std::stoi(datetime_str.substr(11, 2));
				minute = std::stoi(datetime_str.substr(14, 2));

				// Unspecified variables = 0 
				second = 0;
			}
			else {
				std::cout << "Error: please input a format (fmt)\n";
				exit(-1);
			}

		}

		// Return minutes from year start
		double minutes_since_start_of_year() {

			return (second / 60) + (minute)+(hour * 60) + (day * 1440) + (month * 44640);// +(year * 535680);
		}
	};

	class Timer {
	private:
		const char* stackframeName;		// Name of stackframe for timer

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
			end = std::chrono::high_resolution_clock::now();		// end the timer
			duration = end - start;									// Calculate duration
			std::cout << "Time taken for " << stackframeName << ": " << duration.count() << " s\n";
		}
	};

}

