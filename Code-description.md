# mathematical-function-static-library
Numerical function library which can be linked to give functionality similar to Matlab and Python's numerical libraries. 

Program objectives: 

          - To provide any of my projects with a comprehensive list of functions to provide functionality 
            similar to Matlab and Python. 
                    
                    
Notes on code structure: 

          - All functions pass in non-primitive parameters by reference to, if applicable,  preserve size information 
            and/or avoid unneccessary copying. 
            
          - All functions are optimized for maximum execution time. 
          
          - All functions are guarded for invalid inputs. 
          
          - The library has both compile and runtime variants of the functions.
          
          - Compile time variants are template functions with type and non-type parameters which is guarded at compile time. 
          
          - Runtime variants are template functions with type parameters where guards are directly implemented and
            evaluated at runtime. 
          
          - Runtime functions use the std::vector for the 1/2/3D arrays. A function to correct the memory of 
            the std::vector has been created 
            for optimization of functions (e.g. reserve, shave off any excess memory etc. to avoid copying and wasted memory).
            
          - For multi-dimensional arrays, using std::vector< std::vector< prim_type> > is innefficiant since the 
            memory is everywhere; a contiguous array class was created to wrap a 1D vector with functionality to access 
            it in a 2D way. Performance increased in:
              -> Execution time: 31% when using an iterator to loop over data. 
              -> Memory usage:   32% (contigous array class uses same memory as raw dynamic array) 

          
          
Learning objectives:

          - Code safety 
          
          - Use of template functions 
          
          - Creation of interfaces 
          
          - General programming of functions in compile and runtime enviroments.
          
          - Optimization of functions and creating own data structures 
          
          
          
