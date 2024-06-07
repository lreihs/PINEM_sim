					--- General info ---

This project is aimed towards simulating quantum coherent optical phase modulation in an ultrafast 
transmission electron microscope. 

This README file will help with navigating the code.

The code is split into several files, each divided into subsections. A short description of each file
is given below, see further down for an in-depth description. 

_______________________________________________________________________________________________________

					--- File overview ---

01_params		File containing input parameters for all simulation tasks.

02_functions		Code file contiaing all functions involved in all simulation tasks.

03_electron_pulse	Code file used for vizualizing electron energy spectra post interaction with
			the optical field. 

04_propagation  	Code file used for propagation of the elctron pulse post interaction with the
			optical field.

_______________________________________________________________________________________________________

					--- File subsections ---

01_params		1.1	General parameters: 	atomic coordinates: 0 => SI units
							atomic coordinates: 1 => atomic units

			1.2	Parameters defining the initial electron pulse

			1.3	Paramters defining the optical field interacting with the initial  
				electron pulse parametrized in section 1.2. Multiple spectra 
				parameters define a range in optical field strengths interacting 
				with the inital electron pulse to produce field dependent spectra (see
				section 3.2 for visualization).  

			1.4	Parameters defining the propagation of the electron pulse post
				intercation with the optical field.

02_functions		2.1	Importing all parameters from '01_params'. The path to '01_params' 
				in line 10 needs to be set accordingly. 

			2.2	Defining additional parameters from those inputted in '01_params'

			2.3	Functions  

			2.4

03_electron_pulse	3.1

			3.2

04_propagation		4.1

_______________________________________________________________________________________________________

					--- References ---

[1]	

[2]

[3]

_______________________________________________________________________________________________________

Lennard Reihs, 2024
