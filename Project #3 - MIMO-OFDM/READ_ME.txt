This file explains the contents of the files and folders present in this folder.

Note: Some code files appear in multiple folders. With the exception of main they are the exact same code. The main file is different in all three parts.

Folders:
	- Part1_MIMO	: Contains the code and functions for part 1 of this project
		- html 			: Contains published code and results
		- clarkeSpec 		: Generates the clarke spectrum values for a specified number fo points
		- genMIMOdata 		: Generates data for a 2x2 MIMO channel
		- getRayleighFading	: Returns a Rayleigh fading specified by a given duration, Doppler frequency, and number of points
		- main 			: Run this one to see the results
		- mmse_zfMIMO 		: Returns BER for MMSE and ZF equalizers for the given data and channel
		- precodingMIMO 	: Returns BER for SVD channel encoding case
		- rayleighFading2 	: Helps generate Rayleigh fading samples

	- Part2_MIMO 	: Contains the code and functions for part 2 of this project
		- html 			: Same as above
		- bestOFDMconfig	: Finds the best OFDM rate and modulation configurations for the given channels
		- main 			: Same as above
		- OFDMdemod 		: Sends the OFDM symbol through the channel and demodulates it. Returns BER
		- OFDMsymGen 		: Genereates an OFDM symbols given the rate and modulation order

	- Part3_MIMO 	: Contains the code and functions for part 3 of this project
		- clarkeSpec		: Same as above
		- genMIMO_OFDMdata 	: Generate a MIMO-OFDM data stream according to the input parameters
		- getRayleighFading 	: Same as above
		- main 			: Same as above
		- mmse_MIMO_OFDM 	: Sends the data and recovers it using MMSE and ZF equalizers. Returns BER
		- OFDMdemod		: Same as above
		- OFDMsymGen 		: Same as above
		- precodingMIMO_OFDM 	: Performs SVD on channel to send MIMO-OFDM signal	
		- rayleighFading2 	: Same as above

Files:
	- Project #3 - MIMO-OFDM Report : The report for this project
	- READ_ME 	: Reference READ_ME for more information.