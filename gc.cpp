#include <iostream>

int parse_arguments(const int argc, char * argv[], char * time_series, 
	char * interactions, char * prefix, int lookahead, int width)
{
	std::cout << "Granger Causality version 0.1" << std::endl << std::endl;
	std::cout << "Usage: gc <options>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl << std::endl;
	std::cout << "\t--time_series=<file name>" << std::endl;
	std::cout << "\t--interactions=<file name>" << std::endl;
	std::cout << "\t--prefix=<string> : This string is used to name the output files." << std::endl;
	std::cout << "\t--lookahead=<integer value> : Optional, default value is 4" << std::endl;
	std::cout << "\t--width=<integer value> : Optional, default value is 4" << std::endl;
	
	std::cout << std::endl;
	return 0;
}

int main(int argc, char * argv[])
{
	int lookahead, width;
	char * time_series;
	char * interactions;
	char * prefix;

	if(parse_arguments(argc, argv, time_series, interactions, prefix, lookahead, width))
	{
		return 1;
	}

	//print the arguments
	std::cout << "Time series file: " << time_series << ", Interactions file: " << interactions << ", Prefix: " << prefix << ", Lookahead: " << lookahead << ", Width: " << width << std::endl;

	//read the users parameters

	//read the uesrs social network

	return 0;
}
