#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include "user.h"
#include <stdlib.h>
#include <vector>
#include <iomanip>

void print_help()
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
}

int parse_arguments(const int argc, char * argv[], char * &time_series, 
	char * &interactions, char * &prefix, int &lookahead, int &width)
{
	if(argc <= 1)
	{
		print_help();
		return 1; //error
	}
	else
	{
		width = 4;
		lookahead = 4;

		//Start to process the commandline arguments
		for(int i = 1; i < argc; i++)
		{
			//std::cout << argv[i] << std::endl;

			//ensure that the first 2 characters are --
			char * current_arg = argv[i];
			char * option;
			char * option_value;

			if(strncmp(current_arg, "--", 2) == 0)
			{
				//now split the string using '=' as delimiter
				option = strtok(&current_arg[2], "=");
				//std::cout << option << std::endl;

				option_value = strtok(NULL, "=");
				//std::cout << option_value << std::endl;

				if(strcmp(option, "time_series") == 0)
				{
					time_series = option_value;
				}
				else if(strcmp(option, "interactions") == 0)
				{
					interactions = option_value;
				}
				else if(strcmp(option, "prefix") == 0)
				{
					prefix = option_value;
				}
				else if(strcmp(option, "lookahead") == 0)
				{
					lookahead = atoi(option_value);
				}
				else if(strcmp(option, "width") == 0)
				{
					width = atoi(option_value);
				}
			} // if(strncmp(current_arg, "--", 2) == 0)
			else
			{
				std::cout << argv[i] << " needs to begin with '--'." << std::endl;
				return 1; //error
			}
		} // for(int n = 1; n < argc; n++)
	} // if(argc > 1)
	return 0;
}

int split(std::string &tuple, char delim)
{
	int pos = tuple.find_first_of(delim);
	return atoi(tuple.substr(pos + 1).c_str());
}

/*
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
*/

double * get_posterior(std::string &line, const int T, int &K)
{
	int pos = line.find_first_of(':');

	//std::cout << line << std::endl;

	std::string tmp1 = line.substr(pos + 3, line.length() - pos - 5);
	//std::cout << tmp1 << std::endl;

	//estimate the K
	pos = (T == 1) ? tmp1.length() : tmp1.find("),(");
	std::string cur = tmp1.substr(0, pos);

	//std::cout << cur << std::endl;


	std::stringstream ss1(cur);
	std::string tmp2;

	std::vector<double> vec;

	while(std::getline(ss1, tmp2, ','))
	{
		vec.push_back(atof(tmp2.c_str()));
	}
	K = vec.size();

	//std::cout << "K=" << K << std::endl;

	double * posterior = new double[T * K];
	int k=0;
	for(std::vector<double>::iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		posterior[k++] = *iter;
	}

	for(int t=1;t<T;t++)
	{
		tmp1 = tmp1.substr(pos + 3);
		pos = (t == T-1) ? tmp1.length() : tmp1.find("),(");
		cur = tmp1.substr(0, pos);
		//std::cout << cur << std::endl;

		ss1.str(cur);
		ss1.clear();
		ss1.seekg(0, std::ios::beg);
		
		k=0;
		while(std::getline(ss1, tmp2, ','))
		{
			//std::cout << tmp2 << std::endl;
			//std::cout << "T=" << T << ", K=" << K << ", TxK=" << T*K << ", t=" << t 
			//	<< ", k=" << k << ", txK+k=" << t*K + k << std::endl;
			posterior[t*K + k] = atof(tmp2.c_str()) ;
			k++;
		}
	}
	return posterior;
}

void write_gc(char * prefix, std::map<std::string, User *> &users)
{
	std::string out_file(prefix);
	out_file.append("_gc.txt");
	std::ofstream out_stream(out_file.c_str());

	if(out_stream.is_open())
	{
		User * user;
		out_stream << std::fixed << std::setprecision(8);
		for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
		{
			user = iter->second;
			if(user->Tn < 1) continue;

			for(std::map<int, std::map<User *, double> *>::iterator iter2 = user->neighbors_gc.begin(); 
				iter2 != user->neighbors_gc.end(); iter2++)
			{
				int t = iter2->first;
				out_stream << user->id << '\t' << t;
				
				for(std::map<User *, double>::iterator iter3 = iter2->second->begin(); 
					iter3 != iter2->second->end(); iter3++)
				{
					if(iter3 == iter2->second->begin())
						out_stream << '\t';
					else
						out_stream << ' ';
					
					out_stream << iter3->first->id << ':' << iter3->second;
				}
				
				out_stream << std::endl;
			}
		}
		out_stream.close();
	}
	else
	{
		std::cout << "Unable to open " << out_file << std::endl;
	}
}

int main(int argc, char * argv[])
{
	int t;
	std::string n, line, tuple;
	std::string tf, tl, mu, posterior;
	int lookahead, width;
	char * time_series;
	char * interactions;
	char * prefix;
	std::map<std::string, User *> users;
	User * user, * tar;

	if(parse_arguments(argc, argv, time_series, interactions, prefix, lookahead, width))
	{
		return 1;
	}

	//print the arguments
	std::cout << "Time series file: " << time_series << ", Interactions file: " << interactions << ", Prefix: " << prefix << ", Lookahead: " << lookahead << ", Width: " << width << std::endl;

	//read the users parameters
	std::ifstream * file_ptr = new std::ifstream(time_series);
	if (!file_ptr->is_open())
	{
		std::cout << "Unable to open " << time_series << std::endl;
	}

	while(file_ptr->good())
	{
		getline(*file_ptr, line);
		
		if(line.length() < 1) break;

		//std::cout << line << std::endl;

		std::istringstream iss (line, std::istringstream::in);

		iss >> n >> tf >> tl >> mu >> posterior;
		//std::cout << n << ", " << tf << ", " << tl << ", " << mu << ", " << posterior << std::endl;

		int t_first = split(tf, ':');
		int t_last  = split(tl, ':');
		
		//std::cout << t_first << ", " << t_last << std::endl;

		//get posterior
		int K;
		double * vec = get_posterior(posterior, t_last - t_first + 1, K);
		
		if(users.find(n) == users.end())
		{
			//create new user here
			user = new User(n, t_first, t_last, vec, K);
			users[n] = user;
		}
		else
		{
			std::cout << n << std::endl;
			std::cout << "Error: Duplicate entry in the parameters file." << std::endl;
		}
	}

	file_ptr->close();
	delete file_ptr;

	//read the uesrs social network
	
	file_ptr = new std::ifstream(interactions);
	if (!file_ptr->is_open())
	{
		std::cout << "Unable to open " << interactions << std::endl;
	}

	while (file_ptr->good())
	{
		getline(*file_ptr, line);

		if(line.length() < 1) break;
		
		// Open the line as a stream
		std::istringstream iss (line, std::istringstream::in);
		
		// Get first two values 
		iss >> n >> t;
		
		std::map<std::string, User *>::iterator iter;
		iter = users.find(n);

		if(iter == users.end())
		{
			//create new user here
			user = new User(n);
			users[n] = user;
		}
		else
		{
			user = iter->second;
		}

		while(iss.good())
		{
			iss >> tuple;
			int pos = tuple.find_first_of(':');
			
			std::string tar_id = tuple.substr(0, pos);

			std::map<std::string, User *>::iterator iter;
			iter = users.find(tar_id);
			if(iter == users.end())
			{
				tar = new User(tar_id);
				users[tar_id] = tar;
			}
			else
			{
				tar = iter->second;
			}


			int count = atoi(tuple.substr(pos+1).c_str());

			user->add_social(t, tar, count);
		}
	}
	file_ptr->close();
	delete file_ptr;
	
	//calculate the GC

	for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
	{
		user = iter->second;
		user->calculate_gc(width, lookahead);
	}

	//write out the results
	write_gc(prefix, users);

	return 0;
}
