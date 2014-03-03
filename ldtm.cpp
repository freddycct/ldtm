#include <fstream>
#include <iostream>
#include <sstream>
#include <gsl/gsl_rng.h>
#include "user.h"
#include <vector>
#include <map>
#include <string.h>
#include <ctime>
#include <iomanip>
#include <stdlib.h>

double kl_divergence(double * p, double * q, int K)
{
	double sum = 0;
	for(int k=0;k<K;k++)
	{
		sum += p[k] * log(p[k]/q[k]) ;
	}
	return sum;
}

double js_divergence(double * p, double * q, int K)
{
	double sum = 0;
	double * m = new double[K];
	for(int k=0;k<K;k++)
	{
		m[k] = 0.5 * (p[k] + q[k]);
	}
	sum = 0.5 * (kl_divergence(p, m, K) + kl_divergence(q, m, K)) ;
	delete [] m;
	return sum;
}

int sample(const int M, const int K, const double * prior_x, const int * phi, 
			const int w, const int * lambda, const int * sum_lambda, const gsl_rng * r, 
			const double ALPHA, const double BETA, double * p, double * cum_p)
{	
	double sum_p = 0;
	for(int k = 0; k < K; k++)
	{
		p[k] = (prior_x[k] + phi[k] + ALPHA) * ((lambda[k * M + w] + BETA) / (sum_lambda[k] + M * BETA));
		sum_p += p[k];
	}
	double rand = gsl_rng_uniform(r);
	rand = rand * sum_p;
	cum_p[0] = p[0];
	for(int k = 0; k < K - 1; k++)
	{
		if(rand < cum_p[k])
		{
			return k;
		}
		else
		{
			cum_p[k+1] = cum_p[k] + p[k+1];
		}
	}
	
	return K-1;
}

void print_help()
{
	// lda_lds(Yt, M, N, T, K, ITERATIONS, 0.5, false, false);
	// lda_lds(Yt, M, N, T, K, ITERATIONS, 0.5, false, true, false, false, 10.0);

	std::cout << "Linear Dynamical Topic Model version 0.1" << std::endl << std::endl;
	std::cout << "Usage: cat <file_name> | ldtm <options>" << std::endl << std::endl;
	std::cout << "ldtm reads in from standard input." << std::endl << std::endl;
	std::cout << "Options:" << std::endl << std::endl;
	//std::cout << "\t--input=<file name>" << std::endl;
	std::cout << "\t--num_topics=<integer value>" << std::endl;
	std::cout << "\t--iterations=<integer value>" << std::endl;
	std::cout << "\t--init_mu=<real value between 0 and 1>" << std::endl;
	std::cout << "\t--learning_rate=<real value>" << std::endl;
	std::cout << "\t--prefix=<string> : This string is used to name the output files." << std::endl;
	std::cout << "\t--record_logll=<true/false> : Optional, Record the Log Likelihood. Default value is false." << std::endl;
	std::cout << "\t--estimate=<true/false> : Optional, Estimate the dynamics matrix. Default is true." << std::endl;
	std::cout << "\t--single_A=<true/false> : Optional, Dynamics matrix for every time step. Default is false." << std::endl;
	std::cout << "\t--single_mu=<true/false> : Optional, Multiple or single mu in dynamics matrix. Default is false." << std::endl;

	std::cout << std::endl;
}

int parse_arguments(const int argc, char * argv[], int &K, 
	int &iterations, double &mu, char * &prefix, double &learning_rate, bool &estimate, 
	bool &single_A, bool &single_mu, bool &record_logll)
{
	if(argc <= 1)
	{
		print_help();
		return 1; //error
	}
	else
	{
		record_logll = false;
		estimate = true;
		single_A = false;
		single_mu = false;

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
				//std::string to_tokenize(current_arg);
				//int pos = to_tokenize.find_first_of("=");
				//current_arg[pos] = '\0';

				//option = &current_arg[2];

				option = strtok(&current_arg[2], "=");
				//std::cout << option << std::endl;

				option_value = strtok(NULL, "=");
				//option_value = &current_arg[pos+1];
				//std::cout << option_value << std::endl;

				/*
				if(strcmp(option, "input") == 0)
				{
					input_file = option_value;
				}
				*/

				if(strcmp(option, "num_topics") == 0)
				{
					K = atoi(option_value);
				}
				else if(strcmp(option, "iterations") == 0)
				{
					iterations = atoi(option_value);
				}
				else if(strcmp(option, "init_mu") == 0)
				{
					mu = atof(option_value);
				}
				else if(strcmp(option, "prefix") == 0)
				{
					prefix = option_value;
				}
				else if(strcmp(option, "learning_rate") == 0)
				{
					learning_rate = atof(option_value);
				}
				else if(strcmp(option, "record_logll") == 0)
				{
					if(strcasecmp(option_value, "true") == 0)
					{
						record_logll = true;
					}
				}
				else if(strcmp(option, "estimate") == 0)
				{
					if(strcasecmp(option_value, "false") == 0)
					{
						estimate = false;
					}
				}
				else if(strcmp(option, "single_A") == 0)
				{
					if(strcasecmp(option_value, "true") == 0)
					{
						single_A = true;
					}
				}
				else if(strcmp(option, "single_mu") == 0)
				{
					if(strcasecmp(option_value, "true") == 0)
					{
						single_mu = true;
					}
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

void init(std::map<std::string, User *> &users, const int M, const int K, 
	int * lambda, int * sum_lambda, gsl_rng * r)
{
	set_zero<int>(K, sum_lambda);
	set_zero<int>(M*K, lambda);

	User * user;
	//initialize
	for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
	{
		user = iter->second;
		for(std::map<int, std::list<int> *>::iterator iter2 = user->adoption.begin(); iter2 != user->adoption.end(); iter2++)
		{
			int t = iter2->first;
			std::list<int> * adoption_t = iter2->second;

			int adjusted_t = t - user->t_first;
			int * Z_t = user->Z[adjusted_t];

			int i=0;
			for(std::list<int>::iterator iter3 = adoption_t->begin(); iter3 != adoption_t->end(); iter3++)
			{
				int k = gsl_rng_uniform_int(r, K);
				int w = *iter3;

				lambda[k * M + w]++; //increment to the topic word distribution
	 			sum_lambda[k]++; //the denominator
	 			
	 			user->phi[adjusted_t * K + k]++;
	 			Z_t[i++] = k;
			}
		}
	}
}

void gibbs_sampling(const int iterations, const bool record_logll, const bool estimate, const int M, const int K, 
	const gsl_rng * r, const double ALPHA, const double BETA, std::map<std::string, User *> &users, 
	int * lambda, int * sum_lambda, double * logll, double * logll_kl, double * time_elapsed, const double learning_rate)
{
	
	User * user;
	double * prior_x        = new double[K];
	double * prior_dist     = new double[K];
	double * posterior_dist = new double[K];
	double * tmp_p = new double[K];
	double * cum_p = new double[K];
	
	for(int i = 0; i < iterations; i++)
	{
		printf("\nIterations: %d/%d", i+1, iterations);
		std::clock_t begin, end; //for calculating the time taken during each iteration

		if(record_logll)
		{
			logll[i] = 0;
			logll_kl[i] = 0;
			begin = std::clock();
		}

		for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
		{
			user = iter->second;
			int Tn = user->Tn;

			int prev_t = -1;
			for(std::map<int, std::list<int> *>::iterator iter2 = user->adoption.begin(); iter2 != user->adoption.end(); iter2++)
			{
				int adjusted_t = iter2->first - user->t_first;
				
				if(adjusted_t - prev_t > 1)
				{
					//have to handle the case when adjusted_t is not consecutive with previous
					for(int tt = prev_t + 1; tt <adjusted_t; tt++)
					{
						for(int k = 0; k < K; k++)
						{
							user->posterior_x[(tt * K) + k] = user->posterior_x[(tt-1)*K + k];
						}
					}
				}

				std::list<int> * adoption_t = iter2->second;
				int * Z_t = user->Z[adjusted_t];
				int * phi_t = &(user->phi[adjusted_t * K]);

				if(adjusted_t == 0)
				{
					set_zero<double>(K, prior_x);
				}
				else
				{
					user->get_prior_x(adjusted_t, K, prior_x);
				}

				int j = 0;
				for(std::list<int>::iterator iter3 = adoption_t->begin(); iter3 != adoption_t->end(); iter3++)
				{
					int w = *iter3;
					int k = Z_t[j];

					phi_t[k]--;
					lambda[k * M + w]--; //increment to the topic word distribution
	 				sum_lambda[k]--; //the denominator
					
					k = sample(M, K, prior_x, phi_t, w, lambda, sum_lambda, r, ALPHA, BETA, tmp_p, cum_p);

	 				sum_lambda[k]++; //the denominator
	 				lambda[k * M + w]++; //increment to the topic word distribution
					phi_t[k]++;

					Z_t[j++] = k; //store the newly inferred topic
				}

				for(int k = 0; k < K; k++)
				{
					user->posterior_x[(adjusted_t * K) + k] = prior_x[k] + phi_t[k];
				}
				prev_t = adjusted_t;
			}

			//learn dynamics here
			if(Tn > 1 && estimate)
			{
				user->estimate_dynamics(K, ALPHA, learning_rate, tmp_p);
			}

			if(record_logll)
			{
				end = std::clock();
				if(i==0)
					time_elapsed[i] = ((double)(end - begin) / CLOCKS_PER_SEC) ;
				else
					time_elapsed[i] = ((double)(end - begin) / CLOCKS_PER_SEC) + time_elapsed[i-1];

				// Compute Log Likelihood here
				// This part computes the log likelihood based on the observed words

				for(std::map<int, std::list<int> *>::iterator iter2 = user->adoption.begin(); iter2 != user->adoption.end(); iter2++)
				{
					int adjusted_t = iter2->first - user->t_first;
					std::list<int> * adoption_t = iter2->second;

					user->get_prior(adjusted_t, ALPHA, K, prior_dist);
					user->get_posterior(adjusted_t, ALPHA, K, posterior_dist);

					logll_kl[i] -= log(1 + kl_divergence(posterior_dist, prior_dist, K));

					for(std::list<int>::iterator iter3 = adoption_t->begin(); iter3 != adoption_t->end(); iter3++)
					{
						int w = *iter3;
						double p = 0;
						for(int k=0;k<K;k++)
						{
							p += ((lambda[ k*M + w ] + BETA) / (sum_lambda[k] + M*BETA)) *  posterior_dist[k] ;
						}
						logll[i] += log(p) ;
					}
				}
			} //if(record_logll)
		}
		printf(", logll: %f", logll[i]);
	} //for(int i = 0; i < iterations; i++)
	printf("\n");

	delete [] cum_p;
	delete [] tmp_p;
	delete [] posterior_dist;
	delete [] prior_dist;
	delete [] prior_x;
}

void write_logll(const char * prefix, const int iterations, const double * logll, 
	const double * logll_kl, const double * time_elapsed)
{
	std::string out_file(prefix);
	out_file.append("_logll.txt");

	std::ofstream out_stream(out_file.c_str());

	if(out_stream.is_open())
	{
		out_stream << std::fixed << std::setprecision(8);
		for(int i=0;i<iterations;i++)
		{
			out_stream << time_elapsed[i] << '\t' << logll[i] << '\t' << logll_kl[i] << std::endl;
		}
		out_stream.close();
	}
	else
	{
		std::cout << "Unable to open " << out_file << std::endl;
	}
}

void write_settings(const char * prefix, const int M, const int N, const int K, 
	const int iterations, const double mu, const double learning_rate, 
	const bool estimate, const bool single_A, const bool single_mu, const bool record_logll)
{
	std::string out_file(prefix);
	out_file.append("_settings.txt");

	std::ofstream out_stream(out_file.c_str());
	if(out_stream.is_open())
	{
		out_stream << "prefix=" << prefix << std::endl;
		out_stream << "number of items=" << M << std::endl;
		out_stream << "number of users=" << N << std::endl;
		out_stream << "num_topics=" << K << std::endl;
		out_stream << "iterations=" << iterations << std::endl;
		out_stream << "init_mu=" << mu << std::endl;
		out_stream << "learning_rate=" << learning_rate << std::endl;

		out_stream << "estimate=";
		if(estimate)
			out_stream << "true";
		else 
			out_stream << "false";
		out_stream << std::endl;

		out_stream << "single_A=";
		if(single_A)
			out_stream << "true";
		else 
			out_stream << "false";
		out_stream << std::endl;

		out_stream << "single_mu=";
		if(single_mu)
			out_stream << "true";
		else 
			out_stream << "false";
		out_stream << std::endl;

		out_stream << "record_logll=";
		if(record_logll)
			out_stream << "true";
		else 
			out_stream << "false";
		out_stream << std::endl;

		out_stream.close();
	}
	else
	{
		std::cout << "Unable to open " << out_file << std::endl;
	}
}

void write_topic(const char * prefix, const int K, const int M, 
	const int * lambda, const int * sum_lambda, const double BETA)
{
	std::string out_file(prefix);
	out_file.append("_topics.txt");
	std::ofstream out_stream(out_file.c_str());

	if(out_stream.is_open())
	{
		out_stream << std::fixed << std::setprecision(8);

		for(int m=0;m<M;m++)
		{
			out_stream << m + 1;
			for(int k=0;k<K;k++)
			{
				out_stream << '\t' << (lambda[k*M + m] + BETA) / (sum_lambda[k] + M*BETA) ;
			}
			out_stream << std::endl;
		}

		out_stream.close();
	}
	else
	{
		std::cout << "Unable to open " << out_file << std::endl; 
	}
}

void write_users(const char * prefix, const int K, std::map<std::string, User *> &users, 
	const double ALPHA, const bool single_mu, const bool single_A)
{
	std::string out_file(prefix);
	out_file.append("_users.txt");
	std::ofstream out_stream(out_file.c_str());

	double * tmp = new double[K];

	if(out_stream.is_open())
	{
		User * user;
		out_stream << std::fixed << std::setprecision(8);
		for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
		{
			user = iter->second;
			out_stream << iter->first << "\tt_first:" << user->t_first << "\tt_last:" << user->t_last;
			out_stream << "\tmu:[";

			int Tn = user->Tn;

			if(Tn > 1)
			{
				if(single_A && single_mu)
				{
					out_stream << user->mu_t[0];
				}
				else if(!single_A && single_mu)
				{
					for(int t=0;t<Tn-1;t++)
					{
						out_stream << user->mu_t[t];
						if(t < Tn-2)
							out_stream << ",";
					}
				}
				else if(single_A && !single_mu)
				{
					//not implemented!
				}
				else
				{
					for(int t=0;t<Tn-1;t++)
					{
						out_stream << "(";
						for(int k=0;k<K;k++)
						{
							//mu_t[t*K + k] = mu;
							out_stream << user->mu_t[t*K + k];
							if(k<K-1)
								out_stream << ",";
						}
						out_stream << ")";
						if(t<Tn-2)
							out_stream << ",";
					}
				}
			}
			out_stream << "]";

			out_stream << "\tposterior:[" ;
			for(int t=0;t<Tn;t++)
			{
				out_stream << "(";
				user->get_posterior(t, ALPHA, K, tmp);
				for(int k=0;k<K;k++)
				{
					out_stream << tmp[k];
					if(k<K-1)
						out_stream << ",";
				}
				out_stream << ")";
				if(t<Tn-1)
					out_stream << ",";
			}
			out_stream << "]" << std::endl;
		}
		out_stream.close();
	}
	else
	{
		std::cout << "Unable to open " << out_file << std::endl; 
	}

	delete [] tmp;
}

int main(int argc, char * argv[])
{
	//change it to read in from file instead of standard input
	//use map points in users

	int * lambda;
	int * sum_lambda;
	double * logll;
	double * logll_kl;
	double * time_elapsed;
	char * prefix;
	
	int M, N, K, iterations;
	double mu, learning_rate;
	bool estimate, single_A, single_mu, record_logll;
	std::map<std::string, User *> users;

	//parse the commandline arguments
	if(parse_arguments(argc, argv, K, iterations, mu, prefix, learning_rate, estimate, 
		single_A, single_mu, record_logll))
	{
		return 1;
	}

	//print the arguments
	std::cout << "Topics: " << K << ", Iterations: " << iterations << ", Initial Mu: " << mu 
		<< ", Learning Rate: " << learning_rate << ", Record Logll: " << record_logll 
		<< ", Estimate: " << estimate << ", single_A: " << single_A << ", single_Mu: " << single_mu 
		<< ", Prefix: " << prefix << std::endl;

	//Now read in the input file
	//std::ifstream * file_ptr = new std::ifstream(input_file);
	
	if(record_logll)
	{
		logll        = new double[iterations];
		logll_kl     = new double[iterations];
		time_elapsed = new double[iterations];
	}

	int t;
	std::string n, line, tuple;
	User * user;

	//while(std::cin.good())
	M = -1;
	while(!std::cin.eof())
	{
		getline(std::cin, line);
		//std::cout << line << std::endl;

		if(line.length() < 1) break;

		/* Open the line as a stream */
		std::istringstream iss(line, std::istringstream::in);
		
		/* Get first two values */
		iss >> n >> t;
		
		std::map<std::string, User *>::iterator iter;
		iter = users.find(n);

		if(iter == users.end())
		{
			//create new user here
			user = new User(n, single_A, single_mu);
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
			int m = atoi(tuple.substr(0, pos).c_str());
			int count = atoi(tuple.substr(pos+1).c_str());

			if (m > M)
			{
				M = m;
				// if(M > 100000)
				// {
				// 	std::cout << "error has occured " << M << std::endl;
				// }
			}

			user->add(t, m, count);
		}
		// if(M > 100000)
		// {
		// 	std::cout << "error has occured again and again " << M << std::endl;
		// }
	}

	N = users.size();
	//std::cout << N << ", " << M << std::endl;
	//std::cout.flush();

	for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
	{
		iter->second->init_mu(K, mu);
	}

	lambda = new int[M*K];
	sum_lambda = new int[K];

	//Setup the random generator
	gsl_rng * r;
	const gsl_rng_type * rng_T;
	gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	r = gsl_rng_alloc (rng_T);
	//End of random generator

	//initialize the gibbs sampler
	std::cout << "Initialize Gibbs" << std::endl;
	//std::cout.flush();

	init(users, M, K, lambda, sum_lambda, r);

	double ALPHA =  50.0 / K;
	double BETA  = 200.0 / M;

	//gibbs sampling inference
	std::cout << "Gibbs Sampling" << std::endl;
	//std::cout.flush();

	gibbs_sampling(iterations, record_logll, estimate, M, K, r, ALPHA, BETA, users, 
		lambda, sum_lambda, logll, logll_kl, time_elapsed, learning_rate);

	//output the results

	//output the settings
	std::cout << "Writing settings file" << std::endl;
	//std::cout.flush();
	write_settings(prefix, M, N, K, iterations, mu, learning_rate, estimate, 
		single_A, single_mu, record_logll);

	std::cout << "Writing logll file" << std::endl;
	//std::cout.flush();
	//output the log likelihood
	write_logll(prefix, iterations, logll, logll_kl, time_elapsed);

	std::cout << "Writing topic file" << std::endl;
	//std::cout.flush();
	//output the topic file
	write_topic(prefix, K, M, lambda, sum_lambda, BETA);

	std::cout << "Writing users file" << std::endl;
	//std::cout.flush();
	//output the users
	write_users(prefix, K, users, ALPHA, single_A, single_mu);
	
	std::cout << "Free memory" << std::endl;
	//std::cout.flush();
	//Free memory this is so important, omg....
	gsl_rng_free(r);

	delete [] sum_lambda;
	delete [] lambda;

	for(std::map<std::string, User *>::iterator iter = users.begin(); iter != users.end(); iter++)
	{
		delete iter->second;
	}

	if(record_logll)
	{
		delete [] time_elapsed;
		delete [] logll_kl;
		delete [] logll;
	}
	
	return 0;
}
