#include "user.h"

void User::init_mu(const int K, const double mu)
{
	Tn = t_last - t_first + 1;
	Z = new int * [Tn];
	for(int t=0;t<Tn;t++)
	{
		if(adoption.find(t + t_first) != adoption.end())
			Z[t] = new int[adoption[t + t_first]->size()];
		else
			Z[t] = NULL;
	}

	phi = new int[Tn * K];
	for(int i = 0; i < Tn * K; i++)
	{
		phi[i] = 0;
	}

	posterior_x = new double[Tn * K];

	if(Tn > 1)
	{
		if(single_A && single_mu)
		{
			mu_t = new double[1];
			mu_t[0] = mu;
		}
		else if(!single_A && single_mu)
		{
			mu_t = new double [Tn - 1];
			for(int t=0;t<Tn-1;t++)
			{
				mu_t[t] = mu;
			}
		}
		else if(single_A && !single_mu)
		{
			//not implemented!
		}
		else
		{
			mu_t = new double[ (Tn-1) * K ] ; //the arrangement should be K x (Tn-1)
			for(int t=0;t<Tn-1;t++)
			{
				for(int k=0;k<K;k++)
				{
					mu_t[t*K + k] = mu;
				}
			}
		}
	} //if(Tn > 1)
}

User::User(bool single_A, bool single_mu)
{
	this->single_A  = single_A;
	this->single_mu = single_mu;
	
	t_first = -1;
	t_last = -1;
	Tn = 0;
}

void User::add(const int t, const int m)
{
	//int adjusted_index = t - t_first;
	//std::cout << adjusted_index << std::endl;

	if(adoption.find(t) == adoption.end())
	{
		adoption[t] = new std::list<int>();		

		if(t > t_last)
		{
			t_last = t;
		}

		if(t < t_first || t_first < 0)
		{
			t_first = t;
		}
		
	}
	if(m < 1)
	{
		//throw an error
		std::cout << "error occurred here, m must be greater than or equal to one." << std::endl;
	}
	
	adoption[t]->push_back(m - 1);
	//Z[t]->push_back(-1);
}
/*

User::User(int M, int T, int K, 
	double * sr, int * irs, int * jcs, 
	double mu, bool single_A, bool single_mu)
{
	t_first = getFirst(T, jcs);
	t_last  = getLast(T, jcs);
	Tn = t_last - t_first + 1;
		
	adoption = new int * [Tn];
	Z   = new int * [Tn];
	Mn = new int[Tn];
	this->single_A = single_A;
	this->single_mu = single_mu;

	if(Tn > 1)
	{
		if(single_A && single_mu)
		{
			mu_t = new double[1];
			mu_t[0] = mu;
		}
		else if(!single_A && single_mu)
		{
			mu_t = new double [Tn - 1];
			for(int t=0;t<Tn-1;t++)
			{
				mu_t[t] = mu;
			}
		}
		else if(single_A && !single_mu)
		{
			//not implemented!
		}
		else
		{
			mu_t = new double[ (Tn-1) * K ] ; //the arrangement should be K x (Tn-1)
			for(int t=0;t<Tn-1;t++)
			{
				for(int k=0;k<K;k++)
				{
					mu_t[t*K + k] = mu;
				}
			}
		}
	}
	
	//Allocate space for gsl vectors & matrices
	phi         = new gsl_vector * [Tn];
	posterior_x = new gsl_vector * [Tn];
	//End
	
	int sr_ptr = 0;
	for(int t = 0; t < Tn; t++)
	{
		//find out the number of items in this time step
		int nz = jcs[t_first + t] - jcs[t_first + t - 1];
		Mn[t] = 0;
		for(int n = 0; n < nz; n++)
		{
			Mn[t] += sr[sr_ptr + n];
		}
		
		adoption[t] = new int[Mn[t]];
		Z[t] = new int[Mn[t]];
		
		int i = 0;
		for(int n = 0; n < nz; n++)
		{
			for(int k = 0; k < sr[sr_ptr + n]; k++)
			{
				adoption[t][i++] = irs[sr_ptr + n]; //minus one because of the difference in indexing between matlab and C
			}
		}
		sr_ptr += nz;
		
		phi[t]         = gsl_vector_calloc(K);
		posterior_x[t] = gsl_vector_calloc(K);
	}
}
*/

void User::get_prior_x(const int t, const int K, double * out)
{
	if(t == 0)
	{
		for(int k = 0; k < K; k++)
		{
			out[k] = 0;
		}
	}
	else
	{
		double * mu;
		if(single_A && single_mu)
		{
			//mu is a scalar
			mu = &mu_t[0];
		}
		else if(!single_A && single_mu)
		{
			mu = &mu_t[t-1];
		}
		else if(!single_A && !single_mu)
		{
			mu = &mu_t[(t-1)*K];
		}

		for(int k = 0; k < K; k++)
		{
			if(single_mu)
				out[k] = mu[0] * posterior_x[(t-1)*K + k] ;
			else
				out[k] = mu[k] * posterior_x[(t-1)*K + k] ;
		}
	}
}

void User::get_prior(const int t, const double alpha, const int K, double * out)
{
	get_prior_x(t, K, out);
	double sum = sum_vector(K, out);
	sum += K * alpha;

	for(int k = 0; k < K; k++)
	{
		out[k] = (out[k] + alpha) / sum ;
	}
}

void User::get_posterior(int t, const double alpha, const int K, double * out)
{
	double sum = 0;
	if(t == 0)
	{
		sum = sum_vector(K, &phi[t * K]);
		sum += K * alpha;
		for(int k = 0; k < K; k++)
		{
			out[k] = (phi[t*K + k] + alpha) / sum;
		}
	}
	else
	{
		double * mu;
		if(single_A && single_mu)
		{
			//mu is a scalar
			mu = &mu_t[0];
		}
		else if(!single_A && single_mu)
		{
			mu = &mu_t[t-1];
		}
		else if(!single_A && !single_mu)
		{
			mu = &mu_t[(t-1)*K];
		}

		for(int k = 0; k < K; k++)
		{
			if(single_mu)
				sum += mu[0] * posterior_x[(t-1)*K + k] + phi[t*K + k];
			else
				sum += mu[k] * posterior_x[(t-1)*K + k] + phi[t*K + k];
		}
		sum += K * alpha;
		for(int k = 0; k < K; k++)
		{
			if(single_mu)
				out[k] = (mu[0] * posterior_x[(t-1)*K + k] + phi[t*K + k] + alpha) / sum;
			else
				out[k] = (mu[k] * posterior_x[(t-1)*K + k] + phi[t*K + k] + alpha) / sum;
		}
	}
}

double User::get_decay(const int t)
{
	return single_A ? mu_t[0] : mu_t[t] ;
}

double * User::get_decay_vector(const int t, const int K)
{
	return &mu_t[t*K];
}

void User::estimate_dynamics(const int K, const double ALPHA, const double learning_rate)
{
	//Estimate the decay parameter here
	//First Calculate the Objective Function here.

	// double L = 0;
	// for(int t = 1; t < Tn; t++)
	// {
	// 	//obtain the prior
	//  get_prior(t, ALPHA, K, prior_dist);

	// 	//obtain the posterior
	// 	get_posterior(t, ALPHA, K, posterior_dist);

	// 	//get the KL divergence
	// 	L += kl_divergence(posterior_dist, prior_dist, K);
	// }
	// printf("n:%d, Before, L = %f\n", n, L);

	static double * posterior_dist = new double[K];

	double gradient = 0;
	for(int t=1;t<Tn;t++)
	{
		double sum_phi_nt = sum_vector<int>(K, &phi[t * K]);
		double sum_x_nt_1_t_1 = sum_vector<double>(K, &posterior_x[(t - 1) * K]);

		double mu_tmp;
		double * mu_t_vec;

		double mu_sum_x_nt_1_t_1;
		if(single_mu)
		{
			mu_tmp = get_decay(t-1);
			mu_sum_x_nt_1_t_1 = mu_tmp * sum_x_nt_1_t_1;
		}
		else
		{
			mu_t_vec = get_decay_vector(t-1, K);
			mu_sum_x_nt_1_t_1 = sum_vector<double>(K, &posterior_x[(t - 1) * K], mu_t_vec);
		}

		double denom = mu_sum_x_nt_1_t_1 + sum_phi_nt + K * ALPHA ;
		double denom2 = pow(denom, 2);

		double eqn3_2 = sum_x_nt_1_t_1 / denom ;
		double denom_2 = mu_sum_x_nt_1_t_1 + K * ALPHA ;

		get_posterior(t, ALPHA, K, posterior_dist);

		double gradient_t = 0;

		for(int k=0;k<K;k++)
		{
			double mu_k = single_mu ? mu_tmp : mu_t_vec[k] ;

			double x_nt_1_t_1 = posterior_x[(t-1)*K + k];
			double phi_nt_k = phi[t*K + k];

			//Eqn (1)
			double eqn1 = (x_nt_1_t_1 * (sum_phi_nt + K*ALPHA) 
				- sum_x_nt_1_t_1 * (phi_nt_k + ALPHA)) / denom2 ;

			//Eqn (2)
			double eqn2_1 = mu_k * x_nt_1_t_1 + phi_nt_k + ALPHA;
			double eqn2 = log(eqn2_1) - log(denom);

			//Eqn (3)
			double eqn3 = (x_nt_1_t_1/eqn2_1) - eqn3_2;

			//Eqn (4)
			double eqn4_1 = mu_k * x_nt_1_t_1 + ALPHA;
			double eqn4 = log(eqn4_1) - log(denom_2);

			//Eqn (5)
			double eqn5 = (x_nt_1_t_1/eqn4_1) - (sum_x_nt_1_t_1 / denom_2) ;

			double gradient_tk = eqn1 * (eqn2 - eqn4) + posterior_dist[k] * (eqn3 - eqn5);

			if(!single_A && !single_mu)
			{
				int idx = (t-1)*K + k;
				//printf("n:%d, t=%d, g=%f, old_mu=%f", n, t-1 + users[n]->t_first, gradient_tk, users[n]->mu_t[idx]);
				mu_t[idx] -= learning_rate * gradient_tk ;
				mu_t[idx] = std::min(std::max(0.0, mu_t[idx]), 1.0);
				//printf(", new_mu=%f\n", users[n]->mu_t[idx]);
			}

			gradient_t += gradient_tk;
		}

		if(!single_A && single_mu)
		{
			// printf("n:%d, t=%d, g=%f, old_mu=%f", n, t-1 + users[n]->t_first, gradient_t, users[n]->mu_t[t-1]);
			mu_t[t-1] -= learning_rate * gradient_t ;
			mu_t[t-1] = std::min(std::max(0.0, mu_t[t-1]), 1.0);
			// printf(", new_mu=%f\n", users[n]->mu_t[t-1]);
		}

		gradient += gradient_t;
	}

	if(single_A && single_mu)
	{
		mu_t[0] -= learning_rate * gradient;
		mu_t[0] = std::min(std::max(0.0, mu_t[0]), 1.0);
	}

	// L = 0;
	// for(int t=1;t<Tn;t++)
	// {
	// 	//obtain the prior
	// 	users[n]->get_prior(t, ALPHA, K, prior_dist);

	// 	//obtain the posterior
	// 	users[n]->get_posterior(t, ALPHA, K, posterior_dist);

	// 	//get the KL divergence
	// 	L += kl_divergence(posterior_dist, prior_dist, K);
	// }
	// printf("n:%d, After, L = %f\n", n, L);
}

User::~User()
{
	for(std::map<int, std::list<int> *>::iterator iter = adoption.begin(); iter != adoption.end(); iter++)
	{
		delete iter->second;
	}

	for(int t = 0;t < Tn; t++)
	{
		if(Z[t] != NULL)
			delete [] Z[t];
	}

	delete [] Z;

	if(Tn > 1)
		delete [] mu_t;

	delete [] phi;
	delete [] posterior_x;

	/*
	for(int t = Tn-1; t >= 0; t--)
	{
		gsl_vector_free(posterior_x[t]);
		gsl_vector_free(phi[t]);
		delete [] Z[t];
		delete [] adoption[t];
	}

	delete [] Mn;
	delete [] Z;
	delete [] adoption;
	*/
}
