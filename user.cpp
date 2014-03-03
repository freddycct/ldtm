#include "user.h"
#include "linear_regression.h"

double User::calculate_gc(User * tar, const int tau, const int width, const int lookahead,
	double * tar_tseries, double * src_tseries)
{
	int tar_start = tau - width ;
	int tar_end   = tau + lookahead ;

	int src_start = tar_start ;
	int src_end   = tar_end - 1 ;

	if (tau + lookahead < tar->t_first || tau < t_first)
    	return nan("");
	else if (src_end > t_last)
		return nan("");
	else if(tar_end > tar->t_last)
    	return nan("");

	int midpoint = tar_start > tar->t_first ? tar_start : tar->t_first;
	//double * tar_matrix = new double[(width+lookahead+1)*K];
	//x = [zeros(K, midpoint - x_start) + 1/K, vectors{target_id}{midpoint:x_end}] ;

	for(int t = tar_start; t < midpoint; t++)
	{
		int adjusted_t = t - tar_start;
		for(int k = 0; k < tar->K; k++)
		{
			tar_tseries[adjusted_t * tar->K + k] = 1.0 / (tar->K);
		}
	}

	for(int t = midpoint; t <= tar_end; t++)
	{
		int adjusted_t1 = t - tar_start;
		int adjusted_t2 = t - tar->t_first;
		for(int k = 0; k < tar->K; k++)
		{
			tar_tseries[adjusted_t1 * tar->K + k] = tar->posterior_x[adjusted_t2 * tar->K + k];
		}
	}
	
	midpoint = src_start > t_first ? src_start : t_first;
	//y = [zeros(K, midpoint - y_start) + 1/K, vectors{source_id}{midpoint:y_end}] ;

	for(int t = src_start; t < midpoint; t++)
	{
		int adjusted_t = t - src_start;
		for(int k = 0; k < K; k++)
		{
			src_tseries[adjusted_t * K + k] = 1.0 / K;
		}
	}

	for(int t = midpoint; t <= src_end; t++)
	{
		int adjusted_t1 = t - src_start;
		int adjusted_t2 = t - t_first;
		for(int k = 0; k < K; k++)
		{
			src_tseries[adjusted_t1 * K + k] = posterior_x[adjusted_t2 * K + k];
		}
	}
	
	return gc(src_tseries, K, width+lookahead, tar_tseries, tar->K, width+lookahead+1, width);
}

void User::calculate_gc(const int width, const int lookahead)
{
	if(Tn < 1) return;

	double * tar_tseries = new double[(width+lookahead+1) * K];
	double * src_tseries = new double[(width+lookahead+1) * K];
	for(std::map<int, std::map<User *, int> *>::iterator iter = neighbors.begin(); iter != neighbors.end(); iter++)
	{
		int t = iter->first;
		std::map<User *, int> * neighbors_t = iter->second;
		std::map<User *, double> * neighbors_gc_t = neighbors_gc[t];
		
		for(std::map<User *, int>::iterator iter2 = neighbors_t->begin(); iter2 != neighbors_t->end(); iter2++)
		{
			User * tar = iter2->first;
			(*neighbors_gc_t)[tar] = calculate_gc(tar, t, width, lookahead, src_tseries, tar_tseries);
		}
	}
	delete [] src_tseries;
	delete [] tar_tseries;
}

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
	posterior_x = new double[Tn * K];

	for(int i = 0; i < Tn * K; i++)
	{
		phi[i] = 0;
	}

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
			for(int t = 0; t < Tn - 1; t++)
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
			mu_t = new double[ (Tn - 1) * K ] ; //the arrangement should be K x (Tn-1)
			for(int t = 0; t < Tn - 1; t++)
			{
				for(int k = 0; k < K; k++)
				{
					mu_t[t*K + k] = mu;
				}
			}
		}
	} //if(Tn > 1)
}

User::User(std::string id)
{
	this->id = id;
	t_first = -1;
	t_last = -1;
	Tn = 0;
}

User::User(std::string id, const int t_first, const int t_last, double * posterior_x, const int K)
{
	this->id = id;
	this->t_first = t_first;
	this->t_last = t_last;
	Tn = t_last - t_first + 1;
	this->posterior_x = posterior_x;
	this->K = K;
}

User::User(std::string id, const bool single_A, const bool single_mu)
{
	this->id = id;
	t_first = -1;
	t_last = -1;
	Tn = 0;
	this->single_A  = single_A;
	this->single_mu = single_mu;
}

void User::add_social(const int t, User * tar, const int freq)
{
	std::map<User *, int> * neighbors_t;
	std::map<int, std::map<User *, int> *>::iterator iter;
	iter = neighbors.find(t);

	if(iter == neighbors.end())
	{
		neighbors_t = new std::map<User *, int>();
		neighbors[t] = neighbors_t;
		neighbors_gc[t] = new std::map<User *, double>();
	}
	else
	{
		neighbors_t = iter->second;
	}

	std::map<User *, int>::iterator iter2;
	iter2 = neighbors_t->find(tar);

	if(iter2 == neighbors_t->end())
	{
		(*neighbors_t)[tar] = freq;
	}
	else
	{
		iter2->second += freq;
	}
}

void User::add(const int t, const int m, const int freq)
{
	//int adjusted_index = t - t_first;
	//std::cout << adjusted_index << std::endl;

	std::list<int> * adoption_t;
	std::map<int, std::list<int> *>::iterator iter;
	iter = adoption.find(t);

	if(iter == adoption.end())
	{
		adoption_t = new std::list<int>();
		adoption[t] = adoption_t;

		if(t > t_last)
		{
			t_last = t;
		}

		if(t < t_first || t_first < 0)
		{
			t_first = t;
		}
	}
	else
	{
		adoption_t = iter->second;
	}

	if(m < 1)
	{
		//throw an error
		std::cout << "error occurred here, m must be greater than or equal to one." << std::endl;
	}
	
	for(int n=0;n<freq;n++)
	{
		adoption_t->push_back(m - 1);
	}
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

void User::estimate_dynamics(const int K, const double ALPHA, 
	const double learning_rate, double * posterior_dist)
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
