#include <map>
#include <list>
#include <iostream>
#include <math.h>

#ifndef USER
#define USER

class User
{
public:
	int t_first, t_last, Tn;
	double * mu_t;
	bool single_A, single_mu;

	std::map<int, std::list<int> *> adoption;

	int ** Z;
	int * phi;
	double * posterior_x;

	User(bool, bool);
	~User();

	void init_mu(const int, const double);
	void add(const int, const int);

	void estimate_dynamics(const int, const double, const double);
	
	void get_prior_x(const int, const int, double *);
	void get_prior(const int, const double, const int, double *);
	void get_posterior(const int, const double, const int, double *);
	double get_decay(const int);
	double * get_decay_vector(const int, const int);
};

template<typename T>
T sum_vector(const int K, const T * vec, const double * mu_t)
{
	T sum = 0;
	for(int k = 0; k < K; k++)
	{
		sum += mu_t[k] * vec[k];
	}
	return sum;
}

template<typename T>
T sum_vector(const int K, const T * vec)
{
	T sum = 0;
	for(int k = 0; k < K; k++)
	{	
		sum += vec[k];
	}
	return sum;
}

template<typename T>
void set_zero(const int K, T * out)
{
	for(int k = 0; k < K; k++)
	{
		out[k] = 0;
	}
}
#endif
