#include <map>
#include <list>
#include <iostream>
#include <math.h>

#ifndef USER
#define USER

class User
{
public:
	std::string id;
	int t_first, t_last, Tn;
	double * mu_t;
	bool single_A, single_mu;

	std::map<int, std::list<int> *> adoption;
	std::map<int, std::map<User *, int> *> neighbors;
	std::map<int, std::map<User *, double> *> neighbors_gc;

	int K;
	int ** Z;
	int * phi;
	double * posterior_x;

	User(std::string, const int, const int, double *, const int);
	User(std::string, const bool, const bool);
	User(std::string);
	~User();

	void init_mu(const int, const double);
	void add(const int, const int, const int);
	void add_social(const int, User *, const int);

	void estimate_dynamics(const int, const double, const double, double *);
	
	void get_prior_x(const int, const int, double *);
	void get_prior(const int, const double, const int, double *);
	void get_posterior(const int, const double, const int, double *);
	double get_decay(const int);
	double * get_decay_vector(const int, const int);

	void calculate_gc(const int, const int);
	double calculate_gc(User *, const int, const int, const int, double *, double *);
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
