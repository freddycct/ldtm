#ifndef LINEAR_REGRESSION
#define LINEAR_REGRESSION

double dot_product(double * vec1, double * vec2, int K)
{
	double sum = 0;
	for(int k=0;k<K;k++)
	{
		sum += vec1[k] * vec2[k];
	}
	return sum;
}

double sum_array(double * array, int N)
{
	double sum = 0;
	for(int n=0;n<N;n++)
	{
		sum += array[n];
	}
	return sum;
}

void reset_array(double * array, int N)
{
	for(int n=0;n<N;n++)
	{
		array[n] = 0.0;
	}
}

double gc(double * src_tseries, const int src_K, const int src_len,
	double * tar_tseries, const int tar_K, const int tar_len,
	const int W)
{
	int iterations = 100;

	double gamma, sse;
	
	double alpha_0 = 0.0;
	
	double * alpha = new double[W];
	double * beta  = new double[W];
	for(int i = 0; i < W; i++)
	{
		alpha[i] = 0.0;
		beta[i] = 0.0;
	}

	//first linear regression problem
	double * tmp_vector = new double[tar_K];
	for(int iter = 0; iter < iterations; iter++)
	{
		//estimate alpha_0
		double tmp = 0;
		double tmp1 = 0;
		for(int t = W; t < tar_len; t++)
		{
			reset_array(tmp_vector, tar_K);
			for(int k = 0; k < tar_K; k++)
			{
				tmp_vector[k] = tar_tseries[t * tar_K + k];
			}
			for(int i = 0; i < W; i++)
			{
				for(int k = 0; k < tar_K; k++)
				{
					tmp_vector[k] -= alpha[i] * tar_tseries[(t - 1 - i) * tar_K + k];
				}
			}
			tmp += sum_array(tmp_vector, tar_K);
		}
		alpha_0 = tmp / ((tar_len - W) * tar_K);

		//estimate alpha_n
		for(int n = 0; n < W; n++)
		{
			tmp  = 0;
			for(int t = W; t < tar_len; t++)
			{
				reset_array(tmp_vector, tar_K);
				for(int k = 0; k < tar_K; k++)
				{
					tmp_vector[k] = tar_tseries[t * tar_K + k] - alpha_0;
				}

				for(int i = 0; i < W; i++)
				{
					if(n != i)
					{
						for(int k = 0; k < tar_K; k++)
						{
							tmp_vector[k] -= alpha[i] * tar_tseries[(t - 1 - i) * tar_K + k];
						}
					}
				}
				tmp += dot_product( &tar_tseries[(t - 1 - n) * tar_K], tmp_vector, tar_K );
			}

			//sum the denominator
			tmp1 = 0;
			for(int t = W; t < tar_len; t++)
			{
				tmp1 += dot_product(&tar_tseries[(t - 1 - n) * tar_K], 
					&tar_tseries[(t - 1 - n) * tar_K], tar_K);
			}
			alpha[n] = tmp / tmp1;
		}
	}

	//second linear regression
	double * constants = new double[ (tar_len - W) * tar_K ];
	//constants has to be a matrix
	for(int t = W; t < tar_len; t++)
	{
		//t is a index for the column in x
		//t - W is a index normalized to start from 0 for constants
		for(int k = 0; k < tar_K; k++)
		{
			constants[(t - W) * tar_K + k] = tar_tseries[t * tar_K + k] - alpha_0;
		}

		for(int i = 0; i < W; i++)
		{
			for(int k = 0; k < tar_K; k++)
			{
				constants[(t - W) * tar_K + k] -= alpha[i] * tar_tseries[(t - 1 - i) * tar_K + k];
			}
		}
	}
	
	//double * tmp = new double[K];
	
	for(int iter = 0; iter < iterations; iter++)
	{
		double tmp, tmp1;
		
		//estimate beta_n
		for(int n = 0; n < W; n++)
		{
			tmp = 0;
			for(int t = W; t < tar_len; t++)
			{
				for(int k = 0; k < tar_K; k++)
				{
					tmp_vector[k] = constants[(t - W) * tar_K + k];
				}
				for(int i = 0; i < W; i++)
				{
					if(n != i)
					{
						for(int k = 0; k < src_K; k++)
						{
							tmp_vector[k] -= beta[i] * src_tseries[(t - 1 - i) * src_K + k];
						}
					}
				}
				tmp += dot_product(&src_tseries[(t - 1 - n) * src_K], tmp_vector, src_K);
			}

				//sum the denominator
			tmp1 = 0;
			for(int t = W; t < tar_len; t++)
			{
				tmp1 += dot_product(&src_tseries[(t - 1 - n) * src_K], 
					&src_tseries[(t - 1 - n) * src_K], src_K);
			}
			beta[n] = tmp / tmp1;
		}
	}

	int N = tar_len + src_len;
	
	//for computing the Sum-of-Squared errors (SSE)
	double R2 = 0;
	double R1 = 0;
	for(int t = W; t < tar_len; t++)
	{
		for(int k = 0; k < tar_K; k++)
		{
			tmp_vector[k] = tar_tseries[t * tar_K + k] - alpha_0;
		}
		for(int i = 0; i < W; i++)
		{
			for(int k = 0; k < tar_K; k++)
			{
				tmp_vector[k] -= alpha[i] * tar_tseries[(t - 1 - i) * tar_K + k];
			}
		}

		R2 += dot_product(tmp_vector, tmp_vector, tar_K);
	
		for(int i = 0; i < W; i++)
		{
			for(int k = 0; k < src_K; k++)
			{
				tmp_vector[k] -= beta[i] * src_tseries[(t - 1 - i) * src_K + k];
			}
		}
		
		R1 += dot_product(tmp_vector, tmp_vector, src_K);
	}

	delete [] constants;
	delete [] tmp_vector;
	delete [] beta;
	delete [] alpha;

	return ((R2 - R1) / R1) * (N - 2 * W - 1) / W ;
}

#endif
