#include "tkalman_simulation.hpp"

void do_tkalman_simulation(	gsl_vector ** x,
							gsl_vector ** y,
							const gsl_vector * t0,
							const gsl_matrix * sqrt_q0,
							const gsl_matrix * f,
							const gsl_matrix * sqrt_q,
							const unsigned int n,
							gsl_vector * vect_t_1,
							gsl_vector * vect_t_1_view_x,
							gsl_vector * vect_t_1_view_y,
							gsl_vector * vect_t_2,
							gsl_rng * r)
{
	unsigned int size_y, 
				 size_x,
				 size_t;
	size_x = x[0]->size;
	size_y = f->size1 - size_x;
	size_t = size_x + size_y;
	//Tirage de t0
	for ( unsigned int i = 0; i < size_t; ++ i )
	{
		vect_t_2->data[i * vect_t_2->stride] = gsl_ran_ugaussian (r);
	}
	gsl_vector_memcpy(vect_t_1, t0);
	gsl_blas_dgemv(CblasTrans,
				   1.0,
				   sqrt_q0,
				   vect_t_2,
				   1.0,
				   vect_t_1);
	 		   
	gsl_vector_memcpy(x[0],
					  vect_t_1_view_x);
	//Tirage des tn
	for ( unsigned int i = 1; i <= n; ++ i)
	{
		//Calcul de F T
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   f,
					   vect_t_1,
					   0.0,
					   vect_t_2);
		
		//Bruit
		for (unsigned int j = 0; j < size_t; ++ j)
		{
			vect_t_1->data[j * vect_t_1->stride] = gsl_ran_ugaussian (r);	
		}
		gsl_blas_dgemv(CblasTrans,
					   1.0,
					   sqrt_q,
					   vect_t_1,
					   1.0,
					   vect_t_2);
		gsl_vector_memcpy(vect_t_1, vect_t_2);
		gsl_vector_memcpy(x[i], vect_t_1_view_x);
		gsl_vector_memcpy(y[i - 1],  vect_t_1_view_y);
	}
}

void do_no_tkalman_simulation(	gsl_vector ** y,
								const gsl_vector * const * x,
								const gsl_matrix * sqrt_q_yy,
								const unsigned int n,
								gsl_vector * vect,
						        gsl_rng * r)
{
	unsigned int size;
	size = x[0]->size;
	//Tirage des tn
	for ( unsigned int i = 0; i < n; ++ i)
	{
		gsl_vector_memcpy(y[i], x[i]);
		//Bruit
		for (unsigned int j = 0; j < size; ++ j)
		{
			vect->data[j * vect->stride] = gsl_ran_ugaussian (r);	
		}
		gsl_blas_dgemv(CblasTrans,
					   1.0,
					   sqrt_q_yy,
					   vect,
					   1.0,
					   y[i]);
	}
}

void do_tkalman_simulation_markov(	gsl_vector ** x,
									gsl_vector ** y,
									const gsl_vector * t0,
									const gsl_matrix * sqrt_q0,
									const gsl_matrix * const * f,
									const gsl_matrix * const * sqrt_q,
									const unsigned int * process,
									const unsigned int n,
									gsl_vector * vect_t_1,
									gsl_vector * vect_t_1_view_x,
									gsl_vector * vect_t_1_view_y,
									gsl_vector * vect_t_2,
									gsl_rng * r)
{
	unsigned int size_y, 
				 size_x,
				 size_t;
	size_x = x[0]->size;
	size_y = f[0]->size1 - size_x;
	size_t = size_x + size_y;
	//Tirage de t0
	for ( unsigned int i = 0; i < size_t; ++ i )
	{
		vect_t_2->data[i * vect_t_2->stride] = gsl_ran_ugaussian (r);
	}
	gsl_vector_memcpy(vect_t_1, t0);
	gsl_blas_dgemv(CblasTrans,
				   1.0,
				   sqrt_q0,
				   vect_t_2,
				   1.0,
				   vect_t_1);
	 		   
	gsl_vector_memcpy(x[0],
					  vect_t_1_view_x);
	//Tirage des tn
	for ( unsigned int i = 1; i <= n; ++ i)
	{
		//Calcul de F T
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   f[process[i - 1]],
					   vect_t_1,
					   0.0,
					   vect_t_2);
		
		//Bruit
		for (unsigned int j = 0; j < size_t; ++ j)
		{
			vect_t_1->data[j * vect_t_1->stride] = gsl_ran_ugaussian (r);	
		}
		gsl_blas_dgemv(CblasTrans,
					   1.0,
					   sqrt_q[process[i - 1]],
					   vect_t_1,
					   1.0,
					   vect_t_2);
		gsl_vector_memcpy(vect_t_1, vect_t_2);
		gsl_vector_memcpy(x[i], vect_t_1_view_x);
		gsl_vector_memcpy(y[i - 1],  vect_t_1_view_y);
	}
	
	
	
}
