#include "c_pupil_tracking.hpp"
#include "iris_default.hpp"
int c_pupil_tracking :: default_setup ( unsigned int width, 
											unsigned int height,
											ostream * _err_stream)
{
	int q = 0;
	gsl_matrix * tmp;
	TRACKING__T_0(tmp);
	gsl_vector t_0;
		t_0.size = tmp->size1 * tmp->size2;
		t_0.data = tmp->data;
		t_0.stride = 1;
	
	
	gsl_matrix * sqrt_q0;
	TRACKING__SQRT_Q_0(sqrt_q0);
	
	
	gsl_matrix * f;
	TRACKING__F(f);
	
	gsl_matrix * sqrt_q;
	TRACKING__SQRT_Q(sqrt_q);
	
	double n;
	TRACKING__N_SIGMA( n );
	
	q = setup (	width,
				height,
				&t_0,
				sqrt_q0,
				f,
				sqrt_q,
				n,
				_err_stream );
	
	gsl_matrix_free(tmp);
	gsl_matrix_free(sqrt_q0);
	gsl_matrix_free(f);
	gsl_matrix_free(sqrt_q);
	
	return q;
}
