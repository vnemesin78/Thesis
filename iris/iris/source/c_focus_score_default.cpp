#include "c_focus_score.hpp"
#include "iris_default.hpp"
int c_focus_score :: default_setup(	)
{
	int q = 0;
	gsl_matrix * kernel;
	SCORE__KERNEL(kernel);
	
	double c;
	SCORE__C(c);
	
	unsigned int width;
	SCORE__WIDTH(width);
	
	unsigned int height;
	SCORE__HEIGHT(height);
	
	
	
	
	q = setup( width, height, kernel, c );
	
	
	
	gsl_matrix_free( kernel );
	
	return q;
}

