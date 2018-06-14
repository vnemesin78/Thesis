#include "c_get_iris_template.hpp"
#include "iris_default.hpp"
int c_get_iris_template :: default_setup ( ostream * stream )
{
	int q =0;
	
	if ( image_thread->default_setup(stream  ) )
		q = 1;
	if ( pupil_thread->default_setup(stream ) )
		q = 1;
	if ( iris_thread->default_setup(stream ) )
		q = 1;
	return q;
	
}
