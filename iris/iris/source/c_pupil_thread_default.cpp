#include "c_pupil_thread.hpp"
#include "iris_default.hpp"
int c_pupil_thread :: default_setup ( ostream * _err_stream )
{
	int q = 0;
	setup();
	unsigned int d_width, d_height;
	WIDTH(d_width);
	HEIGHT(d_height);
	
	preprocessing = new c_preprocessing;
	pupil_tracking = new c_pupil_tracking;
	pupil_segmentation = new c_pupil_segmentation;
	p_data = new pupil_data;	
	
	if ( preprocessing->default_setup( 	d_width, 
										d_height, 
										_err_stream ) )
		q = 1;
	if ( pupil_tracking->default_setup(	d_width, 
										d_height, 
										_err_stream ) )
		q = 1;
	if ( pupil_segmentation->default_setup( d_width, 
											d_height, 
											_err_stream ) )
		q = 1;
	
	p_data->setup( cvSize ( d_width, 
							d_height ) );
		
	if ( q )
		return 1;
	pupil_tracking->reset_tracking();
	return 0;
}
