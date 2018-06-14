#include "c_iris_thread.hpp"
#include "iris_default.hpp"
int c_iris_thread :: default_setup ( ostream * _err_stream )
{
	int q = 0;
	free();
	initialize();
	unsigned int d_width, d_height;
	err_stream = _err_stream;
	WIDTH(d_width);
	HEIGHT(d_height);
	SCORE__TYPE(_score_type);
	IRIS__VALIDNESS_FACTOR(_validness_factor);
	IRIS__SPOT_THRESHOLD(_spot_threshold);
	iris_segmentation = new c_iris_segmentation;
	polar_iris = new c_polar_iris;
	eyelid_segmentation = new c_eyelids_segmentation;
	iris_code = new c_iris_code;
	i_data = new iris_data;

	if ( iris_segmentation->default_setup(  d_width, d_height ) )
		q = 1;

	if ( polar_iris->default_setup( d_width, d_height ) )
		q = 1;
	if ( eyelid_segmentation->default_setup(	d_width * (1 + 2 * iris_segmentation->r_ratio() ),
												d_height * (1 + 2 * iris_segmentation->r_ratio() ) ) )
		q = 1;

	if ( iris_code->default_setup())
		q = 1;
	if ( _score_type == 0 )
	{
		focus_score = new c_focus_score;
		if ( focus_score->default_setup() )
			q = 1;
	} 

	iris_data_params args;
		args.img_width = d_width;
		args.img_height = d_height;
		args.polar_width = polar_iris->nb_directions();
		args.polar_height = polar_iris->nb_samples();
		args.nb_samples_iris = polar_iris->nb_samples_iris();
		args.iris_code_width = iris_code->nb_directions();
		args.iris_code_height = iris_code->nb_samples();
		args.iris_width = eyelid_segmentation->width();
		args.iris_height = eyelid_segmentation->height();

	if (i_data->setup( args ))
		q = 1;

	tmp_p_img = cvCreateImage( 	cvSize ( 	iris_code->nb_directions(),
											iris_code->nb_samples() ),
								IPL_DEPTH_64F,
								1 );

	tmp_p_mask = cvCreateImage( cvSize ( 	iris_code->nb_directions(),
											iris_code->nb_samples() ),
								IPL_DEPTH_8U,
								1 );



	if ( q )
		return 1;
	return 0;
}
