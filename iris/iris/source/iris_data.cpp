#include "iris_data.hpp"
iris_data :: iris_data()
{
	initialize();
}

iris_data :: iris_data( const iris_data_params & params )
{
	initialize();
	setup ( params );
}

iris_data :: iris_data( const iris_data & data )
{
	initialize();
	setup( data );
}

void iris_data :: initialize()
{
	x_iris = 0;
	y_iris = 0;
	a_iris = 0;
	b_iris = 0;
	theta_iris = 0;

	new_x_iris = 0;
	new_y_iris = 0;
	new_r_iris = 0;


	//Transformée polaire normalisée
	polar_image = 0;
	polar_mask = 0;
	nb_samples = 0;
	nb_directions = 0;
	nb_samples_iris = 0;
	iris_width = 0;
	iris_height = 0;
	iris_image = 0;
	iris_mask = 0;
	nrj_ratio = 0;

	//Iris code
	nb_samples_code = 0;
	nb_directions_code = 0;
	code = 0;
	code_mask = 0;
	a_upper = 0;
	c_upper = 0;
	theta_upper = 0;
	
	a_lower = 0;
	c_lower = 0;
	theta_lower = 0;
	
	p_data.setup();
}

void iris_data :: free()
{
	if ( polar_image )
		cvReleaseImage( & polar_image );
	if ( polar_mask )
		cvReleaseImage( & polar_mask );
	if ( code )
		cvReleaseImage ( &code );
	if ( code_mask )
		cvReleaseImage ( &code_mask );
	if ( iris_image )
		cvReleaseImage ( &iris_image );
	if ( iris_mask )
		cvReleaseImage ( &iris_mask );
}

void iris_data :: alloc()
{

	polar_image = cvCreateImage ( 	cvSize( nb_directions, nb_samples),
									IPL_DEPTH_64F,
									1 );
	polar_mask = cvCreateImage ( 	cvSize( nb_directions, nb_samples),
									IPL_DEPTH_8U,
									1 );


	code =
		cvCreateImage ( cvSize( nb_directions_code, 2 * nb_samples_code),
						IPL_DEPTH_8U,
						1 );
	code_mask =
		cvCreateImage ( cvSize( nb_directions_code, nb_samples_code),
						IPL_DEPTH_8U,
						1 );

	iris_image = cvCreateImage ( 	cvSize( iris_width, iris_height),
									IPL_DEPTH_8U,
									1 );
	iris_mask = cvCreateImage ( 	cvSize( iris_width, iris_height),
									IPL_DEPTH_8U,
									1 );

}

int iris_data :: setup()
{
	free();
	initialize();

	return 0;
}

int iris_data :: setup ( const iris_data_params & params )
{

	setup();
	nb_directions = params.polar_width;
	nb_directions_code = params.iris_code_width;
	nb_samples =  params.polar_height;
	nb_samples_code = params.iris_code_height;
	nb_samples_iris = params.nb_samples_iris;
	iris_width = params.iris_width;
	iris_height = params.iris_height;
	if ( p_data.setup( cvSize( params.img_width, params.img_height) ) )
		return 1;

	alloc();
	return 0;
}

int iris_data :: setup ( const iris_data & data )
{
	setup();
	nb_directions = data.nb_directions;
	nb_directions_code = data.nb_directions_code;
	nb_samples =  data.nb_samples;
	nb_samples_code = data.nb_samples_code;
	nb_samples_iris = data.nb_samples_iris;
	iris_width = data.iris_width;
	iris_height = data.iris_height;

	if ( p_data.setup( data.p_data ) )
		return 1;
	alloc();

	nrj_ratio = data.nrj_ratio;
	p_data = data.p_data;
		cvSetImageROI ( iris_image,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.iris_image),
								GET_IMAGE_HEIGHT(data.iris_image) ) );
		cvCopyImage( data.iris_image, iris_image );


		cvSetImageROI ( iris_mask,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.iris_mask),
								GET_IMAGE_HEIGHT(data.iris_mask) ) );
		cvCopyImage( data.iris_mask, iris_mask );
		cvSetImageROI ( polar_image,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.polar_image),
								GET_IMAGE_HEIGHT(data.polar_image) ) );
		cvCopyImage( data.polar_image, polar_image );

		cvSetImageROI ( polar_mask,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.polar_mask),
								GET_IMAGE_HEIGHT(data.polar_mask) ) );
		cvCopyImage( data.polar_mask, polar_mask );
		cvSetImageROI ( code,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.code),
								GET_IMAGE_HEIGHT(data.code) ) );
		cvCopyImage( data.code, code );

		cvSetImageROI ( code_mask,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.code_mask),
								GET_IMAGE_HEIGHT(data.code_mask) ) );
		cvCopyImage( data.code_mask, code_mask );
		x_iris = data.x_iris;
		y_iris = data.y_iris;
		a_iris = data.a_iris;
		b_iris = data.b_iris;
		theta_iris = data.theta_iris;

		new_x_iris = data.new_x_iris;
		new_y_iris = data.new_y_iris;
		new_r_iris = data.new_r_iris;
		
		a_lower = data.a_lower;
		c_lower = data.c_lower;
		theta_lower = data.theta_lower;
		
		a_upper = data.a_upper;
		c_upper = data.c_upper;
		theta_upper = data.theta_upper;	
		
		seg_ok = data.seg_ok;
		
		


	return 0;
}

iris_data & iris_data :: operator=( const iris_data & data )
{
	if ( !(  data == *this ) )
	{

		setup( data );
	}
	else
	{
		nrj_ratio = data.nrj_ratio;
		p_data = data.p_data;
		cvSetImageROI ( iris_image,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.iris_image),
								GET_IMAGE_HEIGHT(data.iris_image) ) );
		cvCopyImage( data.iris_image, iris_image );


		cvSetImageROI ( iris_mask,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.iris_mask),
								GET_IMAGE_HEIGHT(data.iris_mask) ) );
		cvCopyImage( data.iris_mask, iris_mask );
		cvSetImageROI ( polar_image,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.polar_image),
								GET_IMAGE_HEIGHT(data.polar_image) ) );
		cvCopyImage( data.polar_image, polar_image );

		cvSetImageROI ( polar_mask,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.polar_mask),
								GET_IMAGE_HEIGHT(data.polar_mask) ) );
		cvCopyImage( data.polar_mask, polar_mask );
		cvSetImageROI ( code,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.code),
								GET_IMAGE_HEIGHT(data.code) ) );
		cvCopyImage( data.code, code );

		cvSetImageROI ( code_mask,
						cvRect( 0,
								0,
								GET_IMAGE_WIDTH(data.code_mask),
								GET_IMAGE_HEIGHT(data.code_mask) ) );
		cvCopyImage( data.code_mask, code_mask );
		x_iris = data.x_iris;
		y_iris = data.y_iris;
		a_iris = data.a_iris;
		b_iris = data.b_iris;
		theta_iris = data.theta_iris;

		new_x_iris = data.new_x_iris;
		new_y_iris = data.new_y_iris;
		new_r_iris = data.new_r_iris;
		
		a_lower = data.a_lower;
		c_lower = data.c_lower;
		theta_lower = data.theta_lower;
		
		a_upper = data.a_upper;
		c_upper = data.c_upper;
		theta_upper = data.theta_upper;
		
		seg_ok = data.seg_ok;
	}
	return *this;
}

iris_data :: ~iris_data()
{
	free();
	initialize();
}

bool iris_data :: operator==( const iris_data & data ) const
{
	if ( 	data.nb_directions == nb_directions 			&&
			data.nb_samples == nb_samples					&&
			data.nb_directions_code == nb_directions_code	&&
			data.nb_samples == nb_samples					&&
			p_data == data.p_data							&&
			iris_width == data.iris_width					&&
			iris_height == data.iris_height					)
		return 1;
	else
		return 0;
}

int iris_data :: save ( const char * rep, unsigned int id ) const
{
	int q = 0;
	stringstream oss;

	//Sauvegarde des données de la pupille
	if ( p_data.save(rep, id) )
		return 1;
	if ( seg_ok == false )
		return 1;


	//Données de seg.
	{
		oss << rep << "/" << "iris_data_" << id << ".txt";
		ofstream file(oss.str().c_str());
		oss.str("");

		file << "#Iris seg. data" << endl;
		file << "x_iris = " << x_iris << endl;
		file << "y_iris = " << y_iris << endl;
		file << "a_iris = " << a_iris << endl;
		file << "b_iris = " << b_iris << endl;
		file << "theta_iris = " << theta_iris << endl;
		file << "new_x_iris = " << new_x_iris << endl;
		file << "new_y_iris = " << new_y_iris << endl;
		file << "new_r_iris = " << new_r_iris << endl;
		
		file << "a_lower = " << a_lower << endl;
		file << "c_lower = " << c_lower << endl;
		file << "theta_lower = " << theta_lower << endl;		
		
		file << "a_upper = " << a_upper << endl;
		file << "c_upper = " << c_upper << endl;
		file << "theta_upper = " << theta_upper << endl;	
		
		file << "nrj_ratio = " << nrj_ratio << endl;
		file.close();
	}

	//Sauvegardes des images
	{
		oss << rep << "/" << "resized_image_"<< id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  iris_image );
		oss.str("");
	}

	{
		oss << rep << "/" << "resized_mask_"<< id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  iris_mask );
		oss.str("");
	}

	{
		oss << rep << "/" << "polar_image_"<< id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  polar_image );
		oss.str("");
	}

	{
		oss << rep << "/" << "polar_mask_"<< id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  polar_mask );
		oss.str("");
	}

	{
		oss << rep << "/" << "iris_code_"<< id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  code );
		oss.str("");
	}

	{
		oss << rep << "/" << "iris_code_mask_"<< id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  code_mask );
		oss.str("");
	}

	return q;
}

int iris_data :: load ( const char * rep, unsigned int id )
{
	free();
	initialize();

	int q = 0;
	stringstream oss;

	if ( p_data.load( rep, id ) )
		q = 0;

	//Données de seg.
	{
		oss << rep << "/" << "iris_data_" << id << ".txt";
			api_parameters params;
			if ( params.load( oss.str().c_str() ) )
				q = 1;
			else
			{
				if ( api_get_double( params,
									  "new_x_iris",
									  &new_x_iris ) )
					q = 1;
				if ( api_get_double( params,
									  "new_y_iris",
									  &new_y_iris ) )
					q = 1;
				if ( api_get_double( params,
									  "new_r_iris",
									  &new_r_iris ) )
					q = 1;
				if ( api_get_double ( params,
									  "x_iris",
									  &x_iris ) )
					q = 1;
				if ( api_get_double ( params,
									  "y_iris",
									  &y_iris ) )
					q = 1;
				if ( api_get_double ( params,
									  "a_iris",
									  &a_iris ) )
					q = 1;
				if ( api_get_double ( params,
									  "b_iris",
									  &b_iris ) )
					q = 1;
				if ( api_get_double ( params,
									  "theta_iris",
									  &theta_iris ) )
					q = 1;
				if ( api_get_double ( params,
									  "nrj_ratio",
									  &nrj_ratio ) )
					q = 1;
					
				if ( api_get_double ( params,
									  "a_upper",
									  &a_upper ) )
					q = 1;
				if ( api_get_double ( params,
									  "c_upper",
									  &c_upper ) )
					q = 1;
				if ( api_get_double ( params,
									  "theta_upper",
									  &theta_upper ) )
					q = 1;		
					
				if ( api_get_double ( params,
									  "a_lower",
									  &a_lower ) )
					q = 1;
				if ( api_get_double ( params,
									  "c_lower",
									  &c_lower ) )
					q = 1;
				if ( api_get_double ( params,
									  "theta_lower",
									  &theta_lower ) )
					q = 1;						

			}



		oss.str("");
	}


	{
		oss << rep << "/" << "resized_image_" << id << ".png";
		iris_image = cvLoadImage ( oss.str().c_str(), 0 );
		oss.str("");
	}
	if ( !iris_image )
		q = 0;

	{
		oss << rep << "/" << "resized_mask_" << id << ".png";
		iris_mask = cvLoadImage ( oss.str().c_str(), 0 );
		oss.str("");
	}
	if ( !iris_mask )
		q = 0;

	{
		oss << rep << "/" << "polar_image_" << id << ".png";
		polar_image = cvLoadImage ( oss.str().c_str(),
									0);
		oss.str("");
	}
	if ( ! polar_image )
		q = 0;

	{
		oss << rep << "/" << "polar_mask_" << id << ".png";
		polar_mask = cvLoadImage ( oss.str().c_str(),
									0);
		oss.str("");
	}
	if ( ! polar_mask )
		q = 0;
	
	q = 0;
	{
		oss << rep << "/" << "iris_code_" << id << ".png";
		code = cvLoadImage ( oss.str().c_str(),
							 0);
		oss.str("");
	}
	if (!code)
		q = 1;

	{
		oss << rep << "/" << "iris_code_mask_" << id << ".png";
		code_mask = cvLoadImage (	 oss.str().c_str(),
									0);
		oss.str("");
	}
	if ( ! code_mask )
		q = 1;

	if (!q)
	{
		//~ iris_width = iris_image->width;
		//~ iris_height = iris_image->height;
		//~ nb_samples = polar_image->height;
		//~ nb_directions = polar_image->width;
		//~ nb_samples_iris = polar_image->height;
		nb_samples_code = code->height / 2;
		nb_directions_code = code->width;
		seg_ok = true;

	}
	return q;
}

void * iris_data_alloc( const void * params )
{
	iris_data_params * data = (iris_data_params *) params;
	return (void*) new iris_data(*data);
}

void iris_data_copy( 	void * data_tg,
						const void * data_src )
{
	iris_data 	* in = (iris_data *) data_src,
				* out = (iris_data *) data_tg;

	(*out) = (*in);
}

void iris_data_free( void * p )
{
	iris_data 	* toto = (iris_data *) p;
	delete toto;
}

