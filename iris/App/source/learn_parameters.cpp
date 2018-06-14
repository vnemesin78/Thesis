#include "lib_iris.hpp"





int main( int argc, char ** argv )
{
	
	if ( argc == 1 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}
	
	if ( 	! strcmp( argv[1], "--help" )	||
			! strcmp( argv[1], "-H" ) )
	{
		cout << "HELP" << endl;
		cout << " Arg[1] : Video name " << endl;
		cout << " Arg[2 - n] : Parameter files" << endl;
		cout << "Author: Valérian Némesin." << endl;
		return 0; 
	}	
	
	if ( argc < 2 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}
	
	//Chargement de la vidéo
	CvCapture * video;
	
	video = cvCaptureFromAVI( argv[1] );
	if ( ! video )
	{
		cout << "Error: can't load " << argv[1] << endl;
		return 1;
	}
	
	//Création des objets de seg.
	c_image_thread image_thread;
	c_pupil_thread pupil_thread;
	
	if ( image_thread.setup( argv + 1, argc - 1, NULL ) )
		return 1;
		
	if ( pupil_thread.setup(argv + 1, argc - 1, NULL ) )
		return 1;	
	
	cout << argv[1] << endl;
	int nb_frames = cvGetCaptureProperty( video , CV_CAP_PROP_FRAME_COUNT );
	int width = cvGetCaptureProperty( video , CV_CAP_PROP_FRAME_WIDTH );
	int height =  cvGetCaptureProperty( video , CV_CAP_PROP_FRAME_HEIGHT );
	cout << "Nb Frames : " << nb_frames << endl;
	
	//Allocation
	unsigned int nb_sig = 0;
	gsl_vector *** sig = new gsl_vector**[nb_frames];
	unsigned int * sizes = new unsigned int[nb_frames];
	sig[0] = new gsl_vector*[nb_frames];
	sizes[0] = 0;
	
	
	//Parcours de toutes les images et segmentation
	IplImage * image;
	while ( ( image = cvQueryFrame(video) ) != NULL )
	{
		int q = 0;
		image_thread.process( image, argv[1] );
		if ( image_thread.data().img_ok )
		{
			pupil_thread.segment_pupil( image_thread.data() );
			if ( pupil_thread.pupil_seg_data().seg_ok )
				q = 1;
		}
		
		
		if ( q )
		{
			sig[nb_sig][sizes[nb_sig]] = gsl_vector_alloc(3);
			sig[nb_sig][sizes[nb_sig]]->data[0] = pupil_thread.pupil_seg_data().x_pupil - width;
			sig[nb_sig][sizes[nb_sig]]->data[1] = pupil_thread.pupil_seg_data().y_pupil - height;
			sig[nb_sig][sizes[nb_sig]]->data[2] = 	sqrt( pow(pupil_thread.pupil_seg_data().a_pupil,2) + pow(pupil_thread.pupil_seg_data().b_pupil,2) );
			sizes[nb_sig] ++;
		}
		else
		{
			if ( sizes[nb_sig] > 0 )
			{
				nb_sig ++;
				sig[nb_sig] = new gsl_vector*[nb_frames];
				sizes[nb_sig] = 0;	
			}
		}
	}
	
	if ( sizes[nb_sig] != 0 )
	{
		nb_sig ++;
	}
	
	cout << "Nb signals : " << nb_sig << endl;
	
	api_parameters params;
	for ( int i = 2; i < argc; ++ i )
		params.load(argv[i] );
		
	gsl_vector t0;
	gsl_matrix sqrt_q0, f, sqrt_q;
	unsigned int nb_iter;
	
	int q = 0;
	if ( api_get_vector( params, "EM::t_0", &t0, &cout ) )
		q = 1;
	if ( api_get_matrix( params, "EM::sqrt_Q_0", &sqrt_q0, &cout ) )
		q = 1;	
	if ( api_get_matrix( params, "EM::sqrt_Q", &sqrt_q, &cout ) )
		q = 1;	
	if ( api_get_matrix( params, "EM::F", &f, &cout ) )
		q = 1;
	if ( api_get_positive_integer( params, "EM::nb_iterations", & nb_iter, &cout ) )
		q = 1;
	if ( q )
		return 1;
	
	tkalman_nc_em em_obj(	& t0,
							&sqrt_q0,
							&f,
							&sqrt_q,
							3,
							nb_frames,
							nb_sig );
							
	em_obj.learn_parameters(	sig,
								sizes,
								nb_sig,
								nb_iter);
				
	cout  << "T0 : " << endl;
	
	FILE * file = fopen( "params.txt","w");
	FILE * file_bis = fopen( "params2.txt","w");
	
	gsl_matrix toto;
	toto.size1 = 1;
	toto.size2 = 3;
	toto.data = em_obj.t0()->data;
	toto.tda = em_obj.t0()->stride;				
	gsl_matrix_fprintf_(stdout,
						API_DEFAULT_MATRIX_FORMAT,
						&toto );	
	fprintf( file_bis, "tracking::t_0 = ");
	gsl_matrix_fprintf_(file_bis,
						API_DEFAULT_MATRIX_FORMAT,
						&toto );						
	fprintf( file_bis, "\n");
						
	gsl_matrix_fprintf_(file,
						"\\begin{pmatrix}^\\end{pmatrix}\n^&^\\\\^%.2lf",
						&toto );				
						
									
	cout << endl;
	cout  << "sqrt_Q0 : " << endl;
	gsl_matrix_fprintf_(stdout,
						API_DEFAULT_MATRIX_FORMAT,
						em_obj.sqrt_q0() );	
						
						
						
						
	fprintf( file_bis, "tracking::sqrt_Q_0 = ");
	gsl_matrix_fprintf_(file_bis,
						API_DEFAULT_MATRIX_FORMAT,
						em_obj.sqrt_q0() );
	fprintf( file_bis, "\n");
	
	gsl_matrix * qq = gsl_matrix_alloc(6,6);
	gsl_blas_dgemm (	CblasTrans, 
						CblasNoTrans,
						1.0, 
						em_obj.sqrt_q0(), 
						em_obj.sqrt_q0(),
						0.0, 
						qq);

	gsl_matrix_fprintf_(file,
						"\\begin{pmatrix}^\\end{pmatrix}\n^&^\\\\^%.2lf",
						qq );
	
	
	
	cout << endl;
	cout  << "F : " << endl;
	gsl_matrix_fprintf_(stdout,
						API_DEFAULT_MATRIX_FORMAT,
						em_obj.f() );	
	fprintf( file_bis, "tracking::F = ");
	gsl_matrix_fprintf_(file_bis,
						API_DEFAULT_MATRIX_FORMAT,
						em_obj.f() );
	fprintf( file_bis, "\n");
			
	gsl_matrix_fprintf_(file,
						"\\begin{pmatrix}^\\end{pmatrix}\n^&^\\\\^%.2lf",
						em_obj.f() );
						
	cout << endl;
	cout  << "sqrt_Q : " << endl;
	gsl_matrix_fprintf_(stdout,
						API_DEFAULT_MATRIX_FORMAT,
						em_obj.sqrt_q() );	
						
	fprintf( file_bis, "tracking::sqrt_Q = ");	
	gsl_matrix_fprintf_(file_bis,
						API_DEFAULT_MATRIX_FORMAT,
						em_obj.sqrt_q() );	
	fprintf( file_bis, "\n");
	cout << endl;
	gsl_blas_dgemm (	CblasTrans, 
						CblasNoTrans,
						1.0, 
						em_obj.sqrt_q(), 
						em_obj.sqrt_q(),
						0.0, 
						qq);
	
	
	gsl_matrix_fprintf_(file,
						"\\begin{pmatrix}^\\end{pmatrix}\n^&^\\\\^%.2lf",
						qq );
	
	fclose(file);
	fclose(file_bis);
	return 0;
}
