#include "c_get_iris_template_pthread.hpp"
#include "distances.hpp"
#include <ctime>
#include <signal.h>
void * image_function ( void * params )
{
	c_get_iris_template_pthread * data = (c_get_iris_template_pthread *) params;
	data->image_acquistion();
	return NULL;
}

void * pupil_function ( void * params )
{
	c_get_iris_template_pthread * data = (c_get_iris_template_pthread *) params;
	data->pupil_segmentation();
	return NULL;
}

void * iris_function ( void * params )
{
	c_get_iris_template_pthread * data = (c_get_iris_template_pthread *) params;
	data->iris_segmentation();
	return NULL;
}


c_get_iris_template_pthread :: c_get_iris_template_pthread ( void )
{
	initialize();
}

int c_get_iris_template_pthread :: setup ( unsigned int argc, char ** argv, const char * save_rep,
											ostream * stream, bool display  )
{
	api_parameters params;
	stringstream oss;
	int q = 0;
	free();
	initialize();
	if ( argc < 1 )
	{
		cout << "Error : Missing argument(s) " << endl;
		return 1;
	}

	//Allco
	image_thread = new c_image_thread;
	pupil_thread = new c_pupil_thread;
	iris_thread = new c_iris_thread;
	focus_score_iris = new c_focus_score;
	focus_score_pupil = new c_focus_score;

	if ( image_thread->setup(	argv,
								argc,
								stream) )
		q = 1;
	if ( pupil_thread->setup(	argv,
								argc,
								stream) )
		q = 1;
	if ( iris_thread->setup(	argv,
								argc,
								stream) )
		q = 1;

	for ( unsigned int i = 0; i < argc; ++ i )
		params.load( argv[i] );

	if ( focus_score_iris->setup( 	params ) )
		q = 1;
	if ( focus_score_pupil->setup( params ) )
		q = 1;





	unsigned int size_buffer_image;
	oss << BUFFER_NAMESPACE << "::" << "size_buffer_image";
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&size_buffer_image,
									stream ) )
		q = 1;
	oss.str("");

	unsigned int size_buffer_pupil;
	oss << BUFFER_NAMESPACE << "::" << "size_buffer_pupil";
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&size_buffer_pupil,
									stream ) )
		q = 1;
	oss.str("");

	unsigned int size_buffer_iris;
	oss << BUFFER_NAMESPACE << "::" << "size_buffer_iris";
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&size_buffer_iris,
									stream ) )
		q = 1;
	oss.str("");
	oss << "score" << "::" << "width";
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&r_width,
									stream ) )
		q = 1;
	oss.str("");
	oss << "score" << "::" << "height";
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&r_height,
									stream ) )
		q = 1;
	oss.str("");

	oss << save_rep << "/" << "fps.m";
	fps_file = new ofstream( oss.str().c_str() );
	if ( ! (*fps_file) )
	{
		if ( stream )
		{
			*stream << "Error : Can't load " <<  oss.str().c_str() << endl;
			return 1;
		}
	}
	*fps_file << "Videos = {};" << endl;
	if ( q )
		return 1;


	//Alloc mémoire
	//Buffer image
	{
		CvSize tmp = cvSize( 	image_thread->data().width,
								image_thread->data().height );
		buffer_image
			= new c_buffer ( 	size_buffer_image,
								image_data_alloc,
								image_data_copy,
								image_data_free,
								(void*) &tmp );

		b_data_image = new c_buffer_data ( (void*) &tmp,
											 image_data_alloc,
											 image_data_copy,
											 image_data_free );
	}

	//Buffer pupille
	{
		CvSize tmp = cvSize( 	image_thread->data().width,
								image_thread->data().height );
		buffer_pupil
			= new c_buffer ( 	size_buffer_pupil,
								pupil_data_alloc,
								pupil_data_copy,
								pupil_data_free,
								(void*) &tmp );

		b_data_pupil = new c_buffer_data ( (void*) &tmp,
											 pupil_data_alloc,
											 pupil_data_copy,
											 pupil_data_free );

	}

	//Buffer iris
	{
		iris_data_params tmp;
			tmp.img_width = image_thread->data().width;
			tmp.img_height = image_thread->data().height;
			tmp.polar_width = iris_thread->iris_seg_data().nb_directions;
			tmp.polar_height = iris_thread->iris_seg_data().nb_samples;
			tmp.nb_samples_iris = iris_thread->iris_seg_data().nb_samples_iris;
			tmp.iris_code_width = iris_thread->iris_seg_data().nb_directions_code;
			tmp.iris_code_height = iris_thread->iris_seg_data().nb_samples_code;
			tmp.iris_width = iris_thread->iris_seg_data().iris_width;
			tmp.iris_height = iris_thread->iris_seg_data().iris_height;
		buffer_iris
			= new c_buffer ( 	size_buffer_iris,
								iris_data_alloc,
								iris_data_copy,
								iris_data_free,
								(void*) &tmp );

		b_data_iris = new c_buffer_data ( (void*) &tmp,
											 iris_data_alloc,
											 iris_data_copy,
											 iris_data_free );

	}
	
			 
						 
	display_on = display;
	if ( display_on )
	{
		
		pthread_mutex_init ( &mutex1 , 
							 NULL );
		pthread_mutex_init ( &mutex2 , 
							 NULL );		 
		pthread_mutex_init ( &mutex3 , 
							 NULL );	
			
		
		unsigned int 	width = image_thread->data().width,
						height = image_thread->data().height,
						p_height = iris_thread->iris_seg_data().nb_samples * (image_thread->data().width / ( (double) iris_thread->iris_seg_data().nb_directions ) ),
						c_height = iris_thread->iris_seg_data().nb_samples_code * (image_thread->data().width / ( (double) iris_thread->iris_seg_data().nb_directions_code ) );
		
		
		image_display_obj = new c_display_image_data;
		pupil_display_obj = new c_display_pupil_data;
		iris_display_obj = new c_display_iris_data;
		
		
		iris_display_obj->setup( width,
								 height,
								 width / 2,
								 height + p_height + c_height,
								 10,
								 INT_RGB(255,255,255),
								 3,
								 INT_RGB(0,255,0),
								 INT_RGB(255,0,0),
								 INT_RGB(255,0,0),
								 INT_RGB(0,0,255),
								 INT_RGB(255,0,0),
								 width,
								 c_height,
								 width,
								 p_height );
		d_iris = cvCreateImage( cvSize( width + width / 2, 
										height + p_height + c_height ),
								IPL_DEPTH_8U,
								3 );
								 
		image_display_obj->setup( 	width,
									height, 
									width / 2,
									height,
									10,
									INT_RGB(255,255,255) );
									
		d_image = cvCreateImage( 	cvSize( width + width / 2, 
											height ),
									IPL_DEPTH_8U,
									3 );			
									
		pupil_display_obj->setup ( 	width,
									height, 
									width / 2,
									height,
									10,
									INT_RGB(255,255,255),
									3,
									INT_RGB(0,255,0) );
		d_pupil = cvCreateImage( 	cvSize( width + width / 2, 
											height ),
									IPL_DEPTH_8U,
									3 );	
	}
	
	
						 
	return 0;
}

unsigned int c_get_iris_template_pthread ::  help ( ostream & out, unsigned int id ) const
{
	out << "argv[" << id ++ << "- N] : .cfg files " << endl;
	return id;
	
}


int c_get_iris_template_pthread :: run ( 	const char * filename,
												int & end )
{
	nb_abb = 0;
	nb_frames = 0;
	nb_i_frames = 0;
	nb_p_frames = 0;
	_filename = filename;
	clock_t time_1 = clock();
	if ( video )
		cvReleaseCapture(&video );
	//Chargement de la vidéo
	if ( ! strcmp( filename, "webcam" ) )
	{
		cout << " Error: not yet implemented" << endl;
		return 1;
	}
	else
	{
		video = cvCaptureFromAVI( filename );
		if ( ! video )
		{
			cout << "Error: can't load " << filename << endl;
			return 1;
		}
	}
	*fps_file << "Videos(" << video_id + 1 << ") = \"" << filename <<"\";"  << endl;
	
	//Reset
		buffer_pupil->reset();
		buffer_image->reset();
		buffer_iris->reset();
		image_thread->reset_id();
		pupil_thread->reset();
		iris_thread->reset();

	//end
		p_end = &end;
		img_end = false;
		pupil_end = false;
		iris_end = false;

	//Création des threads
		pthread_create( tab_threads + 0,
						NULL,
						image_function,
						(void*) ( this ) ) ;

		pthread_create( tab_threads + 1,
						NULL,
						pupil_function,
						(void*) ( this ) ) ;
						//~ 
		pthread_create( tab_threads + 2,
						NULL,
						iris_function,
						(void*) ( this ) ) ;
						
						//~ 
						//~ 
						
						
						
		display();
	//Jonction des threads
	for ( unsigned i = 0; i < 3; ++ i )
		pthread_join( tab_threads[i], NULL);
	clock_t time_2 = clock();
	
	*fps_file << "nb_frames(" << video_id + 1 << ") = " << cvGetCaptureProperty( video , CV_CAP_PROP_FRAME_COUNT )  <<";"  << endl;
	*fps_file << "times(" << video_id + 1 << ") = " << (time_2 - time_1) / ( (double) CLOCKS_PER_SEC )  <<";"  << endl;
	*fps_file << "fps(" << video_id + 1 << ") = " << ( cvGetCaptureProperty( video , CV_CAP_PROP_FRAME_COUNT ) ) * ( (double) CLOCKS_PER_SEC ) / (time_2 - time_1)  <<";"  << endl;
	*fps_file << "p_thread_pupil(" << video_id + 1 << ") = " << nb_p_frames / ((double) nb_frames ) <<";"  << endl;
	*fps_file << "p_thread_iris(" << video_id + 1 << ") = " << nb_i_frames / ((double) nb_frames )  <<";"  << endl;
	*fps_file << "p_abberations(" << video_id + 1 << ") = " <<  nb_abb / ((double) nb_i_frames ) <<";"  << endl;
	
	video_id ++;
	return 0;
}

c_get_iris_template_pthread :: ~c_get_iris_template_pthread()
{
	free();
	initialize();
}

void c_get_iris_template_pthread :: initialize()
{
	video = 0;
	buffer_image = 0; // Score = id
	buffer_pupil = 0;
	buffer_iris = 0;
	fps_file = 0;
	image_thread = 0;
	pupil_thread = 0;
	iris_thread = 0;
	focus_score_pupil = 0;
	focus_score_iris = 0;
	p_end = 0;
	img_end = true;
	pupil_end = true;
	iris_end = true;
	b_data_image = 0;
	b_data_iris = 0;
	b_data_pupil = 0;
	r_width = 0;
	r_height = 0;
	video_id = 0;
	nb_abb = 0;
	nb_frames = 0;
	nb_i_frames = 0;
	nb_p_frames = 0;
	d_image = 0;
	d_pupil = 0; 
	d_iris = 0;
	image_display_obj = 0;
	pupil_display_obj = 0;
	iris_display_obj = 0;
	display_on = false;
	
}

void c_get_iris_template_pthread :: free()
{
	if ( video )
		cvReleaseCapture(&video );
	if ( buffer_image )
		delete buffer_image;
	if ( buffer_pupil )
		delete buffer_pupil;
	if ( buffer_iris )
		delete buffer_iris;
	if ( image_thread )
		delete image_thread;
	if ( pupil_thread )
		delete pupil_thread;
	if ( iris_thread )
		delete iris_thread;
	if ( focus_score_pupil )
		delete focus_score_pupil;
	if ( focus_score_iris )
		delete focus_score_iris;
	if ( b_data_image )
		delete b_data_image;
	if ( b_data_iris )
		delete b_data_iris;
	if ( b_data_pupil )
		delete b_data_pupil;
	if ( fps_file)
	{
		fps_file->close();
		delete fps_file;
	}
	if ( d_image )
		cvReleaseImage( &d_image);
	if ( d_pupil )
		cvReleaseImage( &d_pupil);
	if ( d_iris )
		cvReleaseImage( &d_iris );
	if( image_display_obj )
		delete image_display_obj;
	if ( pupil_display_obj )
		delete pupil_display_obj;
	if ( iris_display_obj )	
		delete iris_display_obj;
	if ( display_on )
	{
		pthread_mutex_destroy ( &mutex1 );	
		pthread_mutex_destroy ( &mutex2 );					
		pthread_mutex_destroy ( &mutex3 );								
	}
							
							
							
}

void c_get_iris_template_pthread :: image_acquistion ( )
{
	clock_t time_1, time_2;
	IplImage * image = NULL;
	do
	{
		time_1 = clock();
		image = cvQueryFrame(video);
		if (image)
		{
			image_thread->process( image, _filename );
			if ( image_thread->data().img_ok )
			{
				nb_frames ++;
				image_thread->set_score(
					focus_score_pupil->get_score( 	image_thread->data().image ) );													
				buffer_image->add_object( 	image_thread->data().frame_id,
											image_thread->data().score,
											(void*) & (image_thread->data()) );
				
				if ( display_on )	
				{
					pthread_mutex_lock ( &mutex1 );
					image_display_obj->display(image_thread->data());
					pthread_mutex_unlock ( &mutex1 );
				}		
				
				time_2 = clock();
				int q = 40000 - ( (time_2 - time_1) * 1000000 ) / CLOCKS_PER_SEC; 
				if ( q > 0 )
					usleep( q );

			}

		}

	} while ( !(*p_end ) && image != NULL);
	img_end = true;
}

void c_get_iris_template_pthread :: pupil_segmentation ( )
{

	bool end = false;	
	do
	{
		if ( img_end )
			end = true;
			
		buffer_image->get_object (	b_data_image,
									0 );
		buffer_image->delete_object( b_data_image );	
		if ( b_data_image->score() > 0 )
		{
			pupil_thread->segment_pupil( *( (image_data*) b_data_image->data() ) );

			if ( pupil_thread->pupil_seg_data().seg_ok )
			{

				nb_p_frames ++;
				//Score
				buffer_pupil->add_object( 	pupil_thread->pupil_seg_data()._img_data.frame_id,
											b_data_image->score(),
											(void*) & (pupil_thread->pupil_seg_data() ) );
		
				if ( display_on )	
				{
					pthread_mutex_lock ( &mutex2 );
					pupil_display_obj->display(pupil_thread->pupil_seg_data());
					pthread_mutex_unlock ( &mutex2 );
				}	
		
											
			}
		}
	} while ( !end  );
	pupil_end = true;
}

void c_get_iris_template_pthread :: iris_segmentation ( )
{
	bool end = false;
	do
	{

		if ( pupil_end )
			end = true;

		buffer_pupil->get_object (	b_data_pupil,
									0 );		
		buffer_pupil->delete_object( b_data_pupil );
		if ( b_data_pupil->score() > 0 )
		{
	
			iris_thread->segment_iris( *( (pupil_data*) b_data_pupil->data() ) );

			if ( iris_thread->iris_seg_data().seg_ok )
			{
				
				
				
				nb_i_frames ++;
				//Détection des chgts d'yeux ou des abérations
				buffer_iris->get_object (	b_data_iris,
											0 );
				
				
				//Score
				CvRect rect = cvRect( 	iris_thread->iris_seg_data().new_x_iris - iris_thread->iris_seg_data().new_r_iris,
										iris_thread->iris_seg_data().new_y_iris - iris_thread->iris_seg_data().new_r_iris,
										2 * iris_thread->iris_seg_data().new_r_iris,
										2 * iris_thread->iris_seg_data().new_r_iris );
				if ( ! ( 	rect.x + rect.width > iris_thread->iris_seg_data().iris_image->width 	||
							rect.y + rect.height > iris_thread->iris_seg_data().iris_image->height 	||
							rect.x <= 0 														   	||
							rect.y <= 0) )
				{
					buffer_iris->add_object( 	iris_thread->iris_seg_data().p_data._img_data.frame_id,
												iris_thread->iris_seg_data().nrj_ratio,
												(void*) & (iris_thread->iris_seg_data() ) );
												
					if ( display_on )	
					{
						pthread_mutex_lock ( &mutex3 );
						//~ pupil_display_obj->display(iris_thread->iris_seg_data().p_data);
						iris_display_obj->display(iris_thread->iris_seg_data());
						pthread_mutex_unlock ( &mutex3 );
					}										
				}								
			}
		}
		usleep(1000);

	} while ( !end );
	
	iris_end = true;
}

void c_get_iris_template_pthread :: display ( )
{
	unsigned int 	width = image_thread->data().width,
					height = image_thread->data().height,
					p_height = iris_thread->iris_seg_data().nb_samples * (image_thread->data().width / ( (double) iris_thread->iris_seg_data().nb_directions ) ),
					c_height = iris_thread->iris_seg_data().nb_samples_code * (image_thread->data().width / ( (double) iris_thread->iris_seg_data().nb_directions_code ) );
	
	if ( display_on )
	{
		char k = 0;
		while ( k != 'q' && k != 'Q' )
		{
			pthread_mutex_lock ( &mutex1 );
				cvSetImageROI(	d_image, 
								cvRect( 0, 0, width, height ) );
				cvCopyImage( 	image_display_obj->seg_image(), 
								d_image );
				cvSetImageROI(	d_image, 
								cvRect( width, 0, width / 2, height ) );
				cvCopyImage( 	image_display_obj->data_image(), 
								d_image );
				cvResetImageROI(d_image);
			pthread_mutex_unlock ( &mutex1 );

			pthread_mutex_lock ( &mutex2 );
				cvSetImageROI(	d_pupil, 
								cvRect( 0, 0, width, height ) );
				cvCopyImage( 	pupil_display_obj->pupil_image(), 
								d_pupil );
				
				cvSetImageROI(	d_pupil, 
								cvRect( width, 0, width / 2, height ) );
				cvCopyImage( 	pupil_display_obj->data_image(), 
								d_pupil);
				cvResetImageROI( d_pupil );
			pthread_mutex_unlock ( &mutex2 );

			
			pthread_mutex_lock ( &mutex3 );
				cvSetImageROI(	d_iris, 
								cvRect( 0, 0, width, height ) );
				cvCopyImage( 	iris_display_obj->iris_image(), 
								d_iris );
				cvSetImageROI(	d_iris, 
								cvRect( width, 0, width / 2, height + c_height + p_height ) );
				cvCopyImage( 	iris_display_obj->data_image(), 
								d_iris);
				cvSetImageROI(	d_iris, 
								cvRect(  0, height, width, p_height ) );
				cvCopyImage( 	iris_display_obj->polar_image(), 
								d_iris);
				cvSetImageROI(	d_iris, 
								cvRect( 0,  height + p_height, width, c_height ) );
				cvCopyImage( 	iris_display_obj->iris_code(), 
								d_iris);

				cvResetImageROI(d_iris);
				
			k = cvWaitKey(40);//Pause de 40ms
			
			if( k == 's' || k == 'S' )
			{
				cvSaveImage("pupil.png", iris_display_obj->pupil_image() );
				cvSaveImage("iris.png", iris_display_obj->iris_image() );
				cvSaveImage("p_image.png", iris_display_obj->polar_image() );
				cvSaveImage("i_code.png", iris_display_obj->iris_code() );
				
				
			}
				
				
			pthread_mutex_unlock ( &mutex3 );
			cvShowImage( "Image thread", d_image );
			cvShowImage( "Pupil thread", d_pupil );
			cvShowImage( "Iris thread", d_iris );

			

			
			
			
			
		}
	}
	
}


int c_get_iris_template_pthread :: save ( const char * rep_name )
{
	int q = 0;
	if ( iris_end )
	{

#ifdef __linux
	mkdir ( rep_name, 014777 );
#elif __APPLE__ 
	mkdir ( rep_name, 014777 );
#elif _WIN32
	_mkdir( rep_name );
#elif _WIN64
	_mkdir( rep_name );
#else
	#error
#endif

		for ( unsigned int i = 0; i < buffer_iris->nb_objects(); ++ i )
		{
			buffer_iris->get_object (	b_data_iris,
										i );

			if ( b_data_iris->score() < 0 )
				break;

			if ( ( (iris_data *) b_data_iris->data() )->save( rep_name, i ) )
				q = 1;
		}
	}
	return q;
}
