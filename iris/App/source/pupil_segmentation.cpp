#include "lib_iris.hpp"
#include <sstream>
using namespace std;
int main( int argc, char ** argv )
{
	
	if ( argc < 2 )
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
			sig[nb_sig][sizes[nb_sig]] = gsl_vector_alloc(2);
			sig[nb_sig][sizes[nb_sig]]->data[0] = pupil_thread.pupil_seg_data().x_pupil - 0.5 * width;
			sig[nb_sig][sizes[nb_sig]]->data[1] = pupil_thread.pupil_seg_data().y_pupil - 0.5 * height;
			//sig[nb_sig][sizes[nb_sig]]->data[2] = 	sqrt( pow(pupil_thread.pupil_seg_data().a_pupil,2) + pow(pupil_thread.pupil_seg_data().b_pupil,2) );
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
	
	api_parameters apero( 2 * nb_sig + 1);
	
	
	//Sauvegarde des signaux
	for ( unsigned int i = 0; i < nb_sig; ++ i )
	{
		gsl_matrix * merde = gsl_matrix_alloc( sizes[i], 2 );
		for ( unsigned int j = 0; j < sizes[i]; ++ j )
		{
			gsl_vector_view toto = gsl_matrix_row( merde, j);
			gsl_vector_memcpy( &toto.vector, sig[i][j]);
		}
		stringstream oss;
		oss << "y(" << i + 1 <<")";
		
		api_variable var( oss.str().c_str(), API_MATRIX, (void*) merde );
		apero.add_variable(var);
		oss.str("");
		oss << "x(" << i + 1 <<")";
		var.set_name(oss.str().c_str());
		apero.add_variable(var);
		gsl_matrix_free( merde );
	}
	
	double q = nb_sig;
	api_variable var( "nb_signals", API_FLOAT, (void*) &q );
	apero.add_variable(var);
	
	apero.save(argv[2]);
	
}
