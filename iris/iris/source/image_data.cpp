#include "image_data.hpp"
image_data :: image_data()
{
	initialize();
}

image_data :: image_data( const CvSize & image_size )
{
	initialize();
	setup ( image_size );
}

image_data :: image_data( const image_data & data )
{
	initialize();
	setup ( data );
}

void image_data :: initialize()
{
	img_ok = false;
	frame_id = 0;
	image = 0;
	score = 0;
	width = 0;
	height = 0;
}

void image_data :: free()
{
	if ( image )
	{
		cvReleaseImage( &image );
	}
}

void image_data :: alloc()
{

	image = cvCreateImage (	cvSize( width, height),
							IPL_DEPTH_8U,
							1 );				
}

int image_data :: setup()
{
	free();
	initialize();
	return 0;
}

int image_data :: setup ( const CvSize & image_size )
{
	free();
	initialize();

	width = image_size.width;
	height = image_size.height;

	if ( ! width || ! height )
		return 1;
	alloc();
	return 0;
}

int image_data :: setup ( const image_data & data )
{
	free();
	initialize();

	width = data.width;
	height = data.height;

	if ( ! width || ! height )
		return 1;
	alloc();

	frame_id = data.frame_id;

	cvSetImageROI( image,
				   cvRect( 	0,
							0,
							GET_IMAGE_WIDTH( data.image ),
							GET_IMAGE_HEIGHT( data.image ) ) );
								
	cvCopyImage( data.image, image );
	name = data.name;
		score = data.score;
	img_ok = data.img_ok;

	return 0;
}

image_data & image_data :: operator=( const image_data & data )
{
	if ( ! ( *this == data ) )
		setup ( data );
	else
	{
		
		frame_id = data.frame_id;

		cvSetImageROI( image,
					   cvRect( 	0,
								0,
								GET_IMAGE_WIDTH( data.image ),
								GET_IMAGE_HEIGHT( data.image ) ) );
								
		cvCopyImage( data.image, image );
		name = data.name;
		score = data.score;
		img_ok = data.img_ok;
	}

	return *this;
}

image_data :: ~image_data()
{
	free();
	initialize();
}

bool image_data :: operator==( const image_data & data ) const
{
	if ( 	data.width == width &&
			data.height == height )
		return true;
	return false;
}

int image_data :: save ( const char * rep,
						   unsigned int id ) const
{
	int q = 0;
	stringstream oss;
	if ( img_ok == false )
		return 1;

	//Données de seg.
	{
		oss << rep << "/" << "image_data_" << id << ".txt";
		ofstream file(oss.str().c_str());
		oss.str("");

		file << "#Image seg. data" << endl;
		file << "id = " << frame_id << endl;
		file << "name = \"" << name << "\"" << endl;
		
		file.close();
		
	}

	{
		oss << rep << "/" << "image_" << id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  image );
		oss.str("");
	}
	
	return q;
}

int image_data :: load ( const char * rep,
						   unsigned int id )
{
	free();
	initialize();

	int q = 0;
	stringstream oss;
	//Données de seg.
	{
		oss << rep << "/" << "image_data_" << id << ".txt";
			api_parameters params;
			if ( params.load( oss.str().c_str() ) )
				q = 1;
			else
			{
				if ( api_get_positive_integer( params,
											   "id",
											   &frame_id ) )
					q = 1;
				if ( api_get_string( 	params,
										"name",
										&name ) )
					q = 1;
			}
		oss.str("");
	}

	//Image
	oss << rep << "/" << "image_" << id << ".png";
	image = cvLoadImage ( oss.str().c_str(),
						  0 );
	oss.str("");
	
	if ( image == NULL )
		q = 1;
	else
	{
		width = image->width;
		height = image->height;
	}
	if ( !q )
		img_ok = true;
	return q;
}

int image_data :: set_dim( unsigned int w, unsigned int h)
{
	if ( w > width || h > height )
		setup ( cvSize( w, h));
	cvSetImageROI(image, cvRect(0, 0, w, h ) );
	return 0;
}


void * image_data_alloc( const void * params )
{
	CvSize * data = (CvSize * ) params;
	return (void*) new image_data(*data);
}

void image_data_copy( 	void * data_tg,
						const void * data_src )
{
	image_data 	* in,
					* out;

	in = ( image_data * ) data_src;
	out = ( image_data * ) data_tg;

	(*out) = (*in);
}

void image_data_free( void * p )
{
	image_data * fshgjkkjf = (image_data*) p;
	delete fshgjkkjf;
}

