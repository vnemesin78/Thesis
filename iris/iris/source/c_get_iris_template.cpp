#include "c_get_iris_template.hpp"

c_get_iris_template :: c_get_iris_template()
{
	image_thread = new c_image_thread;
	pupil_thread = new c_pupil_thread;
	iris_thread = new c_iris_thread;
	
	
	
}
c_get_iris_template :: ~c_get_iris_template()
{
	if ( image_thread )
		delete image_thread;
	if ( pupil_thread )
		delete pupil_thread;
	if ( iris_thread )
		delete iris_thread;
	
}

int c_get_iris_template :: setup( 	char ** argv, 
										unsigned int argc,
										ostream * stream)
{
	int q = 0;
	if (argc == 0 )
	{
		cout << "Error: Missing argument(s)" << endl;
		return 1;
	}
	if ( image_thread->setup(	argv,  
								argc,
								stream) )
		q = 1;
		
	if ( pupil_thread->setup(	argv,  
								argc,
								stream ) )
		q = 1;
	if ( iris_thread->setup(	argv,  
								argc,
								stream ) )
		q = 1;
		

		
		
		
	if ( q )
		return 1;

	return 0;	
}


unsigned int c_get_iris_template :: help( ostream & out, unsigned int id ) const
{
	out << "argv[" << id ++ << "- N] : .cfg files" << endl;
	return id;
	
}







int  c_get_iris_template :: segment( const IplImage * image, 
									  const char * name )
{
	pupil_thread->reset();
	iris_thread->reset();	
	
	image_thread->process( image, name );
	pupil_thread->segment_pupil( image_thread->data() );
	return ( iris_thread->segment_iris( pupil_thread->pupil_seg_data() ) );
}

int c_get_iris_template :: safe_save ( 	const char * template_file, 
											const char * output_file,
											const char * img_filename)
{
	//Status
	{
		ofstream out_file ( output_file );
		if ( !out_file )
		{
			cout << "Can't open " << output_file << "!" << endl;
			return 1;
		}		
			
		if ( data().seg_ok )
		{
			out_file << 1;
			out_file.close();
		}
		else
		{
			out_file << 0;
			out_file.close();
			return 0;
		}
	}
	
	FILE * file = fopen( template_file, "w" );
	
	if ( !file )
	{
		cout << "Can't open " << template_file << "!" << endl;
		return 1;
	}
	int type = 0;
	fwrite ( &type, sizeof ( int ), 1, file );
	fwrite ( &(data().nb_directions_code), sizeof ( unsigned int ), 1, file );
	fwrite ( &(data().nb_samples_code), sizeof ( unsigned int ), 1, file );
	
	
	unsigned char	* re = new unsigned char[data().nb_directions_code * data().nb_samples_code], 
					* im = new unsigned char[data().nb_directions_code * data().nb_samples_code], 
					* mask = new unsigned char[data().nb_directions_code * data().nb_samples_code];
	
	for ( unsigned int i = 0; i < data().nb_samples_code; ++ i )
	{
		for ( unsigned int j = 0; j < data().nb_directions_code; ++ j )
		{
			re[i * data().nb_directions_code + j] = ( (unsigned char*) ( data().code->imageData + i * data().code->widthStep ) )[j];
			im[i * data().nb_directions_code + j] = ( (unsigned char*) ( data().code->imageData + ( i + data().nb_samples_code ) * data().code->widthStep ) )[j];
			mask[i * data().nb_directions_code + j] = ( (unsigned char*) ( data().code_mask->imageData + i * data().code_mask->widthStep ) )[j];
		}		
		
	}
	fwrite( re, sizeof(char), data().nb_directions_code * data().nb_samples_code, file );
	fwrite( im, sizeof(char), data().nb_directions_code * data().nb_samples_code, file );
	fwrite( mask, sizeof(char), data().nb_directions_code * data().nb_samples_code, file );
	
	delete[] re;
	delete[] im;
	delete[] mask;
	
	
	
	
	
	
	
	
	
	
	
	fclose(file);
	
	return 0;
}
