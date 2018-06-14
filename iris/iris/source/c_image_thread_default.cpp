#include "c_image_thread.hpp"
#include "iris_default.hpp"
int c_image_thread :: default_setup ( ostream * _err_stream )
{
	unsigned int width;
	unsigned int height;
	
	WIDTH(width);
	HEIGHT(height);
	DATA_TYPE(type);
	
	_data = new image_data;
	_data->setup( cvSize( width, height ) );

	err_stream = _err_stream;
	
	return 0;
}
