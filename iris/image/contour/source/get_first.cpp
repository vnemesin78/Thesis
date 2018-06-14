#include "get_first.hpp"
template <class T> int get_first(unsigned int & x,
							     unsigned int & y,
							     const T * pixels,
							     unsigned int width,
							     unsigned int height,
							     unsigned int width_step,
							     const T & value)
{
	for(y = 0; y < height; y++)
	{	
		for(x = 0; x < width; x++)
		{
			if (pixels[y * width_step + x] == value)
				return 0;
		}
	}
	return 1;
}

GET_FIRST(unsigned char)
GET_FIRST(unsigned int)

GET_FIRST(char)
GET_FIRST(int)
