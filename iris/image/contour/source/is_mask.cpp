#include "is_mask.hpp"

template <class T> int is_mask(const unsigned int & x,
				               const unsigned int & y,
							   const T * pixels,
							   unsigned int width,
							   unsigned int height,
							   unsigned int width_step,
							   const T & value)
{
	if ( (x >= width) || (y >= height) )
		return 0;

	return ( pixels[y * width_step + x] == value);
}

IS_MASK(unsigned char)
IS_MASK(char)
IS_MASK(int)
IS_MASK(unsigned int)
