#include "get_distance.hpp"
#include <cmath>
#define TMP_GET_DISTANCE( type ) template double get_distance( 	const type & x0,\
																const type & y0,\
																const type & x1,\
																const type & y1);
template<class type> double get_distance(  const type & x0,
										   const type & y0,
										   const type & x1,
										   const type & y1)
{
	double dx = ((double) x0 ) - x1;
	double dy = ((double) y0 ) - y1;
	return sqrt(dx * dx + dy * dy);
	
}
TMP_GET_DISTANCE(unsigned char)
TMP_GET_DISTANCE(char)
TMP_GET_DISTANCE(unsigned short int)
TMP_GET_DISTANCE(short int)
TMP_GET_DISTANCE(unsigned int)
TMP_GET_DISTANCE(int)
TMP_GET_DISTANCE(unsigned long int)
TMP_GET_DISTANCE(long int)
TMP_GET_DISTANCE(float)
TMP_GET_DISTANCE(double)
