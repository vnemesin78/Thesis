#include "get_perimeter.hpp"
#include "get_distance.hpp"
#define TMP_GET_PERIMETER(type) template double get_perimeter(	const type * contour,\
																unsigned int nb_points_contour);

template <class type> double get_perimeter(	const type * contour,
											unsigned int nb_points_contour)
{
	double peri = 0;
	for(unsigned int i = 0; i < nb_points_contour; ++i)
	{
			peri += get_distance(contour[2 * i], 
								 contour[2 * i + 1],
								 contour[2 * ((i + 1) % nb_points_contour)],
								 contour[2 * ((i + 1) % nb_points_contour) + 1]);
	}
	return peri;
}
TMP_GET_PERIMETER(unsigned char)
TMP_GET_PERIMETER(char)
TMP_GET_PERIMETER(unsigned short int)
TMP_GET_PERIMETER(short int)
TMP_GET_PERIMETER(unsigned long int)
TMP_GET_PERIMETER(long int)
TMP_GET_PERIMETER(unsigned int)
TMP_GET_PERIMETER(int)


TMP_GET_PERIMETER(float)
TMP_GET_PERIMETER(double)
