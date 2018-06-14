#include "c_contour.hpp"
#include "compute_contour_8c.hpp"
#include "compute_contour_4c.hpp"
#include "compute_contour_4c_nb.hpp"
#include "digitalize_euclide_contour.hpp"
#include <cmath>
#include "get_perimeter.hpp"
c_contour :: c_contour(unsigned int max_nb_points)
{
	initialize();
	setup(max_nb_points);
}

void c_contour :: setup(unsigned int max_nb_points)
{
	free();
	initialize();
	_max_nb_points_contour = max_nb_points;
	max_nb_points_contour_2 = max_nb_points;
	_contour = new unsigned int[2 * max_nb_points + 4];
	contour_2 = new double[2 * max_nb_points];

}

template <class T> void c_contour :: compute_4c(const T * pixels,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												const T & value)
{
	nb_points_contour_2 = 0;
	compute_contour_4c ( _contour,
						 &_nb_points_contour,
						 _max_nb_points_contour + 2,
						 pixels,
						 width,
						 height,
						 width_step,
						 value);
	
}

template <class T> void c_contour :: compute_4c_nb(const T * pixels,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												const T & value)
{
	nb_points_contour_2 = 0;
	compute_contour_4c_nb ( _contour,
							 &_nb_points_contour,
							 _max_nb_points_contour,
							 pixels,
							 width,
							 height,
							 width_step,
							 value);
		
}

template <class T> void c_contour :: compute_8c(const T * pixels,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												const T & value)
{
	nb_points_contour_2 = 0;
	compute_contour_8c ( _contour,
						 &_nb_points_contour,
						 _max_nb_points_contour + 2,
						 pixels,
						 width,
						 height,
						 width_step,
						 value);
	
}

int c_contour :: digitalize(unsigned int nb_points)
{
	if (nb_points > max_nb_points_contour_2)
		return 42;
	else
	{
		nb_points_contour_2 = nb_points;
		digitalize_euclide_contour(contour_2,
								   nb_points,
								   _contour,
								   _nb_points_contour);
		return 0;
	}
}

c_contour :: ~c_contour()
{
	free();
	initialize();
}

void c_contour :: initialize()
{
	_nb_points_contour = 0;
	_max_nb_points_contour = 0;
	nb_points_contour_2 = 0;
	max_nb_points_contour_2 = 0;
	_contour = 0;
	contour_2 = 0;
}

void c_contour :: free()
{
	if (_contour)
		delete[] _contour;
	if (contour_2)
		delete[] contour_2;
}

int c_contour :: perimeter(double & c1,
						   double & c2)
{
	if (_nb_points_contour == 0)
	{
		return 1;
	}
	else
	{
		c1 = get_perimeter(_contour, _nb_points_contour);
	}
	if (nb_points_contour_2 != 0)
	{
		c2 = 0;
		for (unsigned int i = 0; i < nb_points_contour_2; ++ i)
		{
			double dx = contour_2[2 * i] - contour_2[2 * ((i + 1) % nb_points_contour_2) ];
			double dy = contour_2[2 * i + 1] - contour_2[2 * ((i + 1) % nb_points_contour_2) + 1];
			c2 += sqrt(dx * dx + dy * dy);
		}
	}
	return 0;
}

COMPUTE_4C(char)
COMPUTE_4C(unsigned char)
COMPUTE_4C(int)
COMPUTE_4C(unsigned int)

COMPUTE_8C(char)
COMPUTE_8C(unsigned char)
COMPUTE_8C(int)
COMPUTE_8C(unsigned int)

COMPUTE_4C_NB(char)
COMPUTE_4C_NB(unsigned char)

COMPUTE_4C_NB(int)
COMPUTE_4C_NB(unsigned int)
