#include "discrete_convexe_boundary.hpp"
#define DISCRETE_CONVEX_BOUNDARY(type) template discrete_convex_boundary<type> :: discrete_convex_boundary(unsigned int max_nb_points);
template <class type> discrete_convex_boundary<type> :: discrete_convex_boundary(unsigned int max_nb_points)
: convex_boundary<type> :: convex_boundary(max_nb_points)
{
	initialize();
	left_contour = new unsigned int[max_nb_points];
	right_contour = new unsigned int[max_nb_points];	
}

#define DISCRETE_CONVEX_BOUNDARY_SETUP(type) template void discrete_convex_boundary<type> :: setup(unsigned int max_nb_points);
template <class type> void discrete_convex_boundary<type> :: setup(unsigned int max_nb_points)
{
	convex_boundary<type> :: setup(max_nb_points);
	discrete_convex_boundary<type> :: free();
	discrete_convex_boundary<type> :: initialize();
	left_contour = new unsigned int[max_nb_points];
	right_contour = new unsigned int[max_nb_points];	
}

#define DISCRETE_CONVEX_BOUNDARY_COMPUTE(type) template int discrete_convex_boundary<type> :: compute_segments(type dy);
template <class type> int discrete_convex_boundary<type> :: compute_segments(type dy)
{
	//Recherche du min et du max du contour
	compute_extrema();
	unsigned int i;
	unsigned int nb_pts = (y_max - y_min) / ( (double) dy ) + 1;

	//Test
	if (nb_pts  > convex_boundary<type> :: _max_nb_points)
	{

		return 1;
	}
	//Enveloppe gauche de l'objet
	unsigned int j = 0;

	for (i = y_min; i < y_max; i += dy)
	{
		double x1, x2, 
			   y1, y2;
		while ( (unsigned int) (y2 = ( convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j] + 1] ) ) <= i)
		{

			++ j;
		}
		
		
		y1 = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[(j - 1)] + 1];
		x1 = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[(j - 1)] ];
		x2 = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j] ];
		right_contour[(unsigned int) ((i - y_min) / dy) ] = (unsigned int) (x1 + (i - y1) / (y2 - y1) * (x2 - x1) + 0.5);	
	}

	if (i == y_max)
	{
		while ( (unsigned int) (convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j] + 1] ) != y_max)
		{
			++ j;
		}
		right_contour[(unsigned int) ((i - y_min) / dy) ] = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j]];
		
	}
	
	j = 0;

	for (i = y_min; i < y_max; i += dy)
	{
		double x1, x2, 
			   y1, y2;
		while ( (unsigned int) (y2 = ( convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j] + 1] ) ) <= i)
		{
			-- j;
			if (j == (unsigned int) - 1)
				j =  convex_boundary<type> :: _nb_points - 1;
		}
		y1 = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[(j + 1) % convex_boundary<type> :: _nb_points] + 1];
		x1 = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[(j + 1) % convex_boundary<type> :: _nb_points] ] ;
		x2 = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j] ];
		
		left_contour[(unsigned int) ((i - y_min) / dy)] = (unsigned int) (x1 + (i - y1) / (y2 - y1) * (x2 - x1) + 0.5);	
	}

	if (i == y_max)
	{
		while ( (unsigned int) (convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j] + 1] ) != y_max)
		{
			-- j;
			if (j == (unsigned int) - 1)
				j = convex_boundary<type> :: _nb_points - 1;
		}
		left_contour[(unsigned int) ((i - y_min) / dy)] = convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[j]];
	}

	return 0;
}

#define DISCRETE_CONVEX_BOUNDARY_DEST(type) template discrete_convex_boundary<type> :: ~discrete_convex_boundary();
template <class type> discrete_convex_boundary<type> :: ~discrete_convex_boundary()
{
	discrete_convex_boundary<type> :: free();
	discrete_convex_boundary<type> :: initialize();
}

#define DISCRETE_CONVEX_BOUNDARY_FREE(type) template void discrete_convex_boundary<type> :: free();
template <class type> void discrete_convex_boundary<type> :: free()
{
	if (left_contour)
		delete[] left_contour;
	if (right_contour)
		delete[] right_contour;
}

#define DISCRETE_CONVEX_BOUNDARY_INITIALIZE(type) template void discrete_convex_boundary<type> :: initialize();
template <class type> void discrete_convex_boundary<type> :: initialize()
{
	left_contour = 0;
	right_contour = 0;
	y_min = 0;
	y_max = 0;
}

#define DISCRETE_CONVEX_BOUNDARY_COMPUTE_EXTREMA(type) template void discrete_convex_boundary<type> :: compute_extrema();
template <class type> void discrete_convex_boundary<type> :: compute_extrema()
{
	y_min = (unsigned int) convex_boundary<type> :: points_[2 * convex_boundary<type> :: _labels[0] + 1];
	unsigned id = 0;
	for (unsigned int i = 1; i < convex_boundary<type> :: _nb_points; ++ i)
	{
		unsigned k = convex_boundary<type> :: _labels[i];
		if ( convex_boundary<type> :: points_[2 * k + 1] > convex_boundary<type> :: points_[2 * id + 1] )
			id = k;
	}
	y_max = (unsigned int) (convex_boundary<type> :: points_[2 * id + 1]);
}

#define DISCRETE_CONVEX_BOUNDARY_RIGHT(type) template const unsigned int * discrete_convex_boundary<type> :: get_left_contour() const;
template <class type>  const unsigned int * discrete_convex_boundary<type> :: get_left_contour() const
{
	return left_contour;
}
#define DISCRETE_CONVEX_BOUNDARY_LEFT(type) template const unsigned int * discrete_convex_boundary<type> :: get_right_contour() const;
template <class type>  const unsigned int * discrete_convex_boundary<type> :: get_right_contour() const
{
	return right_contour;
}
#define DISCRETE_CONVEX_BOUNDARY_Y_MIN(type) template unsigned int discrete_convex_boundary<type> :: get_y_min() const;
template <class type> unsigned int discrete_convex_boundary<type> :: get_y_min() const
{
	return y_min;
}
#define DISCRETE_CONVEX_BOUNDARY_Y_MAX(type) template unsigned int discrete_convex_boundary<type> :: get_y_max() const;
template <class type> unsigned int discrete_convex_boundary<type> :: get_y_max() const
{
	return y_max;
}



#define DISCRETE_CONVEX_BOUNDARY_CLASS(type) DISCRETE_CONVEX_BOUNDARY(type)\
											 DISCRETE_CONVEX_BOUNDARY_SETUP(type)\
											 DISCRETE_CONVEX_BOUNDARY_COMPUTE(type)\
											 DISCRETE_CONVEX_BOUNDARY_DEST(type)\
											 DISCRETE_CONVEX_BOUNDARY_FREE(type)\
											 DISCRETE_CONVEX_BOUNDARY_INITIALIZE(type)\
											 DISCRETE_CONVEX_BOUNDARY_COMPUTE_EXTREMA(type)\
											 DISCRETE_CONVEX_BOUNDARY_RIGHT(type)\
											 DISCRETE_CONVEX_BOUNDARY_LEFT(type)\
											 DISCRETE_CONVEX_BOUNDARY_Y_MIN(type)\
											 DISCRETE_CONVEX_BOUNDARY_Y_MAX(type)

DISCRETE_CONVEX_BOUNDARY_CLASS(unsigned int)
DISCRETE_CONVEX_BOUNDARY_CLASS(double)
