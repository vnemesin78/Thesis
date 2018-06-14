#include "convexe_boundary.hpp"
#include <algorithm>
#include <iostream>
using namespace std;
using namespace std;
#define CONVEXE_BOUNDARY_C1(type) template convex_boundary<type> :: convex_boundary(unsigned int max_nb_points);
template <class type> convex_boundary<type> :: convex_boundary(unsigned int max_nb_points)
{
	initialize();
	setup(max_nb_points);
}
#define CONVEXE_BOUNDARY_SETUP(type) template void convex_boundary<type> :: setup(unsigned int max_nb_points);
template <class type> void convex_boundary<type> :: setup(unsigned int max_nb_points)
{
	free();
	initialize();
	_max_nb_points = max_nb_points;
	_labels = new unsigned int[max_nb_points];
	string = new m_c_convexe<type>[max_nb_points];
}

#define CONVEXE_BOUNDARY_COMPUTE(type) template int convex_boundary<type> :: compute( const type * points,\
																					  unsigned int nb_points,\
																					  const unsigned char * mask );
template <class type> int convex_boundary<type> :: compute( const type * points,
															unsigned int nb_points,
															const unsigned char * mask )
{

	nb_points_r = 0;
	if (nb_points > _max_nb_points)
		return 1;
	points_ = points;
	mask_ = mask;
	nb_points_ = nb_points;
	_surface = 0;

	//Ini. de la chaine de maillons
	ini_string();
	//Tri de la chaine de maillons
	sort_string();

	//Calcul de l'enveloppe convexe
	compute_boundary();

	//Création de la liste des points
	compute_label_array();

	return 0;
}

#define CONVEXE_BOUNDARY_LABELS(type) template const unsigned int * convex_boundary<type> :: labels() const;
template <class type>  const unsigned int * convex_boundary<type> :: labels() const
{
	return _labels;
}

#define CONVEXE_BOUNDARY_NB_POINTS(type) template unsigned int convex_boundary<type> :: nb_points() const;
template <class type> unsigned int convex_boundary<type> :: nb_points() const
{
	return _nb_points;
}

#define CONVEXE_BOUNDARY_MAX_NB_POINTS(type) template unsigned int convex_boundary<type> :: max_nb_points() const;
template <class type> unsigned int convex_boundary<type> :: max_nb_points() const
{
	return _max_nb_points;
}

#define CONVEXE_BOUNDARY_DES(type) template convex_boundary<type> :: ~convex_boundary();
template <class type> convex_boundary<type> :: ~convex_boundary()
{
	free();
	initialize();
}

#define CONVEXE_BOUNDARY_INI_STRING(type) template void convex_boundary<type> :: ini_string();
template <class type> void convex_boundary<type> :: ini_string()
{
	unsigned int max_id = 0;
	for (unsigned int i = 0; i < nb_points_; ++ i)
	{
		if ( mask_[i] )
		{
			max_id = i;
			break;
		}
	}

	for (unsigned int i = 1; i < nb_points_; ++ i)
	{
		if ( mask_[i] )
		{

			if (points_[2 * i + 1] < points_[2 * max_id + 1] )
			{
				max_id = i;
			}
			else if (points_[2 * i + 1] == points_[2 * max_id + 1] )
			{
				if (points_[2 * i] < points_[2 * max_id])
					max_id = i;
			}
		}
	}
	for (unsigned int i = 0; i < nb_points_; ++ i)
	{
		if ( mask_[i] )
		{

			string[ nb_points_r].setup(	 string,
										 points_,
										 (max_id + i) % nb_points_);	 			 
			nb_points_r ++;
		}		 
	}


}

#define CONVEXE_BOUNDARY_SURFACE(type) template double convex_boundary<type> :: surface() const;
template <class type> double convex_boundary<type> :: surface() const
{
	return _surface;
}
#define CONVEXE_BOUNDARY_SORT_STRING(type) template void convex_boundary<type> :: sort_string();
template <class type> void convex_boundary<type> :: sort_string()
{
	m_c_convexe<type> * begin,
					  * end;
	begin = string + 1;
	end = string + nb_points_r;
	sort( begin, end);
	string[0].previous(end - 1);
	for (m_c_convexe<type> * p = begin; p < end; ++ p)
	{
		p->previous(p - 1);
	}
	
}

#define CONVEXE_BOUNDARY_COMPUTE_BOUNDARY(type) template void convex_boundary<type> :: compute_boundary();
template <class type> void convex_boundary<type> :: compute_boundary()
{
	m_c_convexe<type> * end;
	end = string + nb_points_r;
	for (m_c_convexe<type> * p = string + 1; p < end; ++ p)
	{
		check_point(p);
	}
}

#define CONVEXE_BOUNDARY_CHECK_POINT(type) template void convex_boundary<type> :: check_point(m_c_convexe<type> * p);
template <class type> void convex_boundary<type> :: check_point(m_c_convexe<type> * p)
{
	if (p != string)
	{
		double det;
		det = p->determinant();
		if (det > 0)
		{
			p->erase();
			check_point(p->previous());
			_surface += det / 2;
		}
		else if (det == 0)
		{
			double ps;
			ps = p->scalar_product();
			if (ps >= 0)
			{

				p->erase();
				check_point(p->previous());
				_surface += det / 2;
			}
		}
	}
}

#define CONVEXE_BOUNDARY_COMPUTE_LABEL_ARRAY(type) template void convex_boundary<type> :: compute_label_array();
template <class type> void convex_boundary<type> :: compute_label_array()
{
	_labels[0] = string->label();
	_nb_points = 1;
	
	for (m_c_convexe<type> * p = string->next(); p != string; p = p->next())
	{

		_labels[_nb_points] = p->label();
		++ _nb_points;
	}
}

#define CONVEXE_BOUNDARY_FREE(type) template void convex_boundary<type> :: free();
template <class type> void convex_boundary<type> :: free()
{
	if (string)
		delete[] string;
	if (_labels)
		delete[] _labels;
}
#define CONVEXE_BOUNDARY_INITIALIZE(type) template void convex_boundary<type> :: initialize();
template <class type> void convex_boundary<type> :: initialize()
{
	string = 0;
	_max_nb_points = 0;
	_nb_points = 0;
	_labels = 0;
	_surface = 0;
	points_ = 0;
	nb_points_ = 0;
	nb_points_r = 0;
	mask_ = 0;
}

#define CONVEXE_BOUNDARY_CLASS(type) CONVEXE_BOUNDARY_C1(type)\
									 CONVEXE_BOUNDARY_SETUP(type)\
									 CONVEXE_BOUNDARY_COMPUTE(type)\
									 CONVEXE_BOUNDARY_LABELS(type)\
									 CONVEXE_BOUNDARY_NB_POINTS(type)\
									 CONVEXE_BOUNDARY_MAX_NB_POINTS(type)\
									 CONVEXE_BOUNDARY_DES(type)\
									 CONVEXE_BOUNDARY_INI_STRING(type)\
									 CONVEXE_BOUNDARY_SORT_STRING(type)\
									 CONVEXE_BOUNDARY_COMPUTE_BOUNDARY(type)\
									 CONVEXE_BOUNDARY_CHECK_POINT(type)\
									 CONVEXE_BOUNDARY_COMPUTE_LABEL_ARRAY(type)\
									 CONVEXE_BOUNDARY_FREE(type)\
									 CONVEXE_BOUNDARY_INITIALIZE(type)\
									 CONVEXE_BOUNDARY_SURFACE(type) 
									 
CONVEXE_BOUNDARY_CLASS(double)
CONVEXE_BOUNDARY_CLASS(unsigned int)
