#include "m_c_convexe.hpp"
#include <iostream>
#include <cmath>
using namespace std;
#define M_C_CONVEXE_CONST(type) template m_c_convexe<type> :: m_c_convexe();
template <class type> m_c_convexe<type> :: m_c_convexe()
{
	_point = 0;
	_label = 0;
	_next = 0;
	_previous = 0;
	_first = 0;
}

#define M_C_CONVEXE_CONSTRUC(type) template m_c_convexe<type> :: m_c_convexe( const m_c_convexe<type> * first,\
																			  const type * point,\
																			  unsigned int label);
template <class type> m_c_convexe<type> :: m_c_convexe(const m_c_convexe<type> * first,
													   const type * point,
													   unsigned int label)
{
	_first = first;
	_point = point + 2 * label;
	_label = label;
	_next = 0;
	_previous = 0;
}


#define M_C_CONVEXE_SETUP(type) template void m_c_convexe<type> :: setup( const m_c_convexe<type> * first,\
																		  const type * point,\
																		  unsigned int label); 
template <class type> void m_c_convexe<type> :: setup( const m_c_convexe<type> * first,
													   const type * point,
													   unsigned int label)
{
	_first = first;
	_point = point + 2 * label;
	_label = label;
}

#define M_C_CONVEXE_NEXT(type) template m_c_convexe<type> * m_c_convexe<type> :: next(m_c_convexe<type> * next);
template <class type> m_c_convexe<type> * m_c_convexe<type> :: next(m_c_convexe<type> * next)
{
	if (next)
	{
		_next = next;
		_next->_previous = this;
	}
	return _next;
}

#define M_C_CONVEXE_PREVIOUS(type) template m_c_convexe<type> * m_c_convexe<type> :: previous(m_c_convexe<type> * previous);
template <class type> m_c_convexe<type> * m_c_convexe<type> :: previous(m_c_convexe<type> * previous)
{
	if (previous)
	{
		_previous = previous;
		_previous->_next = this;
	}
	return _previous;
}

#define M_C_CONVEXE_POINT(type) template const type * m_c_convexe<type> :: point() const;
template <class type> const type * m_c_convexe<type> :: point() const
{
	
	return _point;
}
#define M_C_CONVEXE_LABEL(type) template unsigned int m_c_convexe<type> :: label() const;
template <class type> unsigned int m_c_convexe<type> :: label() const
{
	return _label;
}

#define M_C_CONVEXE_COSINUS(type) template double m_c_convexe<type> :: cosinus() const;
template <class type> double m_c_convexe<type> :: cosinus() const
{
	if ( _first == NULL )
		return -1;
		
	double r = radius();
	if (r == 0)
		return 1;
	else
		return ( ( ( (double) _point[0] ) - ( (double) _first->_point[0] ) ) / r );
}

#define M_C_CONVEXE_SCALAR_PRODUCT(type) template double m_c_convexe<type> :: scalar_product() const;
template <class type> double m_c_convexe<type> :: scalar_product() const
{
	double l_dx,
		   l_dy,
		   r_dx,
		   r_dy;
		   
	if (_previous == 0 || _next == 0)
		return 0;
	
	l_dx = ( (double) _point[0] ) - ( (double) _previous->_point[0] );
	l_dy = ( (double) _point[1] ) - ( (double) _previous->_point[1] );
	
	r_dx = ( (double) _point[0] ) - ( (double) _next->_point[0] );
	r_dy = ( (double) _point[1] ) - ( (double) _next->_point[1] );
	
	return (l_dx * r_dx + l_dy * r_dy);
}

#define M_C_CONVEXE_SCALAR_DET(type) template double m_c_convexe<type> :: determinant() const;
template <class type> double m_c_convexe<type> :: determinant()  const
{
	double l_dx,
		   l_dy,
		   r_dx,
		   r_dy;
		   
	if (_previous == 0 || _next == 0)
		return 0;
	
	l_dx = ( (double) _point[0] ) - ( (double) _previous->_point[0] );
	l_dy = ( (double) _point[1] ) - ( (double) _previous->_point[1] );
	
	r_dx = ( (double) _point[0] ) - ( (double) _next->_point[0] );
	r_dy = ( (double) _point[1] ) - ( (double) _next->_point[1] );
	
	return (l_dx * r_dy - r_dx * l_dy);
}
#define M_C_CONVEXE_SCALAR_ERASE(type) template void m_c_convexe<type> :: erase();
template <class type> void m_c_convexe<type> :: erase()
{
	_next->previous(_previous);
}
#define M_C_CONVEXE_RADIUS(type) template double m_c_convexe<type> :: radius() const;
template <class type> double m_c_convexe<type> :: radius() const
{
	if ( _first == NULL )
		return 0;
	double dx, 
		   dy;
	dx = ( (double) _point[0] ) - ( (double) _first->_point[0] );
	dy = ( (double) _point[1] ) - ( (double) _first->_point[1] );
	return sqrt(dx * dx + dy * dy);
}
#define M_C_CONVEXE_SUP(type) template bool m_c_convexe<type> :: operator<(const m_c_convexe<type> & m) const;
template <class type> bool m_c_convexe<type> :: operator<(const m_c_convexe<type> & m) const
{
	double cos1, cos2;

		

	cos1 = cosinus();
	cos2 = m.cosinus();
	if ( cos1 == cos2 ) 
	{
		return ( radius() < m.radius() );
	}
	else
	{
		return (cos1 > cos2);
	}
	
}
#define M_C_CONVEXE_CLASS(type) M_C_CONVEXE_CONST(type)\
								M_C_CONVEXE_CONSTRUC(type)\
								M_C_CONVEXE_SETUP(type)\
								M_C_CONVEXE_NEXT(type)\
								M_C_CONVEXE_PREVIOUS(type)\
								M_C_CONVEXE_POINT(type)\
								M_C_CONVEXE_LABEL(type)\
								M_C_CONVEXE_COSINUS(type)\
								M_C_CONVEXE_SCALAR_PRODUCT(type)\
								M_C_CONVEXE_SCALAR_DET(type)\
								M_C_CONVEXE_SCALAR_ERASE(type)\
								M_C_CONVEXE_RADIUS(type)\
								M_C_CONVEXE_SUP(type)
								
//Définition de la classe pour type == double et unsigned int
M_C_CONVEXE_CLASS(unsigned int)
M_C_CONVEXE_CLASS(double)
