#include "order_function.hpp"
bool less_than(	const double & x,
				const double & y )
{
	if ( x < y )
		return true;
	return false;
}
				
bool more_than(	const double & x,
				const double & y )
{
	if ( x > y )
		return true;
	return false;
}

bool abs_more_than(	const double & x,
					const double & y )
{
	double x_,
		   y_;
		   
	if ( x < 0 )
		x_ = -x;
	else
		x_ = x;
		
	if ( y < 0 )
		y_ = -y;
	else
		y_ = y;
	
	if ( x_ > y_ )
		return true;
	return false;
}

bool abs_less_than(	const double & x,
					const double & y )
{
	double x_,
		   y_;
		   
	if ( x < 0 )
		x_ = -x;
	else
		x_ = x;
		
	if ( y < 0 )
		y_ = -y;
	else
		y_ = y;
	
	if ( x_ < y_ )
		return true;
	return false;
}
