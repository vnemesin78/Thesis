#include "gibbs_pdf.hpp"
#include <cmath>
double gibbs_pdf( 	const double & _r,
					const double & r,
					const double & r_,
					const void * params )
{
	gibbs_params * data = (gibbs_params*) params;
	
	double dr_ = r - r_,
			_dr = r - _r;
	
	return (data->a * exp( - data->delta * ( dr_ * dr_ + _dr * _dr ) ) );
}

double log_gibbs_pdf( 	const double & _r,
						const double & r,
						const double & r_,
						const void * params )
{
	gibbs_params * data = (gibbs_params*) params;
	
	double dr_ = r - r_,
			_dr = r - _r;
	
	return (log( data->a )  - data->delta * ( dr_ * dr_ + _dr * _dr ) );
}
