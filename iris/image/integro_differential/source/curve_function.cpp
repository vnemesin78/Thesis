#include "curve_function.hpp"
/**@fn
 * @brief
 * Ligne
 * y = a x + b
 * 1 -> a
 * 2 -> b
 */
double line(	const double & x,
				const double * params )
{
	return ( params[0] * x + params[1] );
}

/**@fn
 * @brief
 * Ligne
 * y = a x + b
 * 
 */
double line_derivate(	const double & x,
						const double * params )
{
	return ( params[0] );
}

/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 1 -> a
 * 2 -> x0
 * 3 -> b
 */
double parabol(	const double & x,
				const double * params )
{
	double dx = x - params[1];
	return ( params[0] * dx * dx + params[2]);
}

/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 
 */
double parabol_derivate(	const double & x,
							const double * params )
{
	double dx = x - params[1];
	return ( 2 * params[0] * dx );
}

/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 1 -> a
 * 2 -> x0
 * 3 -> b
 */
double parabol_mod(	const double & x,
					const double * params )
{
	double dx = fmod( x, params[3] ) - params[1];
	return ( params[0] * dx * dx + params[2]);
}

/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 
 */
double parabol_mod_derivate(	const double & x,
								const double * params )
{
	double dx = fmod( x, params[3] ) - params[1];
	return ( 2 * params[0] * dx );
}

/**@fn
 * @brief
 * Spécial
 * y = a (x - x0)^c + b
 * 1 -> a
 * 2 -> x0
 * 3 -> b
 * 4 -> c
 */
double special(	const double & x,
					const double * params )
{
	double dx = x - params[1];
	return ( params[0] * pow( abs(dx),  params[3]) + params[2]);
}

/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 
 */
double special_derivate(	const double & x,
							const double * params )
{
	double dx = x - params[1];
	return ( params[3] * params[0] * pow( abs(dx), params[3] - 1 ) );
}

double x_ellipse_coord ( 	const double & x,
							const double * params )
{
	return ( params[0] + params[2] * cos(x) * cos(params[4]) - params[3] * sin(x) * sin ( params[4] ) );	
}

double y_ellipse_coord ( 	const double & x,
							const double * params )
{
	return ( params[1] + params[2] * cos(x) * sin(params[4]) + params[3] * sin(x) * cos ( params[4] ) );	
}

double x_ellipse_derivate ( 	const double & x,
								const double * params )
{
	return ( - params[2] * sin(x) * cos(params[4]) - params[3] * cos(x) * sin ( params[4] ) );	
}

double y_ellipse_derivate ( 	const double & x,
								const double * params )
{
	return ( - params[2] * sin(x) * sin(params[4]) + params[3] * cos(x) * cos ( params[4] ) );	
}

double x_circle_coord (	const double & x,
							const double * params )
{
	return ( params[0] + params[2] * cos(x) );	
}

double y_circle_coord (	const double & x,
							const double * params )
{
	return ( params[1] + params[2] * sin(x) );	
}

double x_circle_derivate (	const double & x,
							const double * params )
{
	return ( - params[2] * cos(x) );	
}

double y_circle_derivate (	const double & x,
							const double * params )
{
	return ( params[2] * cos(x) );	
}
















