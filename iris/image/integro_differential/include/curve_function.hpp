
#ifndef _CURVE_FUNCTION_HPP_
#define _CURVE_FUNCTION_HPP_

#include <cmath>
#include <complex>
using namespace std;
/**@fn
 * @brief
 * Ligne
 * y = a x + b
 * 1 -> a
 * 2 -> b
 */
double line(	const double & x,
				const double * params );
		
/**@fn
 * @brief
 * Ligne
 * y = a x + b
 * 
 */
double line_derivate(	const double & x,
						const double * params );
						
/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 1 -> a
 * 2 -> x0
 * 3 -> b
 */
double parabol(	const double & x,
				const double * params );
						
double parabol_mod(	const double & x,
					const double * params );

/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 
 */
double parabol_derivate(	const double & x,
							const double * params );
double parabol_mod_derivate(	const double & x,
								const double * params );
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
				const double * params );
						
/**@fn
 * @brief
 * Parabole
 * y = a (x - x0)² + b
 * 
 */
double special_derivate(	const double & x,
							const double * params );
	
	
	
	
	
	
	
	
double x_circle_coord (	const double & x,
							const double * params);
	
double y_circle_coord (	const double & x,
							const double * params );	
	
double x_circle_derivate (	const double & x,
							const double * params);
	
double y_circle_derivate (	const double & x,
							const double * params );

	
double x_ellipse_coord (	const double & x,
							const double * params);
	
double y_ellipse_coord (	const double & x,
							const double * params );	
	
double x_ellipse_derivate (	const double & x,
								const double * params);
	
double y_ellipse_derivate (	const double & x,
								const double * params );
	
	
	
	
	
	
	
	
	
	
						
						
#endif
