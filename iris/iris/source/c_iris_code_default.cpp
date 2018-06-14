#include "c_iris_code.hpp"
#include "iris_default.hpp"
int c_iris_code :: default_setup( )
{
	unsigned int nb_directions;
	IRIS_CODE__NB_DIRECTIONS( nb_directions);
	
	unsigned int nb_samples;
	IRIS_CODE__NB_SAMPLES( nb_samples );
	
	double wavelenght;
	IRIS_CODE__WAVELENGHT( wavelenght );
	
	double sigma;
	IRIS_CODE__SIGMA( sigma );
	
	return setup( nb_directions, nb_samples, wavelenght, sigma );
	
}
