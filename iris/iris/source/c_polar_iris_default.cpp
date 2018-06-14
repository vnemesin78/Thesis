#include "c_polar_iris.hpp"
#include "iris_default.hpp"
int c_polar_iris :: default_setup ( unsigned int width,
									unsigned int height )
{
	unsigned int nb_directions;
	POLAR__NB_DIRECTIONS(nb_directions);
	
	unsigned int nb_samples;
	POLAR__NB_SAMPLES(nb_samples);
	
	unsigned int nb_samples_iris;
	POLAR__NB_SAMPLES_IRIS(nb_samples_iris);
	
	
	
	return ( setup ( 	nb_directions,
						nb_samples,
						nb_samples_iris,
						width,
						height ) );
}
