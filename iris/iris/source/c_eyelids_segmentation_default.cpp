#include "c_eyelids_segmentation.hpp"
#include "iris_default.hpp"
int c_eyelids_segmentation :: default_setup(  	unsigned int width,
												unsigned int height )
{
	
	unsigned int nb_samples;
		EYELIDS__NB_SAMPLES( nb_samples );
		
	unsigned int nb_directions;
		EYELIDS__NB_DIRECTIONS( nb_directions );
	
	unsigned int nb_iter_gem;
		EYELIDS__NB_ITER_GEM( nb_iter_gem );
		
	unsigned int nb_iter_icm;
		EYELIDS__NB_ITER_ICM( nb_iter_icm );
	
	unsigned int nb_iter_em;
		EYELIDS__NB_ITER_EM( nb_iter_em );
	
	double delta;
		EYELIDS__DELTA(delta);
	unsigned int closing;
		EYELIDS__CLOSING( closing );
	unsigned int opening;
		EYELIDS__OPENING( opening );

	int q = 0;
		q = setup( 	 width, height, nb_samples, nb_directions,
					nb_iter_gem, nb_iter_icm, nb_iter_em, delta, closing, opening );
	if ( q )
		return q;
		
	//Lib. mémoire
	
	return 0;
}
