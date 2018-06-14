#include "compute_contour_4c_nb.hpp"
#include "get_first.hpp"
#include "is_mask.hpp"
#include "dir4c.hpp"
#include <iostream>
using namespace std;
template <class T> int compute_contour_4c_nb( unsigned int * contour,
											  unsigned int * pNbPointsContour,
											  unsigned int maxPointsContour,
											  const T * pixels,
											  unsigned int width,
											  unsigned int height,
											  unsigned int width_step,
											  const T & value)
{
	//Gestion d'un contour vide
    if (!maxPointsContour)
    {
		*pNbPointsContour = 0;
		return 1;
	}
	unsigned courant[2];
	//Recherche du premier point du masque
*pNbPointsContour = 0;
	if (get_first(courant[0],
				  courant[1],
				  pixels,
				  width,
				  height,
				  width_step,
				  value))
	{
		*pNbPointsContour = 0;
		return 0;		
	}

	unsigned int suivant[2];
	int dir = 3;
	
	while(true)
    {
        for(int i = 0; i < 4; i++)
        {
            suivant[0] = courant[0] + xDir4c[dir];
            suivant[1] = courant[1] + yDir4c[dir];
            //On teste si le pixel est dans l'image
			if (is_mask(suivant[0],
						suivant[1],
						pixels,
						width,
						height,
						width_step,
						value))
				break;
            //Sauvegarde du point
            if ((*pNbPointsContour) == maxPointsContour)
            {
				cerr << "Warning : Le nombre de points maximal du contour a été atteint!" << endl;
				return -1;
            }

			unsigned int nPoint = (*pNbPointsContour) * 2;
            contour[nPoint] = courant[0] + x4c[dir];
            contour[nPoint + 1] = courant[1] + y4c[dir];
            if (contour[nPoint] == contour[0] && contour[nPoint + 1] == contour[1] && (*pNbPointsContour) > 1) //Condition d'arrêt
                return 0;
            (*pNbPointsContour) ++;

            dir += 1;
            dir %= 4;
        }
        //Prochaine direction de recherche
        dir += 3;
        dir %= 4;
        courant[0] = suivant[0];
        courant[1] = suivant[1];
    }
    return 0;
}
COMPUTE_CONTOUR_4C_NB(char)
COMPUTE_CONTOUR_4C_NB(unsigned char)
COMPUTE_CONTOUR_4C_NB(unsigned int)
COMPUTE_CONTOUR_4C_NB(int)

