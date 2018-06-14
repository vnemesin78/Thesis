#include "compute_contour_4c.hpp"
#include "get_first.hpp"
#include "is_mask.hpp"
#include "dir4c.hpp"
#include <iostream>
using namespace std;
template <class T> int compute_contour_4c(unsigned int * contour,
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
		(*pNbPointsContour) = 0;
		return 1;
	}
	//Recherche du premier point du masque
	if (get_first(contour[0],
				  contour[1],
				  pixels,
				  width,
				  height,
				  width_step,
				  value))
	{
		(*pNbPointsContour) = 0;
		return 0;		
	}
	(*pNbPointsContour) = 1;
	//Recherche de la première direction de recherche
    int dir0 = 0;
    unsigned int suivant[2];
    for(; dir0 < 4; dir0++)
    {
		if (! ((xDir4c[dir0] < 0 && contour[0] == 0) || (yDir4c[dir0] < 0 && contour[1] == 0) ) ) //Test de la validité de la direction
		{
			suivant[0] = contour[0] + xDir4c[dir0];
			suivant[1] = contour[1] + yDir4c[dir0];
			if (is_mask(suivant[0],
						suivant[1],
						pixels,
						width,
						height,
						width_step,
						value))
				break;	
		}

    }
    if (dir0 == 4) //Si c'est un point isolé!
        return 0;
    //Sauvegarde du second point
    contour[2] = suivant[0];
    contour[3] = suivant[1];
	(*pNbPointsContour) = 2;
    //Calcul du contour
    int dir = (dir0 + 3) % 4;
    while(true)
    {
		//Parcours de toutes les directions en laissant le demi tour en dernier
		for(int i = 0; i < 4; i++)
		{
			int nPoint = ((*pNbPointsContour) - 1) * 2;
			if (! ((xDir4c[dir] < 0 && contour[nPoint] == 0) || (yDir4c[dir] < 0 && contour[nPoint + 1] == 0) ) ) //Test de la validité de la direction
			{
				suivant[0] = contour[nPoint] + xDir4c[dir];
				suivant[1] = contour[nPoint + 1] + yDir4c[dir];
				if (is_mask(suivant[0],
							suivant[1],
							pixels,
							width,
							height,
							width_step,
							value))
					break;	
			}
			dir = (dir + 1) % 4;
		}
		//Sauvegarde du point
		//Test mémoire
		if ((*pNbPointsContour) == maxPointsContour)
		{
			cerr << dir0 << endl;
			cerr << "Warning : Le nombre de points maximal du contour a été atteint!" << endl;
			return -1;
		}
		//Recopie du point
		{
			unsigned int nPoint = (*pNbPointsContour) * 2;
			contour[nPoint] = suivant[0];
			contour[nPoint + 1] = suivant[1];
			(*pNbPointsContour) ++;
		}
		//Test d'arrêt
		if (contour[2] == suivant[0] && contour[3] == suivant[1] && dir == dir0)
			break;
		//Changement de direction
		dir = (dir + 3) % 4;
	}   
	(*pNbPointsContour) -= 2;//Suppression des deux points en doublon.

    return 0;    
}

COMPUTE_CONTOUR_4C(char)
COMPUTE_CONTOUR_4C(unsigned char)
COMPUTE_CONTOUR_4C(int)
COMPUTE_CONTOUR_4C(unsigned int)
