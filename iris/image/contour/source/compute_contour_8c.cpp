#include "compute_contour_8c.hpp"
#include "get_first.hpp"
#include "is_mask.hpp"
#include "dir8c.hpp"
#include <iostream>
using namespace std;
template <class T> int compute_contour_8c(unsigned int * contour,
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
    for(; dir0 < 8; dir0++)
    {
		if (! ((xDir8c[dir0] < 0 && contour[0] == 0) || (yDir8c[dir0] < 0 && contour[1] == 0) ) ) //Test de la validité de la direction
		{
			suivant[0] = contour[0] + xDir8c[dir0];
			suivant[1] = contour[1] + yDir8c[dir0];
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
    if (dir0 == 8) //Si c'est un point isolé!
        return 0;
    //Sauvegarde du second point
    contour[2] = suivant[0];
    contour[3] = suivant[1];
	(*pNbPointsContour) = 2;
    //Calcul du contour
    int dir = (dir0 + 5) % 8;
    while(true)
    {
		int i;
		//Parcours de toutes les directions en laissant le demi tour en dernier
		for(i = 0; i < 8; i++)
		{
			
			int nPoint = ((*pNbPointsContour) - 1) * 2;
			if (! ((xDir8c[dir] < 0 && contour[nPoint] == 0) || (yDir8c[dir] < 0 && contour[nPoint + 1] == 0) ) ) //Test de la validité de la direction
			{
				suivant[0] = contour[nPoint] + xDir8c[dir];
				suivant[1] = contour[nPoint + 1] + yDir8c[dir];
				if (is_mask(suivant[0],
							suivant[1],
							pixels,
							width,
							height,
							width_step,
							value))
					break;	
			}
			dir = (dir + 1) % 8; 
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
		dir = (dir + 5) % 8;
	}   
	(*pNbPointsContour) -= 2;//Suppression des deux points en doublon.
    return 0;    
}

COMPUTE_CONTOUR_8C(char)
COMPUTE_CONTOUR_8C(unsigned char)
COMPUTE_CONTOUR_8C(int)
COMPUTE_CONTOUR_8C(unsigned int)
