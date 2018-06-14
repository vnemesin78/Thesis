#include "digitalize_euclide_contour.hpp"
#include <iostream>
#include "get_perimeter.hpp"
#include "get_distance.hpp"
using namespace std;
#define TMP_DIGITALIZE_EUCLIDE_CONTOUR(type) \
template int digitalize_euclide_contour(	double * dContour,\
											unsigned int nbPointsDContour,\
											const type * contour,\
											unsigned int nbPointsContour);

template <class type> int digitalize_euclide_contour(	double * dContour,
														unsigned int nbPointsDContour,
														const type * contour,
														unsigned int nbPointsContour)
{
	//Cas problématiques
	if ((!nbPointsContour) || (!nbPointsDContour))
		return 1;
	//Calcul du périmètre
	double perimeter = get_perimeter(contour,
								     nbPointsContour);
	//Sauvegarde du premier point
	dContour[0] = contour[0];
	dContour[1] = contour[1];
	
	//Réchantionnage
	double pas = perimeter / nbPointsDContour;
	double distance = pas;
	
	double lCDistance = 0;
	double cDistance = 0;
	
	unsigned int nPoint = 0;
	unsigned int nDPoint = 2;
	
	while (distance < perimeter)
	{
		while(cDistance < distance)
		{

			nPoint += 1;
			lCDistance = cDistance;
			cDistance += get_distance(contour[((nPoint - 1) % nbPointsContour) * 2],
									  contour[((nPoint - 1) % nbPointsContour) * 2 + 1],
									  contour[((nPoint) % nbPointsContour) * 2],
									  contour[((nPoint) % nbPointsContour) * 2 + 1]);	  
		}
		//Calcul du point intermédiaire (x = x0 + (x1 - x0) * distance(m,m0) / distance(m0,m1)
		{
			double r = (distance - lCDistance) / (cDistance - lCDistance);
			type x0 = contour[((nPoint - 1) % nbPointsContour) * 2];
			type y0 = contour[((nPoint - 1) % nbPointsContour) * 2 + 1];
			type x1 = contour[((nPoint) % nbPointsContour) * 2];
			type y1 = contour[((nPoint) % nbPointsContour) * 2 + 1];
			dContour[nDPoint] = x0 + ((double) x1 - x0) * r;
			dContour[nDPoint  + 1] = y0 + ((double) y1 - y0) * r;
		}
		nDPoint += 2;
		distance += pas;
	}
	return 0;
}
TMP_DIGITALIZE_EUCLIDE_CONTOUR(unsigned char)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(char)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(unsigned short int)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(short int)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(unsigned int)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(int)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(unsigned long int)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(long int)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(float)
TMP_DIGITALIZE_EUCLIDE_CONTOUR(double)


