#======================================================================#
#                                 Options                              #
#======================================================================#

#Description
DESC="Pairwise Kalman filter programs"

#Auteur	
AUTHOR="Valérian Némesin"

#Nom de la bibliothèque
NAME=tkalman


#Dépendances
#Lib. statiques
S_LIBS=../../API ../image ../tkalman/tkalman_nc ../iris ../config

#Lib. dynamiques
D_LIBS=gsl opencv fftw3

#Options principales
M_OPTIONS=../../compilation.mk



