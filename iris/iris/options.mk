#======================================================================#
#                                 Options                              #
#======================================================================#

#Nom de la bibliothèque
NAME=lib_iris

#Description
DESC="API for Pairwise Kalman filter algorithmes"

#Auteur
AUTHOR="Valérian Némesin"

#Dépendances
#Lib. statiques
S_LIBS=../image ../../API ../tkalman/tkalman_nc ../config
#Lib. dynamiques
D_LIBS=opencv gsl fftw3		

#Options principales
M_OPTIONS=../../compilation.mk



