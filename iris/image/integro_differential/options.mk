#======================================================================#
#                                 Options                              #
#======================================================================#

#Nom de la bibliothèque
NAME=lib_int_dif

#Description
DESC="API for Pairwise Kalman filter algorithmes"

#Auteur
AUTHOR="Valérian Némesin"

#Dépendances
#Lib. statiques
S_LIBS=../polar ../filters ../interpol_2d ../utilities
#Lib. dynamiques
D_LIBS=fftw3

#Options principales
M_OPTIONS=../../../compilation.mk




