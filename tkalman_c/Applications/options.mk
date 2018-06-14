#======================================================================#
#                                 Options                              #
#======================================================================#

#Répertoire des applications
APP_DIR=/home/vnemesin/Experiences/Executables/

#Description
DESC="Pairwise Kalman filter programs"

#Auteur	
AUTHOR="Valérian Némesin"

#Nom de la bibliothèque
NAME=tkalman


#Dépendances
#Lib. statiques
S_LIBS=../../API ../PKF ../API ../Simulation

#Lib. dynamiques
D_LIBS=gsl

#Options principales
M_OPTIONS=../../compilation.mk



