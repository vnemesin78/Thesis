#ifndef IRIS_DEFAULT_HPP
#define IRIS_DEFAULT_HPP
#define DATA_TYPE(var)\
var = 1;
#define EYELIDS__CLOSING(var)\
var = 3;
#define EYELIDS__DELTA(var)\
var = 0.3;
#define EYELIDS__NB_DIRECTIONS(var)\
var = 48;
#define EYELIDS__NB_ITER_EM(var)\
var = 20;
#define EYELIDS__NB_ITER_GEM(var)\
var = 10;
#define EYELIDS__NB_ITER_ICM(var)\
var = 10;
#define EYELIDS__NB_SAMPLES(var)\
var = 60;
#define EYELIDS__OPENING(var)\
var = 7;
#define FUSION__NB_IMAGES(var)\
var = 15;
#define FUSION__SAVINGS(var)\
var = gsl_matrix_calloc ( 1,7);\
var->data[0 * var->tda + 0] = 1;var->data[0 * var->tda + 1] = 3;var->data[0 * var->tda + 2] = 5;var->data[0 * var->tda + 3] = 7;var->data[0 * var->tda + 4] = 9;var->data[0 * var->tda + 5] = 11;var->data[0 * var->tda + 6] = 13;
#define FUSION__TYPE(var)\
var = 0;
#define HEIGHT(var)\
var = 480;
#define IRIS__DX(var)\
var = 5;
#define IRIS__DX_ELLIPSE(var)\
var = 1;
#define IRIS__DY(var)\
var = 5;
#define IRIS__DY_ELLIPSE(var)\
var = 1;
#define IRIS__KERNEL_SIZE(var)\
var = 11;
#define IRIS__NB_DIRECTIONS(var)\
var = 64;
#define IRIS__NB_SAMPLES(var)\
var = 48;
#define IRIS__NB_SAMPLES_ELLIPSE(var)\
var = 1;
#define IRIS__NB_THETAS(var)\
var = 1;
#define IRIS__NB_X(var)\
var = 10;
#define IRIS__NB_X_ELLIPSE(var)\
var = 1;
#define IRIS__NB_Y(var)\
var = 10;
#define IRIS__NB_Y_ELLIPSE(var)\
var = 1;
#define IRIS__R_MAX(var)\
var = 150;
#define IRIS__R_MIN(var)\
var = 85;
#define IRIS__R_RATIO(var)\
var = 0.05;
#define IRIS__SPOT_THRESHOLD(var)\
var = 160;
#define IRIS__VALIDNESS_FACTOR(var)\
var = 0.4;
#define IRIS_CODE__NB_DIRECTIONS(var)\
var = 240;
#define IRIS_CODE__NB_SAMPLES(var)\
var = 20;
#define IRIS_CODE__SIGMA(var)\
var = 0.5;
#define IRIS_CODE__WAVELENGHT(var)\
var = 18;
#define MATCHING__DISTANCE(var)\
var = new char[8];\
memcpy( var,"Hamming",8);
#define POLAR__NB_DIRECTIONS(var)\
var = 240;
#define POLAR__NB_SAMPLES(var)\
var = 60;
#define POLAR__NB_SAMPLES_IRIS(var)\
var = 40;
#define PRE_PROCESSING__EROSION(var)\
var = 5;
#define PRE_PROCESSING__MEDIAN(var)\
var = 5;
#define PRE_PROCESSING__OPENING(var)\
var = 9;
#define PUPIL__CRITERIA(var)\
var = 0.25;
#define PUPIL__NAME(var)\
var = new char[6];\
memcpy( var,"pupil",6);
#define PUPIL__NB_ITER_ELLIPSE(var)\
var = 1;
#define PUPIL__NB_MIN_MAX(var)\
var = 5;
#define PUPIL__NB_PTS_PUPIL(var)\
var = 32;
#define PUPIL__NB_THRESHOLDS(var)\
var = 5;
#define PUPIL__RADIUS_MAX(var)\
var = 80;
#define PUPIL__RADIUS_MIN(var)\
var = 13;
#define PUPIL__SIGMA(var)\
var = 5;
#define SCORE__C(var)\
var = 0.2;
#define SCORE__HEIGHT(var)\
var = 150;
#define SCORE__KERNEL(var)\
var = gsl_matrix_calloc ( 5,5);\
var->data[0 * var->tda + 0] = -1;var->data[0 * var->tda + 1] = -1;var->data[0 * var->tda + 2] = -1;var->data[0 * var->tda + 3] = -1;var->data[0 * var->tda + 4] = -1;var->data[1 * var->tda + 0] = -1;var->data[1 * var->tda + 1] = -1;var->data[1 * var->tda + 2] = 4;var->data[1 * var->tda + 3] = -1;var->data[1 * var->tda + 4] = -1;var->data[2 * var->tda + 0] = -1;var->data[2 * var->tda + 1] = 4;var->data[2 * var->tda + 2] = 4;var->data[2 * var->tda + 3] = 4;var->data[2 * var->tda + 4] = -1;var->data[3 * var->tda + 0] = -1;var->data[3 * var->tda + 1] = -1;var->data[3 * var->tda + 2] = 4;var->data[3 * var->tda + 3] = -1;var->data[3 * var->tda + 4] = -1;var->data[4 * var->tda + 0] = -1;var->data[4 * var->tda + 1] = -1;var->data[4 * var->tda + 2] = -1;var->data[4 * var->tda + 3] = -1;var->data[4 * var->tda + 4] = -1;
#define SCORE__TYPE(var)\
var = 1;
#define SCORE__WIDTH(var)\
var = 150;
#define TRACKING(var)\
var = 0;
#define TRACKING__F(var)\
var = gsl_matrix_calloc ( 6,6);\
var->data[0 * var->tda + 0] = 1.6729;var->data[0 * var->tda + 1] = 0.108778;var->data[0 * var->tda + 2] = 0.020075;var->data[0 * var->tda + 3] = -0.705716;var->data[0 * var->tda + 4] = -0.101888;var->data[0 * var->tda + 5] = -0.019215;var->data[1 * var->tda + 0] = 0.05462;var->data[1 * var->tda + 1] = 1.55322;var->data[1 * var->tda + 2] = 0.711466;var->data[1 * var->tda + 3] = -0.055986;var->data[1 * var->tda + 4] = -0.627548;var->data[1 * var->tda + 5] = -0.773314;var->data[2 * var->tda + 0] = 0.015382;var->data[2 * var->tda + 1] = -0.01763;var->data[2 * var->tda + 2] = 0.568909;var->data[2 * var->tda + 3] = -0.013123;var->data[2 * var->tda + 4] = 0.013781;var->data[2 * var->tda + 5] = 0.428347;var->data[3 * var->tda + 0] = 1;var->data[3 * var->tda + 1] = 0;var->data[3 * var->tda + 2] = 0;var->data[3 * var->tda + 3] = 0;var->data[3 * var->tda + 4] = 0;var->data[3 * var->tda + 5] = 0;var->data[4 * var->tda + 0] = 0;var->data[4 * var->tda + 1] = 1;var->data[4 * var->tda + 2] = 0;var->data[4 * var->tda + 3] = 0;var->data[4 * var->tda + 4] = 0;var->data[4 * var->tda + 5] = 0;var->data[5 * var->tda + 0] = 0;var->data[5 * var->tda + 1] = 0;var->data[5 * var->tda + 2] = 1;var->data[5 * var->tda + 3] = 0;var->data[5 * var->tda + 4] = 0;var->data[5 * var->tda + 5] = 0;
#define TRACKING__NAME(var)\
var = new char[9];\
memcpy( var,"Tracking",9);
#define TRACKING__N_SIGMA(var)\
var = 9;
#define TRACKING__SIZE_X(var)\
var = 4;
#define TRACKING__SQRT_Q(var)\
var = gsl_matrix_calloc ( 6,6);\
var->data[0 * var->tda + 0] = -7.13678;var->data[0 * var->tda + 1] = 0.092066;var->data[0 * var->tda + 2] = 0.056345;var->data[0 * var->tda + 3] = 0;var->data[0 * var->tda + 4] = 0;var->data[0 * var->tda + 5] = 0;var->data[1 * var->tda + 0] = 0;var->data[1 * var->tda + 1] = -5.95137;var->data[1 * var->tda + 2] = -0.088354;var->data[1 * var->tda + 3] = 0;var->data[1 * var->tda + 4] = 0;var->data[1 * var->tda + 5] = 0;var->data[2 * var->tda + 0] = 0;var->data[2 * var->tda + 1] = 0;var->data[2 * var->tda + 2] = -1.50178;var->data[2 * var->tda + 3] = 0;var->data[2 * var->tda + 4] = 0;var->data[2 * var->tda + 5] = 0;var->data[3 * var->tda + 0] = 0;var->data[3 * var->tda + 1] = 0;var->data[3 * var->tda + 2] = 0;var->data[3 * var->tda + 3] = 0;var->data[3 * var->tda + 4] = -0;var->data[3 * var->tda + 5] = 0;var->data[4 * var->tda + 0] = 0;var->data[4 * var->tda + 1] = 0;var->data[4 * var->tda + 2] = 0;var->data[4 * var->tda + 3] = 0;var->data[4 * var->tda + 4] = 0;var->data[4 * var->tda + 5] = 0;var->data[5 * var->tda + 0] = 0;var->data[5 * var->tda + 1] = 0;var->data[5 * var->tda + 2] = 0;var->data[5 * var->tda + 3] = 0;var->data[5 * var->tda + 4] = 0;var->data[5 * var->tda + 5] = 0;
#define TRACKING__SQRT_Q_0(var)\
var = gsl_matrix_calloc ( 6,6);\
var->data[0 * var->tda + 0] = 100;var->data[0 * var->tda + 1] = 0;var->data[0 * var->tda + 2] = 0;var->data[0 * var->tda + 3] = 100;var->data[0 * var->tda + 4] = 0;var->data[0 * var->tda + 5] = 0;var->data[1 * var->tda + 0] = 0;var->data[1 * var->tda + 1] = 100;var->data[1 * var->tda + 2] = 0;var->data[1 * var->tda + 3] = 0;var->data[1 * var->tda + 4] = 100;var->data[1 * var->tda + 5] = 0;var->data[2 * var->tda + 0] = 0;var->data[2 * var->tda + 1] = 0;var->data[2 * var->tda + 2] = 10;var->data[2 * var->tda + 3] = 0;var->data[2 * var->tda + 4] = 0;var->data[2 * var->tda + 5] = 10;var->data[3 * var->tda + 0] = 0;var->data[3 * var->tda + 1] = 0;var->data[3 * var->tda + 2] = 0;var->data[3 * var->tda + 3] = 0;var->data[3 * var->tda + 4] = 0;var->data[3 * var->tda + 5] = 0;var->data[4 * var->tda + 0] = 0;var->data[4 * var->tda + 1] = 0;var->data[4 * var->tda + 2] = 0;var->data[4 * var->tda + 3] = 0;var->data[4 * var->tda + 4] = 0;var->data[4 * var->tda + 5] = 0;var->data[5 * var->tda + 0] = 0;var->data[5 * var->tda + 1] = 0;var->data[5 * var->tda + 2] = 0;var->data[5 * var->tda + 3] = 0;var->data[5 * var->tda + 4] = 0;var->data[5 * var->tda + 5] = 0;
#define TRACKING__T_0(var)\
var = gsl_matrix_calloc ( 6,1);\
var->data[0 * var->tda + 0] = 0;var->data[1 * var->tda + 0] = 0;var->data[2 * var->tda + 0] = 40;var->data[3 * var->tda + 0] = 0;var->data[4 * var->tda + 0] = 0;var->data[5 * var->tda + 0] = 40;
#define WIDTH(var)\
var = 640;
#endif
