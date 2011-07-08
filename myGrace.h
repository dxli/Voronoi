#pragma once
#ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#  define EXIT_FAILURE -1
#endif
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <time.h>
#include <grace_np.h>

using namespace std;
void my_error_function(const char *);
double axis_normal(double , double );

