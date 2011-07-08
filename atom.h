#pragma once
#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<iterator>
#include<algorithm>
using namespace std;

class atom {
public:
    double xi,yi;
    atom(double,double);
	atom(void);
};

class rxy
// the 2D complex class
{
public:
    double x,y;
    rxy(double,double);
    rxy(double );
    rxy(vector<atom> ::iterator );
	rxy(void);
    double abs2();
};

double det(rxy ,rxy );
