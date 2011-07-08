#pragma once
#include"layer.h"
using namespace std;

class RayP{
	public:
	int site;
	int nb_L,nb_R;
	RayP *pnb_L,*pnb_R, *p3nb_L,*p3nb_R;
	RayP(int, int, int);
	double x3();
	short int compare_y(int);
};
class BeachLine{
	public:
	RayP *Broot;
	int size;
	BeachLine(void);
	BeachLine(int);
	void add(int);
	void remove(RayP);
	
};
class ThreeSiteTree {
	public:
		RayP *root;
};

