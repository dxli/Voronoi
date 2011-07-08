#include "fortune.h"

RayP::RayP(int i){
	site=i;
	nb_L=-1;
	nb_R=-1;
	pnb_L=NULL;
	pnb_R=NULL;
	p3nb_L=NULL;
	p3nb_R=NULL;
}

short int RayP::compare(int){
	double yl;
	if ( pnb_L != NULL ) {
		
}
}

BeachLine::BeachLine(void){
	Broot=NULL;
	size=0;
}

BeachLine::BeachLine(int i){
	Broot= &(RayP(i));
	size=1;
}
