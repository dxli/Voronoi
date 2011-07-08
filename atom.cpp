#include<iostream>
#include"atom.h"

atom::atom(double x0,double y0) {
    xi=x0;
    yi=y0;
}
atom::atom(void){
}
ostream& operator << (ostream& os, atom a0) {
    // Output the atom position to the stream in the format n/d
    os <<'('<<a0.xi<<','<<a0.yi<<')';
    return os;
}

rxy::rxy(void){
}
rxy::rxy(double x0,double y0) {
        x=x0;
        y=y0;
    }

rxy::rxy(vector<atom>::iterator a0) {
        x=a0->xi;
        y=a0->yi;
    }

rxy::rxy(double a0) {
        x=cos(a0);
        y=sin(a0);
    }

double  operator * (rxy a0, rxy a1) {
    return( a0.x * a1.x + a0.y*a1.y);
}

rxy operator / (rxy a0, int n0) {
    return( rxy(a0.x /n0 , a0.y / n0));
}
rxy operator /= (rxy &a0, int n0) {
    a0.x /= n0;
    a0.y /= n0;
    return( a0);
}
rxy operator * (double q0, rxy a0) {
    return( rxy(a0.x *q0 , a0.y * q0));
}

rxy operator += (rxy &a0, rxy a1) {
    a0.x += a1.x;
    a0.y += a1.y;
    return(a0);
}

double rxy::abs2()
{
    return( this->x* this->x + this->y * this->y);
}

double det(rxy a0, rxy a1) {
//determinant
    return(a0.x*a1.y-a1.x*a0.y);
}
