#pragma once
#include<sstream>
#include"atom.h"
#include"parameters.h"
#include"myGrace.h"

class neighbour
{
public:
    int ni;
    int pos;
    int flag;
    double arg;
    rxy r;
    double r2;
    vector<neighbour>::iterator next;
    neighbour(void); 
    neighbour(int ,int ,double ,rxy ,double );
};
class layer
    // simulating box
{
public:
    int n;
    double lx, ly;
    vector < atom >p, oldc;
    double nrx, nry;
    int bx, by;
    vector < vector < vector < int > > >cells;
    double rvs, drvs;
    //vector < double >vz, esave;
    vector < vector < int > >pv1;
    vector < vector < double > >v_nb;
    //vector < vector < int > >pv_nb;	// position of potential in neighbour's list
    //double temperature, tkpb,q0t,u_average;
    //int op_steps;
    //vector<vector < neighbour >::iterator> pvnb;
    vector<int> nlist;
    string fn_data,fn_op, fn_data_chkpt, fn_data_chkpt_bak;
    //layer () { };
    void init (double);
    void voronoi(), write_data (), read_data ();
    void add_nb (int, int, double),add_vnb(int,int);
    double old_v (int), new_v (int),u_total(),gtop(),goop();
    int read_data_chk (string);
    void w_voronoi ();

private:
    double halflx, halfly,half3lx,half3ly, sigma,sigma01, xi, yi, zi, vznew;
    atom new_atom;
    double vnrx,vnry;
    int vbx,vby;
    vector < vector < vector < int > > >vcells;
    vector<vector < neighbour > > vnb;
    double vrvs;
    //double r2 (int, int);
    double dr2 (vector < atom >::iterator);
    void accept (int);
    rxy displ (rxy ,rxy );
};

