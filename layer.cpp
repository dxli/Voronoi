#include"layer.h"

neighbour::neighbour(int n0,int f0,double a0,rxy r0,double r20) {
        ni=n0;
        pos=f0;
        flag=0;
        arg=a0;
        r=r0;
        r2=r20;
    }
neighbour::neighbour(void){
}

bool operator < (const neighbour & p0,const neighbour p1) {
    return( p0.arg < p1.arg);
}

rxy layer::displ (rxy r0,rxy r1) {
        return(rxy(
                   fmod(r0.x -r1.x +half3lx,lx)-halflx,
                   fmod(r0.y -r1.y +half3ly,ly)-halfly
               )
              );
    }
void
layer::init (double temperature0)
{
    //temperature = temperature0;
    //op_steps=0;
    //oldc.resize (n);
    vnb.resize (n);
    //pv1.resize (n);
    //vz.resize (n);
    //v_nb.resize (n);
    //pv_nb.resize (n);
    halflx = 0.5 * lx;
    halfly = 0.5 * ly;
    half3lx = 1.5 * lx;
    half3ly = 1.5 * ly;
    sigma = sqrt (2. / sqrt (3.) * (lx * ly) / n);
    //q0t=4.*M_PI/(sqrt(3.)*sigma);
    //fn_op=fn_data+"op.txt";
    fn_data_chkpt=fn_data+".chkpt";
    fn_data_chkpt_bak=fn_data+".chkpt.bak";
    //bx = (int) (lx / (rv0 * sigma));
    //by = (int) (ly / (rv0 * sigma));
    //nrx = bx / lx;
    //nry = by / ly;
    vbx = (int) (lx / (rvb0 * sigma));
    vby = (int) (ly / (rvb0 * sigma));
    vnrx = vbx / lx;
    vnry = vby / ly;
    //cells.resize (bx);
    vcells.resize (vbx);
    /*
      for (int ii = 0; ii < bx; ii++)
    {
        cells.at(ii).resize (by);
    }
    */
    for (int ii = 0; ii < vbx; ii++)
    {
        vcells.at(ii).resize (vby);
    }
    vnb.resize(n);
    cout << "Simulating box: number of particles= " << n << " (" << lx << ',' <<
         ly << "), sigma= " << sigma << endl;
    //rvs = rv0 > rc0 ? rv0 * sigma : rc0 * sigma;
    //cout << "Potential cutoff rv= " << rvs << endl;
    //rvs *= rvs;
    //drvs = (rv0 - rc0) * sigma;
    //drvs *= drvs;
    vrvs= rvb0*sigma;
    vrvs *= vrvs;
    //mk_potential (table_length);
}

double
layer::dr2 (vector < atom >::iterator p1c)
{
    double dx = halflx - fabs (halflx - fabs (new_atom.xi - p1c->xi));
    double dy = halfly - fabs (halfly - fabs (new_atom.yi - p1c->yi));
    return (dx*dx + dy * dy);
}

void
layer::read_data ()
// read in the data file
{
    int i = 1;
    string fn_chk;
    fn_chk = fn_data + string (".chkpt");
    if (!read_data_chk (fn_chk))
    {
        cout << "Using checkpoint file " << fn_chk << endl;
        i = 0;
    }
    else
    {
        fn_chk = fn_data + string (".chkpt.bak");
        if (!read_data_chk (fn_chk))
        {
            cout << "Using checkpoint file " << fn_chk << endl;
            i = 0;
        }
        else
        {
            ifstream in1 (fn_data.c_str ());
            if (in1.is_open ())
            {
                string linebuf;
                getline (in1, linebuf);
                istringstream iss (linebuf);
                vector < double >va;
                std::copy (istream_iterator < double >(iss),
                           istream_iterator < double >(), back_inserter (va));
                if (va.size () >= 3)
                {
                    n = (int) (va.front() + 0.5);
                    lx = va.at(1);
                    ly = va.at(2);
                    p.resize (n);
                    in1.read ((char *) (&(p.front())), n * sizeof (atom));
                    if (!in1.bad ())
                    {
                        cout << "Using data file " << fn_data << endl;
                        i = 0;
                    }
                }
                in1.close ();
            }

        }
    }
    if (i)
    {
        //Failed on reading

        cout<<"Creating a new data set\n";
        int n0 = 400;
        int nx = n0 / 2;
        int ny = (int) (sqrt (nx * sqrt (3.)) + 0.5);
        nx = (int) (0.5 + (double) nx / ny);
        n = nx * ny * 2;
        lx = sqrt (3.) * nx;
        ly = (double) ny;
        p.resize (0);
        double dx = sqrt (3.), dy = 1., xi0 = 0.45, yi0 = 0.45*sqrt(3.);
        for (int i0 = 0; i0 < nx; i0++)
        {
            for (int j0 = 0; j0 < ny; j0++)
            {
                p.push_back (atom(i0 * dx + xi0*random()/(RAND_MAX+1.0), j0*dy + yi0*random()/(RAND_MAX+1.0)));
                p.push_back (atom(i0 * dx + xi0*random()/(RAND_MAX+1.0)+0.5*sqrt(3.), j0*dy + yi0*random()/(RAND_MAX+1.0)+0.5));
            }
        }


    }
    cout<<"n= "<<n<<" lx= "<<lx<<" ly= "<<ly<<endl;
}

int
layer::read_data_chk (string fn_chk)
// read in a checkpoint file
{
    ifstream in1 (fn_chk.c_str (), ifstream::binary);
    if (!in1.is_open ())
        return (1);
    in1.read ((char *) (&n), sizeof (int));
    if (in1.bad ())
        goto err1;
    in1.read ((char *) (&lx), sizeof (double));
    in1.read ((char *) (&ly), sizeof (double));
    in1.read ((char *) (&sigma01), sizeof (double) );
    if (in1.bad ())
        goto err1;
    p.resize (n);
    in1.read ((char *) (&(p.front())), n * sizeof (atom) );
    if (!in1.bad ())
        return (0);
err1:				//failed
    in1.close ();
    return (1);
}


void layer::add_vnb(int i0,int j0)
//add neighbour for Voronoi construction
{
    double dx= fmod(p.at(j0).xi - p.at(i0).xi + half3lx, lx)-halflx;
    double dy= fmod(p.at(j0).yi - p.at(i0).yi + half3ly, ly)-halfly;
    double a0=atan2(dx,dy);
    double r2=dx*dx+dy*dy;
    vnb.at(i0).push_back(neighbour(j0,vnb.at(j0).size(),a0,rxy(dx,dy),r2));
    vnb.at(j0).push_back(neighbour(i0,vnb.at(i0).size()-1,fmod(2.*M_PI + a0, 2.*M_PI)- M_PI ,rxy(-dx,-dy),r2));
}

int nb_cmp(const void *p0, const void *p1)
{
    neighbour *pt0=(neighbour *) p0,*pt1=(neighbour *) p1;
    if(pt0->arg < pt1->arg) return(-1);
    if(pt0->arg > pt1->arg) return(1);
    else return(0);
}
void layer::voronoi()
// voronoi construction
{
    //find neighbours
    for (int i0 = 0; i0 < vbx; i0++)
    {
        for (int j0 = 0; j0 < vby; j0++)
        {
            vcells.at(i0).at(j0).resize (0);
        }
    }
    vector < atom >::iterator pt7 = p.begin ();
    for (int i0 = 0; i0 < n; i0++)
    {
        vcells.at((int) (pt7->xi * vnrx) % vbx).at((int) (pt7->yi * vnry) % vby).push_back (i0);
        vnb.at(i0).resize (0);
        pt7++;
    }
    for (int i0 = 0; i0 < vbx; i0++)
    {
        for (int j0 = 0; j0 < vby; j0++)
        {
            vector < int >::iterator ptv1 = vcells.at(i0).at(j0).begin ();
            vector < int >::iterator ptv10 = vcells.at(i0).at(j0).end ();
            while (ptv1 != ptv10)
            {
                int l = *ptv1++;
                new_atom=p.at(l);
                vector < int >::iterator ptv2 = ptv1;
                while (ptv2 != ptv10)
                {
                    double dx=dr2( p.begin() + *ptv2);
                    if ( dx < vrvs)
                    {
                        add_vnb (l, *ptv2 );
                    }
                    ptv2++;
                }
                int i1 = 1, j1 = 0;
                while (i1 >= -1)
                {
                    int i2 = (i0 + i1 + vbx) % vbx;
                    int j2 = (j0 + j1 + vby) % vby;
                    vector < int >::iterator ptv2 = vcells.at(i2).at(j2).begin ();
                    vector < int >::iterator ptv20 = vcells.at(i2).at(j2).end ();
                    while (ptv2 != ptv20)
                    {
                        double dx=dr2( p.begin() + *ptv2);
                        if ( dx < vrvs)
                        {
                            add_vnb (l, *ptv2);
                        }
                        ptv2++;
                    }
                    if (j1) {
                        i1--;
                    } else {
                        j1++;
                    }
                }
            }
        }
    }

    /* start the voronoi construction */
    //pvnb.resize(0);
    nlist.resize(0);
    for(int i=0; i<n; i++) {
        new_atom=p.at(i);
        if(vnb.at(i).size()<=0) continue;
        vector<neighbour> ncp;
        for(unsigned int j=0; j<vnb.at(i).size(); j++) {
            if (vnb.at(i).at(j).flag>=0) ncp.push_back(vnb.at(i).at(j));
        }
        //qsort((void *) ( & (ncp.at(0))),ncp.size(),sizeof(neighbour), nb_cmp);
        sort(ncp.begin(),ncp.end());
        /* Voronoi Construction */
        vector<neighbour>::iterator pnb1=ncp.begin();
        vector<neighbour>::iterator pnb2=ncp.end();
        //forming the circul
        int length=0;
        while (pnb1 != pnb2 ) {
            pnb1->r=rxy(
                        fmod(p.at(pnb1->ni).xi - new_atom.xi +half3lx,lx)-halflx,
                        fmod(p.at(pnb1->ni).yi - new_atom.yi +half3ly,ly)-halfly);
            pnb1->next = ++pnb1;
            length++;
        };
        pnb1 --;
        pnb2=pnb1->next = ncp.begin();
        vector<neighbour>::iterator pnb3=pnb2->next;

        int j=0;
        rxy r1=pnb1->r,r2=pnb2->r,r3=pnb3->r;
        do {
            /* loop for each atom*/
            //    cout<<j<<' '<<length<<endl;

            if (pnb2->flag<1) {
                double index=pnb2->r2*det(pnb1->r,pnb3->r)/(pnb1->r2*det(pnb2->r,pnb3->r)+pnb3->r2*det(pnb1->r,pnb2->r));
                //double index=(pnb1->r2*det(pnb2->r,pnb3->r)+pnb3->r2*det(pnb1->r,pnb2->r))/det(pnb1->r,pnb3->r);
                if(index>1.0) { //remove pnb2
               // if(index > 0. && index<pnb2->r2) { //remove pnb2
                    vnb.at(pnb2->ni).at(pnb2->pos).flag=-1;
                    pnb1->next=pnb3;
                    //  delete &( *pnb2);
                    pnb2=pnb3;
                    pnb3= pnb3->next;
                    j=0;
                    length--;
                    continue;
                }
            }
            j++;
            pnb1=pnb2;
            pnb2=pnb3;
            pnb3= pnb3->next;
        } while (j<=length);
        //if(length<4){
        //        cout<<i<<' '<<length<<endl;
        //}
        vnb.at(i).resize(length+1);
        for(j=0; j<length; j++) {
            //if( pnb2->ni > i) vnb.at(pnb2->ni).at(pnb2->pos).flag=1;
            vnb.at(i).at(j)=*pnb2;
            //vnb.at(i).push_back(*pnb2);
            pnb2=pnb2->next;
        }
        vnb.at(i).at(j)=vnb.at(i).at(0);
        //pvnb.push_back(pnb2);
        //cout<<i<<' '<<length<<endl;
        nlist.push_back(length);
    }
}

void layer::w_voronoi()
//plotting
{
    voronoi();
    GraceRegisterErrorFunction(my_error_function);
    if (GraceOpenVA((char *)"xmgrace", 40960, "-nosafe", "-noask","-nosigcatch","-geometry","1400x980+0+0", NULL)==-1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }
    GracePrintf("view xmin 0.15");
    GracePrintf("view xmax 0.95");
    GracePrintf("view ymin 0.15");
    GracePrintf("view ymax 0.95");
    GracePrintf("world xmin %g",-sigma);
    GracePrintf("world xmax %g",lx+sigma);
    GracePrintf("world ymin %g",-sigma);
    GracePrintf("world ymax %g",ly+sigma);
    GracePrintf("autoscale onread none");
    GracePrintf("yaxis tick major %g",axis_normal(0,ly));
    GracePrintf("xaxis tick major %g",axis_normal(0,lx));
    GracePrintf("xaxis tick minor 0");
    GracePrintf("yaxis tick minor 0");

    vector<atom>::iterator pt1=p.begin();
    int i;


    int js=0;
    for(int i=0; i<n; i++) {
        js += nlist.at(i);
        vector<neighbour>::iterator pt7=vnb.at(i).begin();
        int pcolor=0;
        switch(nlist.at(i)) {
        case 4:
            pcolor=8;
            break;
        case 5:
            pcolor=2;
            break;
            /*
                            case 6:
            pcolor=0;
            break;
            */
        case 7:
            pcolor=3;
            break;
        case 8:
            pcolor=4;
            break;
        default:
            pcolor=0;
        }


        rxy r0(pt1 + i);
        vector<double> xf,yf;
        rxy r1(pt1 + pt7++->ni);
        r1 = displ(r1 ,r0);
        rxy r2(pt1 + pt7->ni);
        r2 = displ(r2 ,r0);
        double r12=0.5*r1.abs2(),r22=0.5*r2.abs2();
        for(int j=0; j<nlist.at(i); j++) {
            rxy r3= rxy(
                        det( rxy(r12,r1.y),rxy(r22,r2.y))/det(r1,r2),
                        det( rxy(r1.x,r12),rxy(r2.x,r22))/det(r1,r2)
                    );
            xf.push_back(r0.x + r3.x);
            yf.push_back(r0.y + r3.y);
            r1=r2;
            r12=r22;
            r2=displ(rxy(pt1 + (++pt7)->ni),r0);
            r22=0.5*r2.abs2();
        }
        //cout<<i<<' '<<nlist.at(i)<<endl;
        /*
        for(int j=0;j<xf.size();j++){
                cout<<xf.at(j)<<' '<<yf.at(j)<<endl;
        }
        */
        GracePrintf("s%d on",i);
        GracePrintf("s%d type xy",i);
        GracePrintf("s%d symbol 0",i);
        if(pcolor) {
            GracePrintf("s%d linestyle 0",i);
            //GracePrintf("s%d color %d",i,pcolor);
            GracePrintf("s%d fill 1",i);
            GracePrintf("s%d fill with color",i);
            GracePrintf("s%d fill pattern 1",i);
            GracePrintf("s%d fill color %d",i,pcolor);
            for(unsigned int j=0; j<xf.size(); j++) {
                GracePrintf("s%d point %g,%g",i,xf[j],yf[j]);
            }
        } else {
            GracePrintf("s%d linestyle 1",i);
            GracePrintf("s%d line color 1",i);
            for(unsigned int j=0; j<xf.size(); j++) {
                //if(i==14063)  cout<<'('<<xf[j]<<','<<yf[j]<<')';
                GracePrintf("s%d point %g,%g",i,xf[j],yf[j]);
            }
        }
        GracePrintf("s%d point %g,%g",i,xf[0],yf[0]);
    }
    cout<<"Average number of neighbours = "<<js/n<<' '<<js%n<<endl;
//draw points
    i=n;
    GracePrintf("s%d on",i);
    GracePrintf("s%d type xy",i);
    GracePrintf("s%d symbol 1",i);
    GracePrintf("s%d symbol size %g",i,10./sqrt(n));
    GracePrintf("s%d linestyle 0",i);
    GracePrintf("s%d symbol color 1",i);
    GracePrintf("s%d symbol pattern 1",i);
    GracePrintf("s%d symbol fill color 1",i);
    GracePrintf("s%d symbol fill pattern 1",i);
    for(int j=0; j<n; j++)
    {
        GracePrintf("s%d point %g,%g",i,p.at(j).xi,p.at(j).yi);
    }
//draw box
    i++;
    GracePrintf("s%d on",i);
    GracePrintf("s%d type xy",i);
    GracePrintf("s%d symbol 0",i);
    GracePrintf("s%d linestyle 1",i);
    GracePrintf("s%d line color 14",i);
    GracePrintf("s%d point %g,%g",i,0.,0.);
    GracePrintf("s%d point %g,%g",i,0.,ly);
    GracePrintf("s%d point %g,%g",i,lx,ly);
    GracePrintf("s%d point %g,%g",i,lx,0.);
    GracePrintf("s%d point %g,%g",i,0.,0.);



    GracePrintf("xaxis tick out");
    GracePrintf("yaxis tick out");
    GracePrintf("xaxis ticklabel char size 1.5");
    GracePrintf("yaxis ticklabel char size 1.5");
    GracePrintf("xaxis label char size 1.5");
    GracePrintf("yaxis label char size 1.5");
    GracePrintf("xaxis label \"\\f{Symbol}s\"");
    GracePrintf("yaxis label \"\\f{Symbol}s\"");
    GracePrintf("print to \"%s.ps\"",fn_data.c_str());
    printf("print to \"%s.ps\"\n",fn_data.c_str());
    GracePrintf("redraw");
    GracePrintf("print");
    sleep(2);
    string cmd=string("ps2epsi ")+fn_data+".ps "+fn_data+".eps";
    cout<<cmd<<endl;
    system(cmd.c_str());

    exit(0);
}

