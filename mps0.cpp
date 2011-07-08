//
// voronoi for z3pt32d.cpp results
// show height fluctuations
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iterator>
#include<vector>
#include<algorithm>
#include<complex>
#include <grace_np.h>
#include <stdlib.h>

#include"myGrace.h"
#include"mkRandom.h"
#include"atom.h"
#include"layer.h"
#include<complex>
using namespace std;

int
main (int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Usage:: <datafile> \n";
        exit (0);
    }
    //initilize the random() buffer
    string r_buffer;
    r_buffer.resize (257);
    mk_random (r_buffer.size (), r_buffer);
    //prepare the layer data
    layer layer0;
    layer0.fn_data = string (argv[1]);
    layer0.read_data ();
    layer0.init (1.);
    layer0.w_voronoi ();
}

