#include"myGrace.h"
void my_error_function(const char *msg)
{
cerr<<"library message: \""<<msg<<"\"\n";
}

double axis_normal(double xmin, double xmax)
{
double a=fabs((double) (xmax -xmin));
a=log10(a);
double b=floor(a);
double c=(int) ( pow(10.,a-b));
if(c<2.) {
        b -=1.0;
        a=5.;
}else{
        if(c>5.) {
                a=2.;
        }else {
                a=1.;
        }
}
//cout<<xmin<<' '<<xmax<<' '<<a*pow(10.,b)<<endl;
return(pow(10.,b)*a);
}
