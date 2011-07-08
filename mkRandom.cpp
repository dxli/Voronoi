#include"mkRandom.h"

void mk_random (int nl, string r_buffer)
{
    struct timeval time1;
    gettimeofday (&time1, NULL);
    ifstream in1 ("/dev/urandom", ifstream::binary);
    r_buffer.resize(257);
    if (!in1.is_open ())
    {
        initstate (time1.tv_usec, (char *) r_buffer.c_str (), 256);
    }
    else
    {
        in1.read ((char *) (r_buffer.c_str ()), sizeof (char) * 256);
        initstate (time1.tv_usec, (char *) r_buffer.c_str (), 256);
        in1.close ();
    }
}

