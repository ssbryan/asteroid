
#include "ObjectMgr.h"

// Author Stephen Bryan
// June 22, 2017

int main(int argc, char* argv[])
{
    // read in initialization file, create Objects
    std::string fname = "init.data";

    // get init file from args (default "init.dat")
    if (argc >= 2)
    {
        fname = argv[1];

        if (!fname.compare("-?") || !fname.compare("/?") || !fname.compare("-help") || !fname.compare("/help"))
        {
            printf("Usage:  asteroid [filename]\n");
            printf("\tFile format:\n");
            printf("\t\tOptions: [momtol = 0.0001, ] [dvnormtol = 0.00001, ] [runs = 1e6, ] [tstep = 1.0, ] [minvforcheck = 0.001, ] [startdatasave = 0.0, ] [nthdatasave = 1]\n");
            printf("\t\t<Object data>*\n");
            printf("\t\t<Object data>: xloc yloc zloc xvel yvel zvel mass [all doubles, as many objects, each on its own line, as desired]\n");
            printf("\tExample:\n");
            printf("\t\tOptions: momtol = 2e-2, dvnormtol = 2e-2, runs = 1.0e9, tstep = 1.0, minvforcheck = 0.01, startdatasave = 2e6, nthdatasave = 1024\n");
            printf("\t\t0.0 0.0 0.0 0.0 0.0 0.0 1.0e11\n");
            printf("\t\t-1000.0 0.0 0.0 0.0 - 0.0817 0.0 1.0e5\n");
            printf("\t\t1000.0 0.0 0.0 0.0 0.0817 0.0 1.0e5\n");
            printf("\t\t0.0 - 1000.0 0.0 0.0817 0.0 0.0 1.0e5\n");
            printf("\t\t0.0 1000.0 0.0 - 0.085 0.0 0.0 1.0e5\n");
            return 1;
        }
    }

    ObjectMgr omgr(fname);
    omgr.Run();

    return 0;
}
