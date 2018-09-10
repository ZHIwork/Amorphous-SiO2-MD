//
// Created by ZHI on 5/23/2018.
//

#include "REFERENCE.h"

CONSTANT REF_const;
GROUP REF_group;

void REFERENCE::shortestcoordinate(const long double center[3], const long double comp[3], long double shortest[3])
{
    //get Si coordinates
    xSi = center[0];
    ySi = center[1];
    zSi = center[2];

    //get O coordinates
    xO = comp[0];
    yO = comp[1];
    zO = comp[2];

    //get O coordinates which is nearest to Si
    //for x
    xcomp = abs(xSi-xO);//xO not change
    shortest[0] = xO;
    abs_xOsubs1 = abs(xSi-(xO-1));
    abs_xOadds1 = abs(xSi-(xO+1));
    if (abs_xOsubs1 < xcomp)//xO-1
    {
        shortest[0] = xO-1;
        xcomp = abs_xOsubs1;
        if(abs_xOadds1 < xcomp)//xO+1
            shortest[0] = xO+1;

    }
    else
    {
        if(abs_xOadds1 < xcomp)//xO+1
            shortest[0] = xO+1;
    }

    //for y
    ycomp = abs(ySi-yO);//yO not change
    shortest[1] = yO;
    abs_yOsubs1 = abs(ySi-(yO-1));
    abs_yOadds1 = abs(ySi-(yO+1));
    if (abs_yOsubs1 < ycomp)//yO-1
    {
        shortest[1] = yO-1;
        ycomp = abs_yOsubs1;
        if(abs_yOadds1 < ycomp)//yO+1
            shortest[1] = yO+1;

    }
    else
    {
        if(abs_yOadds1 < ycomp)//yO+1
            shortest[1] = yO+1;
    }

    //for z
    zcomp = abs(zSi-zO);//z not change
    shortest[2] = zO;
    abs_zOsubs1 = abs(zSi-(zO-1));
    abs_zOadds1 = abs(zSi-(zO+1));
    if (abs_zOsubs1 < zcomp)//zO-1
    {
        shortest[2] = zO-1;
        zcomp = abs_zOsubs1;
        if(abs_zOadds1 < zcomp)//zO+1
            shortest[2] = zO+1;

    }
    else
    {
        if(abs_zOadds1 < zcomp)//zO+1
            shortest[2] = zO+1;
    }
};



long double REFERENCE::real_length(long double center[3], long double comp[3])
{
    distance = REF_const.transfer*sqrt((center[0]-comp[0])*(center[0]-comp[0]) +(center[1]-comp[1])*(center[1]-comp[1])+(center[2]-comp[2])*(center[2]-comp[2]));
    return distance;
};



long double REFERENCE::angle(const long double center[3], const long double O1[3], const long double O2[3], long double cO1, long double cO2)
{
    //cosine
    costheta = REF_const.transfer * REF_const.transfer * ( (center[0]-O1[0])*(center[0]-O2[0]) + (center[1]-O1[1])*(center[1]-O2[1]) + (center[2]-O1[2])*(center[2]-O2[2]) ) / (cO1*cO2);
    return costheta;
};