//
// Created by ZHI on 5/23/2018.
//

#include "GROUP.h"

CONSTANT group_const;

void GROUP::Group(long double Si[][3],long double O[][3], long double Graph[][128])
{
    //graph atoms and get Si coordiantes
    for(int i=0; i<group_const.Si_quantity; ++i)
    {
        //get Si atoms' positions
        for(int j=0; j<group_const.O_quantity; ++j)
        {
            //get shortest distance and store into Graph
            Distance = Shortestdistance(Si[i], O[j]);
            if(Distance < 2.10){
                Graph[i][j] = Distance;
            }
            else
                Graph[i][j] = 0.0;
        }
    }
};


long double GROUP::Shortestdistance(const long double center[3], const long double comp[3])
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
    xO_shortest = xO;
    abs_xOsubs1 = abs(xSi-(xO-1));
    abs_xOadds1 = abs(xSi-(xO+1));
    if (abs_xOsubs1 < xcomp)//xO-1
    {
        xO_shortest = xO-1;
        xcomp = abs_xOsubs1;
        if(abs_xOadds1 < xcomp)//xO+1
            xO_shortest = xO+1;

    }
    else
    {
        if(abs_xOadds1 < xcomp)//xO+1
            xO_shortest = xO+1;
    }

    //for y
    ycomp = abs(ySi-yO);//yO not change
    yO_shortest = yO;
    abs_yOsubs1 = abs(ySi-(yO-1));
    abs_yOadds1 = abs(ySi-(yO+1));
    if (abs_yOsubs1 < ycomp)//yO-1
    {
        yO_shortest = yO-1;
        ycomp = abs_yOsubs1;
        if(abs_yOadds1 < ycomp)//yO+1
            yO_shortest = yO+1;

    }
    else
    {
        if(abs_yOadds1 < ycomp)//yO+1
            yO_shortest = yO+1;
    }

    //for z
    zcomp = abs(zSi-zO);//z not change
    zO_shortest = zO;
    abs_zOsubs1 = abs(zSi-(zO-1));
    abs_zOadds1 = abs(zSi-(zO+1));
    if (abs_zOsubs1 < zcomp)//zO-1
    {
        zO_shortest = zO-1;
        zcomp = abs_zOsubs1;
        if(abs_zOadds1 < zcomp)//zO+1
            zO_shortest = zO+1;

    }
    else
    {
        if(abs_zOadds1 < zcomp)//zO+1
            zO_shortest = zO+1;
    }

    //calculate shorteset distance
    distance = group_const.transfer*sqrt(pow((xSi-xO_shortest),2) +pow((ySi-yO_shortest),2)+pow((zSi-zO_shortest),2));
    return distance;
};



