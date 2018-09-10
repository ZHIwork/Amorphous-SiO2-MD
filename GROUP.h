//
// Created by ZHI on 5/23/2018.
//

#ifndef REVISEDSIO2_GROUP_H
#define REVISEDSIO2_GROUP_H

#include "CONSTANT.h"
#include <cmath>
#include <iostream>

using namespace std;

class GROUP {
private:
    long double xO, yO, zO;
    long double xSi, ySi, zSi;
    long double xO_shortest, yO_shortest, zO_shortest;
    long double xcomp, ycomp, zcomp;
    long double abs_xOsubs1,abs_xOadds1;
    long double abs_yOsubs1,abs_yOadds1;
    long double abs_zOsubs1,abs_zOadds1;
    long double distance, Distance;
public:
    GROUP(){};
    void Group(long double Si[][3], long double O[][3], long double Graph[][128]);
    long double Shortestdistance(const long double center[3], const long double comp[3]);
};


#endif //REVISEDSIO2_GROUP_H
