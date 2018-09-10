//
// Created by ZHI on 5/23/2018.
//

#ifndef REVISEDSIO2_REFERENCE_H
#define REVISEDSIO2_REFERENCE_H

#include <cmath>
#include "GROUP.h"
#include "CONSTANT.h"

using namespace std;

class REFERENCE {
private:
    long double costheta;
    long double xO, yO, zO;
    long double xSi, ySi, zSi;
    long double xcomp, ycomp, zcomp;
    long double abs_xOsubs1,abs_xOadds1;
    long double abs_yOsubs1,abs_yOadds1;
    long double abs_zOsubs1,abs_zOadds1;
    long double distance;
public:
    void shortestcoordinate(const long double center[3], const long double comp[3], long double shortest[3]);
    long double real_length(long double center[3], long double comp[3]);
    long double angle(const long double center[3], const long double O1[3], const long double O2[3], long double cO1, long double cO2);
};


#endif //REVISEDSIO2_REFERENCE_H
