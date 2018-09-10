//
// Created by ZHI on 5/23/2018.
//

#ifndef REVISEDSIO2_ENERGY_FORCE_H
#define REVISEDSIO2_ENERGY_FORCE_H

#include "CONSTANT.h"
#include "REFERENCE.h"
#include <cmath>
#include <chrono>
#include <iostream>

using namespace std;

class ENERGY_FORCE {
private:
    long double angular_energy_Si, angular_energy_O;
    long double radial_energy;
    long double repulsive_energy_Si_Si, repulsive_energy_O_O;
    long double center_Si[3], center_O[3];
    long double bound_O[4][3], bound_shortest_O[4][3];
    long double bound_Si[2][3], bound_shortest_Si[2][3];
    int order_O[4], order_Si[2];
    long double length_Si_O1, length_Si_O2, length_Si_O3, length_Si_O4;
    long double length_O_Si1,length_O_Si2;
    long double cosSi_O1_O2, cosSi_O1_O3, cosSi_O1_O4, cosSi_O2_O3, cosSi_O2_O4, cosSi_O3_O4;
    long double cosO_Si1_Si2;
    long double Rf_O1_on_Si[3], Rf_O2_on_Si[3], Rf_O3_on_Si[3], Rf_O4_on_Si[3];
    long double Rf_on_Si[3];
    long double fSiO1O2_on_O1[3], fSiO1O2_on_O2[3], fSiO1O2_on_Si[3];
    long double fSiO1O3_on_O1[3], fSiO1O3_on_O3[3], fSiO1O3_on_Si[3];
    long double fSiO1O4_on_O1[3], fSiO1O4_on_O4[3], fSiO1O4_on_Si[3];
    long double fSiO2O3_on_O2[3], fSiO2O3_on_O3[3], fSiO2O3_on_Si[3];
    long double fSiO2O4_on_O2[3], fSiO2O4_on_O4[3], fSiO2O4_on_Si[3];
    long double fSiO3O4_on_O3[3], fSiO3O4_on_O4[3], fSiO3O4_on_Si[3];
    long double fOSi1Si2_on_Si1[3],fOSi1Si2_on_Si2[3],fOSi1Si2_on_O[3];
    long double repSiSi[3], repSiSishortest[3], repSiSilength, frepSiSi[3];
    long double repOO[3], repOOshortest[3], repOOlength, frepOO[3];
public:
    ENERGY_FORCE(){};
    void GroupfromGraph(int groupSi[][4], int groupO[][2], long double Graph[][128]);
    void Force_Energy(long double Si[][3], long double O[][3],
                      long double forceSi[][3], long double forceO[][3],
                      long double& total_energy,
                      int groupSi[][4], int groupO[][2]);
};


#endif //REVISEDSIO2_ENERGY_FORCE_H
