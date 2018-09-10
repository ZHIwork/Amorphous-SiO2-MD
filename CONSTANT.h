//
// Created by ZHI on 5/23/2018.
//

#ifndef REVISEDSIO2_CONSTANT_H
#define REVISEDSIO2_CONSTANT_H

#include <iostream>
#include <cmath>

using namespace std;

struct CONSTANT
{
    string filename = "SiO2.vasp";
    string outfilename = "final.vasp";
    long double time_interval = 1.0e-15;
    int minimize_times = 5000;
    int cut_times = 8000;
    long double temperature = 6000.0;
    long double Boltzmann_constant = 8.617330350e-5;
    long double resistance = 0.0;
    const int Si_quantity = 64;
    const int O_quantity = 128;
    int total_quantity = Si_quantity + O_quantity;
    long double mass_Si = 46.82e-27;
    long double mass_O = 26.67e-27;
    const long double transfer_vector[3][3] = {
            {14.332000000000 ,  0.000000000000 ,  0.000000000000},
            {0.0000000000000 , 14.332000000000 ,  0.000000000000},
            {0.0000000000000 ,  0.000000000000 , 14.332000000000}
    };
    long double transfer = transfer_vector[0][0];
    long double b0 = 1.55148;
    long double kb = 26.96;
    long double kthetaSi = 1.685, kthetaO = 0.58;
    long double thetaSi = 109.47, thetaO = 180;
    long double cos0_Si = cos(thetaSi * 3.141592653589793/180);
    long double cos0_O = cos(thetaO * 3.141592653589793/180);
    long double repulsive_coe_for_SiSi = 8.0;
    long double repulsive_coe_for_OO = 2.0;

};

#endif //REVISEDSIO2_CONSTANT_H
