//
// Created by ZHI on 5/24/2018.
//

#ifndef REVISEDSIO2_DYNAMICS_H
#define REVISEDSIO2_DYNAMICS_H

#include <random>
#include "REFERENCE.h"
#include "CONSTANT.h"

class DYNAMICS {
private:
    int orderSi[2], ordernotcenterO[2][3];
public:
    void MD(long double forceSi[][3], long double forceO[][3],
            long double vSi[][3], long double vO[][3],
            long double Si[][3], long double O[][3], long double time);
    int generator_int(int begin, int end);
    long double generator_lb(long double begin, long double end);
    void cutbond(long double Graph[][128]);
    bool MC_probability(long double old_energy, long double new_energy);
};


#endif //REVISEDSIO2_DYNAMICS_H
