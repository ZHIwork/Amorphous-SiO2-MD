//
// Created by ZHI on 5/24/2018.
//

#include "DYNAMICS.h"

REFERENCE Dy_REF;
CONSTANT Dy_const;


void DYNAMICS::MD(long double forceSi[][3], long double forceO[][3],
        long double vSi[][3], long double vO[][3],
        long double Si[][3], long double O[][3], long double time)
{
    //update coordinates
    //for Si coordinate
    for(int i=0; i<Dy_const.Si_quantity; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            Si[i][j] += vSi[i][j]*time + 0.5*forceSi[i][j]/Dy_const.mass_Si*time*time;
        }
    }
    //for O coordinate
    for(int i=0; i<Dy_const.O_quantity; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            O[i][j] += vO[i][j]*time + 0.5*forceO[i][j]/Dy_const.mass_O*time*time;
        }
    }

    //update new velocity
    //for Si velocity and initialize force to 0
    for(int i=0; i<Dy_const.Si_quantity; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            vSi[i][j] += forceSi[i][j]/Dy_const.mass_Si*time*Dy_const.resistance;
            forceSi[i][j] = 0;
        }
    }
    //for O velocity and initialize force to 0
    for(int i=0; i<Dy_const.O_quantity; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            vO[i][j] += forceO[i][j]/Dy_const.mass_O*time*Dy_const.resistance;
            forceO[i][j] = 0;
        }
    }

};



//generate a integer in range (begin, end)
int DYNAMICS::generator_int(int begin, int end)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(begin, end);
    return dis(gen);
};


long double DYNAMICS::generator_lb(long double begin, long double end)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(begin, end);
    return dis(gen);
};


//cut bonds
void DYNAMICS::cutbond(long double Graph[][128])
{
    //choose a O atom as center
    int center_order = generator_int(0,Dy_const.O_quantity-1);
    //find two Si atoms bounded to center O
    int a=0;
    for(int i=0; i<Dy_const.Si_quantity; ++i)
    {
        if(Graph[i][center_order]!=0 && a<2)
        {
            orderSi[a] = i;
            ++a;
        }
    }
    //get other three O atoms bounded to two Si atoms
    int b1=0;
    int b2=0;
    for(int i=0; i<Dy_const.O_quantity; ++i)
    {
        if(i != center_order)
        {
            if(Graph[orderSi[0]][i]!=0 && b1<3)
            {
                ordernotcenterO[0][b1] = i;
                ++b1;
            }
            if(Graph[orderSi[1]][i]!=0 && b2<3)
            {
                ordernotcenterO[1][b2] = i;
                ++b2;
            }

        }
    }
    //choose two O atoms randomly
    //for first Si
    int random1 = generator_int(0,2);
    int orderO1 = ordernotcenterO[0][random1];
    //for second Si
    int random2 = generator_int(0,2);
    int orderO2 = ordernotcenterO[1][random2];
    //exchange bonds
    Graph[orderSi[0]][orderO1] = 0;
    Graph[orderSi[1]][orderO2] = 0;

    Graph[orderSi[0]][orderO2] = 1;
    Graph[orderSi[1]][orderO1] = 1;


    //cout << "center O: "<<center_order+1<<endl;
    //cout << "Si1 and O1: "<< orderSi[0]+1<< "    "<<orderO1+1<<endl;
    //cout << "Si2 and O2: "<< orderSi[1]+1<< "    "<<orderO2+1<<endl;

};


bool DYNAMICS::MC_probability(long double old_energy, long double new_energy)
{
    long double T = Dy_const.temperature;
    long double kB = Dy_const.Boltzmann_constant;
    //difference of energy
    long double deltE = new_energy - old_energy;
    //probability
    long double comp = exp(-deltE/(T*kB));
    if(comp<1.0)
    {
        //get another random number
        long double randomNum = generator_lb(0.0,1.0);
        if(randomNum<comp)
            return 1;
        else
            return 0;
    }
    else
        return 1;
}