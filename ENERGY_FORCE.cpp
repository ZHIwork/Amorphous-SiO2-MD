//
// Created by ZHI on 5/23/2018.
//

#include "ENERGY_FORCE.h"

using namespace std;
using  ns = chrono::nanoseconds;
using get_time = chrono::steady_clock ;


CONSTANT EF_const;
REFERENCE EF_Ref;

void ENERGY_FORCE::GroupfromGraph(int groupSi[][4], int groupO[][2], long double Graph[][128])
{
    //group Si and store O's order in array
    for(int i=0; i<EF_const.Si_quantity; ++i)
    {
        int n = 0;
        for (int j = 0; j < EF_const.O_quantity; ++j) {
            if (Graph[i][j]!=0 && n<4) {
                //store O's order
                groupSi[i][n] = j;
                ++n;
            }
        }
    }

    //group O and store Si's order in array
    for(int i=0;i<EF_const.O_quantity;++i)
    {
        int m = 0;
        for (int j = 0; j < EF_const.Si_quantity; ++j)
        {
            if (Graph[j][i]!=0 && m<2)
            {
                groupO[i][m] = j;
                ++m;
            }
        }
    }
};




void ENERGY_FORCE::Force_Energy(long double Si[][3], long double O[][3],
                                long double forceSi[][3], long double forceO[][3],
                                long double& total_energy,
                                int groupSi[][4], int groupO[][2])
{



    //initialization energy parameters
    angular_energy_Si = 0;
    angular_energy_O = 0;
    radial_energy = 0;
    repulsive_energy_Si_Si = 0;
    repulsive_energy_O_O = 0;

    //calculate length between each two bounded atoms (Si and O)
    int a=1;
    for(int i=0; i<EF_const.Si_quantity; ++i) {



        //store Si position
        center_Si[0] = Si[i][0];
        center_Si[1] = Si[i][1];
        center_Si[2] = Si[i][2];

        for(int j=0; j<4; ++j)
        {
            order_O[j] = groupSi[i][j];
            for(int k=0; k<3; ++k)
                bound_O[j][k] = O[order_O[j]][k];
        }

        for(int j=0; j<4; ++j)
            EF_Ref.shortestcoordinate(center_Si, bound_O[j], bound_shortest_O[j]);

        //--------------------------------------------------------------------------------------------------------------


        //length
        length_Si_O1 = EF_Ref.real_length(center_Si, bound_shortest_O[0]);
        length_Si_O2 = EF_Ref.real_length(center_Si, bound_shortest_O[1]);
        length_Si_O3 = EF_Ref.real_length(center_Si, bound_shortest_O[2]);
        length_Si_O4 = EF_Ref.real_length(center_Si, bound_shortest_O[3]);


        //angle
        cosSi_O1_O2 = EF_Ref.angle(center_Si, bound_shortest_O[0], bound_shortest_O[1], length_Si_O1, length_Si_O2);
        cosSi_O1_O3 = EF_Ref.angle(center_Si, bound_shortest_O[0], bound_shortest_O[2], length_Si_O1, length_Si_O3);
        cosSi_O1_O4 = EF_Ref.angle(center_Si, bound_shortest_O[0], bound_shortest_O[3], length_Si_O1, length_Si_O4);
        cosSi_O2_O3 = EF_Ref.angle(center_Si, bound_shortest_O[1], bound_shortest_O[2], length_Si_O2, length_Si_O3);
        cosSi_O2_O4 = EF_Ref.angle(center_Si, bound_shortest_O[1], bound_shortest_O[3], length_Si_O2, length_Si_O4);
        cosSi_O3_O4 = EF_Ref.angle(center_Si, bound_shortest_O[2], bound_shortest_O[3], length_Si_O3, length_Si_O4);




        //--------------------------------------------------------------------------------------------------------------

        /*calculate radial energy and radial force
         * and store forces to both Si array and O array */

        //coefficient
        long double diff_Si_O1 = length_Si_O1 - EF_const.b0;
        long double diff_Si_O2 = length_Si_O2 - EF_const.b0;
        long double diff_Si_O3 = length_Si_O3 - EF_const.b0;
        long double diff_Si_O4 = length_Si_O4 - EF_const.b0;

        //radial energy
        radial_energy += 0.5 * EF_const.kb *
                         (diff_Si_O1*diff_Si_O1 + diff_Si_O2*diff_Si_O2 +
                          diff_Si_O3*diff_Si_O3 + diff_Si_O4*diff_Si_O4);

        //radial force for Si
        long double fcoe = EF_const.transfer * EF_const.kb * diff_Si_O1 / length_Si_O1;
        Rf_O1_on_Si[0] = fcoe * (bound_shortest_O[0][0] - center_Si[0]);//force on x by O1
        Rf_O1_on_Si[1] = fcoe * (bound_shortest_O[0][1] - center_Si[1]);//force on y by O1
        Rf_O1_on_Si[2] = fcoe * (bound_shortest_O[0][2] - center_Si[2]);//force on z by O1

        fcoe = EF_const.transfer * EF_const.kb * diff_Si_O2 / length_Si_O2;
        Rf_O2_on_Si[0] = fcoe * (bound_shortest_O[1][0] - center_Si[0]);//force on x by O2
        Rf_O2_on_Si[1] = fcoe * (bound_shortest_O[1][1] - center_Si[1]);//force on y by O2
        Rf_O2_on_Si[2] = fcoe * (bound_shortest_O[1][2] - center_Si[2]);//force on z by O2

        fcoe = EF_const.transfer * EF_const.kb * diff_Si_O3 / length_Si_O3;
        Rf_O3_on_Si[0] = fcoe * (bound_shortest_O[2][0] - center_Si[0]);//force on x by O3
        Rf_O3_on_Si[1] = fcoe * (bound_shortest_O[2][1] - center_Si[1]);//force on y by O3
        Rf_O3_on_Si[2] = fcoe * (bound_shortest_O[2][2] - center_Si[2]);//force on z by O3

        fcoe = EF_const.transfer * EF_const.kb * diff_Si_O4 / length_Si_O4;
        Rf_O4_on_Si[0] = fcoe * (bound_shortest_O[3][0] - center_Si[0]);//force on x by O4
        Rf_O4_on_Si[1] = fcoe * (bound_shortest_O[3][1] - center_Si[1]);//force on y by O4
        Rf_O4_on_Si[2] = fcoe * (bound_shortest_O[3][2] - center_Si[2]);//force on z by O4

        //sum all the radial forces applied on Si
        Rf_on_Si[0] = Rf_O1_on_Si[0] + Rf_O2_on_Si[0] + Rf_O3_on_Si[0] + Rf_O4_on_Si[0];
        Rf_on_Si[1] = Rf_O1_on_Si[1] + Rf_O2_on_Si[1] + Rf_O3_on_Si[1] + Rf_O4_on_Si[1];
        Rf_on_Si[2] = Rf_O1_on_Si[2] + Rf_O2_on_Si[2] + Rf_O3_on_Si[2] + Rf_O4_on_Si[2];

        //store radial force into forceSi[][3]
        forceSi[i][0] += Rf_on_Si[0];
        forceSi[i][1] += Rf_on_Si[1];
        forceSi[i][2] += Rf_on_Si[2];

        //store radial force into forceO[][3] which just add minus sign to Si forces
        forceO[order_O[0]][0] -= Rf_O1_on_Si[0];//O1, x
        forceO[order_O[0]][1] -= Rf_O1_on_Si[1];//O1, y
        forceO[order_O[0]][2] -= Rf_O1_on_Si[2];//O1, z

        forceO[order_O[1]][0] -= Rf_O2_on_Si[0];//O2, x
        forceO[order_O[1]][1] -= Rf_O2_on_Si[1];//O2, y
        forceO[order_O[1]][2] -= Rf_O2_on_Si[2];//O2, z

        forceO[order_O[2]][0] -= Rf_O3_on_Si[0];//O3, x
        forceO[order_O[2]][1] -= Rf_O3_on_Si[1];//O3, y
        forceO[order_O[2]][2] -= Rf_O3_on_Si[2];//O3, z

        forceO[order_O[3]][0] -= Rf_O4_on_Si[0];//O4, x
        forceO[order_O[3]][1] -= Rf_O4_on_Si[1];//O4, y
        forceO[order_O[3]][2] -= Rf_O4_on_Si[2];//O4, z



        //--------------------------------------------------------------------------------------------------------------

        /*calculate angular energy and angular force
         * and store forces to both Si array and O array */

        //coefficient
        long double diff_Si_O1_O2 = cosSi_O1_O2 - EF_const.cos0_Si;
        long double diff_Si_O1_O3 = cosSi_O1_O3 - EF_const.cos0_Si;
        long double diff_Si_O1_O4 = cosSi_O1_O4 - EF_const.cos0_Si;
        long double diff_Si_O2_O3 = cosSi_O2_O3 - EF_const.cos0_Si;
        long double diff_Si_O2_O4 = cosSi_O2_O4 - EF_const.cos0_Si;
        long double diff_Si_O3_O4 = cosSi_O3_O4 - EF_const.cos0_Si;

        //cout << cosSi_O1_O2 << "    " << cosSi_O1_O3 << "    " << cosSi_O1_O4 << "    "<<cosSi_O2_O3<<"    "<<cosSi_O2_O4<<"    "<<cosSi_O3_O4<<endl;
        //cout << diff_Si_O1_O2 << "    " << diff_Si_O1_O3 << "    " << diff_Si_O1_O4 << "    "<<diff_Si_O2_O3<<"    "<<diff_Si_O2_O4<<"    "<<diff_Si_O3_O4<<endl;

        //angular energy for Si
        angular_energy_Si += 0.5 * EF_const.kthetaSi * EF_const.b0 * EF_const.b0 *
                             (diff_Si_O1_O2*diff_Si_O1_O2 + diff_Si_O1_O3*diff_Si_O1_O3 +
                              diff_Si_O1_O4*diff_Si_O1_O4 + diff_Si_O2_O3*diff_Si_O2_O3 +
                              diff_Si_O2_O4*diff_Si_O2_O4 + diff_Si_O3_O4*diff_Si_O3_O4);

        //angular force for Si
        fcoe = -1.0 * EF_const.transfer * EF_const.kthetaSi * EF_const.b0 * EF_const.b0;

        /* force calculation (Si, O1 and O2) */
        for (int m = 0; m < 3; m++)//force between Si O1 and O2
            fSiO1O2_on_O1[m] = fcoe *((diff_Si_O1_O2 / length_Si_O1) *
                                      ((bound_shortest_O[1][m] - center_Si[m]) / length_Si_O2 - cosSi_O1_O2 * (bound_shortest_O[0][m] - center_Si[m]) / length_Si_O1));
        for (int m = 0; m < 3; m++)//force between Si O1 and O2
            fSiO1O2_on_O2[m] = fcoe * ((diff_Si_O1_O2 / length_Si_O2) *
                                       ((bound_shortest_O[0][m] - center_Si[m]) / length_Si_O1 - cosSi_O1_O2 * (bound_shortest_O[1][m] - center_Si[m]) / length_Si_O2));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
            fSiO1O2_on_Si[m] = -(fSiO1O2_on_O1[m] + fSiO1O2_on_O2[m]);


        /* force calculation (Si, O1 and O3) */
        for (int m = 0; m < 3; m++)//force between Si O1 and O2
            fSiO1O3_on_O1[m] = fcoe * ((diff_Si_O1_O3 / length_Si_O1) *
                                       ((bound_shortest_O[2][m] - center_Si[m]) / length_Si_O3 - cosSi_O1_O3 * (bound_shortest_O[0][m] - center_Si[m]) / length_Si_O1));
        for (int m = 0; m < 3; m++)//force between Si O1 and O2
            fSiO1O3_on_O3[m] = fcoe * ((diff_Si_O1_O3 / length_Si_O3) *
                                       ((bound_shortest_O[0][m] - center_Si[m]) / length_Si_O1 - cosSi_O1_O3 * (bound_shortest_O[2][m] - center_Si[m]) / length_Si_O3));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
            fSiO1O3_on_Si[m] = -(fSiO1O3_on_O1[m] + fSiO1O3_on_O3[m]);


        /* force calculation (Si, O1 and O4) */
        for (int m = 0; m < 3; m++)//force between Si O1 and O4
            fSiO1O4_on_O1[m] = fcoe * ((diff_Si_O1_O4 / length_Si_O1) *
                                       ((bound_shortest_O[3][m] - center_Si[m]) / length_Si_O4 - cosSi_O1_O4 * (bound_shortest_O[0][m] - center_Si[m]) / length_Si_O1));
        for (int m = 0; m < 3; m++)//force between Si O1 and O4
            fSiO1O4_on_O4[m] = fcoe * ((diff_Si_O1_O4 / length_Si_O4) *
                                       ((bound_shortest_O[0][m] - center_Si[m]) / length_Si_O1 - cosSi_O1_O4 * (bound_shortest_O[3][m] - center_Si[m]) / length_Si_O4));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
            fSiO1O4_on_Si[m] = -(fSiO1O4_on_O1[m] + fSiO1O4_on_O4[m]);


        /* force calculation (Si, O2 and O3) */
        for (int m = 0; m < 3; m++)//force between Si O2 and O3
            fSiO2O3_on_O2[m] = fcoe * ((diff_Si_O2_O3 / length_Si_O2) *
                                       ((bound_shortest_O[2][m] - center_Si[m]) / length_Si_O3 - cosSi_O2_O3 * (bound_shortest_O[1][m] - center_Si[m]) / length_Si_O2));
        for (int m = 0; m < 3; m++)//force between Si O2 and O3
            fSiO2O3_on_O3[m] = fcoe * ((diff_Si_O2_O3 / length_Si_O3) *
                                       ((bound_shortest_O[1][m] - center_Si[m]) / length_Si_O2 - cosSi_O2_O3 * (bound_shortest_O[2][m] - center_Si[m]) / length_Si_O3));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O3 on Si
            fSiO2O3_on_Si[m] = -(fSiO2O3_on_O2[m] + fSiO2O3_on_O3[m]);


        /* force calculation (Si, O2 and O4) */
        for (int m = 0; m < 3; m++)//force between Si O2 and O4
            fSiO2O4_on_O2[m] = fcoe * ((diff_Si_O2_O4 / length_Si_O2) *
                                       ((bound_shortest_O[3][m] - center_Si[m]) / length_Si_O4 - cosSi_O2_O4 * (bound_shortest_O[1][m] - center_Si[m]) / length_Si_O2));
        for (int m = 0; m < 3; m++)//force between Si O2 and O4
            fSiO2O4_on_O4[m] = fcoe * ((diff_Si_O2_O4 / length_Si_O4) *
                                       ((bound_shortest_O[1][m] - center_Si[m]) / length_Si_O2 - cosSi_O2_O4 * (bound_shortest_O[3][m] - center_Si[m]) / length_Si_O4));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O4 on Si
            fSiO2O4_on_Si[m] = -(fSiO2O4_on_O2[m] + fSiO2O4_on_O4[m]);


        /* force calculation (Si, O3 and O4) */
        for (int m = 0; m < 3; m++)//force between Si O3 and O4
            fSiO3O4_on_O3[m] = fcoe * ((diff_Si_O3_O4 / length_Si_O3) *
                                       ((bound_shortest_O[3][m] - center_Si[m]) / length_Si_O4 - cosSi_O3_O4 * (bound_shortest_O[2][m] - center_Si[m]) / length_Si_O3));
        for (int m = 0; m < 3; m++)//force between Si O3 and O4
            fSiO3O4_on_O4[m] = fcoe * ((diff_Si_O3_O4 / length_Si_O4) *
                                       ((bound_shortest_O[2][m] - center_Si[m]) / length_Si_O3 - cosSi_O3_O4 * (bound_shortest_O[3][m] - center_Si[m]) / length_Si_O4));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O3 and O4 on Si
            fSiO3O4_on_Si[m] = -(fSiO3O4_on_O3[m] + fSiO3O4_on_O4[m]);

        //store angular force into forceSi[][3]
        for (int m = 0; m < 3; m++)
            forceSi[i][m] += fSiO1O2_on_Si[m] + fSiO1O3_on_Si[m] + fSiO1O4_on_Si[m] +
                             fSiO2O3_on_Si[m] + fSiO2O4_on_Si[m] + fSiO3O4_on_Si[m];


        //store angular force into forceO[][3]
        for (int m = 0; m < 3; m++) {
            forceO[order_O[0]][m] += fSiO1O2_on_O1[m] + fSiO1O3_on_O1[m] + fSiO1O4_on_O1[m];
            forceO[order_O[1]][m] += fSiO1O2_on_O2[m] + fSiO2O3_on_O2[m] + fSiO2O4_on_O2[m];
            forceO[order_O[2]][m] += fSiO1O3_on_O3[m] + fSiO2O3_on_O3[m] + fSiO3O4_on_O3[m];
            forceO[order_O[3]][m] += fSiO1O4_on_O4[m] + fSiO2O4_on_O4[m] + fSiO3O4_on_O4[m];
        }


        // calculate repulsive energy and force between Si and Si
        for(int k=a; k<EF_const.Si_quantity; ++k)
        {
            repSiSi[0] = Si[k][0];
            repSiSi[1] = Si[k][1];
            repSiSi[2] = Si[k][2];

            EF_Ref.shortestcoordinate(center_Si, repSiSi, repSiSishortest);

            repSiSilength = EF_Ref.real_length(center_Si, repSiSishortest);

            if (repSiSilength <= 3.10)
            {
                //calculate energy
                repulsive_energy_Si_Si += 0.5 * EF_const.repulsive_coe_for_SiSi * (3.10 - repSiSilength)*(3.10 - repSiSilength)*(3.10 - repSiSilength);

                //calculate force
                for (int m=0; m<3; ++m)
                {
                    frepSiSi[m] = -3 * EF_const.transfer * EF_const.repulsive_coe_for_SiSi *
                                   (3.10 - repSiSilength)*(3.10 - repSiSilength) / repSiSilength *
                                   (repSiSishortest[m] - center_Si[m]);
                    forceSi[i][m] += frepSiSi[m];
                    forceSi[k][m] -= frepSiSi[m];
                }
            }
        }
        ++a;


    }


    //--------------------------------------------------------------------------------------------------------------



    //calculate angular energy of O
    int b=1;
    for(int i=0;i<EF_const.O_quantity;++i)
    {
        //store O position
        center_O[0] = O[i][0];
        center_O[1] = O[i][1];
        center_O[2] = O[i][2];

        for(int j=0;j<2;++j)
        {
            order_Si[j] = groupO[i][j];
            for(int k=0; k<3; ++k)
                bound_Si[j][k] = Si[order_Si[j]][k];
        }

        for(int j=0;j<2;++j)
            EF_Ref.shortestcoordinate(center_O, bound_Si[j], bound_shortest_Si[j]);

        //length
        length_O_Si1 = EF_Ref.real_length(center_O, bound_shortest_Si[0]);
        length_O_Si2 = EF_Ref.real_length(center_O, bound_shortest_Si[1]);

        //angle
        cosO_Si1_Si2 = EF_Ref.angle(center_O, bound_shortest_Si[0], bound_shortest_Si[1], length_O_Si1, length_O_Si2);

        //calculate angular energy
        long double diff_O_Si1_Si2 = cosO_Si1_Si2 - EF_const.cos0_O;
        angular_energy_O += 0.5 * EF_const.kthetaO * EF_const.b0 * EF_const.b0 * diff_O_Si1_Si2 * diff_O_Si1_Si2;

        //calculate angular force
        long double fcoe = -1.0 * EF_const.transfer * EF_const.kthetaO * EF_const.b0 * EF_const.b0;
        //force calculation (Si, O2 and O4)
        for (int m = 0; m < 3; m++)//force between O Si1 and Si2
            fOSi1Si2_on_Si1[m] = fcoe * ((diff_O_Si1_Si2 / length_O_Si1) *
                                         ((bound_shortest_Si[1][m] - center_O[m]) / length_O_Si2 - cosO_Si1_Si2 * (bound_shortest_Si[0][m] - center_O[m]) / length_O_Si1));
        for (int m = 0; m < 3; m++)//force between Si O2 and O4
            fOSi1Si2_on_Si2[m] = fcoe * ((diff_O_Si1_Si2 / length_O_Si2) *
                                         ((bound_shortest_Si[0][m] - center_O[m]) / length_O_Si1 - cosO_Si1_Si2 * (bound_shortest_Si[1][m] - center_O[m]) / length_O_Si2));
        for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O4 on Si
            fOSi1Si2_on_O[m] = -(fOSi1Si2_on_Si1[m] + fOSi1Si2_on_Si2[m]);

        //store angular force into forceO[][3]
        for (int m = 0; m < 3; m++)
            forceO[i][m] += fOSi1Si2_on_O[m];

        //store angular force into forceSi[][3]
        for (int m = 0; m < 3; m++)
        {
            forceSi[order_Si[0]][m] += fOSi1Si2_on_Si1[m];
            forceSi[order_Si[1]][m] += fOSi1Si2_on_Si2[m];
        }



        //calculate repulsive energy between O and O
        for(int j=b;j<EF_const.O_quantity;++j)
        {
            repOO[0] = O[j][0];
            repOO[1] = O[j][1];
            repOO[2] = O[j][2];

            EF_Ref.shortestcoordinate(center_O, repOO, repOOshortest);

            repOOlength = EF_Ref.real_length(center_O, repOOshortest);

            if(repOOlength <= 2.533)
            {
                //calculate energy
                repulsive_energy_O_O += 0.5*EF_const.repulsive_coe_for_OO*(2.533-repOOlength)*(2.533-repOOlength)*(2.533-repOOlength);

                //calculate force for Si and O
                for(int m=0;m<3;++m)
                {
                    frepOO[m] = -3 * EF_const.transfer * EF_const.repulsive_coe_for_OO * (2.533 - repOOlength)*(2.533-repOOlength) / repOOlength *
                                (repOOshortest[m] - center_O[m]);
                    forceO[i][m] += frepOO[m];
                    forceO[j][m] -= frepOO[m];
                }
            }
        }
        ++b;

    }

    //--------------------------------------------------------------------------------------------------------------

    /*cout << "radial energy: " << radial_energy << endl;
    cout << "Si angular energy: " << angular_energy_Si << endl;
    cout << "O angular energy: " << angular_energy_O << endl;
    cout << "Si repulsive energy: " << repulsive_energy_Si_Si << endl;
    cout << "O repulsive energy: " << repulsive_energy_O_O << endl;*/

    total_energy = radial_energy + angular_energy_Si + angular_energy_O +
                   repulsive_energy_O_O + repulsive_energy_Si_Si;



};