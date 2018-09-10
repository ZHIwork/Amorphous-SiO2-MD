#include <iostream>
#include "mpi.h"
#include <stdio.h>
#include "IO.h"
#include "CONSTANT.h"
#include "GROUP.h"
#include "ENERGY_FORCE.h"
#include "DYNAMICS.h"

IO main_IO;
CONSTANT main_const;
GROUP main_group;
ENERGY_FORCE main_EF;
DYNAMICS main_dy;

using namespace std;


const int Si_quantity = main_const.Si_quantity;
const int O_quantity = main_const.O_quantity;

void init2darray3(long double a[][3], int n);
void initgraph(long double a[][128], int n);

int main(int argc, char** argv) {

    //set up MPI
    long double oldenergy=0;
    int size, rank;
    MPI_Status status;
    struct{
        long double energy;
        int rank;
    }total_energy, lowenergy;

    lowenergy.energy = 0;

    long double Si[Si_quantity][3], O[O_quantity][3];
    init2darray3(Si,Si_quantity);
    init2darray3(O,O_quantity);

    long double oldSi[Si_quantity][3], oldO[O_quantity][3];
    init2darray3(oldSi,Si_quantity);
    init2darray3(oldO,O_quantity);

    long double Graph[Si_quantity][128];
    initgraph(Graph, Si_quantity);

    long double forceSi[Si_quantity][3], forceO[O_quantity][3];
    init2darray3(forceSi,Si_quantity);
    init2darray3(forceO,O_quantity);

    long double vSi[Si_quantity][3], vO[O_quantity][3];
    init2darray3(vSi,Si_quantity);
    init2darray3(vO,O_quantity);

    int groupSi[Si_quantity][4], groupO[O_quantity][2];
    int accept=0, i, j;

    //---------------------------------------------------------------------------------------------
    MPI_Init(NULL, NULL);
    // Get the number of processes
    size = MPI::COMM_WORLD.Get_size();
    // Get the rank of the process
    rank = MPI::COMM_WORLD.Get_rank();

    //---------------------------------------------------------------------------------------------

    //obtain coordinates of Si and O
    main_IO.Import(Si, O);
    //Graph Si and O;
    main_group.Group(Si, O, Graph);

    //---------------------------------------------------------------------------------------------
    for(i=0; i<main_const.cut_times; ++i)
    {

        //start MPI
        //master rank 0
        if (rank == 0)
        {
            //cout<<"Times: "<<i<<"  Total energy: "<<total_energy.energy<<" lowenergy: "<<lowenergy.energy<<"\n";
			total_energy.energy = 1e7;
            total_energy.rank = rank;
        }
            //slaves other ranks
        else 
        {
            //value rank
            total_energy.rank = rank;
            //cout << "Start" << endl;

            //cut bonds
            main_dy.cutbond(Graph);
            //group
            main_EF.GroupfromGraph(groupSi, groupO, Graph);

            //initilaze force and velocity;
            init2darray3(forceSi,Si_quantity);
            init2darray3(forceO,O_quantity);
            init2darray3(vSi,Si_quantity);
            init2darray3(vO,O_quantity);

            //cycle(minimize energy)
            for (j = 0; j < main_const.minimize_times; ++j) {
                //calculate energy(update energy in total_energy.energy)
                main_EF.Force_Energy(Si, O, forceSi, forceO, total_energy.energy, groupSi, groupO);
                //update coordinate
                main_dy.MD(forceSi, forceO, vSi, vO, Si, O, main_const.time_interval);
            }
	    //cout<<" Total energy:"<<total_energy.energy<<"    Old energy:"<<oldenergy<<"\n";
            bool p = main_dy.MC_probability(oldenergy, total_energy.energy);
            if (p == 0)//use original coordinates
                total_energy.energy = 1e7;
        }
        //obtain all total energy from each node and find the lowest one and its rank
        MPI_Allreduce(&total_energy, &lowenergy, 1, MPI_LONG_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        if (rank == lowenergy.rank && rank != 0) {
              //update all arrays needed
              //cout << "low"<<"\n";
		//tag=0
                MPI::COMM_WORLD.Send(&Graph, Si_quantity * 128, MPI::LONG_DOUBLE, 0, 0);
		//tag=1
                MPI::COMM_WORLD.Send(&Si, Si_quantity * 3, MPI::LONG_DOUBLE, 0, 1);
		//tag=2
                MPI::COMM_WORLD.Send(&O, O_quantity * 3, MPI::LONG_DOUBLE, 0, 2);

                main_IO.Output(Si, O, i);
             
        }

        if (rank == 0 && lowenergy.rank!=0) {
            //receive arrays from lowest energy's rank
            //cout << "Receive start"<<"\n";
            MPI::COMM_WORLD.Recv(&Graph, Si_quantity * 128, MPI::LONG_DOUBLE, MPI::ANY_SOURCE, 0);
            //cout << "Receive 1"<<"\n";	
            MPI::COMM_WORLD.Recv(&Si, Si_quantity * 3, MPI::LONG_DOUBLE, MPI::ANY_SOURCE, 1);
            MPI::COMM_WORLD.Recv(&O, O_quantity * 3, MPI::LONG_DOUBLE, MPI::ANY_SOURCE, 2);
	    	oldenergy = lowenergy.energy;
	    	cout << "Accpeted Energy: "<<oldenergy << "\n";
	    	++accept;
        }

        MPI::COMM_WORLD.Bcast(&Graph, Si_quantity * 128, MPI::LONG_DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&Si, Si_quantity * 3, MPI::LONG_DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&O, O_quantity * 3, MPI::LONG_DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&oldenergy, 1, MPI::LONG_DOUBLE, 0);

    }

    MPI::COMM_WORLD.Barrier();
    if(rank == 0)
    {
        //output file
       // main_IO.Output(Si, O);
        //print accept times
        cout << "Accepted times are: "<<accept<<endl;
    }

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI::Finalize();

    return 0;
}

//----------------------------------------------------------------------------------------------------

void init2darray3(long double a[][3], int n)
{
      for(int i=0; i<n; ++i)
      {
          for(int j=0; j<3; ++j)
              a[i][j] = 0;
      }
};

void initgraph(long double a[][128], int n)
{
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<128; ++j)
            a[i][j] = 0;
    }
};
