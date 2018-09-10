//
// Created by ZHI on 5/23/2018.
//

#include "IO.h"
#include "CONSTANT.h"

CONSTANT IO_const;


IO::IO(){};

/* store all atoms positions into 2D arrays */
void IO::Import(long double Si[][3], long double O[][3])
{
    string filename;
    filename = IO_const.filename;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile.is_open()) // failed to open file
    {
        cout << "Could not open the file " << filename << endl;
    }
    else{

        string temp;//store useless line
        long double file_position;//get position from file and import to the array
        for (int i = 0; i< 8+(IO_const.total_quantity);i++)
        {
            if (i>7 && i <= (7+IO_const.O_quantity) )
            {

                inFile >> file_position;
                O[i-8][0] = file_position;
                inFile >> file_position;
                O[i-8][1] = file_position;
                inFile >> file_position;
                O[i-8][2] = file_position;

            }
            else if (i>(7+IO_const.O_quantity) )
            {
                inFile >> file_position;
                Si[i-8-IO_const.O_quantity][0] = file_position;
                inFile >> file_position;
                Si[i-8-IO_const.O_quantity][1] = file_position;
                inFile >> file_position;
                Si[i-8-IO_const.O_quantity][2] = file_position;

            }
            getline(inFile,temp);

        };
    };

};



void IO::Output(long double Si[][3], long double O[][3], int i)
{
    string filename;
    filename = IO_const.filename;
    /*cout << "Please input the file's name: "<<endl;
    getline(cin,filename);*/
    ifstream inFile;
    inFile.open(filename);
    if (!inFile.is_open()) // failed to open file
    {
        cout << "Could not open the file " << filename << endl;
    }
    else
    {
        string temp;//store line

        ofstream Savefile(IO_const.outfilename+to_string(i));//create a file to output data
        Savefile.setf(ios::fixed,ios::floatfield);
        Savefile.precision(12);
        for (int i = 0; i < 8; i++)
        {
            getline(inFile,temp);
            Savefile << temp;
        }
        Savefile << endl;
        Savefile.clear();
        //output position
        for(int a = 0; a<2; a++)
        {
            if (a == 0)
            {
                for(int b = 0; b < IO_const.O_quantity; b++)
                {

                    Savefile << O[b][0];
                    Savefile << "    ";
                    Savefile << O[b][1];
                    Savefile << "    ";
                    Savefile << O[b][2];
                    Savefile << "    ";
                    Savefile << endl;
                }
            }
            else
            {
                for(int b = 0; b < IO_const.Si_quantity; b++)
                {

                    Savefile << Si[b][0];
                    Savefile << "    ";
                    Savefile << Si[b][1];
                    Savefile << "    ";
                    Savefile << Si[b][2];
                    Savefile << "    ";
                    Savefile << endl;

                }

            };
        };
    }
};
