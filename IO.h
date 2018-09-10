//
// Created by ZHI on 5/23/2018.
//

#ifndef REVISEDSIO2_IO_H
#define REVISEDSIO2_IO_H

#include <fstream>

class IO {
public:
    IO();
    void Import(long double Si[][3], long double O[][3]);
    void Output(long double Si[][3], long double O[][3], int i);

};


#endif //REVISEDSIO2_IO_H
