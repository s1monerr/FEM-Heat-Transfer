/****************** /
// * FEM PROGRAM - 2021/22
// * @Author: Szymon Rewilak
// * input data: https://home.agh.edu.pl/~pkustra/MES/TC_2d.pdf
/******************/


#include <iostream>
#include "FEM_lib.cpp"
#include "file_reader.cpp"
#include "test_no.h"

int main(){
    InitData data = InitData();
    read_file(TEST_BASIC, &data);

    final_alghoritm(INTEGRATION_POINTS, data, PRINTABLE, PRECISION);

    return 0;
}   