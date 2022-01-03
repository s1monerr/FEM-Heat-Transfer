/****************** /
// * MES PROGRAM - 2021
// * input data: https://home.agh.edu.pl/~pkustra/MES/TC_2d.pdf
// * @Author: Szymon Rewilak
/******************/



#include <iostream>
#include "HBC_matrix.cpp"
#include "file_reader.cpp"
#include "test_no.h"

#define PRECISION 6
#define PRINTABLE false

int main(){
    el_9_2d element_9 = el_9_2d();
    //    macierz pochodnych dN/dksi - 9 wezlow
    double **matrix_ksi_9points = new double*[9];
    for(int i = 0; i < 9; i++){
        matrix_ksi_9points[i] = new double[4];
    }

    // macierz pochodnych dN/dn - 9 wezlow
    double **matrix_n_9points = new double*[9];
    for(int i = 0; i < 9; i++){
        matrix_n_9points[i] = new double[4];
    }
    std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dksi ===="<<std::endl;
    element_9.count_matrix_9points_ksi(matrix_ksi_9points);
    print_custom(matrix_ksi_9points, 9, 4, 4);

    std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dksi ===="<<std::endl;
    element_9.count_matrix_9points_n(matrix_n_9points);
    print_custom(matrix_n_9points, 9, 4, 4);


    grid grid_1 = grid(HEIGTH, WIDTH, NODES_HEIGTH, NODES_WIDTH);
    createNodes(grid_1);
    grid_1.printNodes();
    createElements(grid_1, NODES_HEIGTH);
    grid_1.printElements();

    double **matrix_4points_HBC = new double*[8];

    for(int i = 0; i < 8; i++) matrix_4points_HBC[i] = new double[4];

    el_4_2d_HBC element_HBC = el_4_2d_HBC();
    element_HBC.count_matrix_4points(matrix_4points_HBC);


        // macierz pochodnych dN/dksi - 4 wezly
    double **matrix_ksi_4points = new double*[4];
    for(int i = 0; i < 4; i++){
        matrix_ksi_4points[i] = new double[4];
    }


    // macierz pochodnych dN/dn - 4 wezly
    double **matrix_n_4points = new double*[4];
    for(int i = 0; i < 4; i++){
        matrix_n_4points[i] = new double[4];
    }

    el_4_2d element = el_4_2d("H");
    element.count_matrix_4points_ksi(matrix_ksi_4points);
    element.count_matrix_4points_n(matrix_n_4points);

    std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dksi ===="<<std::endl;
    print_custom(matrix_n_4points, 4, 4, 6);
    std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dksi ===="<<std::endl;
    print_custom(matrix_n_4points, 4, 4, 6);

    ELEMENT_2D_ARRAYS ELEMENT_1 = ELEMENT_2D_ARRAYS(matrix_ksi_9points, matrix_n_9points);
    ELEMENT_2D_ARRAYS ELEMENT_2 = ELEMENT_2D_ARRAYS(matrix_ksi_4points, matrix_n_4points);

    node temp_element[4];
    node n0 = node(0, 0, 1);
    node n1 = node(0.025, 0, 1);
    node n2 = node(0.025, 0.025, 1);
    node n3 = node(0, 0.025, 1);

    temp_element[0] = n0;
    temp_element[1] = n1;
    temp_element[2] = n2;
    temp_element[3] = n3;

    double **H_sum = new double*[4];
        for(int k = 0; k < 4; k++){
            H_sum[k] = new double[4];
        }

    // for(int j = 0; j < 9; j++){
    //     std::cout<<"j = "<<j<<std::endl;
    //             double **jacobian = new double*[2];
    //                 for(int k = 0; k < 2; k++){
    //                     jacobian[k] = new double[2];
    //                 }
    //         jacobi(3, j, jacobian, ELEMENT_1, temp_element);
    //         H_sum = sum_matrix(count_H_matrix(3, j, jacobian, ELEMENT_1, 30, true), H_sum, 4);
    //         for(int k = 0; k < 2; k++)
    //             delete jacobian[k];
    //         delete[] jacobian;
    //     }
    // print_custom(H_sum, 4, 4, 4);
        for(int i = 0; i < 9; i++){
        double **jacobian = new double*[2];
        for(int i = 0; i < 2; i++){
               jacobian[i] = new double[2];
            }
        jacobi(3, i, jacobian, ELEMENT_1, temp_element);
        std::cout<<"JACOBI REVERSED: "<<std::endl;
        // const_multi_matrix(1/determinant(jacobian), jacobian, 2);
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                std::cout<<"   "<<jacobian[i][j]<<"     ";
            }
            std::cout<<std::endl;
        }
    }
    // count_HBC_grid(&grid_1, false, matrix_4points_HBC, ALPHA);
    // grid_1.print_HBC_matrices();


    // count_H_matrix_grid(INTEGRATION_POINTS, &grid_1, ELEMENT_1, CONDUCTIVITY, true);

    //BACKUP
    // InitData data = InitData();
    // read_file(TEST_BASIC, &data);
    // data.print();
    // final_alghoritm_2POINTS(INTEGRATION_POINTS, data, PRINTABLE, PRECISION);
    // return 0;
} 