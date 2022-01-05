/****************** /
// * MES PROGRAM - 2021/22
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
    std::cout<<"matrix[00] = "<<matrix_ksi_9points[0][0]<<std::endl;

    std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dN ===="<<std::endl;
    element_9.count_matrix_9points_n(matrix_n_9points);
    print_custom(matrix_n_9points, 9, 4, 4);


    grid grid_1 = grid(HEIGTH, WIDTH, NODES_HEIGTH, NODES_WIDTH);
    createNodes(grid_1);
    // grid_1.printNodes();
    createElements(grid_1, NODES_HEIGTH);
    // grid_1.printElements();

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

    el_4_2d element = el_4_2d();
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

    // double **H_sum = new double*[4];
    //     for(int k = 0; k < 4; k++){
    //         H_sum[k] = new double[4];
    //     }

    // for(int j = 0; j < 9; j++){
    //     std::cout<<"j = "<<j<<std::endl;
    //             double **jacobian = new double*[2];
    //                 for(int k = 0; k < 2; k++){
    //                     jacobian[k] = new double[2];
    //                 }
    //         jacobi(3, j, jacobian, ELEMENT_1, temp_element);
    //         const_multi_matrix(1/determinant(jacobian), jacobian, 2);
    //         std::cout<<"jacobi for pc "<<j<<std::endl;
    //         print_custom(jacobian,2,2,4);
    //         H_sum = sum_matrix(count_H_matrix(3, j, jacobian, ELEMENT_1, 30, true), H_sum, 4);
    //         // std::cout<<"Summed H after iteration\n:";
    //         // print_custom(H_sum, 4, 4, 4);
    //         for(int k = 0; k < 2; k++)
    //             delete jacobian[k];
    //         delete[] jacobian;
    //     }
    // print_custom(H_sum, 4, 4, 4);
    //     for(int i = 0; i < 9; i++){
    //     double **jacobian = new double*[2];
    //     for(int i = 0; i < 2; i++){
    //            jacobian[i] = new double[2];
    //         }
    //     jacobi(3, i, jacobian, ELEMENT_1, temp_element);
    //     std::cout<<"JACOBI REVERSED: "<<std::endl;
    //      const_multi_matrix(1/determinant(jacobian), jacobian, 2);
    //     for(int i = 0; i < 2; i++){
    //         for(int j = 0; j < 2; j++){
    //             std::cout<<"   "<<std::setprecision(4)<<jacobian[i][j]<<"     ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // count_HBC_grid(&grid_1, false, matrix_4points_HBC, ALPHA);
    // grid_1.print_HBC_matrices();
    // count_H_matrix_grid(2, &grid_1, ELEMENT_2, 25, true);

    // count_H_matrix_grid(3, &grid_1, ELEMENT_1, 25, true);


    // count_H_matrix_grid(INTEGRATION_POINTS, &grid_1, ELEMENT_1, CONDUCTIVITY, true);

    // el_4_2d_HBC elhbc = el_4_2d_HBC();

    // double **matrix_9points_HBC = new double*[12];
    // for(int i = 0; i < 12; i++) matrix_9points_HBC[i] = new double[4];

    // elhbc.count_matrix_9points(matrix_9points_HBC);

    // // count_HBC_grid(2, &grid_1, true, matrix_4points_HBC, 25);

    // count_HBC_grid(3, &grid_1, true, matrix_9points_HBC, 25);

    // el_4_2d_C element_C = el_4_2d_C();
    //     double **matrix_4points_C = new double*[4];
    //     for(int i = 0; i < 4; i++) matrix_4points_C[i] = new double[4];

    // el_9_2d_C element_C_3 = el_9_2d_C();
    //     double **matrix_9points_C = new double*[9];
    //     for(int i = 0; i < 9; i++) matrix_9points_C[i] = new double[4];


    // element_C_3.count_matrix_9points(matrix_9points_C);
    // std::cout<<"C shape values 3p: \n";
    // print_custom(matrix_9points_C, 9, 4, 4);


    
    // element_C.count_matrix_4points(matrix_4points_C);
    // double **matrix_test_C = new double*[4];
    //     for(int i = 0; i < 4; i++) matrix_test_C[i] = new double[4];
    
    // // matrix_test_C = count_C_matrix(2, temp_element, matrix_4points_C, 700, 7800, ELEMENT_2, true);
    // matrix_test_C = count_C_matrix(3, temp_element, matrix_9points_C, 700, 7800, ELEMENT_1, true);

    // std::cout<<"C summed: "<<std::endl;
    // print_custom(matrix_test_C, 4, 4, 6);

    double **HBC_3p_values = new double*[12];
    for(int i = 0; i < 12; i++)HBC_3p_values[i] = new double[4];

    el_9_2d_HBC el3 = el_9_2d_HBC();
    el3.count_matrix_9points(HBC_3p_values);

    print_custom(HBC_3p_values, 12, 4, 6);

    // count_bound_HBC(2, temp_element, matrix_4points_HBC, 25, true);
    count_bound_HBC(3, temp_element, HBC_3p_values, 25, true);


    
    //BACKUP
    // InitData data = InitData();
    // read_file(TEST_KWADRAT_31, &data);
    // data.print();
    // final_alghoritm(INTEGRATION_POINTS, data, PRINTABLE, PRECISION);
    return 0;
} 