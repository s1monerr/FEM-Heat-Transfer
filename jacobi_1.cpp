// IMPLEMENTACJA JAKOBIANU W SIATCE
// schemat 2 punktowy calkowania
// PROGRAM OBLICZA JAKOBIAN DLA KAZDEGO ELEMENTU W SIATCE
// NASTEPNIE MACIERZ H DLA KAZDEGO ELEMENTU W SIATCE W KAZDYM PUNKCIE CALKOWANIA
// ORAZ MACIERZ WARUNKU BRZEGOWEGO Hbc

#include <iostream>
#include <iomanip>
#include <cmath>
// #include "grid.h"

#define EPSILON 0.00000001
// PARAMETRY SIATKI
//#ifndef WIDTH
#define WIDTH 0.1 // szerokosc
#define HEIGTH 0.2 // wysokosc
#define NODES_HEIGTH 5 // ilosc wezlow wysokosc
#define NODES_WIDTH 4 // ilosc wezlow szerokosc
//#endif

void print(double **, int);

class node{

    public:
        float x, y;
        int weight;

    node(float X, float Y, int W){
        x = X;
        y = Y;
        weight = W;
    }

    node(float X, float Y){
        x = X;
        y = Y;
        weight = -1;
    }

    node(){
        ;
    }
};

class element{
    public:
        int id[4]; // node's id
        
        double **H; // H matrix for element

        element(int i1, int i2, int i3, int i4){
            id[0] = i1;
            id[1] = i2;
            id[2] = i3;
            id[3] = i4;
            for(int i = 0; i < 4; i++){
                H[i] = new double[4];
            }
        }

        element(){;}

        void print_H_matrix(){
            print(H, 4);
        }
};

// SIATKA
class grid{

    public:
    // wysokosc i szerokosc siatki
        float heigth;
        float width;
        int nH, nW; // ilosc wezlow -> szerokosc i wysokosc
        int nodes_number, elm_number; // liczba wezlow; elementow

        // tablice wezlow i elementow
        node ** nodes;
        element ** elements;

    // constructor
    grid(float high, float widt, int number_H,  int number_W){
        heigth = high;
        width = widt;
        nH = number_H;
        nW = number_W;
        nodes_number = nH*nW;
        elm_number = (nH - 1)*(nW - 1);

        // tablice elementow i wezlow 2d
        nodes = new node*[nW];
        for(int i = 0; i < nW; i++){
            nodes[i] = new node[nH];
        }

        elements = new element *[nW - 1];
        for(int i = 0; i < nW-1; i++){
            elements[i] = new element[nH - 1];
        }
    }

    void printElements(){
        int c = 0;
        std::cout<<std::endl<<"ELEMENTS NUMBER: "<<elm_number<<std::endl;
        std::cout<<"================="<<std::endl;
        for(int i = 0; i < nW - 1; i++){
            for(int j = 0; j < nH - 1; j++){
                std::cout<<"Element "<<c++<<" : ";
                std::cout<<"id1 = "<<elements[i][j].id[0]<<", id2 = "<<elements[i][j].id[1]<<", id3 = "
                <<elements[i][j].id[2]<<", id4 = "<<elements[i][j].id[3]<<std::endl;
            }
        }
        std::cout<<std::endl;
    }

    void printNodes(){
        std::cout<<"NODES NUMBER: "<<nodes_number<<std::endl;
        std::cout<<"================="<<std::endl;
        int counter = 0;
        for(int i = 0; i < nW; i++){
            for(int j = 0; j < nH; j++){
                std::cout<<"Node  "<<std::setprecision(2)<<counter++<<", x = "<<nodes[i][j].x<<", y = "<<nodes[i][j].y<<"   |    weight = "<<nodes[i][j].weight<<std::endl;
            }
        }
    }

    void print_H_matrices(){
        std::cout<<"=================== H MATRICES FOR GRID ==================="<<std::endl;
        int counter = 0;
        for(int i = 0; i < nW-1; i++){
            for(int j = 0; j < nH-1; j++){
                std::cout<<"ELEMENT "<<counter<<std::endl;
                elements[i][j].print_H_matrix();
                counter++;
        }
    }
    }

};

// funkcja do uzupelnienia tablicy z wezlami
void createNodes(grid gr){
    // wyznaczanie dx i dy - szerokosc i wysokosc podzielona na ilosc wezlow - 1
    float dx = gr.width/(gr.nW-1);
    float dy = gr.heigth/(gr.nH-1);

    for(int i = 0; i < gr.nW; i++){
        float y = 0.0f; // w kazdej kolumnie zaczyna od dolu do gory
        float x = 0.0f + i*dx; // kolejna kolumna - zwieksz x

        // interior loop - iteracja po rzedach
        for(int j = 0; j < gr.nH; j++){
            if(i == 0 || j == 0 | i == gr.nW - 1 || j == gr.nH - 1) // jesli wezel jest skrajny w siatce
                gr.nodes[i][j] = node(x, y, 1); // dodaj nowy wezel do tablicy - waga = 1
            else
                gr.nodes[i][j] = node(x, y, 0);
            y += dy; // zwieksz y o dy
        }
    }
}


// funckja do uzupelnienia tablicy elementow
void createElements(grid gr, int node_heigth){
    int n = gr.elm_number;
    int w = node_heigth; // zmienna do zwiekszania id w kolejnej kolumnie
    // MAIN LOOP - columns
    for(int i = 0; i < gr.nW - 1; i++){
        // INTERIOR LOOP - rows
        for(int j = 0; j < gr.nH - 1; j++){
            gr.elements[i][j].id[0] = j + w - node_heigth; // zawsze dodawaj w, ale mmiejsze o jedna wysokosc od id2
            gr.elements[i][j].id[1] = j + w;
            // kolejne id sa odpowiednio wieksze o 1
            gr.elements[i][j].id[2] = gr.elements[i][j].id[1] + 1;
            gr.elements[i][j].id[3] = gr.elements[i][j].id[0] + 1;
        }
        w += node_heigth; // w kolejnej kolumnie zwieksz o kolejna liczbe wezlow w kolumnie - czyli wysokosc wezlow
    }

}

// klasa z wspolrzednymi calkowania roznych schematow
class factors{
    public:
        double w_2p[2]; // wektor wag schematu 2 punkt.
        double w_3p[3]; // wektor wag schematu 3 punkt.
        
        double point_2p[2]; // wektor punktow schematu 2 punkt.
        double point_3p[3]; // wektor punktow schematu 3 punkt.

    // w konstruktorze uzupelnia wektory o wartosci kwadratur Gaussa
    factors(){
        // schemat 2 punktowy
        w_2p[0] = 1; // wagi
        w_2p[1] = 1;
        point_2p[0] = -sqrt(3)/3; // wezly (punkty)
        point_2p[1] = -point_2p[0];
        

        // schemat 3 punktowy
        w_3p[0] = 0.55555555555; // wagi
        w_3p[1] = 0.88888888888;
        w_3p[2] = w_3p[0];
        point_3p[0] = -sqrt(0.6); // wezly
        point_3p[1] = 0;
        point_3p[2] = -point_3p[0];
    }    
};

class point{
    public:
        double n, ksi;

        point(){
            n = 0;
            ksi = 0;
        }
        point(double xx, double yy){
            n = xx;
            ksi = yy;
        }
};

class el_4_2d{

factors fac;

    public:

        point *array;
        // dN/dPsi - funkcja wybierana przez 2 argument
        double integral_ksi(point p, int index){
            if(index == 0) return -0.25*(1-p.n);
            
            else if (index == 1) return 0.25*(1-p.n);
        
            else if(index == 2) return 0.25*(1+p.n);
        
            else if (index == 3) return -0.25*(1+p.n);

            return -1;
        }

        double integral_n(point p, int index){
            if(index == 0) return -0.25*(1-p.ksi);
            
            else if (index == 1) return -0.25*(1+p.ksi);
        
            else if(index == 2) return 0.25*(1+p.ksi);
        
            else if (index == 3) return 0.25*(1-p.ksi);

            return -1;
        }

        el_4_2d(std::string FACTOR_CONST){ // rodzaj macierzy pochodnych - HBC, jesli dla obliczen macierzy obciazen
            fac = factors();
            // uzupelnij tablice punktami calkowania
            if(FACTOR_CONST == "H"){
                array = new point[4];
                array[0] = point(fac.point_2p[0], fac.point_2p[0]);
                array[1] = point(fac.point_2p[0], fac.point_2p[1]);
                array[2] = point(fac.point_2p[1], fac.point_2p[1]);
                array[3] = point(fac.point_2p[1], fac.point_2p[0]);
            }
            else if(FACTOR_CONST == "HBC"){
                array = new point[4];
                array[0] = point(-1, fac.point_2p[0]);
                std::cout<<"point 1 ksi = "<<array[0].ksi<<" n = "<<array[0].n<<std::endl;
                array[1] = point(-1, fac.point_2p[1]);
                std::cout<<"point 2 ksi = "<<array[1].ksi<<" n = "<<array[1].n<<std::endl;
                array[2] = point(fac.point_2p[0], 1);
                std::cout<<"point 3 ksi = "<<array[2].ksi<<" n = "<<array[2].n<<std::endl;
                array[3] = point(fac.point_2p[1], 1);
                std::cout<<"point 4 ksi = "<<array[3].ksi<<" n = "<<array[3].n<<std::endl;
                array[4] = point(1, fac.point_2p[1]);
                std::cout<<"point 5 ksi = "<<array[4].ksi<<" n = "<<array[4].n<<std::endl;
                array[5] = point(1, fac.point_2p[0]);
                std::cout<<"point 6 ksi = "<<array[5].ksi<<" n = "<<array[5].n<<std::endl;
                array[6] = point(fac.point_2p[1], -1);
                std::cout<<"point 7 ksi = "<<array[6].ksi<<" n = "<<array[6].n<<std::endl;
                array[7] = point(fac.point_2p[0], -1);
                std::cout<<"point 8 ksi = "<<array[7].ksi<<" n = "<<array[7].n<<std::endl;
            }
        }

        void count_matrix_4points_ksi(double **matrix){
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    matrix[i][j] = integral_ksi(array[i], j);
                }
            }
        }

        void count_matrix_4points_n(double **matrix){
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    matrix[i][j] = integral_n(array[i], j);
                }
            }
        }

        void count_matrix_4points_ksi_HBC(double **matrix){
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 4; j++){
                    matrix[i][j] = integral_ksi(array[i], j);
                }
            }
        }

        void count_matrix_4points_n_HBC(double **matrix){
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 4; j++){
                    matrix[i][j] = integral_n(array[i], j);
                }
            }
        }

};

// klasa - macierz jakobiego dla elementu 9 wezlowego
class el_9_2d{

factors fac;

    public:

        point array[9];
        // dN/dPsi - funkcja wybierana przez 2 argument
        double integral_ksi(point p, int index){
            if(index == 0) return -0.25*(1-p.n);
            
            else if (index == 1) return 0.25*(1-p.n);
        
            else if(index == 2) return 0.25*(1+p.n);
        
            else if (index == 3) return -0.25*(1+p.n);

            return -1;
        }

        double integral_n(point p, int index){
            if(index == 0) return -0.25*(1-p.ksi);
            
            else if (index == 1) return -0.25*(1+p.ksi);
        
            else if(index == 2) return 0.25*(1+p.ksi);
        
            else if (index == 3) return 0.25*(1-p.ksi);

            return -1;
        }

        el_9_2d(){
            fac = factors();
            int counter = 0;
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    array[counter++] = point(fac.point_3p[i], fac.point_3p[j]);
                }
            }
        }

        void count_matrix_9points_ksi(double **matrix){
            for(int i = 0; i < 9; i++){
                for(int j = 0; j < 4; j++){
                    matrix[i][j] = integral_ksi(array[i], j);
                }
            }
        }

        void count_matrix_9points_n(double **matrix){
            for(int i = 0; i < 9; i++){
                for(int j = 0; j < 4; j++){
                    matrix[i][j] = integral_n(array[i], j);
                }
            }
        }

};

void print(double **matrix, int n){
    std::cout<<std::endl;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            std::cout<<"["<<i<<"]["<<j<<"] "<<matrix[i][j]<<"  ;  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print(double **matrix, int n, int precision){
    std::cout<<std::endl;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            std::cout<<"["<<i<<"]["<<j<<"] "<<std::setprecision(precision)<<matrix[i][j]<<"  ;  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print_9(double **matrix){
    std::cout<<std::endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 9; j++){
            std::cout<<"["<<i<<"]["<<j<<"], "<<matrix[j][i]<<"  ;  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

// final arrays from el_4_2d
class ELEMENT_2D_ARRAYS{
    public:
        double **MATRIX_4P_DKSI;
        double **MATRIX_4P_DN;

        ELEMENT_2D_ARRAYS(double **arr1, double **arr2){
            MATRIX_4P_DKSI = arr1;
            MATRIX_4P_DN = arr2;
        }
};

double determinant(double **matrix) { return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]; }


// COUNT JACOBI MATRIX
// node elmnt[4] - reprezentacja elementu z 4 wezlami w siatce
// iter - iteracja -> iterowany punkt calkowania 
// elm_number -> numer elementu
void jacobi(/*int elm_number,*/ int iter, double** J, double** J_inv, ELEMENT_2D_ARRAYS elmnt_4p, node elmnt[4]){
    // count values in jacobi matrix - punkt calkowania pc[iter]
   // std::cout<<"GAUSS: pc = "<<iter<<std::endl;

    for(int i = 0; i < 4; i++){
        //std::cout<<elmnt_4p.MATRIX_4P_DKSI[iter][i]<<"*"<<elmnt[i].x<<std::endl;
        J[0][0] +=  elmnt_4p.MATRIX_4P_DKSI[iter][i]*elmnt[i].x;
    }
   // std::cout<<std::endl;

    for(int i = 0; i < 4; i++){
     //   std::cout<<elmnt_4p.MATRIX_4P_DKSI[iter][i]<<"*"<<elmnt[i].y<<std::endl;
        J[0][1] +=  elmnt_4p.MATRIX_4P_DKSI[iter][i]*elmnt[i].y;
    }
    //std::cout<<std::endl;

    for(int i = 0; i < 4; i++){
      //  std::cout<<elmnt_4p.MATRIX_4P_DN[iter][i]<<"*"<<elmnt[i].x<<std::endl;
        J[1][0] +=  elmnt_4p.MATRIX_4P_DN[iter][i]*elmnt[i].x;
    }
    //std::cout<<std::endl;

    for(int i = 0; i < 4; i++){
      //  std::cout<<elmnt_4p.MATRIX_4P_DN[iter][i]<<"*"<<elmnt[i].y<<std::endl;
        J[1][1] +=  elmnt_4p.MATRIX_4P_DN[iter][i]*elmnt[i].y; 
    }

    // print results
    // std::cout<<std::endl;
    // std::cout<<"JACOBI MATRIX: pc = "<<iter<<std::endl;
    // for(int i = 0; i < 2; i++){
    //     for(int j = 0; j < 2; j++){
    //         if(J[i][j] < EPSILON)
    //             J[i][j] = 0;
    //         std::cout<<std::setprecision(5)<<"["<<j<<"]["<<i<<"], "<<J[i][j]<<"  ;  ";
    //     }
    //     std::cout<<std::endl;
    // }

    // std::cout<<std::setprecision(5)<<"DETERMINANT: "<<determinant(J)<<std::endl;
    // std::cout<<std::setprecision(5)<<"1/DETERMINANT: "<<1/determinant(J)<<std::endl;
    // std::cout<<std::endl;
}

// COUNT JACOBI MATRIX FOR GRID
// ELEMENT_1 - przechowuje macierze dN/dKsi oraz dN/dn
void jacobi_grid(grid GRID, ELEMENT_2D_ARRAYS ELEMENT_1){
    // tablica 1d wezlow w siatce
    node *nodes = new node[GRID.nodes_number];
    element *elements = new element[GRID.elm_number];
    // przepisywanie wezlow do tablicy 1d
    int counter = 0;
    for(int i = 0; i < GRID.nW; i++){
        for(int j = 0; j < GRID.nH; j++){
            nodes[counter] = GRID.nodes[i][j];
                //std::cout<<"NODES ["<<counter<<"] : x = "<<nodes[counter].x<<" y = "<<nodes[counter].y<<std::endl;
                counter++;
        }
    }
    
    // przepisywanie elementow do tablicy 1d
    counter = 0;
    for(int i = 0; i < GRID.nW-1; i++){
        for(int j = 0; j < GRID.nH-1; j++){
            elements[counter] = GRID.elements[i][j];
                /*std::cout<<"ELEMENTS ["<<counter<<"] : id1 = "<<elements[counter].id[0]
                <<" id2 = "<<elements[counter].id[1]<<" id3 = "<<elements[counter].id[2]<<
                " id4 = "<<elements[counter].id[3]<<std::endl;*/
                counter++;
        }
    }

    // deklaracja macierzy do przechowania Jakobianu
    double **jacobian = new double*[2];
    for(int i = 0; i < 2; i++)
        jacobian[i] = new double[2];
    
    // MAIN ALGHORITM
    // oblicz Jakobian dla kazdego elementu -> elemnt ma id wezlow
    for(int i = 0; i < GRID.elm_number; i++){
        // dla kazdego elementu stworz temp_element z 4 wezlami tego elementu
        
        node temp_element[4];
        temp_element[0] = nodes[elements[i].id[0]];
        temp_element[1] = nodes[elements[i].id[1]];
        temp_element[2] = nodes[elements[i].id[2]];
        temp_element[3] = nodes[elements[i].id[3]];
        std::cout<<"=== COUNTING JACOBI MATRIX FOR ELEMENT "<<i<<" ==="<<std::endl;
        // inner loop -> 4 punkty calkowania dla kazdego 
        for(int j = 0; j < 4; j++){
                double **jacobian = new double*[2];
                    for(int k = 0; k < 2; k++){
                        jacobian[k] = new double[2];
                    }
            jacobi(/*i,*/j, jacobian, jacobian, ELEMENT_1, temp_element);
                    for(int k = 0; k < 2; k++){
                        delete jacobian[k];
                    }
                    delete[] jacobian;
        }
    }
    delete[] nodes;
    delete[] elements;
}

// *** COUNTING H MATRIX *** //
void const_multi_matrix(double const_value, double **matrix, int size){    // function to multiply matrix by const
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            matrix[i][j] *= const_value;
        }
    }
}

// multiplication of matrice by transponed matrice
double **multi_matrix_transponed(int size, double *matrix) {
    double **result = new double*[size];
    for(int i = 0; i < size; i++){
        result[i] = new double[size];
    }
 
//  std::cout<<"counting: result "<<i<<"  "<<j<<" matrix[j] = "<<matrix[j]<<" matrixx[i]"<<matrix[k]<<std::endl;
//                 result[i][k] += matrix[k] * matrix[k];
//                 std::cout<<"RESULT[i][j] = "<<result[i][j]<<std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = matrix[i]*matrix[j];
        }
    }
    return result;
}

// summing square matrices
double **sum_matrix(double **matrix_1, double **matrix_2, int size){
    double **result = new double*[size];
    for(int i = 0; i < size; i++){
        result[i] = new double[size];
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < size; j++){
            result[i][j] = matrix_1[i][j] + matrix_2[i][j];
        }
    }
    return result;
}

// FUNCTION TO ADD ROW TO H MATRIX
// matrix_dx, matrix_dy - temp matrices, iteration - punkt calkowania, element - integral dksi dn arrays
void add_row(double **matrix_dx, double **matrix_dy, int iteration, double **jacobi, ELEMENT_2D_ARRAYS element){

    double row[4]; // final row
    // temp structures for f
    double row_temp_x[4] = {0, 0, 0, 0};
    double row_temp_y[4] = {0, 0, 0, 0};

    int counter = 0; // row iterator
    for(int i = 0; i < 4; i++){ // main loop - 4 integrals, inner loop - multiply dn/dX dn/dY by jacobi elements
        for(int j = 0; j < 2; j++){
          //  std::cout<<"counter = "<<counter<<"j = "<<j<<" jacobi[0][j] = "<<jacobi[0][j]<<" element = "<<element.MATRIX_4P_DKSI[iteration][j]<<std::endl;
            if(j == 0)
                row_temp_x[counter] += jacobi[0][j]*element.MATRIX_4P_DKSI[iteration][i];
            else if(j == 1)
                row_temp_x[counter] += jacobi[0][j]*element.MATRIX_4P_DN[iteration][i];
          //  std::cout<<"counter = "<<counter<<"j = "<<j<<" jacobi[1][j] = "<<jacobi[1][j]<<" element = "<<element.MATRIX_4P_DKSI[iteration][j]<<std::endl;
            if(j == 0)
                row_temp_y[counter] += jacobi[1][j]*element.MATRIX_4P_DKSI[iteration][i];
            else if(j == 1)
                row_temp_y[counter] += jacobi[1][j]*element.MATRIX_4P_DN[iteration][i];
        }
       // std::cout<<"counter = "<<counter<<" x = "<<row_temp_x[counter]<<" y = "<<row_temp_y[counter]<<std::endl;
        ++counter;
    }

    // fill matrices
    for(int i = 0; i < 4; i++){
        matrix_dx[iteration][i] = row_temp_x[i];
        matrix_dy[iteration][i] = row_temp_y[i];
    }
    
}

// FUNCTION TO COUNT H MATRIX FOR ELEMENT
// iteration - punkt calkowania, jacobi - jakobian, element - tablica pochodnych dN i dKsi, print - flaga czy drukowac wyniki
double** count_H_matrix(int iteration, double **jacobi, ELEMENT_2D_ARRAYS element, bool if_print){
    double **DX_matrix = new double*[4];
    double **DY_matrix = new double*[4]; // temp matrices inside integral
    double ***H = new double**[4];
    for(int i = 0; i  < 4; i++){
        H[i] = new double*[4];
        for(int j = 0; j < 4; j++)
            H[i][j] = new double[4];
        DX_matrix[i] = new double[4];
        DY_matrix[i] = new double[4];
    }

    for(int i = 0; i < 4; i++)
        add_row(DX_matrix, DY_matrix, i, jacobi, element);

    // count row by transponed row for each pc (punkt calkowania)
    double **X_mat = new double*[4];
    double **Y_mat = new double*[4];
    for(int i = 0; i < 4; i ++){
        X_mat[i] = new double[4];
        Y_mat[i] = new double[4];
    }

    double k = 30; // k(t) factor
    double dV = 1/determinant(jacobi); // dV value
    // final counting H
    // 1 => X_mat + Y_mat

        X_mat = multi_matrix_transponed(4, DX_matrix[iteration]);
        Y_mat = multi_matrix_transponed(4, DY_matrix[iteration]);
        H[iteration] = sum_matrix(X_mat, Y_mat, 4);
        const_multi_matrix(k, H[iteration], 4);
        const_multi_matrix(dV, H[iteration], 4);

        if(if_print){
            std::cout<<"H FOR PC "<<iteration<<std::endl;
            print(H[iteration], 4, 5);
        }
        for(int i = 0; i < 4; i ++){
            delete X_mat[i];
            delete Y_mat[i];
            delete DX_matrix[i];
            delete DY_matrix[i];
        }
        delete[] X_mat;
        delete[] Y_mat;
        delete[] DY_matrix;
        delete[] DX_matrix;
        return H[iteration];
}

// *** COUNTING H MATRIX FOR ETNIRE GRID
// parameters - grid, dN dKsi integrals element, boolean to print results or not 
void H_matrix_grid(grid *GRID, ELEMENT_2D_ARRAYS ELEMENT_1, bool if_print){
    // tablica 1d wezlow w siatce
    node *nodes = new node[GRID->nodes_number];
    element *elements = new element[GRID->elm_number];
    // przepisywanie wezlow do tablicy 1d
    int counter = 0;
    for(int i = 0; i < GRID->nW; i++){
        for(int j = 0; j < GRID->nH; j++){
            nodes[counter] = GRID->nodes[i][j];
                counter++;
        }
    }
    
    // przepisywanie elementow do tablicy 1d
    counter = 0;
    for(int i = 0; i < GRID->nW-1; i++){
        for(int j = 0; j < GRID->nH-1; j++){
            elements[counter] = GRID->elements[i][j];
            counter++;
        }
    }

    // deklaracja macierzy do przechowania Jakobianu
    double **jacobian = new double*[2];
    for(int i = 0; i < 2; i++)
        jacobian[i] = new double[2];
    
    // MAIN ALHORITM
    for(int i = 0; i < GRID->elm_number; i++){
        // dla kazdego elementu stworz temp_element z 4 wezlami tego elementu
        double **H_sum = new double*[4];
        for(int k = 0; k < 4; k++){
            H_sum[k] = new double[4];
        }
        node temp_element[4];
        temp_element[0] = nodes[elements[i].id[0]];
        temp_element[1] = nodes[elements[i].id[1]];
        temp_element[2] = nodes[elements[i].id[2]];
        temp_element[3] = nodes[elements[i].id[3]];
        // inner loop -> 4 punkty calkowania dla kazdego 
        if(if_print)
            std::cout<<"================= COUNGTING H MATRIX FOR ELEMENT "<<i<<" ========================="<<std::endl;
        for(int j = 0; j < 4; j++){
                double **jacobian = new double*[2];
                    for(int k = 0; k < 2; k++){
                        jacobian[k] = new double[2];
                    }
            jacobi(/*i,*/j, jacobian, jacobian, ELEMENT_1, temp_element);
            H_sum = sum_matrix(count_H_matrix(j, jacobian, ELEMENT_1, if_print), H_sum, 4);
            for(int k = 0; k < 2; k++)
                delete jacobian[k];
            delete[] jacobian;
        }
        if(if_print){
            std::cout<<"SUMMED H MATRIX FOR ELEMENT "<<i<<std::endl;
            print(H_sum, 4);}

        elements[i].H = H_sum;
    }

    // save H matrix for each element in grid
    counter = 0;
    for(int i = 0; i < GRID->nW-1; i++){
        for(int j = 0; j < GRID->nH-1; j++){
            GRID->elements[i][j] = elements[counter];
            counter++;
        }
    }
    delete[] elements;
    delete[] nodes;
}

// bound - bound to count boundary codition on
void count_Hbc_matrix(int bound, node elmnt[4]){
    double **Hpc1 = new double*[4];
    double **Hpc2 = new double*[4];
    double *n1 = new double[4];
    double *n2 = new double[4];
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            Hpc1[i] = new double[4];
            Hpc2[i] = new double[4];
        }
    }
}

int main(){
    el_4_2d element = el_4_2d("H");
    el_9_2d element_9 = el_9_2d();

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

    // macierz pochodnych dN/dksi - 9 wezlow
    // double **matrix_ksi_9points = new double*[9];
    // for(int i = 0; i < 9; i++){
    //     matrix_ksi_9points[i] = new double[4];
    // }

    // // macierz pochodnych dN/dn - 9 wezlow
    // double **matrix_n_9points = new double*[9];
    // for(int i = 0; i < 9; i++){
    //     matrix_n_9points[i] = new double[4];
    // }

    // MATRIX 4X4, dN/dksi
    std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dksi ===="<<std::endl;
    element.count_matrix_4points_ksi(matrix_ksi_4points);
    print(matrix_ksi_4points, 4);

    // MATRIX 4X4, dN/dn
    std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dn ===="<<std::endl;
    element.count_matrix_4points_n(matrix_n_4points);
    print(matrix_n_4points, 4);

  //  std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dksi ===="<<std::endl;
  //  element_9.count_matrix_9points_ksi(matrix_ksi_9points);
  //  print_9(matrix_ksi_9points);

   // std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dn ===="<<std::endl;
  //  element_9.count_matrix_9points_n(matrix_n_9points);
  //  print_9(matrix_n_9points);

    ELEMENT_2D_ARRAYS ELEMENT_1 = ELEMENT_2D_ARRAYS(matrix_ksi_4points, matrix_n_4points);

    // creating grid
    grid grid_1 = grid(HEIGTH, WIDTH, NODES_HEIGTH, NODES_WIDTH);
    createNodes(grid_1);
    grid_1.printNodes();
    createElements(grid_1, NODES_HEIGTH);
    grid_1.printElements();

    // count J for grid
   // jacobi_grid(grid_1, ELEMENT_1);

    node temp_element[4];
    node n0 = node(0, 0, 1);
    node n1 = node(0.025, 0, 1);
    node n2 = node(0.025, 0.025, 1);
    node n3 = node(0, 0.025, 1);

    temp_element[0] = n0;
    temp_element[1] = n1;
    temp_element[2] = n2;
    temp_element[3] = n3;

    // // count Jacobi example for 4 iterations (4 punkty calkowania)
    // for(int i = 0; i < 4; i++){
    //     double **jacobian = new double*[2];
    //     double **jacobian_inv = new double*[2];
    //     for(int i = 0; i < 2; i++){
    //            jacobian[i] = new double[2];
    //            jacobian_inv[i] = new double[2];
    //         }
    //     jacobi(i, jacobian, jacobian_inv, ELEMENT_1, temp_element);
    //     std::cout<<"JACOBI REVERSED: "<<std::endl;
    //     const_multi_matrix(1/determinant(jacobian), jacobian, 2);
    //     for(int i = 0; i < 2; i++){
    //         for(int j = 0; j < 2; j++){
    //             std::cout<<"   "<<jacobian[i][j]<<"     ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
        double **jacobian = new double*[2];
        for(int i = 0; i < 2; i++){
               jacobian[i] = new double[2];
            }
        jacobi(0, jacobian, jacobian, ELEMENT_1, temp_element);

        // count matrix for example element - uncomment printing result in count_H_matrix function
        std::cout<<"======================================================= "<<std::endl;
        std::cout<<"Example element: "<<std::endl;
        for(int i = 0; i < 4; i++)
            count_H_matrix(i, jacobian, ELEMENT_1, true);
        
        std::cout<<"======================================================="<<std::endl;
        
        // count H matrices for grid
        H_matrix_grid(&grid_1, ELEMENT_1, true);
        // grid_1.print_H_matrices();

        // counting Hbc
        el_4_2d element_HBC = el_4_2d("HBC");
        // macierz pochodnych dN/dksi - 4 wezly
        for(int i = 0; i < 4; i++){
            delete matrix_ksi_4points[i];
            delete matrix_n_4points[i];
        }
        delete[] matrix_n_4points;
        delete[] matrix_ksi_4points;
        delete[] jacobian;
        
        double **H = new double*[8];
        // const int s = 8;
        // double **matrix_n_4points_HBC = new double*[8];
        // for(int i = 0; i < 8; i++){
        //     std::cout<<"here i = "<<i<<std::endl;
        //     //matrix_ksi_4points_HBC[i] = new double[4];
        //     matrix_n_4points_HBC[i] = new double[0];
        // }

    return 0;
}