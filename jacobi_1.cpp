// IMPLEMENTACJA JAKOBIANU - part 1 
// schemat 2 punktowy calkowania
// PROGRAM OBLICZA JAKOBIAN DLA KAZDEGO ELEMENTU W SIATCE

#include <iostream>
#include <iomanip>
#include <cmath>
// #include "grid.h"

// PARAMETRY SIATKI
//#ifndef WIDTH
#define WIDTH 0.1 // szerokosc
#define HEIGTH 0.2 // wysokosc
#define NODES_HEIGTH 5 // ilosc wezlow wysokosc
#define NODES_WIDTH 4 // ilosc wezlow szerokosc
//#endif

// #ifndef WIDTH
#define WIDTH 0.1 // szerokosc
#define HEIGTH 0.2 // wysokosc
#define NODES_HEIGTH 5 // ilosc wezlow wysokosc
#define NODES_WIDTH 4 // ilosc wezlow szerokosc
// #endif

class node{

    public:
        float x, y;

    node(float X, float Y){
        x = X;
        y = Y;
    }

    node(){
        ;
    }
};

class element{

    public:
        //int id1, id2, id3, id4; // id wezlow
        int id[4];

        element(int i1, int i2, int i3, int i4){
            id[0] = i1;
            id[1] = i2;
            id[2] = i3;
            id[3] = i4;
        }

        element(){
            ;
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
    }

    void printNodes(){
        std::cout<<"NODES NUMBER: "<<nodes_number<<std::endl;
        std::cout<<"================="<<std::endl;
        int counter = 0;
        for(int i = 0; i < nW; i++){
            for(int j = 0; j < nH; j++){
                std::cout<<"Node  "<<std::setprecision(2)<<counter++<<", x = "<<nodes[i][j].x<<", y = "<<nodes[i][j].y<<std::endl;
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
            gr.nodes[i][j] = node(x, y); // dodaj nowy wezel do tablicy
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

        point array[4];
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

        el_4_2d(){
            fac = factors();
            int counter = 0;
            for(int i = 0; i < 2; i++){
                for(int j = 0; j < 2; j++){
                    array[counter++] = point(fac.point_2p[i], fac.point_2p[j]);
                }
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
            std::cout<<"["<<i<<"]["<<j<<"], "<<matrix[i][j]<<"  ;  ";
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

class element_grid{

    public:
        double x[4]; // id wezlow
        double y[4];

        element_grid(double i1, double i2, double i3, double i4, double y1, double y2, double y3, double y4){
            x[0] = i1;
            x[1] = i2;
            x[2] = i3;
            x[3] = i4;
            y[0] = y1;
            y[1] = y2;
            y[2] = y3;
            y[3] = y4;
        }

        element_grid(){
            ;
        }
};

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

// double multi_matrix_2x2(double **m1,){ 
//     double sum = 0.0;
//     for(int i = 0; i < 2; i++){
//         for(int j = 0; j < 2; j++)
//             sum += matrix
//     }
//  }

// COUNT JACOBI MATRIX
// node elmnt[4] - reprezentacja elementu z 4 wezlami w siatce
// iter - iteracja -> iterowany punkt calkowania 
// elm_number -> numer elementu
void jacobi(/*int elm_number,*/ int iter, double** J, double** J_inv, ELEMENT_2D_ARRAYS elmnt_4p, node elmnt[4]){
    // count values in jacobi matrix - punkt calkowania pc[iter]
    std::cout<<" ; pc = "<<iter<<std::endl;

    for(int i = 0; i < 4; i++){
       // std::cout<<elmnt_4p.MATRIX_4P_DKSI[iter][i]<<"*"<<elmnt[i].x<<std::endl;
        J[0][0] +=  elmnt_4p.MATRIX_4P_DKSI[iter][i]*elmnt[i].x;
    }
    std::cout<<std::endl;

    for(int i = 0; i < 4; i++){
       // std::cout<<elmnt_4p.MATRIX_4P_DKSI[iter][i]<<"*"<<elmnt[i].y<<std::endl;
        J[0][1] +=  elmnt_4p.MATRIX_4P_DKSI[iter][i]*elmnt[i].y;
    }
    std::cout<<std::endl;

    for(int i = 0; i < 4; i++){
     //   std::cout<<elmnt_4p.MATRIX_4P_DN[iter][i]<<"*"<<elmnt[i].x<<std::endl;
        J[1][0] +=  elmnt_4p.MATRIX_4P_DN[iter][i]*elmnt[i].x;
    }
    std::cout<<std::endl;

    for(int i = 0; i < 4; i++){
      //  std::cout<<elmnt_4p.MATRIX_4P_DN[iter][i]<<"*"<<elmnt[i].y<<std::endl;
        J[1][1] +=  elmnt_4p.MATRIX_4P_DN[iter][i]*elmnt[i].y; 
    }

    // print results
    std::cout<<std::endl;
    std::cout<<"JACOBI MATRIX: pc = "<<iter<<std::endl;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            std::cout<<std::setprecision(5)<<"["<<j<<"]["<<i<<"], "<<J[i][j]<<"  ;  ";
        }
        std::cout<<std::endl;
    }

    std::cout<<std::setprecision(5)<<"DETERMINANT: "<<determinant(J)<<std::endl;
    std::cout<<std::setprecision(5)<<"1/DETERMINANT: "<<1/determinant(J)<<std::endl;
    std::cout<<std::endl;
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
                std::cout<<"NODES ["<<counter<<"] : x = "<<nodes[counter].x<<" y = "<<nodes[counter].y<<std::endl;
                counter++;
        }
    }
    
    // przepisywanie elementow do tablicy 1d
    counter = 0;
    for(int i = 0; i < GRID.nW-1; i++){
        for(int j = 0; j < GRID.nH-1; j++){
            elements[counter] = GRID.elements[i][j];
                std::cout<<"ELEMENTS ["<<counter<<"] : id1 = "<<elements[counter].id[0]
                <<" id2 = "<<elements[counter].id[1]<<" id3 = "<<elements[counter].id[2]<<
                " id4 = "<<elements[counter].id[3]<<std::endl;
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
            jacobi(/*i,*/j, jacobian, jacobian, ELEMENT_1, temp_element);
        }
    }
}

int main(){
    el_4_2d element = el_4_2d();
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
    double **matrix_ksi_9points = new double*[9];
    for(int i = 0; i < 9; i++){
        matrix_ksi_9points[i] = new double[4];
    }

    // macierz pochodnych dN/dn - 9 wezlow
    double **matrix_n_9points = new double*[9];
    for(int i = 0; i < 9; i++){
        matrix_n_9points[i] = new double[4];
    }

    // MATRIX 4X4, dN/dksi
    std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dksi ===="<<std::endl;
    element.count_matrix_4points_ksi(matrix_ksi_4points);
    print(matrix_ksi_4points, 4);

    // MATRIX 4X4, dN/dn
    std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dn ===="<<std::endl;
    element.count_matrix_4points_n(matrix_n_4points);
    print(matrix_n_4points, 4);

  //  std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dksi ===="<<std::endl;
    element_9.count_matrix_9points_ksi(matrix_ksi_9points);
  //  print_9(matrix_ksi_9points);

   // std::cout<<"==== RESULT FOR MATRIX[4][9]"<<", FUN = dN/dn ===="<<std::endl;
    element_9.count_matrix_9points_n(matrix_n_9points);
  //  print_9(matrix_n_9points);

    ELEMENT_2D_ARRAYS ELEMENT_1 = ELEMENT_2D_ARRAYS(matrix_ksi_4points, matrix_n_4points);

    double **jacobian = new double*[2];
    for(int i = 0; i < 2; i++)
        jacobian[i] = new double[2];

    grid grid_1 = grid(HEIGTH, WIDTH, NODES_HEIGTH, NODES_WIDTH);
    createNodes(grid_1);
    grid_1.printNodes();
    createElements(grid_1, NODES_HEIGTH);
    grid_1.printElements();

    jacobi_grid(grid_1, ELEMENT_1);

    // element_grid temp_element = element_grid(0, 0.025, 0.025, 0,
    //                                     0, 0, 0.025, 0.025);

    //  node temp_element[4];
    //  int counter = 0;
    //  for(int i = 0; i < 2; i++){
        
    //     for(int j = 0; j < 2; j++){
    //         std::cout<<nodes[grid_1.elements[0][0].id[counter]][]<<std::endl;
    //     //temp_element[counter] = grid_1.nodes[grid_1.elements[0][0].id[i]][j];
    //     counter++;
    //     }
    // }
    
    // counter = 0;
    // for(int i = 0; i < 2; i++){
    //     for(int j = 0; j < 2; j++){
    //         std::cout<<"node "<<counter<<" x = "<<temp_element[counter].x<<" y = "<<temp_element[counter].y<<std::endl;
    //         counter++;
    //     }
    // }

    // node temp_element[4];
    // node n0 = node(0,0);
    // node n1 = node(0.025, 0);
    // node n2 = node(0.025, 0.025);
    // node n3 = node(0, 0.025);

    // temp_element[0] = n0;
    // temp_element[1] = n1;
    // temp_element[2] = n2;
    // temp_element[3] = n3;

    // // count Jacobi example for 4 iterations (4 punkty calkowania)
    // for(int i = 0; i < 4; i++)
    //     jacobi(i,0, jacobian, jacobian, ELEMENT_1, temp_element);

    // for(int i = 0; i < 2; i++){
    //     for(int j = 0; j < 2; j++){
    //         std::cout<<std::setprecision(5)<<"["<<j<<"]["<<i<<"], "<<jacobian[i][j]<<"  ;  ";
    //     }
    //     std::cout<<std::endl;
    // }
}