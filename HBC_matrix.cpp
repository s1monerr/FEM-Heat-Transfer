#ifndef HBCfile
#define HBCfile

#include <iostream>
#include <iomanip>
#include <cmath>

/** INPUT DATA
/ * https://home.agh.edu.pl/~pkustra/MES/TC_2d.pdf
*/

#define INIT_TEMP 100
#define SIMULATION_TIME 500
#define TIME_STEP 50
#define TEMPERATURE 1200 // temp otoczenia
#define ALPHA 300 // wspolcz
#define CONDUCTIVITY 25 // przewodnosc
#define C 700 // cieplo wlasciwe
#define DENSITY 7800 // gestosc

#define EPSILON 1e-8

/******
 * PARAMETRY SIATKI - obecnie nieuzywane
*/

#define WIDTH 0.1 // szerokosc
#define HEIGTH 0.1 // wysokosc
#define NODES_HEIGTH 4 // ilosc wezlow wysokosc
#define NODES_WIDTH 4 // ilosc wezlow szerokosc



void print_custom(double**, int, int, int);

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
        double **HBC; // 4 HBC matricees - for each bound
        double **C_mat; // c matrix for element
        double *P_vector; // 4 P vectors for element sumed in  one vector

        element(int i1, int i2, int i3, int i4){
            id[0] = i1;
            id[1] = i2;
            id[2] = i3;
            id[3] = i4;
            H = new double*[4];
            HBC = new double*[4];
            C_mat = new double*[4];
            // std::cout<<"here"<<std::endl;
            for(int i = 0; i < 4; i++){
                H[i] = new double[4];
                HBC[i] = new double[4];
                C_mat[i] = new double[4];
            }
        }

        element(){;}

        void print_H_matrix(){
            print_custom(H, 4, 4, 4);
        }

        void print_HBC_matrix(){
            if(HBC != NULL)
                    print_custom(HBC, 4, 4, 4);
            else
                std::cout<<"No boundary condition in element."<<std::endl;
        }

        void print_C_matrix(){
            print_custom(C_mat, 4, 4, 4);
        }

        void print_P_vectors_matrix(){
            if(P_vector != NULL){
                std::cout<<std::endl;
                for(int i = 0; i < 4; i++){
                        std::cout<<"["<<i<<"] "<<std::setprecision(6)<<P_vector[i]<<"  ;  ";
                    }
                }
                std::cout<<std::endl<<std::endl;
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
    grid(){;}


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

    grid(int number_H, int number_W){
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

    void print_HBC_matrices(){
        std::cout<<"=================== HBC MATRICES FOR GRID ==================="<<std::endl;
            int counter = 0;
            for(int i = 0; i < nW-1; i++){
                for(int j = 0; j < nH-1; j++){
                    std::cout<<"ELEMENT "<<counter<<std::endl;
                    elements[i][j].print_HBC_matrix();
                    counter++;
            }
        }
    }

    void print_P_vectors(){
        std::cout<<"=================== P VECTORS FOR GRID ==================="<<std::endl;
            int counter = 0;
            for(int i = 0; i < nW-1; i++){
                for(int j = 0; j < nH-1; j++){
                    std::cout<<"ELEMENT "<<counter;
                    elements[i][j].print_P_vectors_matrix();
                    counter++;
            }
        }
    }

    void print_C_matrices(){
            std::cout<<"=================== C MATRICES FOR GRID ==================="<<std::endl;
            int counter = 0;
            for(int i = 0; i < nW-1; i++){
                for(int j = 0; j < nH-1; j++){
                    std::cout<<"ELEMENT "<<counter<<std::endl;
                    elements[i][j].print_C_matrix();
                    counter++;
            }
        }
    }

    // read nodes from input
    void readNodes(node **node_arr, int width, int heigth){
        for(int i = 0; i < width; i++){
            for(int j = 0; j < heigth; j++){
                std::cout<<"heheheh"<<std::endl;
                nodes[i][j] = node_arr[i][j];
            }
        }
    }

};


struct InitData{
    double SimulationTime;
    double SimulationStepTime;
    double Conductivity;
    double Alfa;
    double Tot;
    double InitialTemp;
    double Density;
    double SpecificHeat;
    int Nodes_number;
    int Elements_number;
    int Width;
    int Heigth;

    grid grid_input;

    InitData(){
        grid_input = grid();
    }

    void print(){
        std::cout<<"Init data values: \nSimulation Time = "<<std::setprecision(4)<<SimulationTime
        <<"\nStep: "<<SimulationStepTime
        <<"\nConducticity: "<<Conductivity
        <<"\nAlfa: "<<Alfa
        <<"\nToT: "<<Tot
        <<"\nInitialTemp: "<<InitialTemp
        <<"\nDensity: "<<Density
        <<"\nSpecHeat: "<<SpecificHeat
        <<"\nNodes number: "<<Nodes_number
        <<"\nELements number: "<<Elements_number
        <<"\nGrid width: "<<Width
        <<"\nGrid heigth:"<<Heigth
        <<std::endl;
    }
};

// class node{
//     public:
//         float x, y;
//         int weight;

//     node(float X, float Y, int W){
//         x = X;
//         y = Y;
//         weight = W;
//     }
//     node(float X, float Y){
//         x = X;
//         y = Y;
//         weight = -1;
//     }
//     node(){
//         ;
//     }
// };

// class element{
//     public:
//         int id[4]; // node's id
        
//         double **H; // H matrix for element
//         double **HBC; // 4 HBC matricees - for each bound
//         double **C_mat; // c matrix for element
//         double *P_vector; // 4 P vectors for element sumed in  one vector

//         element(int i1, int i2, int i3, int i4){
//             id[0] = i1;
//             id[1] = i2;
//             id[2] = i3;
//             id[3] = i4;
//             H = new double*[4];
//             HBC = new double*[4];
//             C_mat = new double*[4];
//             // std::cout<<"here"<<std::endl;
//             for(int i = 0; i < 4; i++){
//                 H[i] = new double[4];
//                 HBC[i] = new double[4];
//                 C_mat[i] = new double[4];
//             }
//         }

//         element(){;}

//         void print_H_matrix(){
//             print_custom(H, 4, 4, 4);
//         }

//         void print_HBC_matrix(){
//             if(HBC != NULL)
//                     print_custom(HBC, 4, 4, 4);
//             else
//                 std::cout<<"No boundary condition in element."<<std::endl;
//         }

//         void print_C_matrix(){
//             print_custom(C_mat, 4, 4, 4);
//         }

//         void print_P_vectors_matrix(){
//             if(P_vector != NULL){
//                 std::cout<<std::endl;
//                 for(int i = 0; i < 4; i++){
//                         std::cout<<"["<<i<<"] "<<std::setprecision(6)<<P_vector[i]<<"  ;  ";
//                     }
//                 }
//                 std::cout<<std::endl<<std::endl;
//             }
// };


// // SIATKA
// class grid{

//     public:
//     // wysokosc i szerokosc siatki
//         float heigth;
//         float width;
//         int nH, nW; // ilosc wezlow -> szerokosc i wysokosc
//         int nodes_number, elm_number; // liczba wezlow; elementow

//         // tablice wezlow i elementow
//         node ** nodes;
//         element ** elements;

//     // constructor
//     grid(){;}


//     grid(float high, float widt, int number_H,  int number_W){
//         heigth = high;
//         width = widt;
//         nH = number_H;
//         nW = number_W;
//         nodes_number = nH*nW;
//         elm_number = (nH - 1)*(nW - 1);

//         // tablice elementow i wezlow 2d
//         nodes = new node*[nW];
//         for(int i = 0; i < nW; i++){
//             nodes[i] = new node[nH];
//         }

//         elements = new element *[nW - 1];
//         for(int i = 0; i < nW-1; i++){
//             elements[i] = new element[nH - 1];
//         }
//     }

//     grid(int number_H, int number_W){
//         nH = number_H;
//         nW = number_W;

//         nodes_number = nH*nW;
//         elm_number = (nH - 1)*(nW - 1);

//         // tablice elementow i wezlow 2d
//         nodes = new node*[nW];
//         for(int i = 0; i < nW; i++){
//             nodes[i] = new node[nH];
//         }

//         elements = new element *[nW - 1];
//         for(int i = 0; i < nW-1; i++){
//             elements[i] = new element[nH - 1];
//         }
//     }

//     void printElements(){
//         int c = 0;
//         std::cout<<std::endl<<"ELEMENTS NUMBER: "<<elm_number<<std::endl;
//         std::cout<<"================="<<std::endl;
//         for(int i = 0; i < nW - 1; i++){
//             for(int j = 0; j < nH - 1; j++){
//                 std::cout<<"Element "<<c++<<" : ";
//                 std::cout<<"id1 = "<<elements[i][j].id[0]<<", id2 = "<<elements[i][j].id[1]<<", id3 = "
//                 <<elements[i][j].id[2]<<", id4 = "<<elements[i][j].id[3]<<std::endl;
//             }
//         }
//         std::cout<<std::endl;
//     }

//     void printNodes(){
//         std::cout<<"NODES NUMBER: "<<nodes_number<<std::endl;
//         std::cout<<"================="<<std::endl;
//         int counter = 0;
//         for(int i = 0; i < nW; i++){
//             for(int j = 0; j < nH; j++){
//                 std::cout<<"Node  "<<std::setprecision(2)<<counter++<<", x = "<<nodes[i][j].x<<", y = "<<nodes[i][j].y<<"   |    weight = "<<nodes[i][j].weight<<std::endl;
//             }
//         }
//     }

//     void print_H_matrices(){
//             std::cout<<"=================== H MATRICES FOR GRID ==================="<<std::endl;
//             int counter = 0;
//             for(int i = 0; i < nW-1; i++){
//                 for(int j = 0; j < nH-1; j++){
//                     std::cout<<"ELEMENT "<<counter<<std::endl;
//                     elements[i][j].print_H_matrix();
//                     counter++;
//             }
//         }
//     }

//     void print_HBC_matrices(){
//         std::cout<<"=================== HBC MATRICES FOR GRID ==================="<<std::endl;
//             int counter = 0;
//             for(int i = 0; i < nW-1; i++){
//                 for(int j = 0; j < nH-1; j++){
//                     std::cout<<"ELEMENT "<<counter<<std::endl;
//                     elements[i][j].print_HBC_matrix();
//                     counter++;
//             }
//         }
//     }

//     void print_P_vectors(){
//         std::cout<<"=================== P VECTORS FOR GRID ==================="<<std::endl;
//             int counter = 0;
//             for(int i = 0; i < nW-1; i++){
//                 for(int j = 0; j < nH-1; j++){
//                     std::cout<<"ELEMENT "<<counter;
//                     elements[i][j].print_P_vectors_matrix();
//                     counter++;
//             }
//         }
//     }

//     void print_C_matrices(){
//             std::cout<<"=================== C MATRICES FOR GRID ==================="<<std::endl;
//             int counter = 0;
//             for(int i = 0; i < nW-1; i++){
//                 for(int j = 0; j < nH-1; j++){
//                     std::cout<<"ELEMENT "<<counter<<std::endl;
//                     elements[i][j].print_C_matrix();
//                     counter++;
//             }
//         }
//     }

//     // read nodes from input
//     void readNodes(node **node_arr, int width, int heigth){
//         for(int i = 0; i < width; i++){
//             for(int j = 0; j < heigth; j++){
//                 std::cout<<"heheheh"<<std::endl;
//                 nodes[i][j] = node_arr[i][j];
//             }
//         }
//     }

// };


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

class point{
    public:
        double n, ksi;

        point(){
            n = 0;
            ksi = 0;
        }
        point(double xx, double yy){
            ksi = xx;
            n = yy;
        }
};

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
                // std::cout<<"point 1 ksi = "<<array[0].ksi<<" n = "<<array[0].n<<std::endl;
                array[1] = point(-1, fac.point_2p[1]);
                // std::cout<<"point 2 ksi = "<<array[1].ksi<<" n = "<<array[1].n<<std::endl;
                array[2] = point(fac.point_2p[0], 1);
                // std::cout<<"point 3 ksi = "<<array[2].ksi<<" n = "<<array[2].n<<std::endl;
                array[3] = point(fac.point_2p[1], 1);
                // std::cout<<"point 4 ksi = "<<array[3].ksi<<" n = "<<array[3].n<<std::endl;
                array[4] = point(1, fac.point_2p[1]);
                // std::cout<<"point 5 ksi = "<<array[4].ksi<<" n = "<<array[4].n<<std::endl;
                array[5] = point(1, fac.point_2p[0]);
                // std::cout<<"point 6 ksi = "<<array[5].ksi<<" n = "<<array[5].n<<std::endl;
                array[6] = point(fac.point_2p[1], -1);
                // std::cout<<"point 7 ksi = "<<array[6].ksi<<" n = "<<array[6].n<<std::endl;
                array[7] = point(fac.point_2p[0], -1);
                // std::cout<<"point 8 ksi = "<<array[7].ksi<<" n = "<<array[7].n<<std::endl;
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
                    std::cout<<"here 2 "<<std::endl;
                    matrix[i][j] = integral_ksi(array[i], j);
                }
            }
        }

        void count_matrix_4points_n_HBC(double **matrix){
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 4; j++){
                   std::cout<<"here 2 "<<std::endl;
                    matrix[i][j] = integral_n(array[i], j);
                }
            }
        }

};

class el_9_2d{

factors fac;

    public:

        point array[9];
        double weights[9];
        // dN/dKsi - funkcja wybierana przez 2 argument
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
                    array[counter] = point(fac.point_3p[j], fac.point_3p[i]);
                    // std::cout<<"point: ksi = "<<array[counter].ksi<<" n = "<<array[counter].n<<std::endl;
                    weights[counter] = fac.w_3p[i]*fac.w_3p[j];
                    // std::cout<<"w[i] = "<<fac.w_3p[i]<<" w[j] = "<<fac.w_3p[j]<<"   ";
                    // std::cout<<"w = "<<weights[counter]<<std::endl<<std::endl;
                    counter++;
                }
            }
            // array[0] = point(fac.point_3p[0], fac.point_3p[0]);
            // array[1] = point(fac.point_3p[1], fac.point_3p[0]);
            // array[2] = point(fac.point_3p[2], fac.point_3p[0]);
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

// final arrays from el_4_2d or 3l_9_2d
class ELEMENT_2D_ARRAYS{
    public:
        double **MATRIX_4P_DKSI;
        double **MATRIX_4P_DN;

        ELEMENT_2D_ARRAYS(double **arr1, double **arr2){
            MATRIX_4P_DKSI = arr1;
            MATRIX_4P_DN = arr2;
        }
};

void jacobi(const int POINTS, int iter, double** J, ELEMENT_2D_ARRAYS integrals, node elmnt[4]);
double determinant(double **matrix);

class el_4_2d_HBC{

factors fac;

    public:

        point *array;
        // funkcja ksztaltu < N1 - N4 > 
        double shape_fun(point p, int index){
            if(index == 0) return 0.25*(1-p.ksi)*(1-p.n);
            
            else if (index == 1) return 0.25*(1-p.n)*(1+p.ksi);
        
            else if(index == 2) return 0.25*(1+p.n)*(1+p.ksi);
        
            else if (index == 3) return 0.25*(1+p.n)*(1-p.ksi);

            return -1;
        }

        el_4_2d_HBC(){ // rodzaj macierzy pochodnych - HBC, jesli dla obliczen macierzy obciazen
            fac = factors();
            // uzupelnij tablice punktami calkowania
                array = new point[8];
                array[0] = point(fac.point_2p[0], -1);
                // std::cout<<"point 1 ksi = "<<array[0].ksi<<" n = "<<array[0].n<<std::endl;
                array[1] = point(fac.point_2p[1], -1);
                // std::cout<<"point 2 ksi = "<<array[1].ksi<<" n = "<<array[1].n<<std::endl;
                array[2] = point(1, fac.point_2p[0]);
                // std::cout<<"point 3 ksi = "<<array[2].ksi<<" n = "<<array[2].n<<std::endl;
                array[3] = point(1, fac.point_2p[1]);
                // std::cout<<"point 4 ksi = "<<array[3].ksi<<" n = "<<array[3].n<<std::endl;
                array[4] = point(fac.point_2p[1], 1);
                // std::cout<<"point 5 ksi = "<<array[4].ksi<<" n = "<<array[4].n<<std::endl;
                array[5] = point(fac.point_2p[0], 1);
                // std::cout<<"point 6 ksi = "<<array[5].ksi<<" n = "<<array[5].n<<std::endl;
                array[6] = point(-1, fac.point_2p[1]);
                // std::cout<<"point 7 ksi = "<<array[6].ksi<<" n = "<<array[6].n<<std::endl;
                array[7] = point(-1, fac.point_2p[0]);
                // std::cout<<"point 8 ksi = "<<array[7].ksi<<" n = "<<array[7].n<<std::endl;
        }

        void count_matrix_4points(double **matrix){
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 4; j++){
                    //std::cout<<"i = "<<i<<" j = "<<j<<" array[i].ksi = "<<array[i].ksi<<" array[i].n = "<<array[i].n<<std::endl;
                    matrix[i][j] = shape_fun(array[i], j);
                }
            }
        }
};


class el_4_2d_C{

factors fac;

    public:

        point *array;
        // funkcja ksztaltu < N1 - N4 > 
        double shape_fun(point p, int index){
            if(index == 0) {return 0.25*(1-p.ksi)*(1-p.n);}
            
            else if (index == 1) return 0.25*(1-p.n)*(1+p.ksi);
        
            else if(index == 2) return 0.25*(1+p.n)*(1+p.ksi);
        
            else if (index == 3) return 0.25*(1+p.n)*(1-p.ksi);

            return -1;
        }

        el_4_2d_C(){ // rodzaj macierzy pochodnych - HBC, jesli dla obliczen macierzy obciazen
            fac = factors();
            // uzupelnij tablice punktami calkowania
                array = new point[4];
                array[0] = point(fac.point_2p[0], fac.point_2p[0]);
                // std::cout<<"point 1 ksi = "<<array[0].ksi<<" n = "<<array[0].n<<std::endl;
                array[1] = point(fac.point_2p[1], fac.point_2p[0]);
                // std::cout<<"point 2 ksi = "<<array[1].ksi<<" n = "<<array[1].n<<std::endl;
                array[2] = point(fac.point_2p[1], fac.point_2p[1]);
                // std::cout<<"point 3 ksi = "<<array[2].ksi<<" n = "<<array[2].n<<std::endl;
                array[3] = point(fac.point_2p[0], fac.point_2p[1]);
                // std::cout<<"point 4 ksi = "<<array[3].ksi<<" n = "<<array[3].n<<std::endl;
        }

        void count_matrix_4points(double **matrix){
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    // std::cout<<"i = "<<i<<" j = "<<j<<" array[i].ksi = "<<array[i].ksi<<" array[i].n = "<<array[i].n<<std::endl;
                    matrix[i][j] = shape_fun(array[i], j);
                    // std::cout<<"counted"<<std::endl;
                }
            }
        }
};

void sum_vectors(double *vec_result, double *vec_added, int size){
    for(int i = 0; i < size; i++)
        vec_result[i] += vec_added[i];
}

void const_multi_matrix(double const_value, double **matrix, int size){    // function to multiply matrix by const
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            matrix[i][j] *= const_value;
        }
    }
}

double **sum_matrix(double **matrix_1, double **matrix_2, int size){
    double **result = new double*[size];
    for(int i = 0; i < size; i++){
        result[i] = new double[size];
    }

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            result[i][j] = matrix_1[i][j] + matrix_2[i][j];
        }
    }
    return result;
}


// multiplication of matrice by transponed matrice
double **multi_vectors(int size, double *matrix) {
    double **result = new double*[size];
    for(int i = 0; i < size; i++){
        result[i] = new double[size];
    }
 
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // std::cout<<"counting: i = "<<i<<" j = "<<j<<" matrix[i] = "<<matrix[i]
            // <<" matrix[j]"<<matrix[j]<<std::endl;
            result[i][j] = matrix[i]*matrix[j];
        }
    }
    return result;
}

double* multi_matrix_vector(double **matrix, double *vector, int size){
    double *result = new double[size];

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            result[i] += matrix[i][j]*vector[j];
        }
    }
    return result;
}

void fix_values(double **array, int size){
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if(fabs(array[i][j]) < EPSILON)
                array[i][j] = 0;
        }
    }
}

void fix_values(double *array, int size){
    for(int i = 0; i < size; i++){
            if(fabs(array[i]) < EPSILON)
                array[i] = 0;
    }
}

void print_custom(double **matrix, int rows, int columns, int precision){
    std::cout<<std::endl;
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            // std::cout<<"["<<i<<"]["<<j<<"] "<<std::setprecision(precision)<<matrix[i][j]<<"  ;  ";
            std::cout<<matrix[i][j]<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

// grid - grid to fill element's HBC, id - element id, nw nh - number of nodes in element's array
// nodes - element nodes, values - shape function values, cond - conductivity in fourier formula
double** count_bound_HBC(node* nodes, double **values, double cond, bool if_print){
    double ** result_HBC = new double*[4]; // final HBC for element
    for(int i = 0; i < 4; i++){
        result_HBC[i] = new double[4];
        for(int j = 0; j < 4; j++){
            result_HBC[i][j] = 0;
        }
    }

    bool HBC_flag = false; // check if there is any bound condition in element

    // MAIN ALGHORITM - 4 bounds
    for(int i = 0; i < 4; i++){
        // check bound values
        if(i != 3){
            if(nodes[i].weight == 1&& nodes[i + 1].weight == 1){
                HBC_flag = true;
                double J = 0;
                if(fabs(nodes[i].x - nodes[i+1].x)<EPSILON)
                    J = fabs((nodes[i].y - nodes[i+1].y)/2);
                else if(fabs(nodes[i].y - nodes[i+1].y)<EPSILON)
                    J = fabs((nodes[i].x - nodes[i+1].x)/2);
                else { // no square grid case
                    double x = pow(nodes[i].x - nodes[i+1].x, 2);
                    double y = pow(nodes[i].y - nodes[i+1].y, 2);
                    J = sqrt(x+y)/2;
                }
                // std::cout<<"i = "<<i<<" node 1: x = "<<nodes[i].x<<" y = "<<nodes[i].y<<
                // " node 2 : x = "<<nodes[i+1].x<<" y = "<<nodes[i+1].y<<std::endl;
                // std::cout<<"J = "<<J<<std::endl;
                // temp structures
                double **bound_p1 = new double*[4];
                double **bound_p2 = new double*[4];
                double **bound_HBC = new double*[4];
                for(int j = 0; j < 4; j++){
                    bound_HBC[j] = new double[4];
                    bound_p1[j] = new double[4];
                    bound_p2[j] = new double[4];
                    for(int k = 0; k < 4; k++){ 
                            bound_HBC[j][k] = 0;
                            bound_p1[j][k] = 0;
                            bound_p2[j][k] = 0;
                        }
                }
                bound_p1 = multi_vectors(4, values[2*i]);
                bound_p2 = multi_vectors(4, values[(2*i)+1]);
                bound_HBC = sum_matrix(bound_p1, bound_p2, 4);
                const_multi_matrix(cond, bound_HBC, 4);
                const_multi_matrix(J, bound_HBC, 4);
                result_HBC = sum_matrix(bound_HBC, result_HBC, 4);
                if(if_print){
                    std::cout<<"HBC for bound "<<i<<std::endl;
                    print_custom(bound_HBC, 4, 4, 5);
                }
                for(int j = 0; j < 4; j++){
                    delete bound_HBC[j];
                    delete bound_p1[j];
                    delete bound_p2[j];
                }
                delete[] bound_HBC;
                delete[] bound_p1;
                delete[] bound_p2;
            }
            else if(if_print)
                std::cout<<"No boundary condition on bound "<<i<<std::endl;
        }
        else if(i == 3){
             if(nodes[i].weight == 1&& nodes[0].weight == 1){
                HBC_flag = true;
                double J = 0;
                if(fabs(nodes[i].x - nodes[0].x)<EPSILON)
                    J = fabs((nodes[i].y - nodes[0].y)/2);
                else if(fabs(nodes[i].y - nodes[0].y)<EPSILON)
                    J = fabs((nodes[i].x - nodes[0].x)/2);
                else { // no square grid case
                    double x = pow(nodes[i].x - nodes[0].x, 2);
                    double y = pow(nodes[i].y - nodes[0].y, 2);
                    J = sqrt(x+y)/2;
                }
                // std::cout<<"J = "<<J<<std::endl;
                // std::cout<<"i = "<<i<<" node 1: x = "<<nodes[2*i].x<<" y = "<<nodes[i].y<<
                // " node 2 : x = "<<nodes[0].x<<" y = "<<nodes[0].y<<std::endl;
                // temp structures
                double **bound_p1 = new double*[4];
                double **bound_p2 = new double*[4];
                double **bound_HBC = new double*[4];
                for(int j = 0; j < 4; j++){
                    bound_HBC[j] = new double[4];
                    bound_p1[j] = new double[4];
                    bound_p2[j] = new double[4];
                    for(int k = 0; k < 4; k++){ 
                            bound_HBC[j][k] = 0;
                            bound_p1[j][k] = 0;
                            bound_p2[j][k] = 0;
                        }
                }
                bound_p1 = multi_vectors(4, values[2*i]);
                bound_p2 = multi_vectors(4, values[2*i+1]);
                bound_HBC = sum_matrix(bound_p1, bound_p2, 4);
                const_multi_matrix(cond, bound_HBC, 4);
                const_multi_matrix(J, bound_HBC, 4);
                result_HBC = sum_matrix(bound_HBC, result_HBC, 4);
                if(if_print){
                    std::cout<<"HBC for bound "<<i<<std::endl;
                    print_custom(bound_HBC, 4, 4, 5);
                }
                for(int j = 0; j < 4; j++){
                    delete bound_HBC[j];
                    delete bound_p1[j];
                    delete bound_p2[j];
                }
                delete[] bound_HBC;
                delete[] bound_p1;
                delete[] bound_p2;
            }
            else if(if_print)
                std::cout<<"No boundary condition on bound "<<i<<std::endl<<std::endl;
        }
    }
    if(HBC_flag&&if_print){
        std::cout<<"HBC: "<<std::endl;
        print_custom(result_HBC, 4, 4, 4);
    }
    else if(!HBC_flag&&if_print){
        std::cout<<"No boundary conditions in element!"<<std::endl;
    }

    if(HBC_flag)
        return result_HBC;
    else return NULL;
}

void count_HBC_grid(grid *GRID, bool if_print, double **values, double cond){
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

    for(int i = 0; i < GRID->elm_number; i++){
        // dla kazdego elementu stworz temp_element z 4 wezlami tego elementu
        node temp_element[4];
        temp_element[0] = nodes[elements[i].id[0]];
        temp_element[1] = nodes[elements[i].id[1]];
        temp_element[2] = nodes[elements[i].id[2]];
        temp_element[3] = nodes[elements[i].id[3]];
        // inner loop -> 4 punkty calkowania dla kazdego 

        //      std::cout<<"element = "<<i<<" id[0] = "<<elements[i].id[0]<<" id[1] = "<<elements[i].id[1]<<" id[0] = "<<elements[i].id[2]<<" id[3] = "<<elements[i].id[3]<<std::endl;
        // std::cout<<"element = "<<i<<" weights = "<<temp_element[0].weight<<"  "<<temp_element[1].weight<<"  "<<temp_element[2].weight<<"  "<<temp_element[3].weight<<std::endl;
        
        if(if_print)
            std::cout<<"================= COUNGTING HBC MATRICES FOR ELEMENT "<<i<<" ========================="<<std::endl;
        elements[i].HBC = count_bound_HBC(temp_element, values, cond, if_print);
        if(elements[i].HBC != NULL)
            fix_values(elements[i].HBC, 4);
    }

    // save HBC matrix for each element in grid
    counter = 0;
    for(int i = 0; i < GRID->nW-1; i++){
        for(int j = 0; j < GRID->nH-1; j++){
            if(elements[counter].HBC != NULL)
                GRID->elements[i][j].HBC = elements[counter].HBC;
            counter++;
        }
    }
}


// HBC - schemat 3 punktowy
double** count_bound_HBC_3(node* nodes, double **values, double cond, bool if_print){
    double ** result_HBC = new double*[4]; // final HBC for element
    for(int i = 0; i < 4; i++){
        result_HBC[i] = new double[4];
        for(int j = 0; j < 4; j++){
            result_HBC[i][j] = 0;
        }
    }

    bool HBC_flag = false; // check if there is any bound condition in element

    // MAIN ALGHORITM - 4 bounds
    for(int i = 0; i < 4; i++){
        // check bound values
        if(i != 3){
            if(nodes[i].weight == 1&& nodes[i + 1].weight == 1){
                HBC_flag = true;
                double J = 0;
                if(nodes[i].x == nodes[i+1].x)
                    J = fabs((nodes[i].y - nodes[i+1].y)/2);
                else
                    J = fabs((nodes[i].x - nodes[i+1].x)/2);
                // std::cout<<"i = "<<i<<" node 1: x = "<<nodes[i].x<<" y = "<<nodes[i].y<<
                // " node 2 : x = "<<nodes[i+1].x<<" y = "<<nodes[i+1].y<<std::endl;
                // std::cout<<"J = "<<J<<std::endl;
                // temp structures
                double **bound_p1 = new double*[4];
                double **bound_p2 = new double*[4];
                double **bound_HBC = new double*[4];
                for(int j = 0; j < 4; j++){
                    bound_HBC[j] = new double[4];
                    bound_p1[j] = new double[4];
                    bound_p2[j] = new double[4];
                    for(int k = 0; k < 4; k++){ 
                            bound_HBC[j][k] = 0;
                            bound_p1[j][k] = 0;
                            bound_p2[j][k] = 0;
                        }
                }
                bound_p1 = multi_vectors(4, values[2*i]);
                bound_p2 = multi_vectors(4, values[(2*i)+1]);
                bound_HBC = sum_matrix(bound_p1, bound_p2, 4);
                const_multi_matrix(cond, bound_HBC, 4);
                const_multi_matrix(J, bound_HBC, 4);
                result_HBC = sum_matrix(bound_HBC, result_HBC, 4);
                if(if_print){
                    std::cout<<"HBC for bound "<<i<<std::endl;
                    print_custom(bound_HBC, 4, 4, 5);
                }
                for(int j = 0; j < 4; j++){
                    delete bound_HBC[j];
                    delete bound_p1[j];
                    delete bound_p2[j];
                }
                delete[] bound_HBC;
                delete[] bound_p1;
                delete[] bound_p2;
            }
            else if(if_print)
                std::cout<<"No boundary condition on bound "<<i<<std::endl;
        }
        else if(i == 3){
             if(nodes[i].weight == 1&& nodes[0].weight == 1){
                HBC_flag = true;
                double J = 0;
                if(nodes[i].x == nodes[0].x)
                    J = fabs((nodes[i].y - nodes[0].y)/2);
                else
                    J = fabs((nodes[i].x - nodes[0].x)/2);
                // std::cout<<"J = "<<J<<std::endl;
                // std::cout<<"i = "<<i<<" node 1: x = "<<nodes[2*i].x<<" y = "<<nodes[i].y<<
                // " node 2 : x = "<<nodes[0].x<<" y = "<<nodes[0].y<<std::endl;
                // temp structures
                double **bound_p1 = new double*[4];
                double **bound_p2 = new double*[4];
                double **bound_HBC = new double*[4];
                for(int j = 0; j < 4; j++){
                    bound_HBC[j] = new double[4];
                    bound_p1[j] = new double[4];
                    bound_p2[j] = new double[4];
                    for(int k = 0; k < 4; k++){ 
                            bound_HBC[j][k] = 0;
                            bound_p1[j][k] = 0;
                            bound_p2[j][k] = 0;
                        }
                }
                bound_p1 = multi_vectors(4, values[2*i]);
                bound_p2 = multi_vectors(4, values[2*i+1]);
                bound_HBC = sum_matrix(bound_p1, bound_p2, 4);
                const_multi_matrix(cond, bound_HBC, 4);
                const_multi_matrix(J, bound_HBC, 4);
                result_HBC = sum_matrix(bound_HBC, result_HBC, 4);
                if(if_print){
                    std::cout<<"HBC for bound "<<i<<std::endl;
                    print_custom(bound_HBC, 4, 4, 5);
                }
                for(int j = 0; j < 4; j++){
                    delete bound_HBC[j];
                    delete bound_p1[j];
                    delete bound_p2[j];
                }
                delete[] bound_HBC;
                delete[] bound_p1;
                delete[] bound_p2;
            }
            else if(if_print)
                std::cout<<"No boundary condition on bound "<<i<<std::endl<<std::endl;
        }
    }
    if(HBC_flag&&if_print){
        std::cout<<"HBC: "<<std::endl;
        print_custom(result_HBC, 4, 4, 4);
    }
    else if(!HBC_flag&&if_print){
        std::cout<<"No boundary conditions in element!"<<std::endl;
    }

    if(HBC_flag)
        return result_HBC;
    else return NULL;
}

void count_HBC_grid_3(grid *GRID, bool if_print, double **values, double cond){
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

    for(int i = 0; i < GRID->elm_number; i++){
        // dla kazdego elementu stworz temp_element z 4 wezlami tego elementu
        node temp_element[4];
        temp_element[0] = nodes[elements[i].id[0]];
        temp_element[1] = nodes[elements[i].id[1]];
        temp_element[2] = nodes[elements[i].id[2]];
        temp_element[3] = nodes[elements[i].id[3]];
        // inner loop -> 4 punkty calkowania dla kazdego 

        //      std::cout<<"element = "<<i<<" id[0] = "<<elements[i].id[0]<<" id[1] = "<<elements[i].id[1]<<" id[0] = "<<elements[i].id[2]<<" id[3] = "<<elements[i].id[3]<<std::endl;
        // std::cout<<"element = "<<i<<" weights = "<<temp_element[0].weight<<"  "<<temp_element[1].weight<<"  "<<temp_element[2].weight<<"  "<<temp_element[3].weight<<std::endl;
        
        if(if_print)
            std::cout<<"================= COUNGTING HBC MATRICES FOR ELEMENT "<<i<<" ========================="<<std::endl;
        elements[i].HBC = count_bound_HBC(temp_element, values, cond, if_print);
        if(elements[i].HBC != NULL)
            fix_values(elements[i].HBC, 4);
    }

    // save HBC matrix for each element in grid
    counter = 0;
    for(int i = 0; i < GRID->nW-1; i++){
        for(int j = 0; j < GRID->nH-1; j++){
            if(elements[counter].HBC != NULL)
                GRID->elements[i][j].HBC = elements[counter].HBC;
            counter++;
        }
    }
}

// P VECTOR COUNTING
void print_vector(double *vec, int size, int precision){
    for(int i = 0; i < size; i++)
        std::cout<<"["<<std::setprecision(precision)<<vec[i]<<"] ";
    std::cout<<std::endl;
}

void multi_vector_const(double *vector, double constant, int size){
    for(int i = 0; i < size; i++){
        vector[i] *= constant;
    }
}

// TEMP - temperature
double *multi_vector_shape_values(double *vector, double **shape_values, int index1, int index2, int size, double TEMP){
    for(int i = 0; i < size; i++) vector[i] = 0;
        
    for(int i = 0; i < size; i++){
        vector[i] += shape_values[index1][i]*TEMP;
        vector[i] += shape_values[index2][i]*TEMP;
    }
    return vector;
}

// values - HBC matrix values in element
double *count_P_vectors(node nodes[4], double **values, double temperature, double cond, bool if_print){
    double *result = new double[4];
    for(int i = 0; i < 4; i++){
        result[i] = 0.0;
    }

    // count for every bound
    bool bound_flag = false; // flag if boundary condition exists in element
    for(int i = 0; i < 4; i++){
        bool loop_bound_flag = false; // check if bound vector will be null - save or not save in 4x4 result array
        double *bound_vector = new double[4];  
        if(i != 3){
            if((nodes[i].weight == 1) && (nodes[i+1].weight == 1)){
                bound_flag = true;
                loop_bound_flag = true;
                double J = 0;

                if(nodes[i].x == nodes[i+1].x)
                    J = fabs((nodes[i].y - nodes[i+1].y)/2);
                else
                    J = fabs((nodes[i].x - nodes[i+1].x)/2);
   
                bound_vector = multi_vector_shape_values(bound_vector, values, 2*i, 2*i+1, 4, temperature);
                multi_vector_const(bound_vector, cond, 4);
                multi_vector_const(bound_vector, J, 4);
                if(if_print)
                    print_vector(bound_vector, 4, 4);
            }
            else if(if_print){
                std::cout<<"No boundary condition on bound "<<i<<"!"<<std::endl;
                // bound_vector = NULL;
            }
            // else{
            //     bound_vector = NULL;
            // }
        }
        else if(i == 3){
            if(nodes[i].weight == 1 && nodes[0].weight == 1){
                bound_flag = true;
                loop_bound_flag = true;

                double J = 0;

                if(nodes[i].x == nodes[0].x)
                    J = fabs((nodes[i].y - nodes[0].y)/2);
                else
                    J = fabs((nodes[i].x - nodes[0].x)/2);
  
                bound_vector = multi_vector_shape_values(bound_vector, values, 2*i, 2*i+1, 4, temperature);
                multi_vector_const(bound_vector, cond, 4);
                multi_vector_const(bound_vector, J, 4);
                if(if_print)
                    print_vector(bound_vector, 4, 4);
            }
            else if(if_print){
                std::cout<<"No boundary condition on bound "<<i<<"!"<<std::endl;
                // bound_vector = NULL;
            }
            // else{
            //     bound_vector = NULL;
            // }
        }
        // add boundary condition to result matrix
        if(loop_bound_flag)
        for(int k = 0; k < 4; k++){
            result[k] += bound_vector[k];
        }
        // else
        //     for(int k = 0; k < 4; k++){
        //     result[k] += 0;
        // }
        delete[] bound_vector;
    }
   
   fix_values(result, 4);
    if(bound_flag)
        return result;
    else{
        double *zero = new double[4]{0, 0, 0, 0};
        fix_values(zero, 4);
        return zero;
    }
}

void count_P_vectors_grid(grid *GRID, bool if_print, double **values, double cond, double temperature){
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

    for(int i = 0; i < GRID->elm_number; i++){
        // dla kazdego elementu stworz temp_element z 4 wezlami tego elementu
        node temp_element[4];
        temp_element[0] = nodes[elements[i].id[0]];
        temp_element[1] = nodes[elements[i].id[1]];
        temp_element[2] = nodes[elements[i].id[2]];
        temp_element[3] = nodes[elements[i].id[3]];
        // std::cout<<"element = "<<i<<" id[0] = "<<elements[i].id[0]<<" id[1] = "<<elements[i].id[1]<<" id[0] = "<<elements[i].id[2]<<" id[3] = "<<elements[i].id[3]<<std::endl;
        // std::cout<<"element = "<<i<<" weights = "<<temp_element[0].weight<<"  "<<temp_element[1].weight<<"  "<<temp_element[2].weight<<"  "<<temp_element[3].weight<<std::endl;
        // inner loop -> 4 punkty calkowania dla kazdego 
        if(if_print)
            std::cout<<"================= COUNGTING P VECTORS FOR ELEMENT "<<i<<" ========================="<<std::endl;
        elements[i].P_vector = count_P_vectors(temp_element, values, temperature, cond, if_print);
    }

    // save P vectors matrix for each element in grid
    counter = 0;
    for(int i = 0; i < GRID->nW-1; i++){
        for(int j = 0; j < GRID->nH-1; j++){
            // if(elements[counter].HBC != NULL)
                GRID->elements[i][j].P_vector = elements[counter].P_vector;
            counter++;
        }
    }

}

// p vectors aggregation
void count_P_vectors_2d_aggregation(grid GRID, double *result){
    // ZEROWANIE TABLICY 
    for(int i = 0; i > GRID.elm_number; i++)
        result[i] = 0;

    element *elements = new element[GRID.elm_number];
    
    // przepisywanie elementow do tablicy 1d
    int counter = 0;
    for(int i = 0; i < GRID.nW-1; i++){
        for(int j = 0; j < GRID.nH-1; j++){
            elements[counter] = GRID.elements[i][j];
            counter++;
        }
    }

    for(int i = 0; i < GRID.elm_number; i++){ // main loop = elements loop
        // FOR EACH ELEMENT ALGHORITM
        // std::cout<<"MAIN ITERATOR : i = "<<i<<std::endl;
        if(elements[i].P_vector != NULL){
            for(int j = 0; j < 4; j++){
                    result[elements[i].id[j]] += elements[i].P_vector[j];
                }
            }
        }
}
    
// C MATRIX

// multiplication of matrice by transponed matrice
double **multi_matrix_transponed(int size, double *matrix) {
    double **result = new double*[size];
    for(int i = 0; i < size; i++){
        result[i] = new double[size];
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = matrix[i]*matrix[j];
        }
    }
    return result;
}

// nodes - temp element, C_shape_values - shape functions values
// c - specific heat, density - element density. if_print - printing results flag
// JACOBIAN - @TODO J - jacobian
double **count_C_matrix(node *nodes, double **C_shape_values, double c, double density, ELEMENT_2D_ARRAYS ELEMENT, bool if_print){
    double **result = new double*[4];
    for(int i = 0; i < 4; i++){
        result[i] = new double[4];
    }

    // MAIN LOOP - integration points
    for(int i = 0; i < 4; i++){
        double *vec_bound = new double[4];
        for(int j = 0; j < 4; j++){
            vec_bound[j] = 0;
            vec_bound[j] = C_shape_values[i][j];
        }

            double **C_mat_bound = new double*[4];
            for(int m = 0; m < 4; m ++){
                C_mat_bound[m] = new double[4];
            }

            C_mat_bound = multi_matrix_transponed(4, vec_bound);
            const_multi_matrix(c, C_mat_bound, 4);
            const_multi_matrix(density, C_mat_bound, 4);

            // JACOBIAN TO BE COUNTED
            double **jac = new double*[2];
            for(int i = 0; i < 2; i++) jac[i] = new double[2];
            jacobi(2, i, jac, ELEMENT, nodes);
            double J = determinant(jac);
            const_multi_matrix(J, C_mat_bound, 4);

            if(if_print){
                std::cout<<"=== C MATRIX FOR BOUND "<<i<<" ==="<<std::endl;
                print_custom(C_mat_bound, 4, 4, 4);
            }
            
            // save results in result matrix
            for(int k = 0; k < 4; k++){
                for(int l = 0; l < 4; l++){
                    result[k][l] += C_mat_bound[k][l];
                }
            }

            for(int k = 0; k < 4; k++){
                delete C_mat_bound[k];
            }
            delete[] C_mat_bound;
            delete[] vec_bound;
        }
    return result;
}

void count_C_matrix_grid(grid *GRID, double **C_shape_values, double c, double density, ELEMENT_2D_ARRAYS ELEMENT, bool if_print){
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
        double **C_sum = new double*[4];
        for(int k = 0; k < 4; k++){
            C_sum[k] = new double[4];
        }
        node temp_element[4];
        temp_element[0] = nodes[elements[i].id[0]];
        temp_element[1] = nodes[elements[i].id[1]];
        temp_element[2] = nodes[elements[i].id[2]];
        temp_element[3] = nodes[elements[i].id[3]];
        // inner loop -> 4 punkty calkowania dla kazdego 
        if(if_print)
            std::cout<<"================= COUNGTING C MATRIX FOR ELEMENT "<<i<<" ========================="<<std::endl;
        

        C_sum = sum_matrix(count_C_matrix(temp_element, C_shape_values, c, density, ELEMENT, if_print), C_sum, 4);
            
        if(if_print){
            std::cout<<"SUMMED C MATRIX FOR ELEMENT "<<i<<std::endl;
            print_custom(C_sum, 4, 4, 4);}

        elements[i].C_mat = C_sum;
    }

    // save C matrix for each element in grid
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

double determinant(double **matrix) { return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]; }


// COUNT JACOBI MATRIX
// points - liczba punktow calkowania
// node elmnt[4] - reprezentacja elementu z 4 wezlami w siatce
// iter - iteracja -> iterowany punkt calkowania 
// elm_number -> numer elementu
void jacobi(const int POINTS, int iter, double** J, ELEMENT_2D_ARRAYS integrals, node elmnt[4]){
    // count values in jacobi matrix - punkt calkowania pc[iter]
   // std::cout<<"GAUSS: pc = "<<iter<<std::endl;

    el_9_2d temp = el_9_2d();
    if(POINTS == 2 || POINTS == 3){
        for(int i = 0; i < 4; i++){
            //std::cout<<integrals.MATRIX_4P_DKSI[iter][i]<<"*"<<elmnt[i].x<<std::endl;
            if(POINTS == 2)
                J[0][0] +=  integrals.MATRIX_4P_DKSI[iter][i]*elmnt[i].x;
            else   
                J[0][0] +=  integrals.MATRIX_4P_DKSI[iter][i]*elmnt[i].x*temp.weights[iter];
        }
    // std::cout<<std::endl;

        for(int i = 0; i < 4; i++){
        //   std::cout<<integrals.MATRIX_4P_DKSI[iter][i]<<"*"<<elmnt[i].y<<std::endl;
            if(POINTS == 2)
                J[0][1] +=  integrals.MATRIX_4P_DKSI[iter][i]*elmnt[i].y;
            else
                J[0][1] +=  integrals.MATRIX_4P_DKSI[iter][i]*elmnt[i].y*temp.weights[iter];
        }
        //std::cout<<std::endl;

        for(int i = 0; i < 4; i++){
        //  std::cout<<integrals.MATRIX_4P_DN[iter][i]<<"*"<<elmnt[i].x<<std::endl;
            if(POINTS == 2)
                J[1][0] +=  integrals.MATRIX_4P_DN[iter][i]*elmnt[i].x;
            else
                J[1][0] +=  integrals.MATRIX_4P_DN[iter][i]*elmnt[i].x*temp.weights[iter];
        }
        //std::cout<<std::endl;

        for(int i = 0; i < 4; i++){
        //  std::cout<<integrals.MATRIX_4P_DN[iter][i]<<"*"<<elmnt[i].y<<std::endl;
            if(POINTS == 2)
                J[1][1] +=  integrals.MATRIX_4P_DN[iter][i]*elmnt[i].y;
            else
                J[1][1] +=  integrals.MATRIX_4P_DN[iter][i]*elmnt[i].y*temp.weights[iter];
        }
    }
    else
        std::cout<<"Error: integrationpoints must be 2 or 3!"<<std::endl;
}

// COUNT JACOBI MATRIX FOR GRID
// ELEMENT_1 - przechowuje macierze dN/dKsi oraz dN/dn
// void jacobi_grid(grid GRID, ELEMENT_2D_ARRAYS ELEMENT_1){
//     // tablica 1d wezlow w siatce
//     node *nodes = new node[GRID.nodes_number];
//     element *elements = new element[GRID.elm_number];
//     // przepisywanie wezlow do tablicy 1d
//     int counter = 0;
//     for(int i = 0; i < GRID.nW; i++){
//         for(int j = 0; j < GRID.nH; j++){
//             nodes[counter] = GRID.nodes[i][j];
//                 //std::cout<<"NODES ["<<counter<<"] : x = "<<nodes[counter].x<<" y = "<<nodes[counter].y<<std::endl;
//                 counter++;
//         }
//     }
    
//     // przepisywanie elementow do tablicy 1d
//     counter = 0;
//     for(int i = 0; i < GRID.nW-1; i++){
//         for(int j = 0; j < GRID.nH-1; j++){
//             elements[counter] = GRID.elements[i][j];
//                 /*std::cout<<"ELEMENTS ["<<counter<<"] : id1 = "<<elements[counter].id[0]
//                 <<" id2 = "<<elements[counter].id[1]<<" id3 = "<<elements[counter].id[2]<<
//                 " id4 = "<<elements[counter].id[3]<<std::endl;*/
//                 counter++;
//         }
//     }

//     // deklaracja macierzy do przechowania Jakobianu
//     double **jacobian = new double*[2];
//     for(int i = 0; i < 2; i++)
//         jacobian[i] = new double[2];
    
//     // MAIN ALGHORITM
//     // oblicz Jakobian dla kazdego elementu -> elemnt ma id wezlow
//     for(int i = 0; i < GRID.elm_number; i++){
//         // dla kazdego elementu stworz temp_element z 4 wezlami tego elementu
        
//         node temp_element[4];
//         temp_element[0] = nodes[elements[i].id[0]];
//         temp_element[1] = nodes[elements[i].id[1]];
//         temp_element[2] = nodes[elements[i].id[2]];
//         temp_element[3] = nodes[elements[i].id[3]];
//         std::cout<<"=== COUNTING JACOBI MATRIX FOR ELEMENT "<<i<<" ==="<<std::endl;
//         // inner loop -> 4 punkty calkowania dla kazdego 
//         for(int j = 0; j < 4; j++){
//                 double **jacobian = new double*[2];
//                     for(int k = 0; k < 2; k++){
//                         jacobian[k] = new double[2];
//                     }
//             jacobi(/*i,*/j, jacobian, jacobian, ELEMENT_1, temp_element);
//                     for(int k = 0; k < 2; k++){
//                         delete jacobian[k];
//                     }
//                     delete[] jacobian;
//         }
//     }
//     delete[] nodes;
//     delete[] elements;
// }

// FUNCTION TO ADD ROW TO H MATRIX
// matrix_dx, matrix_dy - temp matrices, iteration - punkt calkowania, element - integral dksi dn arrays
void add_row(const int integration_points, double **matrix_dx, double **matrix_dy, int iteration, double **jacobi, ELEMENT_2D_ARRAYS element){

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
double** count_H_matrix(const int integration_points, int iteration, double **jacobi, ELEMENT_2D_ARRAYS element, double cond, bool if_print){

    int integrals = integration_points*integration_points;

    double **DX_matrix = new double*[integrals];
    double **DY_matrix = new double*[integrals]; // temp matrices inside integral
    double ***H = new double**[4];
    for(int i = 0; i  < integrals; i++){
        H[i] = new double*[4];
        for(int j = 0; j < 4; j++)
            H[i][j] = new double[4];
        DX_matrix[i] = new double[4];
        DY_matrix[i] = new double[4];
    }



    for(int i = 0; i < integrals; i++){
        // std::cout<<"here"<<std::endl;
        add_row(integration_points, DX_matrix, DY_matrix, i, jacobi, element);
    }
    // transponowany rzad dla kazdego punktu calkowania
    double **X_mat = new double*[4];
    double **Y_mat = new double*[4];
    for(int i = 0; i < 4; i ++){
        X_mat[i] = new double[4];
        Y_mat[i] = new double[4];
    }

    double k = cond; // k(t) factor
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
            print_custom(H[iteration], 4, 4, 4);
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
void count_H_matrix_grid(const int integration_points, grid *GRID, ELEMENT_2D_ARRAYS ELEMENT_1, double cond, bool if_print){
    int integrals = integration_points*integration_points;
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
        for(int j = 0; j < integrals; j++){
                double **jacobian = new double*[2];
                    for(int k = 0; k < 2; k++){
                        jacobian[k] = new double[2];
                    }
            jacobi(integration_points, j, jacobian, ELEMENT_1, temp_element);
            H_sum = sum_matrix(count_H_matrix(integration_points, j, jacobian, ELEMENT_1,cond, if_print), H_sum, 4);
            for(int k = 0; k < 2; k++)
                delete jacobian[k];
            delete[] jacobian;
        }
        if(if_print){
            std::cout<<"SUMMED H MATRIX FOR ELEMENT "<<i<<std::endl;
            print_custom(H_sum, 4, 4, 4);}

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

// COUNT H AGGREGATED
// void H_aggregate(double **result, int nodes_number, grid GRID){
//     for(int i = 0; i < nodes_number; i++){
//         for(int j = 0; j < nodes_number; j++){

//         }
//     }
// }

void count_H_2d_aggregation(grid GRID, double **result){
    // tablica 1d wezlow w siatce
    node *nodes = new node[GRID.nodes_number];
    element *elements = new element[GRID.elm_number];
    // przepisywanie wezlow do tablicy 1d
    int counter = 0;
    for(int i = 0; i < GRID.nW; i++){
        for(int j = 0; j < GRID.nH; j++){
            nodes[counter] = GRID.nodes[i][j];
                counter++;
        }
    }
    
    // przepisywanie elementow do tablicy 1d
    counter = 0;
    for(int i = 0; i < GRID.nW-1; i++){
        for(int j = 0; j < GRID.nH-1; j++){
            elements[counter] = GRID.elements[i][j];
            counter++;
        }
    }

    for(int i = 0; i < GRID.elm_number; i++){ // main loop = elements loop
        // FOR EACH ELEMENT ALGHORITM
        // std::cout<<"MAIN ITERATOR : i = "<<i<<std::endl;
        for(int j = 0; j < 4; j++){
            for(int k = 0; k < 4; k++){
            //    std::cout<<"here: j = "<<j<<" id[j] = "<<elements[i].id[j]<<" k = "<<k<<" id[k] = "<<elements[i].id[k]<<std::endl;
                result[elements[i].id[j]][elements[i].id[k]] += elements[i].H[j][k];
                if(elements[i].HBC != NULL){
                    // std::cout<<"adding "<<elements[i].HBC[j][k]<<std::endl;
                    result[elements[i].id[j]][elements[i].id[k]] += elements[i].HBC[j][k];
                }
            }
        }
    }
}

void count_C_2d_aggregation(grid GRID, double **result){
    // tablica 1d wezlow w siatce
    node *nodes = new node[GRID.nodes_number];
    element *elements = new element[GRID.elm_number];
    // przepisywanie wezlow do tablicy 1d
    int counter = 0;
    for(int i = 0; i < GRID.nW; i++){
        for(int j = 0; j < GRID.nH; j++){
            nodes[counter] = GRID.nodes[i][j];
                counter++;
        }
    }
    
    // przepisywanie elementow do tablicy 1d
    counter = 0;
    for(int i = 0; i < GRID.nW-1; i++){
        for(int j = 0; j < GRID.nH-1; j++){
            elements[counter] = GRID.elements[i][j];
            counter++;
        }
    }

    for(int i = 0; i < GRID.elm_number; i++){ // main loop = elements loop
        // FOR EACH ELEMENT ALGHORITM
        // std::cout<<"MAIN ITERATOR : i = "<<i<<std::endl;
        for(int j = 0; j < 4; j++){
            for(int k = 0; k < 4; k++){
            //    std::cout<<"here: j = "<<j<<" id[j] = "<<elements[i].id[j]<<" k = "<<k<<" id[k] = "<<elements[i].id[k]<<std::endl;
                result[elements[i].id[j]][elements[i].id[k]] += elements[i].C_mat[j][k];
            }
        }
    }
}


// *** FINAL COUNTING
// H, C_mat - H, C matrices
// step - simulation step time, iteration - number of simulation iteration, size - result matrix size
double **iteration_H_matrix(double **H, double **C_mat, double step, int size){
    double dT_fixed = 1/step; // count time passed - C matrix multiplied by 1/dT = matrix/dT
    const_multi_matrix(dT_fixed, C_mat, size);
    return sum_matrix(H, C_mat, size);
}


// GAUSS MATRIX SOLVER
void swap(double& a, double& b) {
	double buf = a;
	a = b;
	b = buf;
}

bool gaussMatrix(double **factors, int n, int begin){

    if(begin+1==n){ // recursion break condition - entire matrix is done
       return true; 
    }
    
    int maxIndex = begin;
    double maxValue = fabs(factors[begin][begin]);

    // ALGHORITM - find row with max value
    for(int i = begin; i<n;i++){
        if(fabs(factors[i][begin]>maxValue)){
            maxIndex = i;
            maxValue = factors[i][begin];
        }
    }
    
    if(maxIndex != begin){ // if first component isnt the biggest
        for(int i = begin; i <= n; i++){
            swap(factors[begin][i],factors[maxIndex][i]); // swap entire row
        }
    }


    //GAUSS ALGHORITM - TRIANGLE MATRIX - start from second row (index 1)
    for(int i = begin+1; i<n; i++){ // row loop
        const double multiplier = factors[i][begin]/factors[begin][begin]; // const multiplier for entire row - a[1][0]
            for(int j = begin; j <= n; j++){ // column loop
            factors[i][j] -= multiplier*factors[begin][j]; // ex. a[1][0] = a[1][0]/a[0][0]*a[0][0] (entire row 1) 
            if(factors[i][j]<EPSILON&&factors[i][j]>-EPSILON){ // double mistakes
                factors[i][j] = 0.0;
            }
        }
    }

    gaussMatrix(factors, n, ++begin); // recursion -> with [n-1][n-1] matrix
}


void calcSolution(double **factors, int size, double *results){
    size--; // array counted from 0 

    for(int i = size; i>=0; i--){
        double tempSum = 0.0;
        for(int j = i; j < size; j++){
            tempSum += results[j+1]*factors[i][j+1]; // calculate sum of non-counted factors
        }    

        results[i] = (factors[i][size+1] - tempSum)/factors[i][i]; // calculate result[i]
    }
}

void printResults(double **factors, int size, double *results, bool if_print){
    std::cout<<"-------------------------"<<std::endl;
    std::cout<<"         RESULTS         "<<std::endl;
    std::cout<<"-------------------------"<<std::endl;
    if(if_print){
        std::cout<<"TRANSFORMATED MATRIX: "<<std::endl<<std::endl;
        for(int i = 0; i < size; i++){
            for(int j = 0; j<= size; j++){
                if(j==size){
                    std::cout<<" |";
                }
                std::cout<<std::setprecision(3)<<factors[i][j]<<"    ";
            }
            std::cout<<std::endl;
        }
    }
    std::cout<<std::endl;
    std::cout<<"Temperature solution: "<<std::endl;
    for(int i = 0; i < size; i++){
        std::cout<<std::setprecision(4)<<"T"<<i+1<<" = "<<results[i]<<", ";
    }
    std::cout<<std::endl;
}

void printInitial(double **factors, int size){
        std::cout<<"-------------------------"<<std::endl;
        std::cout<<"      INITIAL MATRIX     "<<std::endl;
        std::cout<<"-------------------------"<<std::endl;
        // print initial matrix
        for(int i = 0; i < size; i++){
            for(int j = 0; j<= size; j++){
                if(j==size){
                    std::cout<<" |";
                }
                std::cout<<std::setprecision(6)<<factors[i][j]<<"    ";
            }
            std::cout<<std::endl;
        }
}

// find extreme values
double* find_extremes(double *vector, int size){
    double *result = new double[2];
    double max = vector[0];
    double min = max;
    for(int i = 0; i < size; i++){
        if(vector[i] < min){
            min = vector[i];
            result[0] = min;
        }
        if(vector[i] >= max){
            max = vector[i];
            result[1] = max;
        }
    }
    return result;
}


/*** 
 * MAIN LOOP COUNTING 
 * 
 * 2 punkty calkowania
 * */

void final_alghoritm_2POINTS(const int integration_points, InitData data, bool print_init, int precision = 6){

    // read input
    double simulation_time = data.SimulationTime;
    double time_step = data.SimulationStepTime;
    double cond = data.Conductivity;
    double density = data.Density;
    double alpha = data.Alfa;
    double tot = data.Tot;
    double init_temp = data.InitialTemp;
    double specific_heat = data.SpecificHeat;

    grid grid_1 = data.grid_input;

    int GSIZE = grid_1.nodes_number; // rozmiar siatki
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
        element.count_matrix_4points_ksi(matrix_n_4points);
        ELEMENT_2D_ARRAYS ELEMENT_1 = ELEMENT_2D_ARRAYS(matrix_ksi_4points, matrix_n_4points);
        //     // MATRIX 4X4, dN/dksi
        // std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dksi ===="<<std::endl;
        element.count_matrix_4points_ksi(matrix_ksi_4points);
        // print_custom(matrix_ksi_4points, 4, 4, 4);

        // MATRIX 4X4, dN/dn
        // std::cout<<"==== RESULT FOR MATRIX[4][4]"<<", FUN = dN/dn ===="<<std::endl;
        element.count_matrix_4points_n(matrix_n_4points);
        // print_custom(matrix_n_4points, 4, 4, 4);
        count_H_matrix_grid(integration_points, &grid_1, ELEMENT_1, cond ,false);

        double **matrix_4points_HBC = new double*[8];

        for(int i = 0; i < 8; i++) matrix_4points_HBC[i] = new double[4];


        el_4_2d_HBC element_HBC = el_4_2d_HBC();
        element_HBC.count_matrix_4points(matrix_4points_HBC);

        count_HBC_grid(&grid_1, false, matrix_4points_HBC, alpha);
        if(print_init){
            std::cout<<"==== RESULT FOR HBC MATRIX - SHAPE FUNCTIONS===="<<std::endl;
            print_custom(matrix_4points_HBC, 8, 4, 4);
            grid_1.print_HBC_matrices();
            // count_H_2d_aggregation(grid_1, aggregation_2d_H_matrix);
        }

        // P vectors
        count_P_vectors_grid(&grid_1, false, matrix_4points_HBC, alpha, tot);
        // grid_1.print_P_vectors();

        // P VECTOR AGGREGATION
        // double *P_aggregated = new double[GSIZE];
        // count_P_vectors_2d_aggregation(grid_1, P_aggregated);
        


        // C MATRIX
        // shape function values
        el_4_2d_C element_C = el_4_2d_C();
        double **matrix_4points_C = new double*[4];
        for(int i = 0; i < 4; i++) matrix_4points_C[i] = new double[4];


        element_C.count_matrix_4points(matrix_4points_C);

        if(print_init){
            std::cout<<"==== RESULT FOR C MATRIX - SHAPE FUNCTIONS===="<<std::endl;
            print_custom(matrix_4points_C, 4, 4, 4);
        }

        double **matrix_test_C = new double*[4];
        for(int i = 0; i < 4; i++) matrix_test_C[i] = new double[4];

        count_C_matrix_grid(&grid_1, matrix_4points_C, specific_heat, density, ELEMENT_1, false);

        double *T0 = new double[GSIZE];
        for(int i = 0; i < GSIZE; i++) T0[i] = init_temp;

        // iteration P vector
        double *result_P_iteration = new double[GSIZE];
        for(int i = 0; i < GSIZE; i++) result_P_iteration[i] = 0;

        int iterations = (int)simulation_time/time_step;

        double **aggregation_2d_H_matrix = new double*[GSIZE];
        double **aggregation_2d_C_matrix = new double*[GSIZE];
        double *P_aggregated = new double[GSIZE];
        double **result_H_iteration = new double*[GSIZE];
        for(int k = 0; k < GSIZE; k++) result_H_iteration[k] = new double[GSIZE];

         for(int k = 0; k < GSIZE; k++){
            aggregation_2d_H_matrix[k] = new double[GSIZE];
            aggregation_2d_C_matrix[k] = new double[GSIZE];
            for(int j = 0; j < GSIZE; j++){ 
                aggregation_2d_H_matrix[k][j] = 0;
                aggregation_2d_C_matrix[k][j] = 0;
            }
        }
        count_H_2d_aggregation(grid_1, aggregation_2d_H_matrix);
        count_C_2d_aggregation(grid_1, aggregation_2d_C_matrix);
        count_P_vectors_2d_aggregation(grid_1, P_aggregated);
        result_H_iteration = iteration_H_matrix(aggregation_2d_H_matrix, aggregation_2d_C_matrix, time_step, GSIZE);
        result_P_iteration = multi_matrix_vector(aggregation_2d_C_matrix, T0, GSIZE);

        fix_values(P_aggregated, GSIZE);
        sum_vectors(result_P_iteration, P_aggregated, GSIZE);
        double **result_iteration = new double*[GSIZE];

        for(int k = 0; k<=GSIZE; k++){ // <= - factors is [x][x+1] matrix - last column is column of const.
            result_iteration[k] = new double[GSIZE+1];
        }
        for(int k = 0; k<GSIZE; k++){ // read factors
            for(int j = 0; j<=GSIZE; j++){
                if(j == GSIZE){
                    result_iteration[k][j] = result_P_iteration[k];
                }
                else{
                    result_iteration[k][j] = result_H_iteration[k][j];
                }
            }
        }

        // MAIN LOOP
    for(int i = 0; i < iterations; i++){
        if(print_init){
            std::cout<<"ITERATION 0 "<<std::endl;
            std::cout<<"============ H + HBC MATRICES AGGREGATION ============"<<std::endl;
            print_custom(aggregation_2d_H_matrix, GSIZE, GSIZE, 4);
            std::cout<<"============ C MATRICES AGGREGATION ============"<<std::endl;
            print_custom(aggregation_2d_C_matrix, GSIZE, GSIZE, 4);
            std::cout<<"============ P VECTORS AGGREGATION ============"<<std::endl;
            print_vector(P_aggregated, GSIZE, 6);
            std::cout<<std::endl;
            std::cout<<"============       [H]+[C]/dT      ============"<<std::endl;
            print_custom(result_H_iteration, GSIZE, GSIZE, 4);
        }

        result_P_iteration = multi_matrix_vector(aggregation_2d_C_matrix, T0, GSIZE);

        fix_values(P_aggregated, GSIZE);
        sum_vectors(result_P_iteration, P_aggregated, GSIZE); 

        // GAUSS SOLVING
        
        for(int k = 0; k<GSIZE; k++){ // read factors
            for(int j = 0; j<=GSIZE; j++){
                if(j == GSIZE){
                    result_iteration[k][j] = result_P_iteration[k];
                }
                else{
                    result_iteration[k][j] = result_H_iteration[k][j];
                }
            }
        }

        // Gauss counting
        if(gaussMatrix(result_iteration, GSIZE, 0)){ // if matrix is solved correctly
            double *results = new double[GSIZE]; // results array
            calcSolution(result_iteration, GSIZE, results);
            // printInitial(result_iteration, GSIZE);
            // printResults(result_iteration, GSIZE, results, true);
            std::cout<<std::endl;
            std::cout<<"----------------------------"<<std::endl;
            std::cout<<"ITERATION "<<i<<" TIME = "<<(i+1)*time_step<<std::endl;
            std::cout<<"----------------------------"<<std::endl;
            // printResults(result_iteration, GSIZE, results, false);
            std::cout<<"MIN TEMPERATURE: "<<std::setprecision(6)<<find_extremes(results, GSIZE)[0]<< " MAX TEMPERATURE = "<<find_extremes(results, GSIZE)[1]<<std::endl;

            std::cout<<"RES:"<<std::endl;
            print_vector(results, GSIZE, 6);
            for(int k = 0; k < GSIZE; k++)
                T0[k] = results[k]; // overwrite temp. vector with new temperature in each step (iteration)
            std::cout<<"T0 after:"<<std::endl;
            print_vector(T0, GSIZE, 6);
            // delete[] results;
        }

        // delete[] aggregation_2d_H_matrix;
        // delete[] aggregation_2d_C_matrix;
        // delete[] P_aggregated;
        // delete[] result_H_iteration;
    }
}
#endif