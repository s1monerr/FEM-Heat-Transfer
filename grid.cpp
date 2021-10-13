// IMPLEMENTACJA SIATKI DO MES


#include <iostream>
#include <iomanip>

#define WIDTH 0.1 // szerokosc
#define HEIGTH 0.2 // wysokosc
#define NODES_HEIGTH 5 // ilosc wezlow wysokosc
#define NODES_WIDTH 4 // ilosc wezlow szerokosc

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
        int id1, id2, id3, id4; // id wezlow

        element(int i1, int i2, int i3, int i4){
            id1 = i1;
            id2 = i2;
            id3 = i3;
            id4 = i4;
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
                std::cout<<"id1 = "<<elements[i][j].id1<<", id2 = "<<elements[i][j].id2<<", id3 = "
                <<elements[i][j].id3<<", id4 = "<<elements[i][j].id4<<std::endl;
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
            gr.elements[i][j].id1 = j + w - node_heigth; // zawsze dodawaj w, ale mmiejsze o jedna wysokosc od id2
            gr.elements[i][j].id2 = j + w;
            // kolejne id sa odpowiednio wieksze o 1
            gr.elements[i][j].id3 = gr.elements[i][j].id2 + 1;
            gr.elements[i][j].id4 = gr.elements[i][j].id1 + 1;
        }
        w += node_heigth; // w kolejnej kolumnie zwieksz o kolejna liczbe wezlow w kolumnie - czyli wysokosc wezlow
    }

}

int main(){
    grid grid_1 = grid(HEIGTH, WIDTH, NODES_HEIGTH, NODES_WIDTH);
    
    createNodes(grid_1);
    grid_1.printNodes();
    createElements(grid_1, NODES_HEIGTH);
    grid_1.printElements();

}
