#ifndef GRID
#define GRID

class node{};
class element{};
class grid{
    public: 
        grid(float high, float widt, int number_H,  int number_W);
        void printNodes();
        void printElements();
};

void createNodes(grid gr);

void createElements(grid gr, int node_heigth);


#endif