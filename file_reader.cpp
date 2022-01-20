#include <iostream>
#include <fstream>
#include <algorithm>
#include "FEM_lib.cpp"
#include "heapsort.cpp"

// Convert string to and double
double string_converter(const std::string& input){
    std::stringstream ss(input);
    double retval;
    ss >> retval;
    return retval;
}

int int_string_converter(const std::string& input){
    std::stringstream ss(input);
    int retval;
    ss >> retval;
    return retval;
}

int read_file(const std::string file, InitData *data/*, grid *GRID*/){
    std::fstream read(file);

    std::string temp;

    if(read.is_open()){
        // READ INITIAL VALUES
        // 1 - read simulation conditions
        int i = 0;
        while(true){
            read >> temp;
            if(temp == "SimulationTime"){
                read >> temp;
                data->SimulationTime = string_converter(temp);
            }

            else if(temp == "SimulationStepTime"){
                read >> temp;
                data->SimulationStepTime = string_converter(temp);
            }

            else if(temp == "Conductivity"){
                read >> temp;
                data->Conductivity = string_converter(temp);
            }

            else if(temp == "Alfa"){
                read >> temp;
                data->Alfa = string_converter(temp);
            }

            else if(temp == "Tot"){
                read >> temp;
                data->Tot = string_converter(temp);
            }

            else if(temp == "InitialTemp"){
                read >> temp;
                data->InitialTemp = string_converter(temp);
            }

            else if(temp == "Density"){
                read >> temp;
                data->Density = string_converter(temp);
            }

            else if(temp == "SpecificHeat"){
                read >> temp;
                data->SpecificHeat = string_converter(temp);
            }

            else if(temp == "Nodes"){
                // double read instruction to avoid blank character
                read >> temp;
                read >> temp;
                data->Nodes_number = int_string_converter(temp);
            }

            else if(temp == "Elements"){
                read >> temp;
                read >> temp;
                data->Elements_number = int_string_converter(temp);
            }

            else if(temp == "Width"){
                read >> temp;
                data->Width = int_string_converter(temp);
            }

            else if(temp == "Heigth"){
                read >> temp;
                data->Heigth = int_string_converter(temp);
                break;
            }

            else{
                std::cout<<"Error 0: incorrent input file!"<<std::endl;
                return -1;
            }
        }

        grid GRID = grid(data->Heigth, data->Width);

        // nodes temp array
        // node 
        node **nodes = new node*[data->Width];
        for(int i = 0; i < data->Width; i++) nodes[i] = new node[data->Heigth];

        // READING NODES
        read >> temp;
        bool end = false; // node loop end flag
        if(temp != "*Node")
            std::cout<<"Error code 1: incorrenct input file!"<<std::endl;
        for(int i = 0; i < data->Width; i++){
            // std::cout<<"i = "<<i<<" wit = "<<data->Width<<std::endl;
            for(int j = 0; j < data->Heigth; j++){
                read >> temp;
                if(temp == "*Element"){ // if element - finish node loop
                    end = true;
                    break;
                }

                // create node
                double x, y;
                read >> temp;
                x = string_converter(temp);
                read >> temp;
                y = string_converter(temp);

                // delete coma from input
                temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
                GRID.nodes[i][j] = node(x, y, 0); // change node value in grid
            }

            if(end){
                end = false; // for next loop
                break;
            }
        }

        // READING ELEMENTS
        read >> temp;
        if(temp == "*Element,"){
            read >> temp; // avoid blank space
            // element reader loop
            for(int i = 0; i < data->Width - 1; i++){
                for(int j = 0; j < data->Heigth - 1; j++){
                    read >> temp;

                    if(temp == "*BC"){ // if element - finish node loop
                        end = true;
                        break;
                    }
                    
                    //read >> temp; // avoid element number
                    
                    int id0, id1, id2, id3;
                    read >> temp;
                    // delete coma from input
                    temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
                    id0 = int_string_converter(temp);
                    id0--;

                    read >> temp;
                    temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
                    id1 = int_string_converter(temp);
                    id1--;

                    read >> temp;
                    temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
                    id2 = int_string_converter(temp);
                    id2--;

                    read >> temp;
                    temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
                    id3 = int_string_converter(temp);
                    id3--;
                    GRID.elements[i][j] = element(id0, id1, id2, id3);

                }
                if(end){
                    end = false;
                    break;
                }
            }
        }

        else{
            std::cout<<"Error 2: incorrect input file!"<<std::endl;
            return -1;
        }

        read >> temp;
        if(temp == "*BC"){
            int rozmiar = data->Width*data->Heigth;
            int *boundary_conditions = new int[rozmiar];
            for(int i = 0; i < rozmiar; i++) boundary_conditions[i] = rozmiar + 100; // elementy > rozmiar: bez warunkow brzegowych
            
            // read >> temp; // iteration 0
            // temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
            // int  temp_int = int_string_converter(temp);
            // int counter = 0; // BC values index from 1
            // for(int i = 0; i < data->Width; i++){
            //     for(int j = 0; j < data->Heigth; j++){
            //         read >> temp;
            //         temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
            //         temp_int = int_string_converter(temp);
            //         std::cout<<"saving: "<<temp_int<<std::endl;
            //         boundary_conditions[counter] += temp_int;
            //         // if(temp_int == counter){ // if iteration is equal to BC 
            //         //     GRID.nodes[i][j].weight = 1;
            //         //     read >> temp; // find next BC node
            //         //     temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
            //         //     temp_int = int_string_converter(temp);
            //         // }
            //         // else do nothing
            //         ++counter;
            //     }
            // }
            int counter = 0;
            while(!read.eof()){
                read >> temp;
                temp.erase(remove(temp.begin(), temp.end(), ','), temp.end());
                int temp_int = int_string_converter(temp);
                temp_int--;
                boundary_conditions[counter] = temp_int;
                counter++;
            }
            
            heapsort(boundary_conditions, rozmiar); // sort boundary conditions

            int counter_conditions = 0;
            counter = 0;
            // save conditions to grid
            for(int i = 0; i < data->Width; i++){
                for(int j = 0; j < data->Heigth; j++){
                    if(boundary_conditions[counter_conditions] == counter){ // if boundary condition in node array
                        GRID.nodes[i][j].weight = 1;
                        ++counter_conditions;
                    }
                    ++counter;
                }
            }
        }
        else{
            std::cout<<"Error code 3: incorrect input file!"<<std::endl;
            return -1;
        }
        data->grid_input = GRID;
    }
    return 0;
}