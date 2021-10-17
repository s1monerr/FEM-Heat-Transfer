// IMPLEMENTACJA CALKOWANIA METODA GAUSSA - warianty 2 i 3 punktowe
#include <iostream>
#include <iomanip>
#include <cmath>

double func_1(double x){ // 5x^2 + 3x + 6
    double sum = 0.0;
    sum += 5*pow(x,2);
    sum += 3*x;
    sum += 6;

    return sum;
}

double func_2(double x, double y){ // 2 zmiennych funkcja
    double sum = 0.0;
    sum += 5*pow(x,2)*pow(y,2);
    sum += 3*x*y;
    sum += 6;

    return sum;
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

// calka Gaussa 1d
double gaussIntegral_1d(int points, factors fac){

    double sum = 0.0;

    if(points == 2){
        for(int i = 0; i < points; i++){
            sum += func_1(fac.point_2p[i])*fac.w_2p[i];
        }
        return sum;
    }

    else if(points == 3){
        for(int i = 0; i < points; i++){
            sum += func_1(fac.point_3p[i])*fac.w_3p[i];
        }
        return sum;
    }

    else{
        std::cout<<"error: wrong points number."<<std::endl;
        return -1;
    }
}

// calka Gaussa 2d
double gaussIntegral_2d(int points, factors fac){

    double sum =0.0;

    if(points == 2){
        for(int i = 0; i < points; i++){
            for(int j = 0; j < points; j++){
                sum += func_2(fac.point_2p[i],fac.point_2p[j])*fac.w_2p[i];
            }
        }
    }

    else if(points == 3){
        for(int i = 0; i < points; i++){
                for(int j = 0; j < points; j++){
                    sum += func_2(fac.point_3p[i], fac.point_3p[j])*fac.w_3p[i]*fac.w_3p[j];
                }
            }
        }

    else   
        return -1;
        
    return sum;
}


int main(){
    factors factors_1 = factors();
    std::cout<<"Gauss 1d, 2 points, result: "<<gaussIntegral_1d(2, factors_1)<<std::endl;
    std::cout<<"Gauss 1d, 3 points, result: "<<gaussIntegral_1d(3, factors_1)<<std::endl;
    std::cout<<"Gauss 2d, 2 points, result: "<<gaussIntegral_2d(2, factors_1)<<std::endl;
    std::cout<<"Gauss 2d, 3 points, result: "<<gaussIntegral_2d(3, factors_1)<<std::endl;
}