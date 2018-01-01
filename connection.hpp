//
//  connection.hpp
//  GSLTest
//
//  Created by Rui Ma on 2/6/17.
//  Copyright © 2017 Yale University. All rights reserved.
//

#ifndef connection_hpp
#define connection_hpp

#include <iostream>
#include <armadillo>
#include <list>

using namespace std;
using namespace arma;

class Filament;

class Connection
{
public:
    Filament* f1;
    Filament* f2;
    list<int>::iterator n1;
    list<int>::iterator n2;
    // Here we use a pointer that is linked to an element of filament's OccupiedList. If the element in OccupiedList is updated,
    // n1 will update simaltaneously. 
    double s1;
    double s2;
    int index1;
    int index2;
    double phi1_1;
    double theta1_1;
    double phi1_2;
    double theta1_2;
    double phi2_1;
    double theta2_1;
    double phi2_2;
    double theta2_2;
    double len;
    double lifetime;
    double transEnergy;
    double torsEnergy;
    Connection();
    Connection(Filament *, Filament *, int, int, double const);
    ~Connection();
    bool GetForceAndTorque(const double kappa, const double r0, const double kappaR1, const double kappaR2, const double Cross_angle, const double modular_angle, const double Max_distance, const double Max_angle, const double dt);
    
    friend bool operator == (Connection &a, Connection &b);
};

#endif /* connection_hpp */
