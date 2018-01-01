//
//  filament.hpp
//  GSLTest
//
//  Created by Rui Ma on 2/3/17.
//  Copyright © 2017 Yale University. All rights reserved.
//

#ifndef filament_hpp
#define filament_hpp


#include <iostream>
#include <armadillo>
#include <list>
#include "connection.hpp"

using namespace std;
using namespace arma;

class Connection;

class Filament
{ public:
    // geometry
    vec3 center_of_mass;
    vec3 orientation;
    vec3 rotation;
    double length;            // length of filament
    double delta;             // length of a single monomer
    double radius;            // radius of filament
    double pitchnumber;       // pitch number of filament, the angle between adjacent monomer is M_PI*(1+1.0/pitch);
    
    
    // mechanics
    vec3 force;
    vec3 torque;
    
    mat33 xi_t;        // tranlsaitonal friction
    mat33 xi_r;        // rotational friction
    mat33 xi_tr;       // translation-rotation coupling
    
    //connection
    list<int> OccupiedList; // store the position of crosslinkers
                            // to have a convenient access to the link information, it's better to use a list<pair<>> class;
    //  list<pair<int, list<Connection>::iterator>> OccupiedList;
    //Connection con;
    
    //topology
    int Branched[200]; // Branched[i]=1 means that a branching filament grows at the (i+1)th monomer
    
    //state
    int BarbedCapping;
    int PointedCapping;
    bool isMobile;
    
    //index
    int index;
    
    // methods
    Filament();        // constructor of a filament with random position in a box of [-50,50]^3, and random orientation and rotation, and 10 monomers long
    Filament(vec3 com);        // constructor of a filament with specified center of mass, but random orientation and rotation, and 10 monomers long
    Filament(int NumMonomers, const double InitialSize, const int InitialAlignment);
    Filament(int NumMonomers, const double InitialSize, const int InitialAlignment, const double pitch);
    Filament(int NumMonomers, const double cylinderradius);
    Filament(vec3 com, vec3 ort, double len);        // constructor of a filament with specified center of mass, orientation, and 10 monomers long,
    ~Filament();
    void GetFricMatrix(double xi_t_para, double xi_t_perp, double xi_r_para, double xi_r_perp);
    mat33 GetCouplingMatrix(vec3);      // return the coupling matrix relative to a specified position
    mat33 GetRotationMatrix(vec3);      // return the rotation matrix relative to a specified position
    vec3 GetRotationMonomer(int n);     // return the rotation vector of monomer n
    vec3 GetPositionMonomer(int n);     // return the position of monomer n
    bool IsOccupiedMomomer(int n);      // return true if monomer n is occupied by a crosslinker
    void Move(const double kT, const double dt);
    void Grow(const int n);             // Grow the filament length by n monomers at the plus end
    void Shrink(const int n);           // Shrink the filament length by n monomers at the minus end
    
    friend bool operator == (Filament &f1, Filament &f2);
    friend bool operator != (Filament &f1, Filament &f2);
    friend bool operator < (Filament &f1, Filament &f2);
};



#endif /* filament_hpp */
