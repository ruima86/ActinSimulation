//
//  filament.h
//  Filament
//
//  Created by Rui Ma on 1/6/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//
#ifndef FILAMENT_H
#define FILAMENT_H

class Filament
{ public:
    double center_of_mass[3];
    double orientation[3];
    double force[3];
    double torque[3];
    double fric_f[3*3];
    double fric_T[3*3];
    double length;
    int Occupied[200];
    int index;
    int N_Connect_Left;
    int N_Connect_Right;
    void PrintCenterOfMass();
    void SetCenterOfMass();
};

#endif

