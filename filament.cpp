//
//  filament.cpp
//  Filament
//
//  Created by Rui Ma on 1/6/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//

#include "filament.h"
#include <iostream>
using namespace std;

void Filament::SetCenterOfMass()
{
    cout << "Please type the center of mass of the filament.\n" << endl;
    cin >> center_of_mass[0] >> center_of_mass[1] >> center_of_mass[2];
}

void Filament::PrintCenterOfMass()
{
    cout << center_of_mass[0] << ' ' << center_of_mass[1] << ' ' << center_of_mass[2] << endl;
    
}