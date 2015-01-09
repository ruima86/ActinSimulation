//
//  BiList.h
//  Filament
//
//  Created by Rui Ma on 1/7/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//
#ifndef NODE_H
#define NODE_H


#include "filament.h"

class node
{ public:
    Filament* f1;
    Filament* f2;
    double s1;
    double s2;
    int index1;
    int index2;
    node* prior;
    node* next;
    node()
    {
        s1=0;
        s2=0;
        next=0;
    };
};

typedef node *BiList;

#endif


