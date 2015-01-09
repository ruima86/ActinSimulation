//
//  BiListConnection.cpp
//  Filament
//
//  Created by Rui Ma on 1/7/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//
#include <iostream>
#include "bilistconnection.h"
using namespace std;

node::node()
{
    s1=0;
    s2=0;
    next=0;
}

void Node_Insert(BiList p, BiList q)
{
    if(p->next!=0)
        p->next->prior=q;
    q->next=p->next;
    p->next=q;
    q->prior=p;
}

void Node_Delete(BiList q)
{
    q->prior->next=q->next;
    if (q->next!=0) {
        q->next->prior=q->prior;
    }
    delete q;
}

void Node_Print(BiList p)
{
    while (p!=0) {
        cout << p->s1 << endl;
        p=p->next;
    }
}
