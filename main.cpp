//
//  main.cpp
//  GSLTest
//
//  Created by Rui Ma on 2/3/17.
//  Copyright © 2017 Yale University. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <list>
#include <random>
#include <armadillo>
#include "filament.hpp"
#include "connection.hpp"

using namespace std;
using namespace arma;

int main(int argc, const char * argv[]) {
    
    void CreateLink(list<Connection*> &ConList, Filament *f1, Filament *f2, const double R_cross, const double dt, const double Max_distance, const double Max_angle, const double Cross_angle);
    void PairwiseTraverse(list<Filament*>& mylist, list<Connection *> &conlist, void (*fun)(list<Connection*> &ConList, Filament *, Filament *,const double, const double, const double, const double, const double), const double R_cross, const double dt, const double Max_distance, const double Max_angle, const double Cross_angle);
    void OutputGeometry(list<Filament*> &FList, ofstream& outfile1, double tt);
    void OutputConnection(list<Connection*> &CList, ofstream& outfile2, double tt);
    void OutputForceTorque(list<Filament*> &FList, ofstream& outfile3, double tt);
    void Boundaryforce(list<Filament *> &FList, const string geometry, const double boxSize, const double height_c, const double radius_c, const double Bkappa);
    void StericForce(list<Filament*>& mylist, const double RepulsiveDistance, const double RepulsiveForce);
    
    
 
    
    /******************************/
    /* Properties of filaments    */
    /******************************/
    const int L_init= 50;                       // length of filaments (number of monomers)
    const int N_init= 2;                        // number of filaments
    
    
    /******************************/
    /* Properties of crosslinkers */
    /******************************/
    
    const double kappa= 0.0;                // Unit pN/nm          spring constant of crosslinkers
    const double kappa_R1= 0.0;           // Unit pN*nm          torsional stiffness of crosslinkers
    const double kappa_R2= 0.0;             // Unit pN*nm
    // kappa_R1 is to align the orientation of the monomer towards the crosslinker
    // kappa_R2 is to align the orientation of the filaments
    const double rest_length=10.0;          // Unit nm             rest length of crosslinkers
    const double Max_distance= 20.0;       // Unit nm             only filament pairs that are less than cric_raidus apart can form links
    const double Max_angle= M_PI;
    const double Cross_angle= 0;
    const double modular_angle = M_PI;     // This is for phase correction
    
    
    const double R_cross=0.0;
    // Unit s^(-1)         crosslinking rate R_cross*N1*N2 is the rate to create a link between two filaments of N1 and N2 monomers
    const double R_decross=0.0;    // de-crosslinking rate
    const double E_decross= 10.0;      // decrosslinking energy, in unit of kT
    
    
    /*************************************/
    /*  Properties of steric interaction */
    /*************************************/
    
    const double RepulsiveDistance = 6;
    const double RepulsiveForce = 100;
    
    
    
    /******************************/
    /* Properties of the medium   */
    /******************************/
    const double kT = 4.1;                     // pN nm
    const double viscosity = 1e-7;             // pN s /nm^2 viscosity of honey
    const double xi_para = (2*M_PI*viscosity*L_init*2.7)/(log(L_init*2.7/2.0/3.0)-0.2);
    const double xi_perp = (4*M_PI*viscosity*L_init*2.7)/(log(L_init*2.7/2.0/3.0)+0.84);
    const double xir_para = M_PI*viscosity*6.0*6.0*L_init*2.7;
    const double xir_perp = M_PI*viscosity*pow(L_init*2.7, 3.0)/3.0/(log(L_init*2.7/2.0/3.0)-0.66);
    

    
    /*************************************/
    /*  Geometry of the central cylinder */
    /*************************************/
    const double height_cylinder=1000;
    const double radius_cylinder=25;
    const double boxSize_cube= 500;
    const string geometry= "cube";       // four options : free, cube, cylinder, cubecylinder
    const double Bkappa= 10;      // boundary stiffness
    
    /*************************************/
    /*  Geometry of Initial conditions   */
    /*************************************/
    const double InitialSize=100;        // filaments are initiall
    const int InitialAlignment=0;        // 1 --> filaments are aligned in parallel
                                         // 0 --> filaments are aligned randomly
    
    
    
    /******************************/
    /* Simulation parameters      */
    /******************************/
    const double TT=10;                            // Unit s              Total simulation time
    vec timescales = {xi_para/(kappa+1),xi_perp/(kappa+1),xir_para/(kappa_R1+1),xir_perp/(kappa_R1+1)};
    const double dt=min(0.01*timescales.min(),1e-5);                          // Unit s              simulation time step
    const double dT=.1;                          // Unit s              recording time
    const int ensemble= 0;                         // the order number of simulations
    
    const double BranchAngle = M_PI*70.0/180.0;      // 70 degree angle of daughter filament relative to mother filament
    
    
    time_t ti = time(0);
    tm *now = localtime(&ti);
    
    ostringstream ostr1, ostr2, ostr3;
    
    ostr1 << "/home/fas/berro/rm796/scratch/cube_withoutBranch5/data/Trajectory_" << ensemble << "_Ninit_" << N_init << "_Linit_" << L_init << "_Rcross_" << R_cross << "_Rdecross_" << R_decross << "_Edecross_" << E_decross << "_kappa_" << kappa << "_kappaR1_" << kappa_R1 << "_kappaR2_" << kappa_R2 << "_Crossangle_" << Cross_angle << "_geometry_" << geometry << "_InitialSize_" << InitialSize << "_InitialAlign_" << InitialAlignment << "_MaxDistance_" << Max_distance << "_RadiusCylinder_" << radius_cylinder << "_BoxSize_" << boxSize_cube << "_Bkappa_" << Bkappa << "_date_" << now->tm_mon+1 << "-" << now->tm_mday << "-" << now->tm_year+1900 << ".txt";
    string ss = ostr1.str();
    ofstream outfile1(ss);
    
    ostr2 << "/home/fas/berro/rm796/scratch/cube_withoutBranch5/data/Connection_" << ensemble << "_Ninit_" << N_init << "_Linit_" << L_init << "_Rcross_" << R_cross << "_Rdecross_" << R_decross << "_Edecross_" << E_decross << "_kappa_" << kappa << "_kappaR1_" << kappa_R1 << "_kappaR2_" << kappa_R2 << "_Crossangle_" << Cross_angle << "_geometry_" << geometry << "_InitialSize_" << InitialSize << "_InitialAlign_" << InitialAlignment << "_MaxDistance_" << Max_distance << "_RadiusCylinder_" << radius_cylinder << "_BoxSize_" << boxSize_cube << "_Bkappa_" << Bkappa << "_date_" << now->tm_mon+1 << "-" << now->tm_mday << "-" << now->tm_year+1900 << ".txt";
    ss = ostr2.str();
    ofstream outfile2(ss);
    
    ostr3 << "/home/fas/berro/rm796/scratch/cube_withoutBranch5/data/Lifetime_" << ensemble << "_Ninit_" << N_init << "_Linit_" << L_init << "_Rcross_" << R_cross << "_Rdecross_" << R_decross << "_Edecross_" << E_decross << "_kappa_" << kappa << "_kappaR1_" << kappa_R1 << "_kappaR2_" << kappa_R2 << "_Crossangle_" << Cross_angle << "_geometry_" << geometry << "_InitialSize_" << InitialSize << "_InitialAlign_" << InitialAlignment << "_MaxDistance_" << Max_distance << "_RadiusCylinder_" << radius_cylinder << "_BoxSize_" << boxSize_cube << "_Bkappa_" << Bkappa << "_date_" << now->tm_mon+1 << "-" << now->tm_mday << "-" << now->tm_year+1900 << ".txt";
    ss = ostr3.str();
    ofstream outfile3(ss);
    
    outfile1 << setprecision(4);
    outfile1 << setiosflags(ios::fixed);
    outfile1 << setiosflags(ios::left);
    
    outfile2 << setprecision(4);
    outfile2 << setiosflags(ios::fixed);
    outfile2 << setiosflags(ios::left);
    
    outfile3 << setprecision(4);
    outfile3 << setiosflags(ios::fixed);
    outfile3 << setiosflags(ios::left);

    
    
    
    list<Filament*> FilamentList;              // Store collection of filaments
    list<Filament*>::iterator it_f;
    list<Connection*> ConnectionList;          // Store collection of crosslinkers
    list<Connection*>::iterator it_c;
    
    
    arma::arma_rng::set_seed_random();         // initialize the seeds of random number of generators;
//    vec3 c1 = {-5,0,0};
//    vec3 c2 = {5,0,0};
//    vec3 c3 = {8,0,0};
//    vec3 nn = {0,0,1};
//    double l = L_init*2.7;
    for (int i=0; i<0; i++) {
        if (geometry == "cube") {
            Filament *f = new Filament(30,50,1);
            f->index = -(i+1);
            f->center_of_mass(2) = -0.5*f->length;
            f->isMobile = false;
            FilamentList.push_back(f);
        }
    }
    
    for (int i=0; i<N_init; i++) {
        if (geometry == "cubecylinder" || geometry == "cylinder") {
            Filament *f = new Filament(L_init,radius_cylinder);
            f->index = i+1;
            FilamentList.push_back(f);
        }
        else
        {
            Filament *f = new Filament(L_init, InitialSize, InitialAlignment);
            f->index = i+1;
            FilamentList.push_back(f);
        }
    }
//    Filament *f1 = new Filament(c1,nn,l);
//    Filament *f2 = new Filament(c2,nn,l);
////    Filament *f3 = new Filament(c3,nn,l);
//    f1->index = 1;
//    f2->index = 2;
////    f3->index = 3;
//    FilamentList.push_back(f1);
//    FilamentList.push_back(f2);
////    FilamentList.push_back(f3);
//    
//    Connection *q = new Connection(f1, f2, 20, 20, 0);
//    ConnectionList.push_front(q);
//    q = new Connection(f1,f2,30,30,0);
//    ConnectionList.push_front(q);
//    q = new Connection(f1,f2,12,12,0,0,0,0);
//    ConnectionList.push_front(q); 
    
    

    for (double tt=0.0; tt<TT; tt+=dt) {
        
        // Go through all the filaments
        // Get Friction Matrix and reset force and torque to zero
        for (it_f = FilamentList.begin(); it_f!=FilamentList.end(); it_f++) {
            (*it_f)->GetFricMatrix(xi_para, xi_perp, xir_para, xir_perp);
            (*it_f)->force.zeros();
            (*it_f)->torque.zeros();
        }
        
        
        // Go through all the connections
        // Either calculate the force and torque or break the connection
        it_c=ConnectionList.begin();
        while (it_c!=ConnectionList.end()) {
            
            (*it_c)->torsEnergy = .5*kappa_R1*((*it_c)->theta1_1)*((*it_c)->theta1_1) + \
            .5*kappa_R1*((*it_c)->theta2_1)*((*it_c)->theta2_1);
            (*it_c)->transEnergy = .5*kappa*((*it_c)->len-rest_length)*((*it_c)->len-rest_length);
            
            double Energy = (*it_c)->torsEnergy + (*it_c)->transEnergy;
            
            if (randu()<(1-exp(-R_decross*exp(Energy/E_decross/kT)*dt))) {
                // break the connection
                Connection* p = (*it_c);
                //outfile3 << setw(20) << tt << setw(20) << p->lifetime << endl;
                ConnectionList.erase(it_c++);
                delete p;
                // !!!!!!!!!!
                // Note that this is to erase an element of a list while going through the list
                // Should be careful when we erase an element
            }
            else
            {
                (*it_c)->GetForceAndTorque(kappa, rest_length, kappa_R1, kappa_R2, Cross_angle, modular_angle, Max_distance, Max_angle, dt);
                ++it_c;
            }
        }
        if (geometry != "free") {
            Boundaryforce(FilamentList, geometry, boxSize_cube, height_cylinder, radius_cylinder, Bkappa);
        }
        
        
        StericForce(FilamentList, RepulsiveDistance, RepulsiveForce);
        
        
        if (fmod(tt+0.01*dt, dT)<dt) {
            OutputGeometry(FilamentList, outfile1, tt);
            OutputConnection(ConnectionList, outfile2, tt);
          //  OutputForceTorque(FilamentList, outfile3, tt);
        }
        
        // Go through all the filaments
        // move the filament

        for (it_f = FilamentList.begin(); it_f != FilamentList.end(); it_f++) {
            if ((*it_f)->isMobile) {
                (*it_f)->Move(kT, dt);
            }
            
        }
        
        
        // Go through each pair of filaments to build connection
        PairwiseTraverse(FilamentList, ConnectionList, &CreateLink, R_cross, dt, Max_distance, Max_angle, Cross_angle);
        
        
       
        
    }

    return 0;
}

void CreateLink(list<Connection*> &ConList, Filament *f1, Filament *f2, const double R_cross, const double dt, const double Max_distance, const double Max_angle, const double Cross_angle)
{
    int N1 = round(f1->length/f1->delta);
    int N2 = round(f2->length/f2->delta);
    double MaxRatio = 0.25;
    if (f1->OccupiedList.size()<MaxRatio*N1 && f2->OccupiedList.size()<MaxRatio*N2) {
        double AvgNumPairs = N1*N2*R_cross*dt;
        
        std::random_device rd;
        default_random_engine generator(rd());
        poisson_distribution<int>   distribution(AvgNumPairs); // Generate a Poisson distribution Pois(AvgNumPairs)
        int numTrials = distribution(generator);
        
        if (numTrials>0) {
            
            vec3 rot1,rot2,p1,p2,p12;
            //double dis, theta1, theta2, phi1, phi2;
            uniform_int_distribution<int> randint1(1,N1);
            uniform_int_distribution<int> randint2(1,N2);
            for (int i=0; i<numTrials; i++) {
                int m = randint1(generator);
                int n = randint2(generator);
                // !!!!!!!
                // Note that this is an approximation
                // The generated pairs (m, n) is likely to have repetitive elements which are supposed to be prevented from occuring.
                // The likelihood of having two repetitive pairs is C(numTrails,2)/(N1*N2), so this approximation is good if filaments
                // are long while dt is small, so is AvgNumPairs.
                if (!f1->IsOccupiedMomomer(m) && !f2->IsOccupiedMomomer(n) && f1->OccupiedList.size()<MaxRatio*N1 && f2->OccupiedList.size()<MaxRatio*N2) {
                    
                    p1 = f1->GetPositionMonomer(m);
                    p2 = f2->GetPositionMonomer(n);
                    double dis = norm(p1-p2);
                    if (dis < Max_distance) {
                        // !!!!!!!!!!
                        // Note that here only the distance match is verified, but not angular match.
                        Connection *u = new Connection(f1,f2,m,n,Cross_angle);
                        ConList.push_front(u);
                    }
                    
                }
                
            }
            
        }

    }
    
}

void PairwiseTraverse(list<Filament*>& mylist, list<Connection *> &conlist, void (*fun)(list<Connection*> &ConList, Filament *, Filament *,const double, const double, const double, const double, const double), const double R_cross, const double dt, const double Max_distance, const double Max_angle, const double Cross_angle)
{
    list<Filament*>::iterator first = mylist.begin();
    list<Filament*>::iterator last = mylist.end();
    list<Filament*>::iterator nextone;
    for(; first != last; ++first)
        for(nextone = std::next(first); nextone != last; ++nextone)
            fun(conlist, *first,*nextone, R_cross, dt, Max_distance, Max_angle, Cross_angle);
    //    mylist.push_front(6);
    //    list<int>::iterator first = mylist.begin();
    //    cout << *first << endl;
    //    fun(*first, 2);
}

void StericForce(list<Filament*>& mylist, const double RepulsiveDistance, const double RepulsiveForce)
{
    vec3 GetDistance(vec3 p1, vec3 p2, vec3 p3, vec3 p4);
    vec3 GetDistance(vec3 p1, vec3 p2, vec3 p3, vec3 p4, double *s, double *t);
    list<Filament*>::iterator first = mylist.begin();
    list<Filament*>::iterator last = mylist.end();
    list<Filament*>::iterator nextone;
    double s=0;
    double t=0;
    for(; first != last; ++first){
        for(nextone = std::next(first); nextone != last; ++nextone){
            
            vec3 p1 = (*first)->center_of_mass-0.5*(*first)->length*(*first)->orientation;
            vec3 p2 = (*first)->center_of_mass+0.5*(*first)->length*(*first)->orientation;
            vec3 p3 = (*nextone)->center_of_mass-0.5*(*nextone)->length*(*nextone)->orientation;
            vec3 p4 = (*nextone)->center_of_mass+0.5*(*nextone)->length*(*nextone)->orientation;
            vec3 dis = GetDistance(p1, p2, p3, p4,&s,&t);
            double len = norm(dis);
            if (len<RepulsiveDistance*pow(2.0,0.0)) {
              //  double RepulsiveForceNew = RepulsiveForce*RepulsiveDistance/6.0*(12.0/len*pow(RepulsiveDistance/len, 12)-6.0/len*pow(RepulsiveDistance/len, 6));
              //  RepulsiveForceNew = min(RepulsiveForceNew,RepulsiveForce*5);
                double RepulsiveForceNew = RepulsiveForce;
                (*first)->force = (*first)->force + RepulsiveForceNew*normalise(dis);
                (*nextone)->force = (*nextone)->force - RepulsiveForceNew*normalise(dis);
                vec3 q1 = p2 + s*(p1-p2);
                vec3 q2 = p4 + t*(p3-p4);
                (*first)->torque = (*first)->torque + cross(q1-(*first)->center_of_mass,RepulsiveForceNew*normalise(dis));
                (*nextone)->torque = (*nextone)->torque - cross(q2-(*nextone)->center_of_mass,RepulsiveForceNew*normalise(dis));
            }
        }
    }
            
}

void OutputGeometry(list<Filament*> &FList, ofstream& outfile1, double tt)
{
    list<Filament*>::iterator it_f;
    Filament *fil;
    for (it_f = FList.begin(); it_f != FList.end(); it_f++) {
        
        fil = *it_f;
        
        outfile1 << setw(20) << tt;
        
        outfile1 << setw(20) << fil->length;
        
        for (int j=0; j<3; j++)
            outfile1 << setw(20) << fil->center_of_mass[j];
        
        for (int j=0; j<3; j++)
            outfile1 << setw(20) << fil->orientation[j];
        
        for (int j=0; j<3; j++)
            outfile1 << setw(20) << fil->rotation[j];
        
        outfile1 << setw(20) << fil->isMobile;
        
        outfile1 << setw(20) << fil->index << endl;
    }
}

void OutputForceTorque(list<Filament*> &FList, ofstream& outfile1, double tt)
{
    list<Filament*>::iterator it_f;
    Filament *fil;
    for (it_f = FList.begin(); it_f != FList.end(); it_f++) {
        
        fil = *it_f;
        
        outfile1 << setw(20) << tt;
        
        for (int j=0; j<3; j++)
            outfile1 << setw(20) << fil->force[j];
        
        for (int j=0; j<3; j++)
            outfile1 << setw(20) << fil->torque[j];
        
        outfile1 << setw(20) << fil->isMobile;
        
        outfile1 << setw(20) << fil->index << endl;
    }
}


void OutputConnection(list<Connection*> &CList, ofstream& outfile2, double tt)
{
    list<Connection*>::iterator it_c;
    Connection *p;
    for (it_c = CList.begin(); it_c != CList.end(); it_c++) {
        p = *it_c;
        outfile2 << setw(20) << tt;
        outfile2 << setw(20) << p->f1->index << setw(20) << p->s1 << setw(20) << p->f2->index << setw(20) << p->s2 << setw(20) << p->transEnergy << setw(20) << p->torsEnergy << setw(20) << p->lifetime << endl;
     //   outfile2 << setw(20) << p->f1->index << setw(20) << p->s1 << setw(20) << p->f2->index << setw(20) << p->s2 << setw(20) << p->phi1_1 << setw(20) << p->theta1_1 << setw(20) << p->phi1_2 << setw(20) << p->theta1_2 << setw(20) << p->phi2_1 << setw(20) << p->theta2_1 << setw(20) << p->phi2_2 << setw(20) << p->theta2_2 << setw(20) << p->len << setw(20) << p->lifetime << endl;
    }
   
    
        // this is wrong
        //outfile2 << setw(20) << p->data.index1 << setw(20) << p->data.s1 << setw(20) << p->data.index2 << setw(20) << p->data.s2 << endl;
    
    
}

void Boundaryforce(list<Filament *> &FList, const string geometry, const double boxSize, const double height_c, const double radius_c, const double BBkappa)
{
    // The cube is [-.5*boxSize,.5*boxSize]*[-.5*boxSize,.5*boxSize]*[-boxSize,0]
    vec3 GetDistance(vec3 p1, vec3 p2, vec3 p3, vec3 p4);
    const double Bkappa = 10;  // boundary stiffness
    list<Filament *>::iterator it_f;
    if (geometry == "cube" || geometry == "cubecylinder") {
        
        for (it_f = FList.begin(); it_f!=FList.end(); it_f++){
            if ((*it_f)->isMobile) {
                Filament *p = *it_f;
                vec3 plus  = p->center_of_mass + 0.5*p->length*p->orientation;
                vec3 minus = p->center_of_mass - 0.5*p->length*p->orientation;
                vec3 x = {1.0, 0.0, 0.0};
                vec3 y = {0.0, 1.0, 0.0};
                vec3 z = {0.0, 0.0, 1.0};
                p->force = p->force + (double)(plus(0)>.5*boxSize)*Bkappa*(.5*boxSize-plus(0))*x     \
                + (double)(plus(0)<-.5*boxSize)*Bkappa*(-.5*boxSize-plus(0))*x   \
                + (double)(plus(1)>.5*boxSize)*Bkappa*(.5*boxSize-plus(1))*y     \
                + (double)(plus(1)<-.5*boxSize)*Bkappa*(-.5*boxSize-plus(1))*y   \
                + (double)(plus(2)>0)*Bkappa*(0.0-plus(2))*z     \
                + (double)(plus(2)<-boxSize)*Bkappa*(-boxSize-plus(2))*z;
                
                p->force = p->force + (double)(minus(0)>.5*boxSize)*Bkappa*(.5*boxSize-minus(0))*x     \
                + (double)(minus(0)<-.5*boxSize)*Bkappa*(-.5*boxSize-minus(0))*x   \
                + (double)(minus(1)>.5*boxSize)*Bkappa*(.5*boxSize-minus(1))*y     \
                + (double)(minus(1)<-.5*boxSize)*Bkappa*(-.5*boxSize-minus(1))*y   \
                + (double)(minus(2)>0)*Bkappa*(0.0-minus(2))*z     \
                + (double)(minus(2)<-boxSize)*Bkappa*(-boxSize-minus(2))*z;
            }
            
        }
    }
    if (geometry == "cylinder" || geometry == "cubecylinder") {
        vec p3 = {0.0,0.0,0.0};
        vec p4 = {0.0,0.0,-height_c};
        
        for (it_f = FList.begin(); it_f!=FList.end(); it_f++){
            if ((*it_f)->isMobile) {
                Filament *p = *it_f;
                vec3 p1  = p->center_of_mass + 0.5*p->length*p->orientation;
                vec3 p2 = p->center_of_mass - 0.5*p->length*p->orientation;
                vec3 dis = GetDistance(p1, p2, p3, p4); // pointing from the central cylinder to the filament
                double len = norm(dis);
                if (len<radius_c) {
                    p->force = p->force + BBkappa*(radius_c-len)*normalise(dis);
                }

            }
        }

    }
    
}

vec3 GetDistance(vec3 p1, vec3 p2, vec3 p3, vec3 p4)
{   // vec returns a vector which is the nereast distance between segment 1 and segment 2
    // segment 1 is given by its end points p1 and p2
    // segment 2 is given by its end points p3 and p3
    // vec points from 2 to 1
    vec3 u = p1 - p2;
    vec3 v = p3 - p4;
    vec3 w = p2 - p4;
    
    double a = dot(u, u);
    double b = dot(u, v);
    double c = dot(v, v);
    double d = dot(u, w);
    double e = dot(v, w);
    double D = a*c-b*b;
    double sD = D;
    double tD = D;
    
    double SMALL_NUM=0.00000001;
    
    double sN,tN;
    if (D < SMALL_NUM)
    {
        sN=0.0;sD=1.0;tN=e;tD=c;
    }
    else
    {
        sN=b*e-c*d;
        tN=a*e-b*d;
        if(sN < 0.0)
        {
            sN=0.0;tN=e;tD=c;
        }
        else if (sN>sD)
        {
            sN=sD;tN=e+b;tD=c;
        }
    }
    
    if (tN<0.0)
    {
        tN = 0.0;
        if (-d < 0.0)
            sN=0.0;
        else if (-d>a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
    {
        tN = tD;
        if (-d+b < 0.0)
            sN=0;
        else if (-d+b>a)
            sN = sD;
        else
        {
            sN = -d+b;
            sD=a;
        }
    }
    
    double sc, tc;
    if (fabs(sN) < SMALL_NUM)
        sc= 0.0;
    else
        sc=sN/sD;
    
    if (fabs(tN)< SMALL_NUM)
        tc = 0.0;
    else
        tc = tN/tD;
    
    return w + sc*u - tc*v;
    
}


vec3 GetDistance(vec3 p1, vec3 p2, vec3 p3, vec3 p4,double *sc, double *tc)
{   // vec returns a vector which is the nereast distance between segment 1 and segment 2
    // segment 1 is given by its end points p1 and p2
    // segment 2 is given by its end points p3 and p3
    // vec points from 2 to 1
    vec3 u = p1 - p2;
    vec3 v = p3 - p4;
    vec3 w = p2 - p4;
    
    double a = dot(u, u);
    double b = dot(u, v);
    double c = dot(v, v);
    double d = dot(u, w);
    double e = dot(v, w);
    double D = a*c-b*b;
    double sD = D;
    double tD = D;
    
    double SMALL_NUM=0.00000001;
    
    double sN,tN;
    if (D < SMALL_NUM)
    {
        sN=0.0;sD=1.0;tN=e;tD=c;
    }
    else
    {
        sN=b*e-c*d;
        tN=a*e-b*d;
        if(sN < 0.0)
        {
            sN=0.0;tN=e;tD=c;
        }
        else if (sN>sD)
        {
            sN=sD;tN=e+b;tD=c;
        }
    }
    
    if (tN<0.0)
    {
        tN = 0.0;
        if (-d < 0.0)
            sN=0.0;
        else if (-d>a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
    {
        tN = tD;
        if (-d+b < 0.0)
            sN=0;
        else if (-d+b>a)
            sN = sD;
        else
        {
            sN = -d+b;
            sD=a;
        }
    }
    
   // double sc, tc;
    if (fabs(sN) < SMALL_NUM)
       *sc = 0.0;
    else
       *sc = sN/sD;
    
    if (fabs(tN)< SMALL_NUM)
        *tc = 0.0;
    else
        *tc = tN/tD;
    
    return w + (*sc)*u - (*tc)*v;
    
}




