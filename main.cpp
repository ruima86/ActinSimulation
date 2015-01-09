//
//  main.cpp
//  Filament
//
//  Created by Rui Ma on 1/5/15.
//  Copyright (c) 2015 Rui Ma. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "filament.h"
#include "BasicLinearAlgebra.h"
#include "mvnrnd.h"
#include "bilistconnection.h"
#define EPS 1e-10

using namespace std;

int main(int argc, const char * argv[])
{
    void Get_Distance(Filament f1, Filament f2, double vec[], double *s, double *t);
    void Node_Insert(BiList p, BiList q);
    void Node_Delete(BiList q);
    void Node_Print(BiList p);
    
    /********************************************************/
    /******************%%%% Parameters **********************/
    /********************************************************/
    const int N_fil=100;
    const int L_fil=30;
    int ensemble=1;
    
    //kinetic rates
    const double R_cross=0.01;
    const double R_decross=2.5e-3;
    
    //mechanical properties of the crosslinker
    const double kappa=1.0;       // pN/nm
    const double fc=3.0;          // pN
    // Simulation
    const double Time=100;      // s       simulation time
    const double dt=1.0e-3;       // s       simulation time step
    const int Nstep=floor(Time/dt/100);
    // Geometry
    const double radius=25.0;    // nm      radius of the cylinder
    const double height=150.0;   // nm      height of the radius
    const double DeltaM=2.7;     // nm      size of actin monomer
    const double Dradius=5.0;    // nm      radius of filament cross-section
    const double cric_radius=7.0;// nm      crictical length of the crosslinker
    
    // Thermodynamics
    const double kT=4.1;         // pN*nm   thermal energy
    const double viscosity=0.00001; // pN*s/nm^2
    
    /********************************************************/
    /*{
     char str[10];
     
     //Creates an instance of ofstream, and opens example.txt
     ofstream a_file ( "example.txt" );
     // Outputs to example.txt through a_file
     a_file<<"This text will now be inside of example.txt";
     // Close the file stream explicitly
     a_file.close();
     //Opens for reading the file
     ifstream b_file ( "example.txt" );
     //Reads one string from the file
     b_file>> str;
     //Should output 'this'
     cout<< str <<"\n";
     cin.get();    // wait for a keypress
     // b_file is closed implicitly here

    }*/
    /********************************************************/
    time_t ti = time(0);
    tm *now = localtime(&ti);
    ostringstream ostr1, ostr2;
    ostr1 << "Trajectory_" << ensemble << "_crossrate_" << R_cross << "_Ndendral_" << N_fil << "_Length_" << L_fil << "_radius_" << radius << "_date_" << now->tm_mon+1 << "-" << now->tm_mday << "-" << now->tm_year+1900 << ".txt";
    string ss = ostr1.str();
    ofstream outfile1(ss);
    
    ostr2 << "Connection_" << ensemble << "_crossrate_" << R_cross << "_Ndendral_" << N_fil << "_Length_" << L_fil << "_radius_" << radius << "_date_" << now->tm_mon+1 << "-" << now->tm_mday << "-" << now->tm_year+1900 << ".txt";
    ss = ostr2.str();
    ofstream outfile2(ss);
    
    outfile1 << setprecision(4);
    outfile1 << setiosflags(ios::fixed);
    outfile1 << setiosflags(ios::left);
    
    outfile2 << setprecision(4);
    outfile2 << setiosflags(ios::fixed);
    outfile2 << setiosflags(ios::left);
    
    /********************************************************/
    /****************** Initialization **********************/
    /********************************************************/
    
    // Initilize random number generator
    const gsl_rng_type *T;
    const gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    // Filaments are equally spaced along the circumrence
    // of the cylinder, randomly orientated, plus end orienting
    // towards the flat surface
    Filament fil[N_fil];
    for(int i=0;i<N_fil;i++)
    {   // center_of_mass
        double RR=radius+5.0;
        VecSet(fil[i].center_of_mass,RR*cos(i*2.0*M_PI/N_fil),RR*sin(i*2.0*M_PI/N_fil),-0.5*L_fil*DeltaM);
        // orientation
        double phi=2.0*M_PI*gsl_rng_uniform(r);
        double theta=M_PI_2*gsl_rng_uniform(r);
        VecSet(fil[i].orientation,sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        // length
        fil[i].length=L_fil*DeltaM;
        // connection
        fil[i].N_Connect_Left=0;
        fil[i].N_Connect_Right=0;
        fil[i].index=i+1;
    }
    
    /********************************************************/
    /************ friction coefficient **********************/
    /********************************************************/
    double xi_para,xi_perp,xi_rpara,xi_rperp;
    {
        double len = L_fil*DeltaM;
        double pp = len/2.0/Dradius;
        double f0 = 6.0*M_PI*viscosity*len*cbrt(3.0/16.0*pp*pp);
        xi_para=f0*(1.009+1.395*1.0e-2*log(pp)+7.880*1e-2*pow(log(pp), 2)+6.040*1e-3*pow(log(pp),3));
        xi_perp=2.0*xi_para;
        double tau0=M_PI*pow(len,3)*viscosity/4.0/pow(pp,2);
        double taua;
        if(pp<0.75)
            taua=tau0*(1.18+0.1744*pow(log(pp)+0.2877,2)-0.2417*pow(log(pp)+0.2877,3)-3.882*1e-2*pow(log(pp)+0.2877,4));
        else
            taua=tau0*(1.18+1.1160*pow(log(pp)+0.2877,2)-0.9729*pow(log(pp)+0.2877,3)+0.4954*pow(log(pp)+0.2877,4));
        double taub=tau0*(1.183+0.2902*log(pp)+0.4406*pow(log(pp),2)-5.850*1e-2*pow(log(pp),3)-9.544*1e-3*pow(log(pp),4));
        xi_rperp=6*taua;
        xi_rpara=1.0/(1.0/taub-5.0/xi_rperp);
    }
    
    BiList head = new node; // the head of the connection list
    
    //dynamics starts from here
    long counter=0;
    for (double tt=0; tt<Time; tt+=dt)
    {
        /********************************************************/
        /***************** stochastic force *********************/
        /********************************************************/
        for(int i=0;i<N_fil;i++)
        {
            double Nmatrix[3*3];
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++)
                {
                    Nmatrix[j*3+k] = fil[i].orientation[j]*fil[i].orientation[k];
                    fil[i].fric_f[j*3+k] = xi_para*Nmatrix[j*3+k]+xi_perp*((j==k)-Nmatrix[j*3+k]);
                    fil[i].fric_T[j*3+k] = xi_rpara*Nmatrix[j*3+k]+xi_rperp*((j==k)-Nmatrix[j*3+k]);            }
            //MatrixPrint(m,6,6);
            double m[3*3];
            VecScalarProd(2.0*kT/dt, fil[i].fric_f, 3*3,m);
            mvnrnd(m,3,r,fil[i].force);
            VecScalarProd(2.0*kT/dt, fil[i].fric_T, 3*3,m);
            mvnrnd(m,3,r,fil[i].torque);
        }
        
        /********************************************************/
        /***************** deterministic force ******************/
        /********************************************************/
        
        // break the connection, otherwise build the force
        { // go through each connection
            BiList q = new node;
            q=head->next;
            BiList p = new node;
            double q1[3],q2[3],vec[3];
            double len;
            while (q!=0)
            {
                p=q->next;
                VecSumScalarProd(q->s1*q->f1->length, q->f1->orientation, q->f1->center_of_mass, q1);
                VecSumScalarProd(q->s2*q->f2->length, q->f2->orientation, q->f2->center_of_mass, q2);
                VecSumScalarProd(-1, q2, q1, vec);
                len=norm2(vec);
                if(len>2*cric_radius || gsl_rng_uniform(r)<R_decross*exp(kappa*len/fc)*dt)
                {//break the connection
                    int n1 = ceil((q->s1+0.5+EPS)*(q->f1->length)/DeltaM);
                    int n2 = ceil((q->s2+0.5+EPS)*(q->f2->length)/DeltaM);
                    if(n1>0) // this is always true without filament shrinkage
                        q->f1->Occupied[n1] = 0;
                    if(n2>0) // this is always true without filament shrinkage
                        q->f2->Occupied[n2] = 0;
                    Node_Delete(q);
                }
                else
                {//build the force
                    double f12[3]={0};
                    double T12[3]={0};
                    double dis[3]={0};
                    // for filament 1
                    VecScalarProd(-kappa, vec, 3, f12);
                    VecSumScalarProd(1, f12, q->f1->force, q->f1->force);
                    VecSumScalarProd(-1, q->f1->center_of_mass, q1, dis);
                    VecCross(dis, f12, T12);
                    VecSumScalarProd(1, T12, q->f1->torque, q->f1->torque);
                    // for filament 2
                    VecScalarProd(kappa, vec, 3, f12);
                    VecSumScalarProd(1, f12, q->f2->force, q->f2->force);
                    VecSumScalarProd(-1, q->f2->center_of_mass, q2, dis);
                    VecCross(dis, f12, T12);
                    VecSumScalarProd(1, T12, q->f2->torque, q->f2->torque);
                }
                q=p;
            }
        }
        
        /********************************************************/
        /********************* boundary force *******************/
        /********************************************************/
        for(int i=0;i<N_fil;i++)
        {
            //upper boundary
            double plus[3],minus[3];
            VecSumScalarProd(0.5*fil[i].length,fil[i].orientation, fil[i].center_of_mass, plus);
            VecSumScalarProd(-0.5*fil[i].length, fil[i].orientation, fil[i].center_of_mass, minus);
            double zdir[3]={0.0,0.0,1.0};
            if(plus[2]>minus[2] && plus[2]>0)
            {
                VecSumScalarProd(-plus[2], zdir, fil[i].force, fil[i].force);
                double temp1[3],temp2[3];
                VecSumScalarProd(-1, fil[i].center_of_mass, plus, temp1);
                VecCross(temp1, zdir, temp2);
                VecSumScalarProd(-plus[2], temp2, fil[i].torque, fil[i].torque);
            }
            
            if(plus[2]<minus[2] && minus[2]>0)
            {
                VecSumScalarProd(-minus[2], zdir, fil[i].force, fil[i].force);
                double temp1[3],temp2[3];
                VecSumScalarProd(-1, fil[i].center_of_mass, minus, temp1);
                VecCross(temp1, zdir, temp2);
                VecSumScalarProd(-minus[2], temp2, fil[i].torque, fil[i].torque);
            }
            
            //cylinder
            Filament fcyl;
            VecSet(fcyl.center_of_mass, 0.0, 0.0, -height/2);
            VecSet(fcyl.orientation, 0.0, 0.0, 1.0);
            fcyl.length=height;
            double vec[3];
            double s,t;
            Get_Distance(fcyl, fil[i], vec , &s , &t);
            double len=norm2(vec);
            if(len<radius)
            {
                VecSumScalarProd(-10/len, vec , fil[i].force, fil[i].force);
                double pos[3];
                VecSumScalarProd(t*fil[i].length, fil[i].orientation, fil[i].center_of_mass, pos);
                double temp1[3],temp2[3];
                VecSumScalarProd(-1, fil[i].center_of_mass, pos, temp1);
                VecCross(temp1, vec , temp2);
                VecSumScalarProd(-10/len, temp2, fil[i].torque, fil[i].torque);
            }
            
            
        }
        /********************************************************/
        /******************** position update  ******************/
        /********************************************************/
        for(int i=0;i<N_fil;i++)
        {
            double Vc[3];
            double Omega[3];
            invMB(fil[i].fric_f, fil[i].force, 3, Vc);
            invMB(fil[i].fric_T, fil[i].torque, 3, Omega);
            double amp=norm2(Omega);
            double u[3];
            VecScalarProd(1.0/amp, Omega, 3, u);
            double ucross[3][3] = {{0,-u[2],u[1]},{u[2],0,-u[0]},{-u[1],u[0],0}};
            double ucrossu[3][3]= {{u[0]*u[0],u[0]*u[1],u[0]*u[2]},{u[1]*u[0],u[1]*u[1],u[1]*u[2]},{u[2]*u[0],u[2]*u[1],u[2]*u[2]}};
            double Rot[3][3];
            for(int i=0;i<3;i++)
                for (int j=0; j<3; j++) {
                    Rot[i][j]=cos(amp*dt)*(i==j)+sin(amp*dt)*ucross[i][j]+(1-cos(amp*dt))*ucrossu[i][j];
                }
            VecSumScalarProd(dt, Vc, fil[i].center_of_mass, fil[i].center_of_mass);
            double temp[3];
            MatrixVecProd(Rot, fil[i].orientation, temp);
            VecSumScalarProd(0, temp, temp, fil[i].orientation);
        }
        
        /********************************************************/
        /******************** Connection build  *****************/
        /********************************************************/
        for(int i=1;i<N_fil;i++)
            for (int j=0; j<i; j++)
            {
                double dis[3];
                VecSumScalarProd(-1, fil[j].center_of_mass, fil[i].center_of_mass, dis);
                if (norm2(dis)<0.5*(fil[j].length+fil[i].length))
                {
                    long nmax1=round(fil[i].length/DeltaM);
                    long n1=randi(r , nmax1);
                    long nmax2=round(fil[j].length/DeltaM);
                    long n2=randi(r , nmax2);
                    if(fil[i].Occupied[n1]==0 && fil[j].Occupied[n2]==0)
                    {
                        double P1[3],P2[3],vec[3];
                        double s1,s2;
                        s1=(n1-0.5)*DeltaM/fil[i].length-0.5;
                        VecSumScalarProd(s1*fil[i].length, fil[i].orientation, fil[i].center_of_mass, P1);
                        s2=(n2-0.5)*DeltaM/fil[j].length-0.5;
                        VecSumScalarProd(s2*fil[j].length, fil[j].orientation, fil[j].center_of_mass, P2);
                        VecSumScalarProd(-1, P2, P1, vec );
                        if (gsl_rng_uniform(r)<nmax1*nmax2*R_cross*dt && norm2(vec)<cric_radius)
                        {
                            BiList x = new node;
                            x->index1 = i+1;
                            x->index2 = j+1;
                            x->s1 = s1;
                            x->s2 = s2;
                            x->f1 = &fil[i];
                            x->f2 = &fil[j];
                            fil[i].Occupied[n1]=1;
                            fil[j].Occupied[n2]=1;
                            Node_Insert(head,x);
                           // VecPrint(fil[i].center_of_mass, 3);
                           // VecPrint(x->f1->center_of_mass, 3);
                           // fil[i].center_of_mass[0]=1.0;
                           // VecPrint(fil[i].center_of_mass, 3);
                           // VecPrint(x->f1->center_of_mass, 3);
                           // x=head->next;
                           // double vec[3];
                           // VecSumScalarProd(-1,x->f1.center_of_mass,fil[x->index1-1].center_of_mass,vec);
                           // cout << norm2(vec) << endl;
                        }
                    }
                }
            }
        
        /********************************************************/
        /************************* Output  **********************/
        /********************************************************/
        if ((counter % Nstep) == 0)
        {
            //output geometry
            for (int i=0; i<N_fil; i++)
            {
                outfile1 << setw(10) << tt;
                
                outfile1 << setw(10) << fil[i].length;
                
                for (int j=0; j<3; j++)
                    outfile1 << setw(10) << fil[i].center_of_mass[j];
                
                for (int j=0; j<3; j++)
                    outfile1 << setw(10) << fil[i].orientation[j];
                
                outfile1 << i+1 << endl;
                
            }
            //output connection
            BiList p = new node;
            p = head->next;
            while (p!=0) {
               // double p1[3],p2[3];
               // double vec[3];
               // cout << p->f1->index-fil[p->index1-1].index << endl;
               // VecSumScalarProd(p->s1*fil[p->index1-1].length, fil[p->index1-1].orientation, fil[p->index1-1].center_of_mass, p1);
               // VecSumScalarProd(p->s2*fil[p->index2-1].length, fil[p->index2-1].orientation, fil[p->index2-1].center_of_mass, p2);
               // VecSumScalarProd(-1, p2, p1, vec);
                outfile2 << setw(10) << tt;
                outfile2 << setw(10) << p->index1 << setw(10) << p->s1 << setw(10) << p->index2 << setw(10) << p->s2 << endl;
                p=p->next;
            }
            
            cout << "t=" << tt << endl;
        }
        counter+=1;
        
    }
    
    outfile1.close();
    outfile2.close();
}

void Get_Distance(Filament f1, Filament f2, double vec[], double *s, double *t)
{   // vec returns a vector which is the nereast distance between filament 1 and filament 2
    // vec points from
    double p1[3],p2[3],p3[3],p4[3];
    VecSumScalarProd(f1.length/2, f1.orientation, f1.center_of_mass, p1);
    VecSumScalarProd(-f1.length/2, f1.orientation, f1.center_of_mass, p2);
    VecSumScalarProd(f2.length/2, f2.orientation, f2.center_of_mass, p3);
    VecSumScalarProd(-f2.length/2, f2.orientation, f2.center_of_mass, p4);
    
    double u[3],v[3],w[3];
    VecSumScalarProd(-1, p2, p1, u);
    VecSumScalarProd(-1, p4, p3, v);
    VecSumScalarProd(-1, p4, p2, w);
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
    
    double dP[3];
    VecSumScalarProd(sc, u, w, dP);
    VecSumScalarProd(-tc, v, dP, dP);
    VecSet(vec, dP[0], dP[1], dP[2]);
    *s=sc-0.5;
    *t=tc-0.5;
}

    