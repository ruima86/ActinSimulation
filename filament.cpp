//
//  filament.cpp
//  GSLTest
//
//  Created by Rui Ma on 2/3/17.
//  Copyright © 2017 Yale University. All rights reserved.
//

#include "filament.hpp"

// constructor
Filament::Filament()
{
    
    center_of_mass = 100*randu(3)-50;      // random position
    double theta=acos(1.0-2.0*randu());
    double phi=2*M_PI*randu();
    vec3 v1 = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
    vec3 v2 = {-sin(phi),cos(phi),0.0};
    orientation = v1;                      // random orientation
    rotation = v2;                         // perpendicular to orientation
    
    // these values are fixed
    delta = 2.7;
    radius = 3.0;
    pitchnumber = 13.0;
    
    length = 10*delta; // minimum length is 10 monomers long
    
    force = zeros(3);
    torque = zeros(3);
    
    xi_t = zeros(3,3);
    xi_r = zeros(3,3);
    xi_tr = zeros(3,3);
    
    for (int i=0; i<200; i++) {
        Branched[i]=0;
    }
    
    //connection
    
    //state
    BarbedCapping=0;
    PointedCapping=0;
    isMobile=true;
    
    // index
    index = 0;
}

Filament::Filament(int NumMonomers,const double InitialSize, const int InitialAlignment)
{
    vec3 r1 = {.5,.5,1.0};
    center_of_mass = InitialSize*(randu(3)-r1);      // random position
    double theta=acos(1.0-2.0*randu());
    double phi=2*M_PI*randu();
    if (InitialAlignment == 1) {
        vec3 v1 = {0.0, 0.0, 1.0};
        orientation = v1;
    }
    else {
        vec3 v1 = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
        orientation = v1;
    }
    
    vec3 v2 = {-sin(phi),cos(phi),0.0};
    // random orientation
    rotation = v2;                         // perpendicular to orientation
    
    // these values are fixed
    delta = 2.7;
    radius = 3.0;
    pitchnumber = 13;
    
    if (NumMonomers>10) {
        length = NumMonomers*delta;
    }
    else
    {
        length = 10*delta; // minimum length is 10 monomers long
        
    }
    
    force = zeros(3);
    torque = zeros(3);
    
    xi_t = zeros(3,3);
    xi_r = zeros(3,3);
    xi_tr = zeros(3,3);
    
    for (int i=0; i<200; i++) {
        Branched[i]=0;
    }
    
    //connection
    
    //state
    BarbedCapping=0;
    PointedCapping=0;
    isMobile = true;
    
    // index
    index = 0;
}



Filament::Filament(int NumMonomers,const double InitialSize, const int InitialAlignment, const double pitch)
{
    vec3 r1 = {.5,.5,1.0}; 
    center_of_mass = InitialSize*(randu(3)-r1);      // random position
    double theta=acos(1.0-2.0*randu());
    double phi=2*M_PI*randu();
    if (InitialAlignment == 1) {
        vec3 v1 = {0.0, 0.0, 1.0};
        orientation = v1;
    }
    else {
        vec3 v1 = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
        orientation = v1;   
    }
   
    vec3 v2 = {-sin(phi),cos(phi),0.0};
                       // random orientation
    rotation = v2;                         // perpendicular to orientation
    
    // these values are fixed
    delta = 2.7;
    radius = 3.0;
    pitchnumber = pitch;
    
    if (NumMonomers>10) {
        length = NumMonomers*delta;
    }
    else
    {
        length = 10*delta; // minimum length is 10 monomers long

    }
    
    force = zeros(3);
    torque = zeros(3);
    
    xi_t = zeros(3,3);
    xi_r = zeros(3,3);
    xi_tr = zeros(3,3);
    
    for (int i=0; i<200; i++) {
        Branched[i]=0;
    }
    
    //connection
    
    //state
    BarbedCapping=0;
    PointedCapping=0;
    isMobile = true;
    
    // index
    index = 0;
}

Filament::Filament(int NumMonomers, const double cylinderradius)
{
    
    double theta=acos(1.0-2.0*randu());
    double phi=2*M_PI*randu();
    vec3 v1 = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
    orientation = v1;
    
    
    vec3 v2 = {-sin(phi),cos(phi),0.0};
    // random orientation
    rotation = v2;                         // perpendicular to orientation
    
    // these values are fixed
    delta = 2.7;
    radius = 3.0;
    pitchnumber = 13.0;
    
    if (NumMonomers>10) {
        length = NumMonomers*delta;
    }
    else
    {
        length = 10*delta; // minimum length is 10 monomers long
        
    }
    
    double r = cylinderradius;
    phi = 2*M_PI*randu();
    vec3 cen = {r*cos(phi),r*sin(phi),-length/2.0};
    center_of_mass = cen;
    
    force = zeros(3);
    torque = zeros(3);
    
    xi_t = zeros(3,3);
    xi_r = zeros(3,3);
    xi_tr = zeros(3,3);
    
    for (int i=0; i<200; i++) {
        Branched[i]=0;
    }
    
    //connection
    
    //state
    BarbedCapping=0;
    PointedCapping=0;
    isMobile = true;
    
    // index
    index = 0;

}


Filament::Filament(vec3 com)        // constructor with center_of_mass as parameter
{
    
    center_of_mass = com;
    double theta=acos(1.0-2.0*randu());
    double phi=2*M_PI*randu();
    vec3 v1 = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
    vec3 v2 = {-sin(phi),cos(phi),0.0};
    orientation = v1;               // random orientation
    rotation = v2;                  // perpendicular to orientation
    
    // these values are fixed
    delta = 2.7;
    radius = 3.0;
    pitchnumber = 13.0;
    
    length = 10*delta;             // minimum length is 10 monomers long
    
    force = zeros(3);
    torque = zeros(3);
    
    xi_t = zeros(3,3);
    xi_r = zeros(3,3);
    xi_tr = zeros(3,3);
    
    for (int i=0; i<200; i++) {
        Branched[i]=0;
    }
    
    //connection
    
    //state
    BarbedCapping=0;
    PointedCapping=0;
    isMobile = true;
    
    // index
    index = 0;
}

Filament::Filament(vec3 com, vec3 ort, double len)        // overload constructor with center_of_mass as parameter
{
    
    center_of_mass = com;
    orientation = ort;               // random orientation
    double theta=acos(1.0-2.0*randu());
    double phi=2*M_PI*randu();
    vec3 v1 = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)}; // generate a random orientation
    rotation = normalise(cross(ort,v1));                  // perpendicular to orientation
    
    // these values are fixed
    delta = 2.7;
    radius = 3.0;
    pitchnumber = 13.0;
    
    length = len;        
    
    force = zeros(3);
    torque = zeros(3);
    
    xi_t = zeros(3,3);
    xi_r = zeros(3,3);
    xi_tr = zeros(3,3);
    
    for (int i=0; i<200; i++) {
        Branched[i]=0;
    }
    
    //connection
    
    //state
    BarbedCapping=0;
    PointedCapping=0;
    isMobile = true;
    
    // index
    index = 0;
    
}


// destructor
Filament::~Filament()
{
   // cout << "Filament's destructor called" << endl;
    OccupiedList.clear();
}


/***********************************************************/
/**************** Operator overload  ***********************/
/***********************************************************/

bool operator == (Filament &f1, Filament &f2)
{
    if (norm(f1.center_of_mass - f2.center_of_mass)<1e-10 && norm(f1.orientation-f2.orientation)<1e-10) {
        return true;
    }
    else
    {
        return false;
    }
    
}

bool operator != (Filament &f1, Filament &f2)
{
    return !(f1==f2);
}

bool operator < (Filament &f1, Filament &f2)
{
    if (f1.length<f2.length) {
        return true;
    }
    else
    {
        return false;
    }
}


void Filament::GetFricMatrix(double xi_t_para, double xi_t_perp, double xi_r_para, double xi_r_perp)
{   // Get friction matrix calculated relative to the center of mass
    
    vec3 n = this->orientation;
    mat33 nn = n*n.t();
    xi_t = xi_t_para*nn+xi_t_perp*(eye(3, 3)-nn);
    xi_r = xi_r_para*nn+xi_r_perp*(eye(3, 3)-nn);
    xi_tr = zeros(3,3);
}

mat33 Filament::GetCouplingMatrix(vec3 ro)
{
    // Transform the coupling matrix from the center-of-mass frame to a new frame specified by ro
    vec3 rc = this->center_of_mass;
    vec3 dr = ro - rc;
    mat33 sigma = {{0, -dr(2), dr(1)},{dr(2),0,-dr(0)},{-dr(1),dr(0),0}};
    return this->xi_tr+xi_t*sigma;
}

mat33 Filament::GetRotationMatrix(vec3 ro)
{
    // Transform the Rotation matrix from the center-of-mass frame to a new frame specified by ro
    vec rc = this->center_of_mass;
    vec3 dr = ro - rc;
    mat33 sigma = {{0, -dr(2), dr(1)},{dr(2),0,-dr(0)},{-dr(1),dr(0),0}};
    return this->xi_r-sigma*this->xi_t*sigma+this->xi_tr.t()*sigma-sigma*xi_tr;
}

vec3 Filament::GetRotationMonomer(int n)
{
    vec3 v = this->orientation;
    // calculate the oritation of the binding surface of the crosslink u
    mat33 vxv=v*v.t();
    mat33 vx={{0,-v(2),v(1)},{v(2),0,-v(0)},{-v(1),v(0),0}};
    // the rotation angle
    double Alpha = M_PI*(1+1.0/this->pitchnumber)*(n-1);
    mat33 Rot;
    Rot = cos(Alpha)*eye(3, 3)+sin(Alpha)*vx+(1-cos(Alpha))*vxv;
    return Rot*this->rotation;
}

vec3 Filament::GetPositionMonomer(int n)
{
    vec3 v = this->orientation;
    vec3 r = this->center_of_mass;
    double s = (n-0.5)*this->delta/this->length-0.5;
    return r+s*this->length*v;
}

bool Filament::IsOccupiedMomomer(int n)
{
    if (std::find(OccupiedList.begin(), OccupiedList.end(), n)!=OccupiedList.end()) {
        return true;
    }
    else
    {
        return false;
    }
    
}


void Filament::Move(const double kT, const double dt)
{
    // Generate Brownian forces
    mat M1 = join_rows(xi_t, xi_tr);
    mat M2 = join_rows(xi_tr.t(), xi_r);
    mat66 Fric = join_cols(M1,M2);
    vec6 v; // v are stochastic force and torque
    if (kT<1e-10) {
        v.zeros();
    }
    else
    {
        mat66 R = chol(2*kT/dt*Fric); // R = chol(X), R is upper triangular, such that R.t()*R = X
        v = R.t()*arma::randn(6);    // v*v.t() = R.t()*R = Fric
    }
   
    vec6 w = join_cols(force, torque);  // w are elastic force and torque
    vec6 sol = arma::solve(Fric, v+w);
    vec3 vc = sol.head(3);
    vec3 omegac = sol.tail(3);
    
    
    mat33 Rot;
    
    double amp = norm(omegac);
    if (amp<1e-10) {
        Rot = eye(3, 3);
    }
    else
    {
        vec3 u = normalise(omegac);
        mat33 ux = {{0,-u[2],u[1]},{u[2],0,-u[0]},{-u[1],u[0],0}};
        mat33 uxu = u*u.t();
        Rot = cos(amp*dt)*eye(3, 3)+sin(amp*dt)*ux + (1-cos(amp*dt))*uxu;
    }
    
    center_of_mass = center_of_mass + vc*dt;
    orientation = Rot*orientation;
    rotation = Rot*rotation;
}

/***********************************************************/
/**************** Binary Tree Threading  *******************/
/***********************************************************/



/***********************************************************/
/******************** Structure change  ********************/
/***********************************************************/



/***********************************************************/
/***************** Elementary processes  *******************/
/***********************************************************/




//void Filament::GetFricMatrix(double viscosity, double DCenterOfMass[])
//{
//    VecZero(this->fric_Matrix, 6*6);
//    double xi_para,xi_perp,xi_rpara,xi_rperp;
//    {// results from Journal of chemical physics vol 119 (2003)
//        // hydrodynamic properties of rodlike and disklike particles in dilute solution
//        double len = this->length;
//        double pp = len/2.0/this->radius;
//        double f0 = 6.0*M_PI*viscosity*len*cbrt(3.0/16.0/pp/pp);
//        xi_para=f0*(1.009+1.395*1.0e-2*log(pp)+7.880*1e-2*pow(log(pp), 2)+6.040*1e-3*pow(log(pp),3));
//        xi_perp=2.0*xi_para;
//        double tau0=M_PI*pow(len,3)*viscosity/4.0/pow(pp,2);
//        double taua;
//        if(pp<0.75)
//            taua=tau0*(1.18+0.1744*pow(log(pp)+0.2877,2)-0.2417*pow(log(pp)+0.2877,3)-3.882*1e-2*pow(log(pp)+0.2877,4));
//        else
//            taua=tau0*(1.18+1.1160*pow(log(pp)+0.2877,2)-0.9729*pow(log(pp)+0.2877,3)+0.4954*pow(log(pp)+0.2877,4));
//        double taub=tau0*(1.183+0.2902*log(pp)+0.4406*pow(log(pp),2)-5.850*1e-2*pow(log(pp),3)-9.544*1e-3*pow(log(pp),4));
//        xi_rperp=6*taua;
//        xi_rpara=1.0/(1.0/taub-5.0/xi_rperp);
//    }
//    
//    double A[3*3];
//    double B[3*3];
//    double N[3*3];
//    
//    for(int j=0;j<3;j++)
//        for(int k=0;k<3;k++)
//        {
//            N[j*3+k] = this->orientation[j]*this->orientation[k];
//            A[j*3+k] = xi_para*N[j*3+k]+xi_perp*((j==k)-N[j*3+k]);  // anisotropic
//            //   A[j*3+k] = xi_para*N[j*3+k]+xi_para*((j==k)-N[j*3+k]);  // isotropic
//            //   B[j*3+k] = xi_rpara*N[j*3+k]+xi_rperp*((j==k)-N[j*3+k]); // anisotropic
//            B[j*3+k] = xi_rpara*N[j*3+k]+xi_rpara*((j==k)-N[j*3+k]); // isotropic
//        }
//    
//    double R[3];
//    VecSumScalarProd(-1, DCenterOfMass, this->center_of_mass, R);
//    double C[3*3]={0,-R[2],R[1],R[2],0,-R[0],-R[1],R[0],0};
//    double AC[3*3];
//    double CAC[3*3];
//    MatrixMatrixProd(A, C, 3, AC);
//    MatrixMatrixProd(C, AC, 3, CAC);
//    for (int j=0; j<3; j++)
//    {
//        for (int k=0; k<3; k++)
//        {
//            this->fric_Matrix[j*6+k]=this->fric_Matrix[j*6+k]+A[j*3+k];
//            this->fric_Matrix[j*6+k+3]=this->fric_Matrix[j*6+k+3]-AC[j*3+k];
//            this->fric_Matrix[(j+3)*6+k]=this->fric_Matrix[(j+3)*6+k]-AC[k*3+j];
//            this->fric_Matrix[(j+3)*6+k+3]=this->fric_Matrix[(j+3)*6+k+3]+B[j*3+k]-CAC[j*3+k];
//        }
//    }
//    
//    
//}


