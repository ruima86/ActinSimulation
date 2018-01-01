//
//  connection.cpp
//  GSLTest
//
//  Created by Rui Ma on 2/6/17.
//  Copyright © 2017 Yale University. All rights reserved.
//

#include "connection.hpp"
#include "filament.hpp"


Connection::Connection()
{
    s1 = 0.0;   // label the position
    s2 = 0.0;
    index1 = 0;   // label the filament bound
    index2 = 0;
    phi1_1 = 0;     // record the angle in previous step for phase correction in torque calculation
    theta1_1 = 0;
    phi1_2 = 0;
    theta1_2 = 0;
    phi2_1 = 0;
    theta2_1 = 0;
    phi2_2 = 0;
    theta2_2 = 0;
    len = 0.0;
    lifetime = 0;
    transEnergy = 0.0;
    torsEnergy = 0.0;
    
}

Connection::Connection(Filament *p1, Filament *p2, int m1, int m2, double const Cross_angle)
{
    // Creater a linker between the m1th monomer in filament p1 and m2th monmer in filament p2
    f1 = p1;
    f2 = p2;
    p1->OccupiedList.push_front(m1);
    p2->OccupiedList.push_front(m2);
    n1 = p1->OccupiedList.begin();
    n2 = p2->OccupiedList.begin();
    s1=(m1-0.5)*p1->delta/p1->length-0.5;
    s2=(m2-0.5)*p2->delta/p2->length-0.5;
    index1 = p1->index;
    index2 = p2->index;
    vec3 q1 = f1->GetPositionMonomer(m1);
    vec3 q2 = f2->GetPositionMonomer(m2);
    vec3 O1 = f1->GetRotationMonomer(m1);
    vec3 O2 = f2->GetRotationMonomer(m2);
    vec3 T1 = cross(O1, f1->orientation);
    vec3 T2 = cross(O2, f2->orientation);
    
    
    len = norm(q1-q2);
    vec3 q12 = normalise(q1 - q2);
    
    // The first subscript indicates filament1 or filament2, the second subscript indicates the bond1 or bond2 in each filament, bond 1 aligns O1 towards -q12, bond 2 aligns f1->orientation towards n1_ref
    double x1_1 = dot(-q12,f1->orientation);
    double y1_1 = dot(-q12,T1);
    double z1_1 = dot(-q12,O1);
    theta1_1 = acos(z1_1);
    phi1_1 = atan2(y1_1, x1_1);
    
    vec3 n1_ref = cos(Cross_angle)*f2->orientation + sin(Cross_angle)*normalise(cross(-q12, f2->orientation));
    double x1_2 = dot(n1_ref,T1);
    double y1_2 = dot(n1_ref,O1);
    double z1_2 = dot(n1_ref,f1->orientation);
    theta1_2 = acos(z1_2);
    phi1_2 = atan2(y1_2, x1_2);
    
    double x2_1 = dot(q12,f2->orientation);
    double y2_1 = dot(q12,T2);
    double z2_1 = dot(q12,O2);
    theta2_1 = acos(z2_1);
    phi2_1 = atan2(y2_1, x2_1);
    
    vec3 n2_ref = cos(Cross_angle)*f1->orientation + sin(Cross_angle)*normalise(cross(q12, f1->orientation));
    double x2_2 = dot(n2_ref,T2);
    double y2_2 = dot(n2_ref,O2);
    double z2_2 = dot(n2_ref,f2->orientation);
    theta2_2 = acos(z2_2);
    phi2_2 = atan2(y2_2, x2_2);
    
    lifetime = 0;
    transEnergy = 0.0;
    torsEnergy = 0.0;
}

Connection::~Connection()
{
   // cout << "Connection's destructor called" << endl;
    f1->OccupiedList.erase(n1);
    f2->OccupiedList.erase(n2);
    
}



bool Connection::GetForceAndTorque(const double kappa, const double r0, const double kappaR1, const double kappaR2, const double Cross_angle, const double modular_angle, const double Max_distance, const double Max_angle, const double dt)
{
    vec3 p1 = f1->GetPositionMonomer(*n1);
    vec3 p2 = f2->GetPositionMonomer(*n2);
    vec3 O1 = f1->GetRotationMonomer(*n1);
    vec3 O2 = f2->GetRotationMonomer(*n2);
    vec3 T1 = cross(O1, f1->orientation);
    vec3 T2 = cross(O2, f2->orientation);
    
    
    double dis = norm(p1-p2);
    vec3 p12 = normalise(p1-p2);
    
    
    // keeping track of the angles
    
    double x1_1 = dot(-p12,f1->orientation);
    double y1_1 = dot(-p12,T1);
    double z1_1 = dot(-p12,O1);
    double phi1_1_new = atan2(y1_1,x1_1); // phi1 ranges from -pi to pi
    // phase correction
    phi1_1 = phi1_1+remainder(phi1_1_new-phi1_1,modular_angle);// this guarantes that the absolute increase of phi is less than pi/2
    
    vec3 S1_1=cos(phi1_1)*f1->orientation+sin(phi1_1)*T1;
    double a1_1 = dot(-p12,S1_1);
    double theta1_1_new = atan2(a1_1,z1_1);
    theta1_1 = theta1_1+remainder(theta1_1_new-theta1_1,2*modular_angle);
    
    
    double x2_1 = dot(p12,f2->orientation);
    double y2_1 = dot(p12,T2);
    double z2_1 = dot(p12,O2);
    double phi2_1_new = atan2(y2_1,x2_1); // phi1 ranges from -pi to pi
    // phase correction
    phi2_1 = phi2_1+remainder(phi2_1_new-phi2_1,modular_angle);// this guarantes that the absolute increase of phi is less than pi/2
    
    vec3 S2_1=cos(phi2_1)*f2->orientation+sin(phi2_1)*T2;
    double a2_1 = dot(p12,S2_1);
    double theta2_1_new = atan2(a2_1,z2_1);
    
    theta2_1 = theta2_1+remainder(theta2_1_new-theta2_1,2*modular_angle);

    
    
    // force
    len = dis;
    vec3 f12 = -kappa * (dis - r0) * p12;
    f1->force = f1->force + f12;
    f2->force = f2->force - f12;
    
    // force-induced torque
    f1->torque = f1->torque + cross(p1-f1->center_of_mass, f12);
    f2->torque = f2->torque - cross(p2-f2->center_of_mass, f12);
    
//    // torsional energy E = (O1-O2)*p12;
//    // torque-induced force
//    vec3 O12 = O1-O2;
//    f12 = -kappaR1/dis*(O12-dot(O12, p12)*p12);
//    f1->force = f1->force + f12;
//    f2->force = f2->force - f12;
//    
//    // torque
//    f1->torque = f1->torque-kappaR1*(cross(O1, p12)+cross(p1-f1->center_of_mass, O1-O2)/dis-dot(O1-O2, p12)*cross(p1-f1->center_of_mass, p12)/dis);
//    f2->torque = f2->torque-kappaR1*(cross(O2, -p12)+cross(p2-f2->center_of_mass, O2-O1)/dis-dot(O2-O1, -p12)*cross(p2-f2->center_of_mass, -p12)/dis);
    
    // another way the energy is E=1/2 kappaR1 *(theta1_1^2+theta1_2^2)
    double sign1 = (sin(theta1_1)>0)-(sin(theta1_1)<0);
    double sign2 = (sin(theta2_1)>0)-(sin(theta2_1)<0);
    vec3 O1perp = O1-dot(O1, p12)*p12;
    vec3 O2perp = O2-dot(O2, p12)*p12;
    
    f1->force = f1->force - kappaR1*theta1_1*sign1/dis*normalise(O1perp);
    f2->force = f2->force + kappaR1*theta1_1*sign1/dis*normalise(O1perp);
    
    
    f1->force = f1->force + kappaR1*theta2_1*sign2/dis*normalise(O2perp);
    f2->force = f2->force - kappaR1*theta2_1*sign2/dis*normalise(O2perp);
    
    
    
    f1->torque = f1->torque + kappaR1*theta1_1*sign1*(normalise(cross(O1, -p12))-cross(p1-f1->center_of_mass,normalise(O1perp))/dis);
    f2->torque = f2->torque + kappaR1*theta1_1*sign1*cross(p2-f2->center_of_mass, normalise(O1perp))/dis;
    
    
    f1->torque = f1->torque + kappaR1*theta2_1*sign2*cross(p1-f1->center_of_mass, normalise(O2perp))/dis;
    f2->torque = f2->torque + kappaR1*theta2_1*sign2*(normalise(cross(O2, p12))-cross(p2-f2->center_of_mass,normalise(O2perp))/dis);
    
    lifetime = lifetime + dt;
    return true;
    
    
    
    
}

