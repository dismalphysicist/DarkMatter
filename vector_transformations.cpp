#include <iostream>
#include <cmath>
#include <array>

using namespace std;

const float pi = acos(-1);

void print (array<float,3> v);

array<float,3> rotate_x (array<float,3> v, float angle) {
    /* Rotates the vector v by an angle 'angle' anticlockwise about the x axis. */
    array<float,3> v_;
    v_[0] = v[0]; //x is unchanged 
    v_[1] = cos(angle)*v[1] - sin(angle)*v[2];
    v_[2] = sin(angle)*v[1] + cos(angle)*v[2];
    return v_;
}

array<float,3> rotate_y (array<float,3> v, float angle) {
    /* Rotates the vector v by an angle 'angle' anticlockwise about the y axis. */
    array<float,3> v_;
    v_[0] = cos(angle)*v[0] + sin(angle)*v[2];
    v_[1] = v[1]; //y is unchanged
    v_[2] = -sin(angle)*v[0] + cos(angle)*v[2];
    return v_;
}

array<float,3> rotate_z (array<float,3> v, float angle) {
    /* Rotates the vector v by an angle 'angle' anticlockwise about the z axis. */
    array<float,3> v_;
    v_[0] = cos(angle)*v[0] - sin(angle)*v[1];
    v_[1] = sin(angle)*v[0] + cos(angle)*v[1];
    v_[2] = v[2]; //z is unchanged 
    return v_;
}

float dot_product (array<float,3> v1, array<float,3> v2) {
    /* Returns the scalar (dot) product between the vectors v1 and v2. */
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; 
}

array<float,3> cross_product (array<float,3> v1, array<float,3> v2) {
    /* Computes the vector (cross) product of vectors v1 and v2, using the right hand rule. */ 
    array<float,3> xp;
    xp[0] = v1[1]*v2[2] - v1[2]*v2[1];
    xp[1] = v1[2]*v2[0] - v1[0]*v2[2];
    xp[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return xp;
}

float mag (array<float,3> v) {
    /* Returns the magnitude/length/absolute value of the vector v. */ 
    //sqrt(x^2 + y^2 + z^2)
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

float theta (array<float,3> v) {
    /* Returns the polar coordinate theta (declination from z axis) of the vector v. */ 
    //tan(theta) = sqrt(x^2+y^2)/z
    return atan(sqrt(v[0]*v[0]+v[1]*v[1])/v[2]);  
}

float phi (array<float,3> v) {
    /* Returns the polar coordinate phi (azimuthal angle in x-y plane) of the vector v. */
    //tan(phi) = y/x 
    return atan(v[1]/v[0]);  
}

float angle_between (array<float,3> v1, array<float,3> v2) {
    /* Returns the angle between vector v1 and vector v2. */
    return acos(dot_product(v1,v2)/(mag(v1)*mag(v2)));
}

float psi (array<float,3> v1, array<float,3> v2) {
    /* Returns the coordinate system rotation angle needed to make the new y-components of v1 and v2 equal. */
    float alpha = angle_between(v2,v1);
    float beta = angle_between({0,1,0},v1);
    return beta + atan((mag(v2)*cos(alpha) - mag(v1)) / sin(alpha));
}

array<float,3> unit_normal (array<float,3> v1, array<float,3> v2) {
    /* Returns a unit vector normal to the plane defined by the vectors v1 and v2. */ 
    array<float,3> xp = cross_product(v1,v2);
    array<float,3> n;
    n[0] = xp[0] / mag(xp);
    n[1] = xp[1] / mag(xp);
    n[2] = xp[2] / mag(xp);
    return n;
}

array<float,3> transformation (array<float,3> v, array<float,3> a1, array<float,3> a2) {
    /* Returns the coordinates of the vector v in a new coordinate system defined by the vectors a1 and a2. 
    The new z-axis will be perpendicular to a1 and a2. */

    array<float,3> n = unit_normal(a1,a2);
    array<float,3> v_ = rotate_z(v,-phi(n));
    //print(v_); //debugging 
    array<float,3> v__ = rotate_y(v_,-theta(n));

    array<float,3> v___ = rotate_z(v__,psi(a1,a2));

    return v___;
}

void print (array<float,3> v) {
    /* A convenient print method for arrays. */ 
    cout << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
}


int main() {
    array<float,3> i = {1,0,0};
    array<float,3> j = {0,1,0};
    array<float,3> k = {0,0,1};

    ////////////////////// tests //////////////////////

    array<float,3> test = {sin(pi/6)*cos(pi/4), sin(pi/6)*sin(pi/4), cos(pi/6)};

    /*
    print(rotate_x(k,-pi/2));
    print(cross_product(i,j));
    print(cross_product(i,k));

    cout << mag(test) << endl;
    cout << theta(test) << endl;
    cout << phi(test) << endl;
    */

    print(unit_normal({1,0,1},{-1,0,1}));
    print(transformation({0,-1,0}, {1,0,1}, {-1,0,1}));
    print(transformation(i, {1,0,1}, {-1,0,1}));
    print(transformation(j, {1,0,1}, {-1,0,1}));
    print(transformation(k, {1,0,1}, {-1,0,1}));

    return 0;
}