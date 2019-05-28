class Vector {
  public:
    double  x, y, z;
  public:
    void init(double x_, double y_, double z_){
        this->x = x_; this->y = y_; this->z = z_;
    }
    Vector(double x_, double y_, double z_){
        init(x_, y_, z_);
    }
    Vector(){}
    Vector operator + (Vector opn){
        return Vector(x + opn.x, y + opn.y, z + opn.z);
    }
    void operator += (Vector opn){
        x += opn.x; y += opn.y; z += opn.z;
    }
    Vector operator - (Vector opn){
        return Vector(x - opn.x, y - opn.y, z - opn.z);
    }
    void operator -= (Vector opn){
        x -= opn.x; y -= opn.y; z -= opn.z;
    }
    Vector operator * (double opn){
        return Vector(x * opn, y * opn, z * opn);
    }
    void operator *= (double opn){
        x *= opn; y *= opn; z *= opn;
    }
    double operator * (Vector opn){     //dot product
        return x * opn.x + y * opn.y + z * opn.z;
    }
    Vector operator & (Vector o){       //cross product
        return Vector(y*o.z-z*o.y, z*o.x-x*o.z, x*o.y-y*o.x);
    }
    double operator ^ (Vector o){
        double cn = *this * o / this->mod() / o.mod(); if (cn > 1) cn = 1; else if (cn < -1) cn = -1;
        double ret = acos(cn);
        if (ret < 0) ret += 3.14159265358979323846264338327950288;
        else if (ret > 3.14159265358979323846264338327950288) ret = 2*3.14159265358979323846264338327950288 - ret;
        return ret;
    }
    Vector operator / (double opn){
        return Vector(x / opn, y / opn, z / opn);
    }
    void operator /= (double opn){
        x /= opn; y /= opn; z /= opn;
    }
    Vector Renorm(){
        double r = mod(); if (r != 0) x /= r; y /= r; z /= r; return *this;
    }
    double pow2(){
        return (x * x + y * y + z * z);
    }
    double mod(){
        return sqrt(x * x + y * y + z * z);
    }
};
Vector RotateTF(Vector vec, double theta,  double phi){
    double ct = cos(theta);
    double cf = cos(phi);
    double st = sin(theta);
    double sf = sin(phi);
    double a[3][3];
    a[0][0] =  ct * cf;
    a[0][1] = -ct * sf;
    a[0][2] =  st;
    a[1][0] =  sf;
    a[1][1] =  cf;
    a[1][2] =  0;
    a[2][0] = -st * cf;
    a[2][1] =  st * sf;
    a[2][2] =  ct;
    double x = a[0][0] * vec.x + a[0][1] * vec.y + a[0][2] * vec.z;
    double y = a[1][0] * vec.x + a[1][1] * vec.y + a[1][2] * vec.z;
    double z = a[2][0] * vec.x + a[2][1] * vec.y + a[2][2] * vec.z;
    return Vector(x, y, z);
}
/*Vector Rotate1(Vector vec, double theta, double psi, double phi){
    double c1 = cos(theta);
    double c2 = cos(psi);
    double c3 = cos(phi);
    double s1 = sin(theta);
    double s2 = sin(psi);
    double s3 = sin(phi);
    double l1 = c2 * c3 - c1 * s2 * s3;
    double l2 = - c2 * s3 - c1 * s2 * c3;
    double l3 = s1 * s2;
    double m1 = s2 * c3 + c1 * c2 * s3;
    double m2 = - c2 * s3 + c1 * c2 * c3;
    double m3 = - s1 * c2;
    double n1 = s1 * s3;
    double n2 = s1 * c3;
    double n3 = c1;
    double z = l1 * vec.x + l2 * vec.y + l3 * vec.z;
    double x = m1 * vec.x + m2 * vec.y + m3 * vec.z;
    double y = n1 * vec.x + n2 * vec.y + n3 * vec.z;
    return Vector(x, y, z);
}
Vector Rotate2(Vector vec, double theta, double psi, double phi){
    double t, p;
    double rxy = sqrt(vec.x * vec.x + vec.y * vec.y);
    double r = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    if (rxy == 0) p = 0; else p = acos(vec.x / rxy); if (vec.y < 0) p = p + 3.1416;
    if (r == 0) t = 0; else t = acos(vec.z / r);
    t += theta; p += phi;
    return Vector(r * sin(t) * cos(p), r * sin(t) * sin(p), r * cos(t));
}*/
Vector RotateZX(Vector vec, double theta, double phi){
    double ct = cos(theta);
    double st = sin(theta);
    double cp = cos(phi);
    double sp = sin(phi);
    double x = vec.x * cp      - vec.y * sp;
    double y = vec.x * ct * sp + vec.y * ct * cp - vec.z * st;
    double z = vec.x * st * sp + vec.y * st * sp + vec.z * ct;
    return Vector(x, y, z);
}
Vector RotateYX(Vector vec, double theta, double phi){
    double ct = cos(theta);
    double st = sin(theta);
    double cp = cos(- phi);
    double sp = sin(- phi);
    double x = vec.x * cp - vec.z * sp;
    double y = - vec.x * st * sp + vec.y * ct - vec.z * st * cp;
    double z = vec.x * ct * sp + vec.y * st + vec.z * ct * cp;
    return Vector(x, y, z);
}

#define Rotate RotateTF
