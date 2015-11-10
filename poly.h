#ifndef POLY_H
#define POLY_H
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <Windows.h>
#include <algorithm>

using namespace std;


class poly
{
    public:
        poly();
        ~poly();
        void getpoints();
        void inter();
        void printvec(std::vector<double>);
        void printans(std::vector< std::vector<double> > mat);
        double power(double,  int);
        double det(std::vector< std::vector<double> > );
        double dot(std::vector< double >, std::vector< double >);
        std::vector< std::vector<double> > mult(std::vector< std::vector<double> >, std::vector< std::vector<double> >);
        std::vector< std::vector<double> > invert(std::vector< std::vector<double> >);
        std::vector< std::vector<double> > round(std::vector< std::vector<double> >);
        std::vector< std::vector<double> > interpolate();
    private:
        int degree;
        std::vector< std::vector<double> > coord;
};

#endif // POLY_H
