#ifndef Orbital_H // Guard to prevent multiple inclusion
#define Orbital_H
#include <vector> // Used for container to hold Body objects
#include <iostream> // Used for input/output operations
#include <fstream> // Used for file operations
#include <iomanip> // Used for formatting output
#include <string> // Used for string operations
#include <limits> // Used for increasing limit for ignoring lines
#include <exception> // Used for error handling
#include <sstream> // Used for string operations and conversions
#include <regex> // Used for data validation
#include <algorithm> // Used for lamda expressions
#include <cmath> // Used for mathematical functions

using namespace std; // To avoid repetetive typing

namespace RK { // To predeclare namespace for use of derived class in base class method
    class Constants; // To predeclare derived class for use in base class method
}
using namespace RK; // To avoid repetetive typing

namespace Celestial{ // Separate namespace for created class to avoid any future name collision
    class Body{ // Body class that defines a celestial body
        protected: // protected variables to avoid redelcaration in derived class RK Constants
            mutable double X; // X-coordinate of the body
            mutable double Y; // Y-coordinate of the body
            mutable double U; // X-velocity of the body
            mutable double V; // Y-velocity of the body
            mutable int planetID; // planet ID / number of body
            mutable double Mass ; // Mass of the body
        public: // Public for use of class methods in main  and source file with ease
            Body(){ setX(0.0).setY(0.0).setU(0.0).setV(0.0).setMass(0.0).setplanetID(0); } // Default constructor to zero values
            Body(double x, double y, double u,double v){ setX(x).setY(y).setU(u).setV(v); }
            Body(int i); // Create body with initial data from "parameters.txt" file
            Body& setX(double x){ X = x; return *this; } // Setter for Body's x position
            Body& setY(double y){ Y = y; return *this; } // Setter for Body's y position
            Body& setU(double u){ U = u; return *this; } // Setter for Body's x velocity
            Body& setV(double v){ V = v; return *this; } // Setter for Body's y velocity
            Body& setMass(double mass){ Mass = mass; return *this; } // Setter for Body's mass
            Body& setplanetID(int i){ planetID = i; return *this; } // Setter for Body's planetID
            double getX()const{ return X; } // Setter for Body's x position
            double getY()const{ return Y; } // Setter for Body's y position
            double getU()const{ return U; } // Setter for Body's x velocity
            double getV()const{ return V; } // Setter for Body's y velocity
            double getMass()const{ return Mass; } // Setter for Body's mass
            int getNum()const{ return planetID; } // Setter for Body's planetID
            void printOut(int t, double h); // Printing function for Body at a given time to "output.txt"
            Body& operator+=(const Constants& K); // Overloading += operator for Constant onto a Body
            ~Body(){} // destructor for when objects go out of scope
            Body operator+( Constants &other) const; // asignment operator overloaded to ease addition of Body and Constant object
            Body operator+( Body &other) const; // asignment operator overloaded to ease addition of two bodies
            Body operator-( Body &other) const; // asignment operator overloaded to ease subtraction of two bodies
            Body operator+( Constants &&other) const {return operator+(other);} // asignment operator overloaded allow for reurn and use of temporary objects and values
            Body operator+( Body &&other) const {return operator+(other);} // asignment operator overloaded allow for reurn and use of temporary objects and values
            Body operator-( Body &&other) const {return operator-(other);} // asignment operator overloaded allow for reurn and use of temporary objects and values
            double lamdaCalc(const Body& dydx, double G);
    };
}
using namespace Celestial; // To avoid repetetive typing

// Namespace RK contains functions related to the Runge-Kutta method for solving the simulation
namespace RK {
    // Constants class that holds temporary values for Runge-Kutta method
    class Constants:public Body{
        public:
        Constants():Body(){}; // Default constructor
        Constants(double x, double y, double u,double v):Body(x, y, u, v){}
        Constants& operator=(Body&Planet); // asignment operator overloaded to ease storing a temp array of planets at each timestep
        Constants operator*( double &k)const; // asignment operator overloaded to ease scalar multiplying of object
        Constants operator/( double &k)const; // asignment operator overloaded to ease scalar division of object
        Constants& operator+=(const Constants&K); // asignment operator overloaded to ease addition assignment of two objects
        Constants operator-( Constants&K)const; // asignment operator overloaded to ease subtraction of two objects
        Constants operator+( Constants&K)const; // asignment operator overloaded to ease addition of two objects
        Constants operator*( double &&k)const{return operator*(k);} // asignment operator overloaded allow for reurn and use of temporary objects and values
        Constants operator/( double &&k)const{return operator/(k);} // asignment operator overloaded allow for reurn and use of temporary objects and values
        Constants operator+( Constants&&K)const{return operator+(K);} // asignment operator overloaded allow for reurn and use of temporary objects and values
        Constants operator-( Constants&&K)const{return operator-(K);} // asignment operator overloaded allow for reurn and use of temporary objects and values
    };
    void zeroConstants(vector<Constants>&); // Zeroes out temporary values of all Constant objects in vector array
}

#endif