#include "Orbital.h" // include the header file "Orbital.h"

// define the Body class custom constructor where initial data is set
Body::Body(int i) { // arguement is equal to planetID
  // Variable manipulated with in order to read planetary data
  string valTemp; // Temporary Variable to use in getline() and stod()
  regex reg("[^0-9\\.\\+-eE]"); // Regular Expression in order to validate arguements passed into stod()

  try { // try block to catch and handle errors from code
      planetID = i; // Set the planet ID based of arguement of constructor
    ifstream Body_Data("parameters.txt"); // open the file "parameters.txt" into Body_Data stream in object
    if (Body_Data.is_open()) { // check file opened correctly to run code
      Body_Data.seekg(0, ios::beg); // move the file pointer to the beginning of the file
      for (int lineCount = 0; lineCount <= i; ++lineCount) { // Iterate through the lines of the "parameters.txt" file
        Body_Data.ignore(numeric_limits<streamsize>::max(), '\n'); // Skip all lines up to the relevant body's data
      }
      for (auto member:{&X, &Y, &U, &V, &Mass}){// iterate through the position, velocity, and Mass of the i-th body
        Body_Data >> valTemp; // pass in each value one by one as a string
        if (regex_search(valTemp, reg)) { // use regular expression to check only numbers, optionally in scientific and decimal form, are passed to stod()
          throw runtime_error("Error: There was an invalid arguement passed to stod().\nSolution: Please check your paramters.txt file for incorrect data entry"); // Thrown error for ease of user and error handling
        }
        *member = stod(valTemp); // dereference and convert each string to doubles and store
      }
      Body_Data.close(); // close the file as no longer needed at this moment in time
    }
    else { // if the file cannot be opened
        throw runtime_error("Error: Could not open parameters file.\nSolution: Please load file correctly and try again."); // throw an error message for user
    }
  }
  catch (const runtime_error &msg){ // handle any other runtime errors that may occur
    cerr << msg .what() << endl; // display error message to user
    exit(1); // quit program to avoid inaccurate results or output
  }
  catch (const invalid_argument &ia) { // handle any other invalid arguement errors that may occur
    cerr << "Error: Invalid arguement passed to " << ia.what() << "()\nSolution: Please double check data being passed to function" << endl; // display error message to user
    exit(1); // quit program to avoid inaccurate results or output
  }
}



// Function to zero all members of every RK Constant in vector
void RK::zeroConstants(vector<Constants> &K) { // arguement is vector array of k1, k2, k3, k4
  for (int count = 0; count < 4; count++) { // iterate through k1, k2, k3, k4
    K[count].setX(0); // using setter to assign x position as 0
    K[count].setY(0); // using setter to assign y position as 0
    K[count].setU(0); // using setter to assign x velocity as 0
    K[count].setV(0); // using setter to assign y velocity as 0
  }
}

// Operator Overload for when Constants = Body
Constants& Constants::operator=(Body&Planet){ // Passing Body by reference to avoid large memory usage
    X = Planet.getX(); // using getter to assign x position of planet to constant
    Y = Planet.getY(); // using getter to assign y position of planet to constant
    U = Planet.getU(); // using getter to assign x velocity of planet to constant
    V = Planet.getV(); // using getter to assign y velocity of planet to constant
    Mass = Planet.getMass(); // using getter to assign mass of planet to constant
    planetID = Planet.getNum(); // using getter to assign planetID of planet to constant
    return *this; // returns current object by dereferncing pointer to current object
}

// Operator Overload for when Constants / double
Constants Constants::operator/( double& k)const { // const function as no paramters changed
    return Constants(X/k, Y/k, U/k, V/k); // return copy to avoid overwiting input
}

// Operator Overload for when Constants * double
Constants Constants::operator*( double& k)const{ // const function as no paramters changed
    return Constants(X*k, Y*k, U*k, V*k); // return copy to avoid overwiting input
}

// Operator Overload for when Constants + Constants
Constants Constants::operator+( Constants& other)const{ // const function as no paramters changed
    return Constants(X + other.getX(), Y + other.getY(), U + other.getU(), V + other.getV()); // return copy to avoid overwiting input
}

// Operator Overload for when Body += Constant
Body& Body::operator+=(const Constants& K) {
    X += K.getX(); // using getter to add Constant object's  and Body object's x position to the Body
    Y += K.getY(); // using getter to add Constant object's  and Body object's y position to the Body
    U += K.getU(); // using getter to add Constant object's  and Body object's x velocity to the Body
    V += K.getV(); // using getter to add Constant object's  and Body object's y velocity to the Body
    return *this; // returns current object by dereferncing pointer to current object
}

// Operator Overload for when Body += Constant
Constants& Constants::operator+=(const Constants& K) {
    X += K.getX(); // using getter to add Constant object's  and Body object's x position to the Body
    Y += K.getY(); // using getter to add Constant object's  and Body object's y position to the Body
    U += K.getU(); // using getter to add Constant object's  and Body object's x velocity to the Body
    V += K.getV(); // using getter to add Constant object's  and Body object's y velocity to the Body
    return *this; // returns current object by dereferncing pointer to current object
}

// Body Class Method to print out a given planet's position and velocity data at a chosen time
void Body::printOut(int t, double h) { // passing current time to function via the size of timestep and number of steps from 0
  ofstream Fileout("output.txt", ios::out | ios::app); // Appends to already created "output.txt" file to add each line of data
  Fileout << scientific << setprecision(2); // set the precision of the output to 2 decimal places with scientific notation 
    // write the body's number, the time, X-coordinate, Y-coordinate, X-velocity, Y-velocity to the file
    Fileout << planetID + 1 << "\t" << t * h << "\t" << X << "\t" << Y << "\t" << U << "\t" << V << "\n";
    Fileout.close(); // Close the file as no longer needed
}

// Operator Overload for when Body + Constant
Body Body::operator+( Constants &other) const { // const function as no paramters changed
        return Body(X + other.getX(), Y + other.getY(), U + other.getU(), V + other.getV()); // return copy to avoid overwiting input
    }

// Operator Overload for when Body + Body
Body Body::operator+( Body &other) const { // const function as no paramters changed
        return Body(X + other.getX(), Y + other.getY(), U + other.getU(), V + other.getV()); // return copy to avoid overwiting input
}

// Operator Overload for when Body - Body
Body Body::operator-( Body &other) const { // const function as no paramters changed
        return Body(X - other.getX(), Y - other.getY(), U - other.getU(), V - other.getV()); // return copy to avoid overwiting input
}

// Function to calculate the values of lamda for use in RK4 scheme
double Body::lamdaCalc(const Body& dydx, double G){ // Passing the objects by reference to avoid large memory usage and current object is jth planet
      double dij = sqrt(dydx.getX()*dydx.getX() + dydx.getY()*dydx.getY()); // calculate the distance between the current planet and the other planet with for a given time and RK constant
      return G * (Mass / pow(dij, 3)); // calculate lamda based off equation from theory
}

// Operator Overload for when Constant - Constant
Constants Constants::operator-( Constants &other) const { // const function as no paramters changed
        return Constants(X - other.getX(), Y - other.getY(), U - other.getU(), V - other.getV()); // return copy to avoid overwiting input
}