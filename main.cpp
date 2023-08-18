#include "Orbital.h" // Includes the Orbital header file

int main() {
  // Initialising parameters to be read and later used in RK4 Sheme
  double G; // Gravitational Constant
  double T; // Time Period
  double h; // Time Step
  int N; // Total Number of Planets
  int numPoints; // Number of Steps in Time Domain
  int NumLines = 0; // Number of lines in "parameters.txt"

  // Initialising variables to be manipulated to help read afformentioned paramters in the code below
  string stringTemp; // Temporary Variable to use in getline() and stod()
  ostringstream oss; // String Stream Object in order to print vector containg blank line numbers
  regex reg("[^0-9\\.\\+-eE]"); // Regular Expression in order to validate arguements passed into stod()
  vector<int>emptyline{}; // Initialising Vector to contain list of empty line numbers for easy user experience
  
  cout << "\n- Checking 'parameters.txt' file..." << endl;

  // Check "parameters.txt" for errors and read parameters
  try { // try block to catch and handle errors from code
    ifstream Parameters("parameters.txt"); // open the file "parameters.txt" into Parameters stream in object
    if (Parameters.is_open()) { // check file opened correctly to run code
    Parameters.seekg(0, ios::beg); // manually move cursor to beginning of file
      while (!Parameters.eof()) { // whilst end of "parameters.txt" file hasnt been reached
        ++NumLines; // add up number of lines to find out how many Planets are in file
        getline(Parameters, stringTemp, '\n'); // store each line of file into stringTemp one by one
        if(stringTemp.find_first_not_of(" \t\r\n") == string::npos){ // if whole line is whitespace then blank line detected
          emptyline.push_back(NumLines); // add blank line number to sotring vector
        }
      }
      if(!emptyline.empty()){ // emptyline was initially empty then if any blank line detected then error thrown
        for(auto it = emptyline.begin(); it != emptyline.end(); it++){ // iterating through emptyline vector to print it in error message
          oss << *it << (it == emptyline.end() - 1 ? "" : ","); // convert vector into comma separated string as long as not last element with string stream object
        }
        throw runtime_error("Error: There are empty lines at lines: " + oss.str()+"\nSolution: Please remove or fix them before running again"); // Thrown error for ease of user and error handling
      }
      N = --NumLines; // if no blank lines detected then provided format of "parameters.txt" means total number of planets is number of lines - 1
      Parameters.seekg(0, ios::beg); // move the file pointer to the beginning of the file to read G,T and h
      for (auto Variables:{&G, &T, &h}){ // read the values of G, T, and h from the first three lines of the file, passed as pointers (Variables) to original memory
        Parameters >> stringTemp; // pass in each value one by one as a string
        if (regex_search(stringTemp, reg)) { // use regular expression to check only numbers, optionally in scientific and decimal form, are passed to stod()
          throw runtime_error("Error: There was an invalid arguement passed to stod().\nSolution: Please check your paramters.txt file for non-numerical data"); // Thrown error for ease of user and error handling
        }
        *Variables = stod(stringTemp); // dereference and convert each string to doubles and store
      }
      numPoints = (T / h + 1); // Calculate the number of points in the mesh of our numerical scheme e.g. 0 -> T
      Parameters.close(); // close the file once paramters have been read
    }
    else { // if the file cannot be opened
      throw runtime_error("Error: Could not open parameters file.\nSolution: Please load file correctly and try again."); // throw an error message for user
    }
  }
  catch (const runtime_error &msg){ // handle any other runtime errors that may occur
    cerr << msg .what() << endl; // display error message to user
    return 1; // quit program to avoid inaccurate results or output
  }
  catch (const invalid_argument &ia) { // handle any other invalid arguement errors that may occur
    cerr << "Error: Invalid arguement passed to " << ia.what() << "()\nSolution: Please double check data being passed to function" << endl; // display error message to user
    return 1; // quit program to avoid inaccurate results or output
  }

  // Create Output file
  ofstream Fileout("output.txt", ios::out | ios::trunc); // open the file "output.txt" into Fileout stream out object to creat/overwrite file for later appending
  Fileout.close(); // close file as no longer needed at this time

  // Create Planet objects with initial data
  vector<Body> Planet(N); // now N is known , create a vector array of planet/body objects for ease of automation, this calls default constructor and zeroes values
  for (int i = 0; i < N; ++i) { // iterate through each object within the vector
    Planet[i] = Body(i); // assign each body object using the Body constructor which will set all initial data
    Planet[i].printOut(0, h); // print each objects initial data to the "output.txt" file
  }

  // Initialising variables to be manipulated to help solve the ODE using 4th order Runge Kutta (RK4) in the code below
  //vector<Constants> K(4) ; // initialise vector of Runge Kutta constants objects to represent k1,k2,k3,k4 in theory, this calls default constructor and zeroes values
  vector<Constants> planetTemp1(N); // initialise vector array of constants objects to keep planetary data accurate to a given timestep, this calls default constructor and zeroes values
  vector<Constants> planetTemp2(N);
  vector<std::vector<Constants>> K(N,vector<Constants>(4)); // N by 4 
  double lamdaij; // variable equal to G*m/(d12^3)
  double kFactor; // multiplication factor used in calculation for k1, k2, k3, k4
  Constants zeroTemp = Constants(); // zeroed out Constants object for use in ternary operators
  
  cout << "\n- Please wait. Only " << numPoints - 1 << " time steps to go..." << endl;

  // Numerically solve for position and velocity data of bodies at each timestep within domain
  for (int t = 1; t < numPoints ; ++t) { // iterate from 0 to T-1 in order to calculate data up to and including T 
    for_each(K.begin(),K.end(),[&](vector<Constants> &temp){zeroConstants(temp);}); // Call function that zeros every member of each RK Consatnt object
    for (int count = 0; count < 4; count++) { // loop through RK Consatnt vector to condense long k1, k2, k3, k4 calculations
      kFactor = ((count == 1 || count == 2) / 2.0 + (count==0 || count==3)); // multiplication factor is 1 when k1 and k4, and 0.5 when k2 and k3
      for (int i = 0; i < N; i++){ // Loop through all bodies in order to solve for data of planet i at time = t + h
        for (int j = 0; j < N; j++){ // Loop through all OTHER bodies (planet j) to calculate their effect upon the planet i
          if (j!=i) { // dont include planet i effect on planet i as not needed / invalid
            auto k_i = (count > 0 ? K[i][count - 1] : zeroTemp); // equal to the previous RK constant for jth body
            auto k_j = (count > 0 ? K[j][count - 1] : zeroTemp); // equal to the previous RK constant for ith body
            auto dydx = Planet[j] + k_j * kFactor - Planet[i] - k_i * kFactor; // equal to the relative distance between i and j body at time = t + h
            auto dudv = Planet[i] + k_i * kFactor; // equal to the u and v velocity of i body at time = t+ h
            lamdaij = Planet[j].lamdaCalc(dydx,G); // Calculate lamdaij using class method by passing planet j
            auto prevKcount = Constants(0,0,K[i][count].getU(),K[i][count].getV()); // Used to add previous constants to account for previous j bodies effect (0 for positions since j doesnt affect it  in ODE)
            auto derivative = Constants(dudv.getU(), dudv.getV(), lamdaij * dydx.getX(), lamdaij*dydx.getY()); // devative calculated at time = t + h
            K[i][count] = prevKcount + derivative*h; //Contributions of every j body at intermediate time steps add up to each RK constant for ith body
          }
        }
      } 
    }
    for(int number = 0; number < N; number++){ // iterate through each planets
      Planet[number] += (K[number][0] + (K[number][1]) * 2.0 + (K[number][2]) * 2.0 + K[number][3]) / 6.0; // Adds total effect of jth bodies to ith body with accurate RK constants
      Planet[number].printOut(t, h); // Print new data for each body at time t
    }
  }

  cout << "\n- Success! Please check output file and plot data to view results." << endl; // Print Success message if no error has occurred and output created successfully
  return 0;
}