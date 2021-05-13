/*************************************************************
* @author   Triston Ruiseco
* @file     PHypoSim.C
* @date     05/09/2021
* @brief    Simulates hypothetical Poisson distributed data
            according to user-input data volume, parameters,
            and uncertainties.
*************************************************************/

// Std Includes
#include <iostream>
#include <fstream>

// ROOT Includes
#include "TRandom.h"

// Directives and declarations for namespaces and namespace members
using std::string, std::stod, std::cout, std::stoi, std::ofstream;

// Program-specific helper function declarations
/**
 * Identical string comparison function
 * param a: first string to be compared
 * param b: second string to be compared
 * return 1 if a & b are identical character-wise and in length, 0 otherwise
 */
bool strsame(string a, string b);

// Begin primary program function
int main(int argc, char** argv){

  // Command line option parsing and tracking variables
  bool argexists = 0;
  bool printhelp = 0;
  bool output = 0;

  // Command line option storage variables
  string filename = "data.txt";
  double exposure = 118*85.3;
  double br_mean = 0.0026;
  double br_dev = 0.0006;
  double sr_min = 0.;
  double sr_max = 0.010;
  double sr_step = 0.00005;
  int epr = 1000000;

  // Parse and process command line options
  for(int i = 1; i < argc; ++i){
    argexists = 0;
    if(strsame(argv[i],"--help")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"-h")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"--output")){
      argexists = 1;
      filename = string(argv[++i]);
      output = 1;
    }
    if(strsame(argv[i],"--exposure")){
      argexists = 1;
      double temp = stod(argv[++i]);
      if(temp > 0.0){
        exposure = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(strsame(argv[i],"--br_mean")){
      argexists = 1;
      double temp = stod(argv[++i]);
      if(temp > 0.0){
        br_mean = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(strsame(argv[i],"--br_dev")){
      argexists = 1;
      double temp = stod(argv[++i]);
      if(temp > 0.0){
        br_dev = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(strsame(argv[i],"--sr_min")){
      argexists = 1;
      double temp = stod(argv[++i]);
      if(temp > 0.0){
        sr_min = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(strsame(argv[i],"--sr_max")){
      argexists = 1;
      double temp = stod(argv[++i]);
      if(temp > 0.0){
        sr_max = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(strsame(argv[i],"--sr_step")){
      argexists = 1;
      double temp = stod(argv[++i]);
      if(temp > 0.0){
        sr_step = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(strsame(argv[i],"--epr")){
      argexists = 1;
      int temp = stoi(argv[++i]);
      if(temp > 0.0){
        epr = temp;
      } else {
        cout << "\n" << argv << " must be greater than 0.\n\n";
        return 0;
      }
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
    }
  }

  /* Print the executable usage instructions if the user adds -h or --help,
     doesn't provide required input, or provides an undefined option */
  if(printhelp || !output || sr_min > sr_max){
    cout << "\nUsage: " << argv[0] << " --output [file] [options]\n"
         << "  options and descriptions:\n"
         << "   --help(-h)            print options\n"
         << "   --output [file]       name of output file\n"
         << "   --exposure [number]   experiment exposure in kg days (118*85.3)\n"
         << "   --br_mean [number]    mean of expected background (0.0026)\n"
         << "   --br_dev [number]     stdev of expected background (0.0006)\n"
         << "   --sr_min [number]     minimum signal rate (0.)\n"
         << "   --sr_max [number]     maximum signal rate (0.01)\n"
         << "   --sr_step [number]    signal rate iteration step size (0.00005)\n"
         << "   --epr [int]           experiments per signal rate (1000000)\n\n";

    return 0;
  }
  cout << "\n";

  // Open output stream
  ofstream outFile(filename);
  if(outFile.is_open()){ // If stream has associated file

    // Simulation and data export helper variables
    double br = 0.;

    // Output formatted simulation parameters
    outFile << exposure << " " << br_mean << " " << br_dev << " " << epr << "\n";

    // Simulate experiments and print results to file
    for(double sr = sr_min; sr <= sr_max; sr += sr_step){ // Loop over true signal rates

      // Print true signal rate preceding experimental results
      outFile << sr << " ";
      for(double e = 0; e < epr; e++){ // Loop over experiments within hypothesis

        // Move background rate by some uncertainty for this experiment
        do {
          br = gRandom->Gaus(br_mean,br_dev);
        } while(br <= 0); // Assert that the background rate be greater than 0

        // Simulate counting experiment and send result to output stream
        outFile << gRandom->Poisson(exposure*(br + sr)) << " ";
      }
      outFile << "\n";
    }
  } else { // If stream does not have associated file
    cout << "\nFailure generating output file.\n\n";
  }

  // Close output stream
  outFile.close();

  return 0;
}

// Program-specific helper function definitions
bool strsame(string a, string b){
  if(a.length()==b.length()){
    int n = a.length();
    for(int i = 0; i < n; i++){
      if(a.at(i)!=b.at(i)){
        return 0;
      }
    }
    return 1;
  }
  return 0;
}
