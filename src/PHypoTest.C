/*************************************************************
* @author   Triston Ruiseco
* @file     PHypoTest.C
* @date     05/09/2021
* @brief    Analyzes and visualizes the significance-signal
						relationship between appropriate Poisson
						distributed data samples.
*************************************************************/

// Std Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

// ROOT Includes
#include "TRandom.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH1D.h"
#include "TLegend.h"
#include "Math/PdfFuncMathCore.h"

// Directives and declarations for namespaces and namespace members
using std::string, std::vector, std::stod, std::cout, std::ifstream,
			std::stringstream, std::sort, std::to_string;

// Global variables and objects
string prepend;

// Program-specific helper function declarations
/**
 * Identical string comparison function
 * param a: first string to be compared
 * param b: second string to be compared
 * return 1 if a & b are identical character-wise and in length, 0 otherwise
 */
bool strsame(string a, string b);

/**
 * return GaussianPDF(x,0,1)
 */
double sigma(int x);

/**
 * return first index of arr such that arr[index] is of lesser value than y
 */
int FirstIndexLess(const vector<double>& arr, double y);

/**
 * Equal type I and type II error probability finding function
 * param arr0: vector of some distribution for some hypothesis 0
 * param arr1: vector of some distribution for some hypothesis 1
 * return significance level of test such that alpha = beta
 */
double FindABSame(const vector<double> &arr0, const vector<double> &arr1);

/**
 * Confidence vs measurements plot generation and exportation function
 */
void PlotConfidence(vector<double>& arr0, vector<double>& arr1, const string& title);

// Begin primary program function
int main(int argc, char** argv){

  // Command line option parsing and tracking variables
  bool argexists = 0;
  bool printhelp = 0;
  bool input = 0;

  // Command line option storage variables
  string filename = "data.txt";

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
    if(strsame(argv[i],"--input")){
      argexists = 1;
      filename = string(argv[++i]);
      input = 1;
    }
    if(strsame(argv[i],"--prepend")){
      argexists = 1;
      prepend = string(argv[++i]);
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
      return 0;
    }

  }

  /* Print the executable usage instructions if the user adds -h or --help,
     doesn't provide required input, or provides an undefined option */
  if(printhelp || !input){
    cout << "\nUsage: " << argv[0] << " --input [file] [options]\n"
         << "  options and descriptions:\n"
         << "   --help(-h)            print options\n"
         << "   --input [file]        name of input file\n"
				 << "		--prepend [string]    string to prepend all output filenames with\n\n";

    return 0;
  }
  cout << "\n";

  // Simulation data storage objects
  vector<vector<double>> hypo_sims;
	double exposure = 0.;
	double br_mean = 0.;
	double br_dev = 0.;
	double epr = 0.;

  // Data import and storage
  ifstream inFile(filename);

  if(inFile.is_open()){ // If stream has associated file

    // Simulation data import helper variables
    string line;
    double tempd;

		// Extract formatted simulation parameters
		getline(inFile, line);
		stringstream params(line);
		params >> exposure >> br_mean >> br_dev >> epr;

    // Import the data from each experiment into a separate vector
    while(getline(inFile,line)){ // Place entire hypothesis sim in one string

      // Place all experiments from one hypothesis into a formatted stream
      stringstream ss(line);

      // Read all experiments from stream into single vector
      vector<double> tempv;
      while(ss >> tempd){
        tempv.push_back(tempd);
      }

      // Store all hypothesis simulations in 2D vector sorted by hypothesis
      hypo_sims.push_back(tempv);
    }

  } else { // If stream does not have associated file
    cout << "\nFailure accessing input file.\n\n";
    return 0;
  }

  // Data analysis helper variables
  double tempLLR = 0.;

  // Data analysis storage objects
  vector<double> altr_LLR;       // LLR of data from alternate hypotheses
  vector<double> null_LLR;       // LLR of data from null hypothesis
  vector<double> ab_sigs;        // significance level of each generated test
  vector<double> signal_rates;   // keeps index-wise track of signals for ab_sigs

  int H = hypo_sims.size(); // Number of simulated hypotheses
  for(int h = 1; h < H; h++){ // Loop across all simulated hypotheses

    // Retrieve signal rate used for simulation of hypothesis number h
    double signal_rate = hypo_sims[h][0];

    // Reset test statistic distribution storage
    altr_LLR.clear();
    null_LLR.clear();

    int E = hypo_sims[h].size(); // Number of experiments simulated
    for(int e = 1; e < E; e++){ // Loop across all simulated experiments

      // Calculate and store log-likelihood ratio under alternate hypothesis
      tempLLR = 0.;
      tempLLR += log(ROOT::Math::poisson_pdf(hypo_sims[h][e], exposure*(br_mean + signal_rate)));
      tempLLR -= log(ROOT::Math::poisson_pdf(hypo_sims[h][e], exposure*br_mean));
      altr_LLR.push_back(tempLLR);

      // Calculate and store log-likelihood ratio under null hypothesis
      tempLLR = 0.;
      tempLLR += log(ROOT::Math::poisson_pdf(hypo_sims[0][e], exposure*(br_mean + signal_rate)));
      tempLLR -= log(ROOT::Math::poisson_pdf(hypo_sims[0][e], exposure*br_mean));
      null_LLR.push_back(tempLLR);
    }

    // Sort test distributions
    sort(null_LLR.begin(),null_LLR.end());
    sort(altr_LLR.begin(),altr_LLR.end());

		// Find and record alpha=beta
		double ab = FindABSame(null_LLR, altr_LLR);
		ab_sigs.push_back(ab);
		signal_rates.push_back(signal_rate);
  }

  // Plot and save results
  PlotConfidence(signal_rates, ab_sigs,
		Form("#lambda_{background} = (%.1f #pm %.1f) #times 10^{-3} keV_{ee}^{-1} kg^{-1} day^{-1}, %.1f kg day exposure", (br_mean*1000), (br_dev*1000), exposure));

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

double sigma(int x){
  return erf(double(x)/sqrt(2));
}

int FirstIndexLess(const vector<double>& arr, double y){
  int n = arr.size();
  for(int i = 0; i < n; ++i){
    if(arr[i] < y){
      return i;
    }
  }
  return n;
}

double FindABSame(const vector<double>& arr0, const vector<double>& arr1){
  int n0 = arr0.size();
  int n1 = arr1.size();

  double a = n0 - 1;
  double b = 0;

  bool stop = 0;

  while(!stop){
    if(abs(arr0[a] - arr1[b]) >= abs(arr0[a-1] - arr1[b+1])){
      a--;
      b++;
    } else {
      stop = 1;
    }
  }

  b /= double(n1);

  return b;
}

void PlotConfidence(vector<double>& arr0, vector<double>& arr1, const string& title){
  TCanvas* canvas = (TCanvas*) new TCanvas("canvas", "Canvas_Title",200,10,500,400);
  // canvas->SetGrid();

	// Configure Canvas
  double lm = 0.15;
  double rm = 0.04;
  double bm = 0.15;
  double tm = 0.07;
  canvas->SetLeftMargin(lm);
  canvas->SetRightMargin(rm);
  canvas->SetBottomMargin(bm);
  canvas->SetTopMargin(tm);
  canvas->SetLogy();
  canvas->Draw();
  canvas->Update();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// Configure Graph
  int N = arr0.size();
  TGraph* graph = new TGraph(N, &(arr0[0]), &(arr1[0]));
  graph->SetLineColor(kAzure+1);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(0);

	// Configure X
	graph->GetXaxis()->SetTitleFont(42);
	graph->GetXaxis()->SetTitleSize(0.05);
	graph->GetXaxis()->SetTitleOffset(1.1);
	graph->GetXaxis()->SetLabelFont(42);
	graph->GetXaxis()->SetLabelSize(0.04);
	graph->GetXaxis()->SetTickSize(0.);
  graph->GetXaxis()->SetTitle("Simulated event rate (day^{-1})");
	graph->GetXaxis()->CenterTitle();

	// Configure Y
	graph->GetYaxis()->SetTitleFont(42);
	graph->GetYaxis()->SetTitleSize(0.05);
	graph->GetYaxis()->SetTitleOffset(1.1);
	graph->GetYaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetTitle("Test Significance #alpha (= #beta)");
	graph->GetYaxis()->CenterTitle();

	// Draw configured graph on canvas
  graph->Draw();
  canvas->Update();

	// Format default latex text drawing settings
  TLatex text;
	text.SetTextFont(42);
  text.SetTextAlign(21);
  text.SetTextSize(0.04);
  text.SetNDC();

	// Draw title
  text.DrawLatex((1.-rm+lm)/2., 1.-tm+0.012, title.c_str());

	// Draw lines at sigma levels of significance and label them
  TLine* line = new TLine();
  line->SetLineWidth(2);

  int n = arr1.size();
  int i = 1;
  int k = FirstIndexLess(arr1, 1 - sigma(i));

	text.SetTextSize(0.035);
	text.SetTextAlign(33);
	text.SetTextAngle(90);

  while(k < n && i <= 5){
    double xndc = (1.-rm-lm)*((arr0[k]-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+lm;
    cout << i << "sigma at " << arr0[k] << "\n";
    line->SetLineColor(kRed+5-i);
    line->DrawLineNDC(xndc,bm,xndc,1.-tm);
    text.DrawLatex(xndc+0.005, 1-tm-0.01, Form("x = %3.1f #times 10^{-3}, #alpha = %d #sigma", (arr0[k]*1000.) , i));
    ++i;
    k = FirstIndexLess(arr1, 1 - sigma(i));
  }

	// Export in log scale as well as normal
  canvas->SaveAs((prepend + "confidenceLogY.png").c_str());
  canvas->SetLogy(0);
  canvas->SaveAs((prepend + "confidence.png").c_str() );
}
