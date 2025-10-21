#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGClient.h"
#include <TRandom3.h>
#include "TStyle.h"
#include <TMatrixD.h> 
#include <TVectorD.h> 

#include <iostream>
using namespace std;

using TMath::Log;

TRandom3* gRand=nullptr; 

//parms
const double xmin=1;
const double xmax=20;
const int npoints=12;
const double sigma=0.2;

double f(double x){
  const double a=0.5;
  const double b=1.3;
  const double c=0.5;
  return a+b*Log(x)+c*Log(x)*Log(x);
}

//_________________________________________________________________________________________
void getX(double *x){
  double step=(xmax-xmin)/npoints;
  for (int i=0; i<npoints; i++){
    x[i]=xmin+i*step;
  }
}

//_________________________________________________________________________________________
void getY(const double *x, double *y, double *ey){
  static TRandom2 tr(0);
  for (int i=0; i<npoints; i++){
    y[i] = f(x[i]) + gRand->Gaus(0.,sigma);
    ey[i]=sigma;
  }
}

// Here, the vector of std::function<double(double)> objects 
// represents one for each term. see below for implementation. 
vector<double> fit_fcn_to_data(
    const vector<function<double(double)>>& fcn_terms, 
    const vector<double>& X, 
    const vector<double>& Y,
    const vector<double>& Y_error
);


//_________________________________________________________________________________________
void leastsq(){
  double x[npoints];
  double y[npoints];
  double ey[npoints];
  getX(x);
  getY(x,y,ey);
  auto tg = new TGraphErrors(npoints,x,y,0,ey);
  tg->Draw("alp");
}

//_________________________________________________________________________________________


//_________________________________________________________________________________________
int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  //initialize the TRandom3 object
  gRand = new TRandom3; 

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  gStyle->SetOptStat(0); // turn off histogram stats box

  TCanvas *tc = new TCanvas("c1","Sample dataset",dw,dh);

  double lx[npoints];
  double ly[npoints];
  double ley[npoints];

  getX(lx);
  getY(lx,ly,ley);
  auto tgl = new TGraphErrors(npoints,lx,ly,0,ley);
  tgl->SetTitle("Pseudoexperiment;x;y");
  
  // An example of one pseudo experiment
  tgl->Draw("alp");
  tc->Draw();
  
  // *** modify and add your code here ***
  auto coeffs = 

  TH2F *h1 = new TH2F("h1","Parameter b vs a;a;b",100,0,1,100,0,1);
  TH2F *h2 = new TH2F("h2","Parameter c vs a;a;c",100,0,1,100,0,1);
  TH2F *h3 = new TH2F("h3","Parameter c vs b;b;c",100,0,1,100,0,1);
  TH1F *h4 = new TH1F("h4","reduced chi^2;;frequency",100,0,1);

  // perform many least squares fits on different pseudo experiments here
  // fill histograms w/ required data
  
  TCanvas *tc2 = new TCanvas("c2","my study results",200,200,dw,dh);
  tc2->Divide(2,2);
  tc2->cd(1); h1->Draw("colz");
  tc2->cd(2); h2->Draw("colz");
  tc2->cd(3); h3->Draw("colz");
  tc2->cd(4); h4->Draw();
  
  tc2->Draw();

  // **************************************
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

//_________________________________________________________________________________________
vector<double> fit_fcn_to_data(
    const vector<function<double(double)>>& fcn_terms,
    const vector<double>& X, 
    const vector<double>& Y,
    const vector<double>& Y_error
)
{
    const int DoF = fcn_terms.size(); 

    if (DoF < 1) return {}; 

    double elems[DoF*DoF] = {0.}; 
    TMatrixD A(DoF, DoF, elems); 
    TVectorD B(DoF, elems); 

    for (int i=0; i<X.size(); i++) {

        double x = X.at(i); 
        double y = Y.at(i); 
        double sigma = (Y_error.at(i) * Y_error.at(i)); 

        vector<double> Xi; 
        for (int j=0; j<DoF; j++) Xi.push_back( (fcn_terms.at(j))(x) ); 

        for (int j=0; j<DoF; j++) {

            B(j) += y * Xi.at(j) / sigma; 

            for (int k=0; k<DoF; k++) A(j, k) += Xi.at(j) * Xi.at(k) / sigma; 
        }
    }

    auto coeffs = A.Invert() * B;
    
    const double* coeff_data = coeffs.GetMatrixArray(); 

    return vector<double>( coeff_data, coeff_data + DoF );  
}
