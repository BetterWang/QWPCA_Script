#include "label.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include "complex"
using namespace std;

typedef std::complex<double> Complex;
const int NETA = 48;

void process(int s1 = 0, int n = 2)
{
	cout << " s1 = " << s1 << "\tn = " << n << endl;
	//addchain(s1);

	chV->Add("test/pca.root/QWPCA/trV");
	chV->SetMakeClass(1);


	int Cent;
	int Mult;

	Int_t           pQeta2_;
	Int_t           pQeta3_;
	Int_t           pQeta4_;
	Double_t	rQeta2[NETA];
	Double_t	iQeta2[NETA];

	Double_t	wQeta[NETA];

	vector<double>  *pQetaW = 0;

	chV->SetBranchStatus("*", 0);

	chV->SetBranchStatus("cent", 1);
	chV->SetBranchStatus("mult", 1);
	chV->SetBranchStatus("pQetaW", 1);

	chV->SetBranchAddress("cent", &Cent);
	chV->SetBranchAddress("mult", &Mult);
	chV->SetBranchAddress("pQetaW", &pQetaW);

	chV->SetBranchStatus(Form("pQeta%i",n), 1);
	chV->SetBranchStatus(Form("pQeta%i._real",n), 1);
	chV->SetBranchStatus(Form("pQeta%i._imag",n), 1);

	chV->SetBranchAddress(Form("pQeta%i",n), &pQeta2_);
	chV->SetBranchAddress(Form("pQeta%i._real",n), rQeta2);
	chV->SetBranchAddress(Form("pQeta%i._imag",n), iQeta2);

	Complex QetaM2[200][48][48];
	Complex QetaM2w2[200][48][48];

	TH1D * hCent = new TH1D("hCent", "hCent", 200, 0, 200);
	TH1D * hMult = new TH1D("hMult", "hMult", 800, 0, 8000);

	int ievt = 0;
	while (chV->GetEntry(ievt++)) {
		if ( ievt % 100 == 0 ) cout << " ! ievt = " << ievt << endl;
		Complex Qeta2[48];
		hCent->Fill(Cent);
		hMult->Fill(Mult);

		for ( int i = 0; i < NETA; i++ ) {
			Qeta2[i] = Complex(rQeta2[i], iQeta2[i]);
		}
	}
}
