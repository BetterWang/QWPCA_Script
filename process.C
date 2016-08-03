#include "label.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include "complex"
using namespace std;

typedef std::complex<double> Complex;
const int NETA = 48;
const int NETA2 = 24;

void process(int s1 = 0, int N = 2)
{
	cout << " s1 = " << s1 << "\tN = " << N << endl;
	addchain(s1);

	//chV->Add("../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v4/160721_182208/0000/pca_99.root/QWPCA/trV");
	chV->SetMakeClass(1);
	TH1::SetDefaultSumw2();


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

	chV->SetBranchStatus(Form("pQeta%i",N), 1);
	chV->SetBranchStatus(Form("pQeta%i._real",N), 1);
	chV->SetBranchStatus(Form("pQeta%i._imag",N), 1);

	chV->SetBranchAddress(Form("pQeta%i",N), &pQeta2_);
	chV->SetBranchAddress(Form("pQeta%i._real",N), rQeta2);
	chV->SetBranchAddress(Form("pQeta%i._imag",N), iQeta2);


	TH1D * hCent = new TH1D("hCent", "hCent", 200, 0, 200);
	TH1D * hMult = new TH1D("hMult", "hMult", 800, 0, 8000);

	TH1D * hQ1Dr[NETA2];
	TH1D * hQ1Di[NETA2];
	TH1D * hQ1Dw[NETA2];
	for ( int i = 0; i < NETA2; i++ ) {
		hQ1Dr[i] = new TH1D(Form("hQ1Dr_%i", i), "", 200, 0 , 200);
		hQ1Di[i] = new TH1D(Form("hQ1Di_%i", i), "", 200, 0 , 200);
		hQ1Dw[i] = new TH1D(Form("hQ1Dw_%i", i), "", 200, 0 , 200);
	}

	TH1D * hQ2Dr[NETA2][NETA2] = {};
	TH1D * hQ2Di[NETA2][NETA2] = {};
	TH1D * hQ2Dw[NETA2][NETA2] = {};
	for ( int i = 0; i < NETA2; i++ ) {
		for ( int j = i; j < NETA2; j++ ) {
			hQ2Dr[i][j] = new TH1D(Form("hQ2Dr_%i_%i", i, j), "", 200, 0, 200);
			hQ2Di[i][j] = new TH1D(Form("hQ2Di_%i_%i", i, j), "", 200, 0, 200);
			hQ2Dw[i][j] = new TH1D(Form("hQ2Dw_%i_%i", i, j), "", 200, 0, 200);
		}
	}

	int ievt = 0;
	while (chV->GetEntry(ievt++)) {
		if ( ievt % 10000 == 0 ) cout << " ! ievt = " << ievt << endl;
		Complex Qeta2[NETA2] = {};
		double weight[NETA2] = {};
		hCent->Fill(Cent);
		hMult->Fill(Mult);

		for ( int i = 0; i < NETA2; i++ ) {
			Qeta2[i] = Complex(rQeta2[2*i] + rQeta2[2*i+1], iQeta2[2*i] + iQeta2[2*i+1]);
			weight[i] = (*pQetaW)[2*i] + (*pQetaW)[2*i+1];

			hQ1Dr[i]->Fill(Cent+1, Qeta2[i].real());
			hQ1Di[i]->Fill(Cent+1, Qeta2[i].imag());
			hQ1Dw[i]->Fill(Cent+1, weight[i]);
		}

		for ( int i = 0; i < NETA2; i++ ) {
			for ( int j = i; j < NETA2; j++ ) {
				Complex c = Qeta2[i] * conj(Qeta2[j]);
				double w = weight[i] * weight[j];
				hQ2Dr[i][j]->Fill(Cent+1, c.real());
				hQ2Di[i][j]->Fill(Cent+1, c.imag());
				hQ2Dw[i][j]->Fill(Cent+1, w);
			}
		}
	}


	TFile * fsave = new TFile(Form("%s/output_%i.root", ftxt[s1], N), "recreate");
	for ( int i = 0; i < NETA2; i++ ) {
		hQ1Dr[i]->Write();
		hQ1Di[i]->Write();
		hQ1Dw[i]->Write();
		for ( int j = i; j < NETA2; j++ ) {
			hQ2Dr[i][j]->Write();
			hQ2Di[i][j]->Write();
			hQ2Dw[i][j]->Write();
		}
	}
	hCent->Write();
	hMult->Write();
}
