#include "label.h"
#include "noff.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "complex"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"

typedef std::complex<double> Complex;
const int NETA = 48;
const int NETA2 = 24;
const int N_GAP = 10;



class GlobalChi2 {
public:
	GlobalChi2(TH2D *h, int gap = 0, int o = 1, bool s=false) :
		order(o), sym(s), d_gap(gap) {

		for ( int i = 0; i < NETA2; i++ ) {
			for ( int j = i+d_gap; j < NETA2; j++ ) {
				covMat[i][j] = h->GetBinContent(i+1, j+1);
				d_covMat[i][j] = h->GetBinError(i+1, j+1);
			}
		}
	}
	double operator() ( const double * );

private:
	double covMat[NETA2][NETA2];
	double d_covMat[NETA2][NETA2];
	int order;
	bool sym;
	int d_gap;
};

double GlobalChi2::operator() ( const double * x ) {
	double x1[NETA2] = {};
	double x2[NETA2] = {};
	double x3[NETA2] = {};
	if ( sym ) {
		for ( int i = 0; i < NETA2/2; i++ ) {
			x1[i] = x[i];
			x1[NETA2 - i -1] = x[i];
			if ( order > 1 ) {
				x2[i] = x[ NETA2/2 + i ];
				x2[NETA2 - i -1] = - x[ NETA2/2 + i ];
			}
			if ( order > 2 ) {
				x3[i] = x[ NETA2 + i ];
				x3[NETA2 - i -1] = x[ NETA2 + i ];
			}
		}
	} else {
		for ( int i = 0; i < NETA2; i++ ) {
			x1[i] = x[i];
			if ( order > 1 ) {
				x2[i] = x[NETA2+i];
			}
			if ( order > 2 ) {
				x3[i] = x[2*NETA2+i];
			}
		}
	}

//	cout << " -> order = " << order << " sym = " << sym << " d_gap = " << d_gap << endl;
//	for ( int i = 0; i < NETA2; i++ ) {
//		cout << "   ->  i = " << i << "\tx1 = " << x1[i] << "\tx2 = " << x2[i] << "\tx3 = " << x3[i] << endl;
//	}

	int bins = 0;
	double chi2 = 0;
	for ( int i = 0; i < NETA2; i++ ) {
		for ( int j = i + d_gap; j < NETA2; j++ ) {
			chi2 += pow( ( x1[i]*x1[j] +x2[i]*x2[j] + x3[i]*x3[j] - covMat[i][j] ) / (d_covMat[i][j]), 2 );
			bins++;
		}
	}
	return chi2 / bins;
}

   // algorithm
   // possible choices are:
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic

void bFit(int s1 = 1, int N = 2, const char * minName = "Minuit2", const char *algoName = "")
{
	cout << " s1 = " << s1 << " N = " << N << endl;

	TH1::SetDefaultSumw2();
	TFile * f = new TFile(Form("%s/outputMatrix_%i.root", ftxt[s1], N));

	TH2D * hConvR[NCent] = {};
	TH2D * hConvI[NCent] = {};
	TH2D * hConvW[NCent] = {};

	for ( int c = 0; c < NCent; c++ ) {
		hConvR[c] = (TH2D*) f->Get(Form("hConvR_%i", c));
		hConvI[c] = (TH2D*) f->Get(Form("hConvI_%i", c));
		hConvW[c] = (TH2D*) f->Get(Form("hConvW_%i", c));
	}

	double vn_start[NCent][NETA2] = {};
	for ( int c = 0; c < NCent; c++ ) {
		for ( int i = 0; i < NETA2; i++ ) {
			double V = hConvR[c]->GetBinContent(i+1, i+1);
			double W = hConvW[c]->GetBinContent(i+1, i+1);
			//vn_start[c][i] = sqrt(V/W);
			vn_start[c][i] = sqrt(V);
		}
	}

//	TH1D * hVn[3][3][NCent] = {};
//	TH1D * hVnSym[3][3][NCent] = {};
//	for ( int c = 0; c < NCent; c++ ) {
//		hVn[c] = new TH1D(Form("hVn_%i", c), "", NETA2, -2.4, 2.4);
//		for ( int i = 0; i < NETA2; i++ ) {
//			hVn[c]->SetBinContent(i+1, vn_start[c][i]);
//		}
//	}

	// prepare fit
	ROOT::Math::Minimizer* pmin;
	pmin = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
	pmin->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
	pmin->SetMaxIterations(10000);  // for GSL
	pmin->SetTolerance(0.001);
	pmin->SetPrintLevel(0);

	for ( int c = 0; c < NCent; c++ ) {
		cout << " ! Fitting for cent = " << c << endl;

		GlobalChi2 gChi2_1(hConvR[c], N_GAP, 1, false);
		GlobalChi2 gChi2_2(hConvR[c], N_GAP, 2, false);
		GlobalChi2 gChi2_3(hConvR[c], N_GAP, 3, false);
		GlobalChi2 gChi2_1sym(hConvR[c], N_GAP, 1, true);
		GlobalChi2 gChi2_2sym(hConvR[c], N_GAP, 2, true);
		GlobalChi2 gChi2_3sym(hConvR[c], N_GAP, 3, true);

		ROOT::Math::WrappedMultiFunction<GlobalChi2 &> wf_1(gChi2_1, NETA2);
		ROOT::Math::WrappedMultiFunction<GlobalChi2 &> wf_2(gChi2_2, NETA2*2);
		ROOT::Math::WrappedMultiFunction<GlobalChi2 &> wf_3(gChi2_3, NETA2*3);
		ROOT::Math::WrappedMultiFunction<GlobalChi2 &> wf_1sym(gChi2_1sym, NETA2/2);
		ROOT::Math::WrappedMultiFunction<GlobalChi2 &> wf_2sym(gChi2_2sym, NETA2);
		ROOT::Math::WrappedMultiFunction<GlobalChi2 &> wf_3sym(gChi2_3sym, NETA2*3/2);

		/////
		pmin->Clear();
		pmin->SetFunction(wf_1);
		for ( int i = 0; i < NETA2; i++ ) {
			pmin->SetVariable(i, Form("vn_%i", i), vn_start[c][i], 0.0001);
		}
		pmin->Minimize();
		cout << "    ! Norder = 1, no sym" << endl;
		pmin->PrintResults();

		/*
		/////
		pmin->Clear();
		pmin->SetFunction(wf_1sym);
		for ( int i = 0; i < NETA2/2; i++ ) {
			pmin->SetVariable(i, Form("vn_%i", i), vn_start[c][i], 0.0001);
		}
		pmin->Minimize();
		cout << "    ! Norder = 1, sym" << endl;
		pmin->PrintResults();

		/////
		pmin->Clear();
		pmin->SetFunction(wf_2);
		int idx = 0;
		for ( int i = 0; i < NETA2; i++ ) {
			pmin->SetVariable(idx, Form("vn1_%i", i), vn_start[c][i], 0.0001);
			idx++;
		}
		for ( int i = 0; i < NETA2/2; i++ ) {
			pmin->SetVariable(idx, Form("vn2_%i", i), -vn_start[c][i]/10, 0.00001);
			idx++;
		}
		for ( int i = NETA2/2; i < NETA2; i++ ) {
			pmin->SetVariable(idx, Form("vn2_%i", i), vn_start[c][i]/10, 0.00001);
			idx++;
		}
		pmin->Minimize();
		cout << "    ! Norder = 2, no sym" << endl;
		pmin->PrintResults();

		/////
		pmin->Clear();
		pmin->SetFunction(wf_2sym);
		idx = 0;
		for ( int i = 0; i < NETA2/2; i++ ) {
			pmin->SetVariable(idx, Form("vn1_%i", i), vn_start[c][i], 0.0001);
			idx++;
		}
		for ( int i = 0; i < NETA2/2; i++ ) {
			pmin->SetVariable(idx, Form("vn2_%i", i), -vn_start[c][i]/10, 0.00001);
			idx++;
		}
		pmin->Minimize();
		cout << "    ! Norder = 2, sym" << endl;
		pmin->PrintResults();
		*/
	}



	// save
//	TFile * fsave = new TFile(Form("%s/outputFit_%i.root", ftxt[s1], N), "recreate");
//	for ( int c = 0; c < NCent; c++ ) {
//		hVn[c]->Write();
//	}
}
