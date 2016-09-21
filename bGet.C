#include "label.h"
#include "noff.h"

typedef std::complex<double> Complex;
const int NETA = 48;
const int NETA2 = 24;
void bGet(int s1 = 1, int N = 2)
{
	cout << " s1 = " << s1 << " N = " << N << endl;
	TH1::SetDefaultSumw2();
	TFile * f = new TFile(Form("%s/output_%i.root", ftxt[s1], N));

	TH1D * hCent = (TH1D*) f->Get("hCent");
	TH1D * hMult = (TH1D*) f->Get("hMult");

	TH1D * hQ1Dr[NETA2];
	TH1D * hQ1Di[NETA2];
	TH1D * hQ1Dw[NETA2];

	TH1D * hQ2Dr[NETA2][NETA2];
	TH1D * hQ2Di[NETA2][NETA2];
	TH1D * hQ2Dw[NETA2][NETA2];

	for ( int i = 0; i < NETA2; i++ ) {
		hQ1Dr[i] = (TH1D*) f->Get(Form("hQ1Dr_%i", i));
		hQ1Di[i] = (TH1D*) f->Get(Form("hQ1Di_%i", i));
		hQ1Dw[i] = (TH1D*) f->Get(Form("hQ1Dw_%i", i));
		for ( int j = i; j < NETA2; j++ ) {
			hQ2Dr[i][j] = (TH1D*) f->Get(Form("hQ2Dr_%i_%i", i, j));
			hQ2Di[i][j] = (TH1D*) f->Get(Form("hQ2Di_%i_%i", i, j));
			hQ2Dw[i][j] = (TH1D*) f->Get(Form("hQ2Dw_%i_%i", i, j));
		}
	}

	// rebin
	TH1D * hCentRebin = (TH1D*) hCent->Rebin(NCent, "hCentRebin", CentBin);

	TH1D * hQ1DrRebin[NETA2];
	TH1D * hQ1DiRebin[NETA2];
	TH1D * hQ1DwRebin[NETA2];

	TH1D * hQ2DrRebin[NETA2][NETA2];
	TH1D * hQ2DiRebin[NETA2][NETA2];
	TH1D * hQ2DwRebin[NETA2][NETA2];

	for ( int i = 0; i < NETA2; i++ ) {
		hQ1DrRebin[i] = (TH1D*) hQ1Dr[i]->Rebin(NCent, Form("hQ1DrRebin_%i", i), CentBin);
		hQ1DiRebin[i] = (TH1D*) hQ1Di[i]->Rebin(NCent, Form("hQ1DiRebin_%i", i), CentBin);
		hQ1DwRebin[i] = (TH1D*) hQ1Dw[i]->Rebin(NCent, Form("hQ1DwRebin_%i", i), CentBin);
		for ( int j = i; j < NETA2; j++ ) {
			hQ2DrRebin[i][j] = (TH1D*) hQ2Dr[i][j]->Rebin(NCent, Form("hQ2DrRebin_%i_%i", i, j), CentBin);
			hQ2DiRebin[i][j] = (TH1D*) hQ2Di[i][j]->Rebin(NCent, Form("hQ2DiRebin_%i_%i", i, j), CentBin);
			hQ2DwRebin[i][j] = (TH1D*) hQ2Dw[i][j]->Rebin(NCent, Form("hQ2DwRebin_%i_%i", i, j), CentBin);
		}
	}

	// average
	for ( int i = 0; i < NETA2; i++ ) {
		hQ1DrRebin[i]->Divide(hCentRebin);
		hQ1DiRebin[i]->Divide(hCentRebin);
		hQ1DwRebin[i]->Divide(hCentRebin);
		for ( int j = i; j < NETA2; j++ ) {
			hQ2DrRebin[i][j]->Divide(hCentRebin);
			hQ2DiRebin[i][j]->Divide(hCentRebin);
			hQ2DwRebin[i][j]->Divide(hCentRebin);

		}
	}

	for ( int i = 0; i < NETA2; i++ ) {
		for ( int j = i; j < NETA2; j++ ) {
			if ( j == i ) {
				hQ2DrRebin[i][j]->Add( hQ1DwRebin[i], -1 );
				hQ2DiRebin[i][j]->Add( hQ1DwRebin[i], -1 );
				hQ2DwRebin[i][j]->Add( hQ1DwRebin[i], -1 );
			}
			for ( int c = 0; c < NCent; c++ ) {
				Complex qi( hQ1DrRebin[i]->GetBinContent(c+1), hQ1DiRebin[i]->GetBinContent(c+1) );
				Complex qj( hQ1DrRebin[j]->GetBinContent(c+1), hQ1DiRebin[j]->GetBinContent(c+1) );
				//cout << " q[" << i << "] = " << qi << endl;
				//cout << " q[" << j << "] = " << qj << endl;
				qi = qi * conj(qj);

				double real = hQ2DrRebin[i][j]->GetBinContent(c+1);
				double imag = hQ2DiRebin[i][j]->GetBinContent(c+1);
				hQ2DrRebin[i][j]->SetBinContent(c+1, real - qi.real());
				hQ2DiRebin[i][j]->SetBinContent(c+1, imag - qi.imag());
			}
		}
	}

	TH2D * hConvR[NCent] = {};
	TH2D * hConvI[NCent] = {};
	TH2D * hConvW[NCent] = {};
	for ( int c = 0; c < NCent; c++ ) {
		hConvR[c] = new TH2D(Form("hConvR_%i", c), "", NETA2, -2.4, 2.4, NETA2, -2.4, 2.4);
		hConvI[c] = new TH2D(Form("hConvI_%i", c), "", NETA2, -2.4, 2.4, NETA2, -2.4, 2.4);
		hConvW[c] = new TH2D(Form("hConvW_%i", c), "", NETA2, -2.4, 2.4, NETA2, -2.4, 2.4);
		for ( int i = 0; i < NETA2; i++ ) {
			for ( int j = i; j < NETA2; j++ ) {
				hConvR[c]->SetBinContent(i+1, j+1, hQ2DrRebin[i][j]->GetBinContent(c+1));
				hConvI[c]->SetBinContent(i+1, j+1, hQ2DiRebin[i][j]->GetBinContent(c+1));
				hConvW[c]->SetBinContent(i+1, j+1, hQ2DwRebin[i][j]->GetBinContent(c+1));

				hConvR[c]->SetBinError(i+1, j+1, hQ2DrRebin[i][j]->GetBinError(c+1));
				hConvI[c]->SetBinError(i+1, j+1, hQ2DiRebin[i][j]->GetBinError(c+1));
				hConvW[c]->SetBinError(i+1, j+1, hQ2DwRebin[i][j]->GetBinError(c+1));
			}

		}
	}

	// save
	TFile * fsave = new TFile(Form("%s/outputMatrix_%i.root", ftxt[s1], N), "recreate");
	for ( int i = 0; i < NETA2; i++ ) {
		hQ1Dr[i]->Write();
		hQ1Di[i]->Write();
		hQ1Dw[i]->Write();
		hQ1DrRebin[i]->Write();
		hQ1DiRebin[i]->Write();
		hQ1DwRebin[i]->Write();
		for ( int j = i; j < NETA2; j++ ) {
			hQ2Dr[i][j]->Write();
			hQ2Di[i][j]->Write();
			hQ2Dw[i][j]->Write();
			hQ2DrRebin[i][j]->Write();
			hQ2DiRebin[i][j]->Write();
			hQ2DwRebin[i][j]->Write();
		}
	}
	for ( int c = 0; c < NCent; c++ ) {
		hConvR[c]->Write();
		hConvI[c]->Write();
		hConvW[c]->Write();
	}
	hCent->Write();
	hCentRebin->Write();
	hMult->Write();
	fsave->Close();
}
