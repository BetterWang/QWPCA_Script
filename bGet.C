#include "label.h"
#include "noff.h"

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

}
