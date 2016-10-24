{
	TFile *f0 = new TFile("txt/HIMB2_PCA_pixel_noeff_v6/outputFit_Minuit2__2.root");
	TFile *f1 = new TFile("txt/HIMB2_PCA_pixel_noeff_v7/outputFit_Minuit2__2.root");

	TH1D * hMult0[11];
	TH1D * hMult1[11];
	for ( int i = 0; i < 11; i++ ) {
		hMult0[i] = (TH1D*) f0->Get(Form("hMult_%i", i));
		hMult1[i] = (TH1D*) f1->Get(Form("hMult_%i", i));
		hMult1[i]->SetLineColor(kRed);
	}

	TCanvas *c = new TCanvas("c", "c", 800, 600);
//	c->Print("comp.pdf(", "pdf");
	for ( int i = 0; i < 11; i++ ) {
		hMult0[i]->SetMinimum(0);
		hMult0[i]->Draw();
		hMult1[i]->Draw("same");
		c->Print(Form("comp_%i.pdf", i));
	}

}

