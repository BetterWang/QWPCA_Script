#include <TChain.h>
#include <iostream>
TChain * chV = new TChain("trV");

char const * fname[] = {
	"test/",					// 0
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v4/160721_182208/0000/",				// 1 PbPb15 HIMB2 pixel private 0.3 < pT < 3.0, |eta|<2.4
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_eff_v2/161006_190443/0000/",				// 2 HIMB5 pp reco 0.3 < pT < 3.0
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v6/160922_142548/0000/",				// 3 PbPb15 HIMB2 official pix 0.3 < pT < 3.0 eff
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v7/161006_184849/0000/",				// 4 PbPb15 HIMB2 official pix 0.3 < pT < 3.0 noeff
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v8/161006_185132/0000/",				// 5 PbPb15 HIMB2 official pix 1.0 < pT < 5.0 noeff
};

char const * ftxt[] = {
	"txt/test/",						// 0
	"txt/HIMB2_PCA_pixel_noeff_v4/",			// 1
	"txt/HIMB5_PCA_ppReco_eff_v1/",				// 2
	"txt/HIMB2_PCA_pixel_noeff_v6/",			// 3
	"txt/HIMB2_PCA_pixel_noeff_v7/",			// 4
	"txt/HIMB2_PCA_pixel_noeff_v8/",			// 5
};

void addchain(int s1)
{
//	chV->SetMakeClass(1);
	std::cout << fname[s1] << std::endl;
	chV->Add(Form("%s/*.root/QWPCA/trV", fname[s1]));
}
