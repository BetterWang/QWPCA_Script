#include <TChain.h>
#include <iostream>
TChain * chV = new TChain("trV");

char const * fname[] = {
	"test/",					// 0
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v4/160721_182208/0000/",				// 1 PbPb15 HIMB2 pixel private 0.3 < pT < 3.0, |eta|<2.4
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_eff_v3/161008_132915/0000/",				// 2 HIMB5 pp reco 0.3 < pT < 3.0
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_eff_v9/161017_094117/0000/",				// 3 PbPb15 HIMB2 official pix 0.3 < pT < 3.0 eff
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v7/161006_184849/0000/",				// 4 PbPb15 HIMB2 official pix 0.3 < pT < 3.0 noeff
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v8/161006_185132/0000/",				// 5 PbPb15 HIMB2 official pix 1.0 < pT < 5.0 noeff
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_noeff_v1/161010_144827/0000/",				// 6 HIMB5 pp reco 0.3 < pT < 3.0 noeff
	"../PCA/Hydjet_Quenched_MinBias_5020GeV_750/crab_Hydjet_PCA_GEN_v1/161012_131811/0000/",		// 7 Hydjet GEN 0.3 < pT < 3.0 noeff
	"../PCA/Hydjet_Quenched_MinBias_5020GeV_750/crab_Hydjet_PCA_pixel_eff_v2/161014_160213/0000/",		// 8 Hydjet pixel 0.3 < pT < 3.0 eff
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_eff_noff_v5/161023_162336/0000/",				// 9 PbPb15 HIMB5 noff 0.3 < pT < 3.0 eff
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_eff_noff_pT0550_v5/161023_162358/0000/",			// 10 PbPb15 HIMB5 noff 0.5 < pT < 5.0 eff
	"../PCA/PAHighPt/crab_pPb5_PCA_eff_noff_v1/161020_073626/0000/",					// 11 pPb5 HM noff 0.3 < pT < 3.0 eff
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_eff_cent_v9/161020_142532/0000/",				// 12 PbPb15 HIMB2 pixel 0.3 < pT < 3.0 eff
	"../PCA/HIMinimumBias2/crab_HIMB2_PCA_pixel_eff_cent_pT0550_v9/161020_174950/0000/",			// 13 PbPb15 HIMB2 pixel 0.5 < pT < 5.0 eff
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_noeff_noff_v5/161023_162438/0000/",			// 14 PbPb15 HIMB5 0.3 < pT < 3.0 noeff
	"../PCA/HIMinimumBias5/crab_HIMB5_PCA_ppReco_noeff_noff_pT0550_v5/161023_162454/0000/",			// 15 PbPb15 HIMB5 0.5 < pT < 5.0 noeff
};

char const * ftxt[] = {
	"txt/test/",						// 0
	"txt/HIMB2_PCA_pixel_noeff_v4/",			// 1
	"txt/HIMB5_PCA_ppReco_eff_v1/",				// 2
	"txt/HIMB2_PCA_pixel_noeff_v6/",			// 3
	"txt/HIMB2_PCA_pixel_noeff_v7/",			// 4
	"txt/HIMB2_PCA_pixel_noeff_v8/",			// 5
	"txt/HIMB5_PCA_ppReco_noeff_v1/",			// 6
	"txt/Hydjet_PCA_GEN_noeff_v1/",				// 7
	"txt/Hydjet_PCA_pixel_eff_v1/",				// 8
	"txt/HIMB5_PCA_ppReco_eff_noff_v4/",			// 9    new code here
	"txt/HIMB5_PCA_ppReco_eff_noff_pT0550_v4/",		// 10
	"txt/HM_pPb5_eff_noff/",				// 11
	"txt/HIMB2_PCA_pixel_eff_cent_v9/",			// 12
	"txt/HIMB2_PCA_pixel_eff_cent_pT0550_v9/",		// 13
	"txt/HIMB5_PCA_ppReco_noeff_noff/",			// 14
	"txt/HIMB5_PCA_ppReco_noeff_noff_pT0550/",		// 15
};

void addchain(int s1)
{
//	chV->SetMakeClass(1);
	std::cout << fname[s1] << std::endl;
	chV->Add(Form("%s/*.root/QWPCA/trV", fname[s1]));
}
