#include <TChain.h>
#include <iostream>
TChain * chV = new TChain("trV");

char const * fname[] = {
	"test/",					// 0
	"../PbPb2015/HIMinimumBias2/crab_HIMB2_PCA_pixel_noeff_v1/160705_122619/0000/",				// 1 PbPb15 HIMB2 pixel private 0.3 < pT < 3.0, |eta|<2.4

};

char const * ftxt[] = {
	"txt/test/",						// 0
	"txt/HIMB2_PCA_pixel_noeff_v1/",			// 1
};

void addchain(int s1)
{
//	chV->SetMakeClass(1);
	std::cout << fname[s1] << std::endl;
	chV->Add(Form("%s/*.root/QWPCA/trV", fname[s1]));
}
