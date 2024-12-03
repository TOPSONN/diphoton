#include<cstdio>
#include<cstdlib>
#include<iostream>
#include <vector>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TEnv.h"
#include "TClass.h"
#include "TMath.h"
#include "THashList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TF3.h"
#include <algorithm>
#include <functional>
#include <fstream>
#include "Math/Vector4D.h"
#include "TClass.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;
#include <math.h>
#include <map>
#include <cstdio>
#include <stdio.h>
#include <string>
#include <cctype>
#include <algorithm>
#include <iterator>
#include <pthread.h>
#include <thread>
#include "TStopwatch.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TLegend.h>
#include <time.h>
#include <TStyle.h>
#include <TColor.h>
#include <cmath> // For log()

const int nHists = 30;
// Create an array of TH1F histograms
TH1F *photonpt1_nocut[nHists] , *deltaetaj1j2_nocut[nHists];
TH1F *photonpt2_nocut[nHists] , *invmassj1j2_nocut[nHists]  ;

TH1F *photoneta1_cut[nHists] , *photoneta2_cut[nHists] ,  *jeteta1_cut[nHists] , *jeteta2_cut[nHists] ;
TH1F *photonpt1_cut[nHists] , *photonpt2_cut[nHists] ,  *jetpt1_cut[nHists] , *jetpt2_cut[nHists] ;
TH1F *photondeltaeta[nHists] ,   *jetdeltaeta[nHists]  ;
TH1F *photondeltaphi[nHists] , *photonphi2_cut[nHists] ,  *jetdeltaphi[nHists]  ;
TH1F *dRp1p2[nHists] , *dRj1j2[nHists] , *dRp1j1[nHists] , *dRp1j2[nHists]  , *dRp2j1[nHists] , *dRp2j2[nHists] ;
TH1F *mpp_cut[nHists] , *mjj_cut[nHists] ;
TH1F *Zi_origin[nHists] , *DPHI_origin[nHists] ;

TH1F *photonIDMVA1_origin[nHists] , *photonIDMVA2_origin[nHists];

TH1F *photonpt1_cut_all[nHists] , *photonpt2_cut_all[nHists] ,  *jetpt1_cut_all[nHists] , *jetpt2_cut_all[nHists] ;
TH1F *photondeltaeta_cut_all[nHists] , *photoneta2_cut_all[nHists] ,  *jetdeltaeta_cut_all[nHists] , *jeteta2_cut_all[nHists] ;
TH1F *photondeltaphi_cut_all[nHists] , *photonphi2_cut_all[nHists] ,  *jetdeltaphi_cut_all[nHists] , *jetphi2_cut_all[nHists] ;
TH1F *mpp_cut_all[nHists] , *mjj_cut_all[nHists] ;

TH1F *data[nHists] , *hRatio[nHists] , *hback_total[nHists] , *hSignificance[nHists] ;

THStack *DPHI_origin_stack = new THStack("DPHI_origin_stack", "DPHI_origin_stack");
THStack *Zi_origin_stack = new THStack("Zi_origin_stack", "Zi_origin_stack");
THStack *photonIDMVA1_origin_stack = new THStack("photonIDMVA1_origin_stack", "photonIDMVA1_origin_stack");
THStack *photonIDMVA2_origin_stack = new THStack("photonIDMVA2_origin_stack", "photonIDMVA2_origin_stack");

THStack *stack_plot[nHists];


TPad *topPad[nHists];
TPad *bottomPad[nHists];

//////////////Global variable/////////////
float array_eff[25]= { 0 };
float array_eff_photonpt[25] = {0};
const double qcd50to100_exs = 11062410000000 , qcd100to200_exs = 1408323000000  , qcd200to300_exs =  92594700000 , qcd300to500_exs = 19306980000 ;
const double qcd500to700_exs = 1799358000 , qcd700to1000_exs = 378736800 , qcd1000to1500_exs = 65192400 , qcd1500to2000_exs = 5955672, qcd2000toInf_exs = 1214895 ;
const double vbs_exs = 5320.464  ,    inter_exs = 868.038  ,  qcddiphoton_exs = 2327106 ;
const double gj40to100_exs = 1242954000 , gj100to200_exs = 552165300, gj200to400_exs = 137608500 ;
const double gj400to600_exs = 16435410 , gj600toInf_exs = 5576577 ;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Signi_log(float array[] , double &sig) {
	double sunofqcd[30] = {0} ;
	double sunofgj[30] = {0} ;
	double signal =  0   ;
	double back_qcd =  0 ;
	double back_qcddiphotonjet =  0 ;
	double back_gjet =  0  ;
	double back_inter = 0 ;
	double back_totalsum = 0 ;

	//loopqcd
	sunofqcd[0] = array[16]*qcd50to100_exs  ;//50to100
	sunofqcd[1] = array[15]*qcd100to200_exs ;//100
	sunofqcd[2] = array[14]*qcd200to300_exs ;//200
	sunofqcd[3] = array[13]*qcd300to500_exs ;//300
	sunofqcd[4] = array[12]*qcd500to700_exs ;//500
	sunofqcd[5] = array[11]*qcd700to1000_exs ;//700
	sunofqcd[6] = array[10]*qcd1000to1500_exs ;//1000
	sunofqcd[7] = array[9]*qcd1500to2000_exs ;//1500
	sunofqcd[8] = array[8]*qcd2000toInf_exs ;//2000toInf

	for(int i = 0 ; i<10 ; i++  ) {
		back_qcd = back_qcd + sunofqcd[i] ;
	}

	//loopgj
	sunofgj[0] = array[7]*gj40to100_exs  ;//40to100
	sunofgj[1] = array[6]*gj100to200_exs ;//100
	sunofgj[2] = array[5]*gj200to400_exs ;//200
	sunofgj[3] = array[4]*gj400to600_exs ;//400
	sunofgj[4] = array[3]*gj600toInf_exs ;//600

	for(int i = 0 ; i<10 ; i++  ) {
		back_gjet = back_gjet + sunofgj[i] ;
	}

	back_inter=		array[1]*inter_exs 				;//inter
	back_qcddiphotonjet= array[2]*qcddiphoton_exs 			;//qcd diphotonjets
	signal=	array[0]*vbs_exs						;//vbs
	back_totalsum = back_qcd + back_qcddiphotonjet + back_gjet + back_inter ;
	sig = sqrt(2*((signal+back_totalsum)*log(1+(signal/back_totalsum))-signal));

	return sig ;
}




float Zi(float jet_eta1, float jet_eta2, float photon_eta1, float photon_eta2) { //Z = [(eta_pho1+eta_pho2) -(eta_j1+eta_j2) ] /2.
	float Z = ((photon_eta1+photon_eta2)-((jet_eta1+jet_eta2)/2));
	return abs(Z) ;
}

float DPHI(float jet_phi1, float jet_phi2, float photon_phi1, float photon_phi2) {
	float photon_phi = photon_phi1 + photon_phi2 ;
	float jet_phi= jet_phi1 + jet_phi2 ;
	float dphi = (photon_phi-jet_phi)/2 ;

	if (dphi > M_PI) dphi -= 2 * M_PI;
	if (dphi < -M_PI) dphi += 2 * M_PI;
	return abs(dphi) ;
}


TCanvas* Frame(TCanvas* canvas , TH1F* frame , THStack* cut_stack  ,float Min_h ,  float Max_h) {
	canvas->cd(1);
	frame->SetMinimum(Min_h);
	frame->SetMaximum(Max_h);
	frame->SetStats(0);
	frame->Draw();
	cut_stack->Draw("hist,same");
	cut_stack->GetYaxis()->SetRangeUser(Min_h , Max_h );
	gPad->SetLogy();
	return canvas;
}

TCanvas* Frame_with_pad(TCanvas* canvas , TH1F* frame , THStack* cut_stack  ,float Min_h ,  float Max_h , TPad* top) {
	top->cd();
	frame->SetMinimum(Min_h);
	frame->SetMaximum(Max_h);
	frame->SetStats(0);
	frame->Draw();
	cut_stack->Draw("hist,same");
	cut_stack->GetYaxis()->SetRangeUser(Min_h , Max_h );
	gPad->SetLogy();
	return canvas;
}


TCanvas* Add_legend( TLegend* legend,  TCanvas* canvas , TH1F* QCD , TH1F* GJETS , TH1F* QCD_Diphoton , TH1F* VBS , TH1F* INTER ) {
	legend = new TLegend(0.7, 0.7, 0.9, 0.9);
	legend->SetTextFont(42);                             // Helvetica font (standard for CMS)
	legend->SetTextSize(0.04);                           // CMS preferred text size
	legend->SetBorderSize(0);                            // No border
	legend->SetFillStyle(0);                             // Transparent background
	legend->SetMargin(0.2);
	legend->SetEntrySeparation(0.02);
	legend->AddEntry(QCD, "QCD", "F");
	legend->AddEntry(GJETS, "G+JETS", "F");
	legend->AddEntry(QCD_Diphoton, "QCD_DiphotonJets", "F");
	legend->AddEntry(INTER, "Interference", "F");
	legend->AddEntry(VBS, "VBS", "F");
	legend->Draw();
	return canvas ;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {
	float dphi = phi1 - phi2;
	if (dphi > M_PI) dphi -= 2 * M_PI;
	if (dphi < -M_PI) dphi += 2 * M_PI;
	float deta = eta1 - eta2;
	return std::sqrt(deta * deta + dphi * dphi);
}

// Struct to store gen-matched photon info
struct MatchedPhoton {
	int recoIndex;  // Index in the reconstructed Photon collection
	float pt;       // Transverse momentum of the reco photon
	float eta;      // Eta of the reco photon
	float phi;      // Phi of the reco photon
	float m ;
	float r9 ;
	float hoe ;
	float phiso;
	float chiso;
	float sieie;
	float mvaID  ;
	int phocutBased  ;
};


struct MatchedObject {
	int index;  // Index in the original collection (reco or gen)
	float pt;   // Transverse momentum
	float eta;  // Eta
	float phi;  // Phi
	float m ;
	int jetid ;
	int jetpuid ;
};


int Filter_first(int cut_check , int &check_condition  , double QCD_phopt1_flt , double QCD_phopt2_flt , double QCD_jetpt1_flt , double QCD_jetpt2_flt  ) {
	if( cut_check ==1 ) { //photonpt
		if( QCD_jetpt1_flt > 100 && QCD_jetpt2_flt > 30 && QCD_phopt1_flt > 10 &&  QCD_phopt2_flt > 10 ) {
			check_condition = check_condition +1;
		}
	}
	if( cut_check ==2 ) { //jetpt
		if( QCD_jetpt1_flt > 10 && QCD_jetpt2_flt > 10 && QCD_phopt1_flt > 40 &&  QCD_phopt2_flt > 25 ) {
			check_condition = check_condition +1;
		}
	}
	if( cut_check >=3 ) { //
		if( QCD_jetpt1_flt > 100 && QCD_jetpt2_flt > 30 && QCD_phopt1_flt > 40 &&  QCD_phopt2_flt > 25 ) {
			check_condition = check_condition +1;
		}
	}

	return check_condition ;
}

int Filter_second(int cut_check , int &check_condition , double etaj1 ,double etaj2 ,double Mjj , double Mpp , double deta_jj ,double etap1 , double etap2) {

	if( cut_check < 3 ) { //pt
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7 &&  Mjj > 500  && Mpp > 10 && abs(deta_jj) > 2.5 &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 3 ) { //abs(etaj1)
		if(  abs(etaj2) < 4.7 &&  Mjj > 500  && Mpp > 10 && abs(deta_jj) > 2.5 &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 4 ) { //abs(etaj2)
		if(  abs(etaj1) < 4.7  &&  Mjj > 500  && Mpp > 10 && abs(deta_jj) > 2.5 &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 5 ) { //Mjj
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7  && Mpp > 10 && abs(deta_jj) > 2.5 &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 6 ) { //Mpp
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7 &&  Mjj > 500   && abs(deta_jj) > 2.5 &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 7 ) { // deta_jj
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7 &&  Mjj > 500  && Mpp > 10  &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 8 ) { //etap1
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7 &&  Mjj > 500  && Mpp > 10 && abs(deta_jj) > 2.5  && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 9 ) { //etap2
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7 &&  Mjj > 500  && Mpp > 10 && abs(deta_jj)> 2.5 &&  abs(etap1) < 2.5  &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check >= 10 ) { //cut all
		if( abs(etaj1) < 4.7 && abs(etaj2) < 4.7 &&  Mjj > 500  && Mpp > 10 && abs(deta_jj)> 2.5 &&  abs(etap1) < 2.5 && abs(etap2) < 2.5 &&  etaj1*etaj2 < 0 ) {
			check_condition = check_condition +1;
		}
	}



	return check_condition ;
}

int Filter_third(int cut_check , int &check_condition  ,  double dRj1j2 , double dRp1j1 , double dRp1j2 , double dRp2j1  , double dRp2j2  ) {

	if( cut_check <10 ) { //cut all
		if( dRj1j2 > 0.5 && dRp1j1 > 0.5 && dRp1j2 > 0.5 && dRp2j1 > 0.5 && dRp2j2>0.5   ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check ==10 ) { //dRj1j2
		if(   dRp1j1 > 0.5 && dRp1j2 > 0.5 && dRp2j1 > 0.5 && dRp2j2>0.5   ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check ==11 ) { //dRp1j1
		if( dRj1j2 > 0.5 &&  dRp1j2 > 0.5 && dRp2j1 > 0.5 && dRp2j2>0.5  ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check ==12 ) { //dRp1j2
		if( dRj1j2 > 0.5 && dRp1j1 > 0.5  && dRp2j1 > 0.5 && dRp2j2>0.5  ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check ==13 ) { //dRp2j1
		if( dRj1j2 > 0.5 && dRp1j1 > 0.5 && dRp1j2 > 0.5 &&  dRp2j2>0.5   ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check ==14 ) { //dRp2j2
		if( dRj1j2 > 0.5 && dRp1j1 > 0.5 && dRp1j2 > 0.5 && dRp2j1 > 0.5   ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check >14 ) { //cut_all
		if( dRj1j2 > 0.5 && dRp1j1 > 0.5 && dRp1j2 > 0.5 && dRp2j1 > 0.5 && dRp2j2>0.5  ) {
			check_condition = check_condition +1;
		}
	}


	return check_condition ;

}


int Filter_fourth(int cut_check , int &check_condition  , double Z_ , double DPHI_ ) {
	if( cut_check < 15 ) { //cut_all
		if( Z_ < 4 && DPHI_ > 1  ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 15 ) { //Z_
		if(  DPHI_ > 1  ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 16 ) { //DPHI_
		if( Z_ < 4 ) {
			check_condition = check_condition +1;
		}
	}
	return check_condition;
}


int ID_check(int &ID,int QCD_jetid1_flt ,int QCD_jetid2_flt ,int QCD_jetpuid1_flt ,int QCD_jetpuid2_flt , float QCD_jetbtagDeepB1_flt , float QCD_jetbtagDeepB2_flt , int QCD_phocutbased1_flt , int QCD_phocutbased2_flt , int &pass_entry) {
//QCD_jetpuid1_flt == 7 && QCD_jetpuid2_flt == 7
	if(  QCD_jetid1_flt == 6 && QCD_jetid2_flt ==6 && QCD_jetpuid1_flt == 7 &&  QCD_jetpuid2_flt == 7 && QCD_phocutbased1_flt ==3 && QCD_phocutbased2_flt ==3   ) {
		ID = 1 ;
		pass_entry = pass_entry + 1 ;
	} else {
		ID = 0 ;
	}

	return ID;
}

TH1F* plot_setting(TH1F* plot , double w , string color ) {

	if( color == "kRed+1") {
		plot->Sumw2();
		plot->SetStats(0);
		plot->SetLineColorAlpha(kRed+1 , 1);
		plot->SetFillColorAlpha(kRed+1 , 1);
		plot->SetFillStyle(4000);
		plot->Scale(w);
	}

	if( color == "kGreen+3") {
		plot->Sumw2();
		plot->SetStats(0);
		plot->SetLineColorAlpha(kGreen+3 , 0.9);
		plot->SetFillColorAlpha(kGreen+3 , 0.9);
		plot->SetFillStyle(4000);
		plot->Scale(w);
	}

	if( color == "kAzure") {
		plot->Sumw2();
		plot->SetStats(0);
		plot->SetLineColorAlpha(kAzure , 0.9);
		plot->SetFillColorAlpha(kAzure , 0.9);
		plot->SetFillStyle(4000);
		plot->Scale(w);
	}

	if( color == "kCyan") {
		plot->Sumw2();
		plot->SetStats(0);
		plot->SetLineColorAlpha(kCyan , 1);
		plot->SetFillColorAlpha(kCyan , 1);
		plot->SetFillStyle(4000);
		plot->Scale(w);
	}

	if( color == "kYellow") {
		plot->Sumw2();
		plot->SetStats(0);
		plot->SetLineColorAlpha(kYellow , 1);
		plot->SetFillColorAlpha(kYellow , 1);
		plot->SetFillStyle(4000);
		plot->Scale(w);
	}

	return plot ;

}

void Add_to_stack(TH1F* variable_plot , THStack* stack_qcd_plot , TCanvas* canvas_plot ,int i ) {
	float Max_h = 1E5 ;
	float Min_h = 1E-5 ;
	variable_plot->GetYaxis()->SetRangeUser(Min_h,Max_h);
	stack_qcd_plot->Add(variable_plot);
	gPad->SetLogy();

	//return canvas_plot ;
}

void Set_TPad(TPad* top , TPad* bot ,TCanvas* canvas ) {

	top->SetBottomMargin(0.05);
	bot->SetTopMargin(0.03);
	bot->SetBottomMargin(0.15);
	canvas->cd(1);
	top->Draw();
	bot->Draw();

}

TCanvas* Set_Ratio(TCanvas* canvas  , TH1F* ratio , TPad* bot ){

bot->cd();
ratio->SetMarkerStyle(20);
ratio->SetMarkerColor(kBlack);
ratio->SetTitle("Significance");
ratio->GetYaxis()->SetTitle("Significance");
ratio->GetYaxis()->SetTitleSize(1.2);
ratio->GetYaxis()->SetTitleOffset(1);

//ratio->GetXaxis()->SetTitle("Observable");
ratio->GetXaxis()->SetTitleSize(0.1);
ratio->GetXaxis()->SetTitleOffset(1.0);

//ratio->GetYaxis()->SetRangeUser(0.00000221632 , 1);
gPad->SetLogy();
double min = ratio->GetMinimum();
double max = ratio->GetMaximum();

// Optionally, expand the range slightly for better visualization
double rangeMargin = 0.1 * (max - min);
ratio->SetMinimum(0.1*min);
ratio->SetMaximum(0.0001*max);

ratio->Draw();

return canvas ; 


}



void Cal_Sig(TH1F* VBS , TH1F* BACK , TH1F* SIG){

// Calculate significance for each bin
    for (int bin = 1; bin <= VBS->GetNbinsX(); ++bin) {
        double dataVal = VBS->GetBinContent(bin);
        double mcVal = BACK->GetBinContent(bin);

        if (mcVal > 0) { // Avoid division by zero
            double significance = (dataVal) / sqrt(mcVal+dataVal);
            SIG->SetBinContent(bin, significance);
            cout << "************************"<<endl;
            cout << "sig=  " << significance   << endl;
            cout<< "******************************" <<endl;
        } else {
            SIG->SetBinContent(bin, 0); // Default for MC = 0
        }
    }


}










void process_HLT(string string_inputroot_process , string string_inputtree_process , double weighted_process , int indexHT_process , int cut_check , string color) {


	TFile *file = new TFile(string_inputroot_process.c_str());
	//TFile *file = TFile::Open(string_inputroot_process);
	TTree *tree;
	file->GetObject(string_inputtree_process.c_str(),tree);
	//tree->Print();
	cout << string_inputroot_process << "  ||  " << string_inputtree_process <<endl;
// std::vector<float> *vec = nullptr;
	std::vector<int> *QCD_jetid= nullptr, *QCD_jetpuid= nullptr ,  *QCD_phoid= nullptr , *QCD_phocutBased= nullptr ;
	std::vector<float> *QCD_phopt= nullptr, *QCD_phoeta= nullptr, *QCD_phophi= nullptr, *QCD_phom= nullptr;
	std::vector<float> *QCD_jetpt= nullptr, *QCD_jeteta= nullptr, *QCD_jetphi= nullptr, *QCD_jetm= nullptr;
	std::vector<float> *QCD_Photon_mvaID= nullptr, *QCD_jetbtagDeepB = nullptr;
	std::vector<float> *QCD_r9= nullptr, *QCD_hoe= nullptr, *QCD_phiso= nullptr, *QCD_chiso= nullptr, *QCD_sieie= nullptr ;


	tree->SetBranchAddress("Photon_pt",   &QCD_phopt  );
	tree->SetBranchAddress("Photon_eta",   &QCD_phoeta  );
	tree->SetBranchAddress("Photon_phi",   &QCD_phophi  );
	tree->SetBranchAddress("Photon_mass",   &QCD_phom  );

	tree->SetBranchAddress("Photon_r9",   &QCD_r9  );
	tree->SetBranchAddress("Photon_hoe",   &QCD_hoe  );
	tree->SetBranchAddress("Photon_pfRelIso03_all",   &QCD_phiso  );
	tree->SetBranchAddress("Photon_pfRelIso03_chg",   &QCD_chiso  );
	tree->SetBranchAddress("Photon_sieie",   &QCD_sieie  );

	tree->SetBranchAddress("Jet_pt",   &QCD_jetpt  );
	tree->SetBranchAddress("Jet_eta",   &QCD_jeteta  );
	tree->SetBranchAddress("Jet_phi",   &QCD_jetphi  );
	tree->SetBranchAddress("Jet_mass",   &QCD_jetm  );

	tree->SetBranchAddress("Jet_jetId",   &QCD_jetid  );
	tree->SetBranchAddress("Photon_pdgId",   &QCD_phoid  );
	tree->SetBranchAddress("Photon_mvaID",   &QCD_Photon_mvaID  );
	
	if (strcmp(string_inputroot_process.c_str(), "2qcd50to100_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd100to200_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd200to300_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd300to500_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd500to700_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd700to1000_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd1000to1500_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd1500to2000_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "2qcd2000toInf_filter.root") == 0) tree->SetBranchAddress("Photon_cutBased",   &QCD_phocutBased  );
	
	if (strcmp(string_inputroot_process.c_str(), "GJET_40to100_.root") == 0) tree->SetBranchAddress("Photon_cutBasedBitmap",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "GJET_100to200_.root") == 0) tree->SetBranchAddress("Photon_cutBasedBitmap",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "GJET_200to400_.root") == 0) tree->SetBranchAddress("Photon_cutBasedBitmap",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "GJET_400to600_.root") == 0) tree->SetBranchAddress("Photon_cutBasedBitmap",   &QCD_phocutBased  );
	if (strcmp(string_inputroot_process.c_str(), "GJET_600toInf_.root") == 0) tree->SetBranchAddress("Photon_cutBasedBitmap",   &QCD_phocutBased  );
	
	tree->SetBranchAddress("Jet_puId",   &QCD_jetpuid  );
	tree->SetBranchAddress("Jet_btagDeepB",   &QCD_jetbtagDeepB  );

	TLorentzVector photon1_t , photon2_t , jet1_t , jet2_t ;

	int total_events = tree->GetEntries() ;
	cout << " total entry=  "   <<  total_events  << endl;
	int pass_entry  = 0 ;
	int pass_entry_photonpt = 0 ;

//begin event loop
	for(int evt= 0   ; evt < total_events  ; evt++ ) {
		int check_condition  = 0 ;
		int ID = 0 ;

		int pho_size =0 , jet_size =0 , QCD_jetpuid1_flt =0 ,  QCD_jetpuid2_flt =0 ;
		int QCD_phoid1_flt =0 ,  QCD_phoid2_flt =0 , QCD_jetid1_flt =0 ,  QCD_jetid2_flt =0 ;
		int QCD_phocutbased1_flt =0 ,  QCD_phocutbased2_flt =0 ;

		float QCD_phopt1_flt = 0 , QCD_phopt2_flt = 0 , QCD_phoeta1_flt = 0 , QCD_phoeta2_flt = 0;
		float QCD_phophi1_flt = 0 , QCD_phophi2_flt = 0 , QCD_phom1_flt = 0 , QCD_phom2_flt = 0;
		float QCD_jetpt1_flt = 0 , QCD_jetpt2_flt = 0 , QCD_jeteta1_flt = 0 , QCD_jeteta2_flt = 0;
		float QCD_jetphi1_flt = 0 , QCD_jetphi2_flt = 0 , QCD_jetm1_flt = 0 , QCD_jetm2_flt = 0 , QCD_jetbtagDeepB1_flt = 0 , QCD_jetbtagDeepB2_flt = 0 ;

		float QCD_phor91_flt = 0 , QCD_phor92_flt = 0 , QCD_phohoe1_flt = 0 , QCD_phohoe2_flt = 0;
		float QCD_phophiso1_flt = 0 , QCD_phophiso2_flt = 0 , QCD_phochiso1_flt = 0 , QCD_phochiso2_flt = 0;
		float QCD_phosieie1_flt = 0 , QCD_phosieie2_flt = 0 , QCD_phoidmva1_flt = 0 , QCD_phoidmva2_flt = 0;

		photon1_t.SetPtEtaPhiM(0, 0, 0, 0) ,photon2_t.SetPtEtaPhiM(0, 0, 0, 0) , jet1_t.SetPtEtaPhiM(0, 0, 0, 0) , jet2_t.SetPtEtaPhiM(0, 0, 0, 0)   ;

		QCD_jetpt->clear();
		QCD_jeteta->clear();
		QCD_jetphi->clear();
		QCD_jetm->clear();
		QCD_phopt->clear();
		QCD_phoeta->clear();
		QCD_phophi->clear();
		QCD_phom->clear();

		QCD_phoid->clear();
		QCD_jetid->clear() ;
		QCD_r9->clear();
		QCD_hoe->clear() ;
		QCD_chiso->clear() ;
		QCD_sieie->clear();
		QCD_Photon_mvaID->clear() ;
		QCD_phiso->clear() ;
		QCD_phocutBased->clear();
		QCD_jetpuid->clear();
		QCD_jetbtagDeepB->clear();

		tree->GetEntry(evt); //read event
		//	cout << "-----entry---- "   <<  evt  << endl;

		if(QCD_phopt->size() >=2 && QCD_jetpt->size() >=2 ) {
			QCD_phopt1_flt = QCD_phopt->at(0);
			QCD_phopt2_flt = QCD_phopt->at(1);
			QCD_phoeta1_flt = QCD_phoeta->at(0);
			QCD_phoeta2_flt = QCD_phoeta->at(1);
			QCD_phophi1_flt = QCD_phophi->at(0);
			QCD_phophi2_flt = QCD_phophi->at(1);
			QCD_phom1_flt = QCD_phom->at(0);
			QCD_phom2_flt = QCD_phom->at(1);
			QCD_phoid1_flt = QCD_phoid->at(0);
			QCD_phoid2_flt = QCD_phoid->at(1);
			QCD_phor91_flt = QCD_r9->at(0);
			QCD_phor92_flt = QCD_r9->at(1);
			QCD_phohoe1_flt = QCD_hoe->at(0);
			QCD_phohoe2_flt = QCD_hoe->at(1);
			QCD_phochiso1_flt = QCD_chiso->at(0);
			QCD_phochiso2_flt = QCD_chiso->at(1) ;
			QCD_phosieie1_flt = QCD_sieie->at(0);
			QCD_phosieie2_flt = QCD_sieie->at(1);
			QCD_phoidmva1_flt = QCD_Photon_mvaID->at(0);
			QCD_phoidmva2_flt = QCD_Photon_mvaID->at(1);
			QCD_phophiso1_flt = QCD_phiso->at(0);
			QCD_phophiso2_flt = QCD_phiso->at(1);
			QCD_jetpt1_flt = QCD_jetpt->at(0);
			QCD_jetpt2_flt = QCD_jetpt->at(1);
			QCD_jeteta1_flt = QCD_jeteta->at(0);
			QCD_jeteta2_flt = QCD_jeteta->at(1);
			QCD_jetphi1_flt = QCD_jetphi->at(0);
			QCD_jetphi2_flt = QCD_jetphi->at(1);
			QCD_jetm1_flt = QCD_jetm->at(0);
			QCD_jetm2_flt = QCD_jetm->at(1);
			QCD_jetid1_flt = QCD_jetid->at(0);
			QCD_jetid2_flt = QCD_jetid->at(1);
			QCD_phocutbased1_flt = QCD_phocutBased->at(0);
			QCD_phocutbased2_flt = QCD_phocutBased->at(1);
			QCD_jetpuid1_flt = QCD_jetpuid->at(0);
			QCD_jetpuid2_flt = QCD_jetpuid->at(1);
			QCD_jetbtagDeepB1_flt = QCD_jetbtagDeepB->at(0);
			QCD_jetbtagDeepB2_flt = QCD_jetbtagDeepB->at(1);

			photon1_t.SetPtEtaPhiM(QCD_phopt1_flt, QCD_phoeta1_flt , QCD_phophi1_flt , QCD_phom1_flt) ;
			photon2_t.SetPtEtaPhiM(QCD_phopt2_flt, QCD_phoeta2_flt , QCD_phophi2_flt , QCD_phom2_flt) ;
			jet1_t.SetPtEtaPhiM(QCD_jetpt1_flt, QCD_jeteta1_flt , QCD_jetphi1_flt , QCD_jetm1_flt) ;
			jet2_t.SetPtEtaPhiM(QCD_jetpt2_flt, QCD_jeteta2_flt , QCD_jetphi2_flt , QCD_jetm2_flt)  ;
			TLorentzVector dijet = jet1_t + jet2_t ;
			TLorentzVector diphoton = photon1_t + photon2_t ;
			double Mjj = dijet.M();
			double Mpp= diphoton.M() ;
			double deltaetajj = QCD_jeteta1_flt - QCD_jeteta2_flt ;
			double deltaetapp = QCD_phoeta1_flt - QCD_phoeta2_flt ;
			double deltaphijj = jet1_t.DeltaPhi(jet2_t);
			double deltaphipp = photon1_t.DeltaPhi(photon2_t);
			double Z = Zi(QCD_jeteta1_flt,QCD_jeteta2_flt ,  QCD_phoeta1_flt ,  QCD_phoeta2_flt) ;
			double dphi = DPHI(QCD_jetphi1_flt , QCD_jetphi2_flt ,  QCD_phophi1_flt , QCD_phophi2_flt ) ;
			double dRp1p2_ = photon1_t.DeltaR(photon2_t);
			double dRj1j2_ = jet1_t.DeltaR(jet2_t);
			double dRp1j1_ = photon1_t.DeltaR(jet1_t);
			double dRp1j2_ = photon1_t.DeltaR(jet2_t);
			double dRp2j1_ = photon2_t.DeltaR(jet1_t);
			double dRp2j2_ = photon2_t.DeltaR(jet2_t);
			
			if(cut_check == -100  ) { //photonpt
				Zi_origin[indexHT_process]->Fill(Z);
				DPHI_origin[indexHT_process]->Fill(dphi);
				photonIDMVA1_origin[indexHT_process]->Fill(QCD_phoidmva1_flt);
				photonIDMVA2_origin[indexHT_process]->Fill(QCD_phoidmva2_flt);

			}
			// pt_cut
			ID_check(ID, QCD_jetid1_flt,QCD_jetid2_flt , QCD_jetpuid1_flt , QCD_jetpuid2_flt  , QCD_jetbtagDeepB1_flt , QCD_jetbtagDeepB2_flt , QCD_phocutbased1_flt ,QCD_phocutbased2_flt , pass_entry );
			Filter_first(cut_check , check_condition  , QCD_phopt1_flt , QCD_phopt2_flt , QCD_jetpt1_flt , QCD_jetpt2_flt  );
			Filter_second(cut_check ,check_condition , QCD_jeteta1_flt , QCD_jeteta2_flt , Mjj ,  Mpp , deltaetajj , QCD_phoeta1_flt ,  QCD_phoeta2_flt );
			Filter_third(cut_check , check_condition ,  dRj1j2_ , dRp1j1_ , dRp1j2_ ,  dRp2j1_ , dRp2j2_);
			
			if(check_condition  ==3 && ID ==1 ) {


				if(cut_check == 1) { //photonpt
					cout << "++++++++++++++++++=  " <<endl;
					photonpt1_cut[indexHT_process]->Fill(abs(QCD_phopt1_flt));
					photonpt2_cut[indexHT_process]->Fill(abs(QCD_phopt2_flt));
				}

				if(cut_check == 2) { //jetpt
					jetpt1_cut[indexHT_process]->Fill(abs(QCD_jetpt1_flt));
					jetpt2_cut[indexHT_process]->Fill(abs(QCD_jetpt2_flt));
				}

				if(cut_check == 3) { //etaj1
					jeteta1_cut[indexHT_process]->Fill(QCD_jeteta1_flt);
				}

				if(cut_check == 4) { //etaj2
					jeteta2_cut[indexHT_process]->Fill(QCD_jeteta2_flt);
				}

				if(cut_check == 5) { //Mjj
					mjj_cut[indexHT_process]->Fill(Mjj);
				}

				if(cut_check == 6) { //Mpp
					mpp_cut[indexHT_process]->Fill(Mpp);
				}

				if(cut_check == 7) { //deltaetajj
					jetdeltaeta[indexHT_process]->Fill(deltaetajj);
				}

				if(cut_check == 8) { //|etap1|
					photoneta1_cut[indexHT_process]->Fill(QCD_phoeta1_flt);
				}

				if(cut_check == 9) { //|etap2|
					photoneta2_cut[indexHT_process]->Fill(QCD_phoeta2_flt);
				}
				/*
								if(cut_check == 10) { //|dRp1p2|
									dRp1p2[indexHT_process]->Fill(dRp1p2_);
								}

								if(cut_check == 11) { //|dRj1j2|
									dRj1j2[indexHT_process]->Fill(dRj1j2_);
								}

								if(cut_check == 12) { //|dRp1j2|
									dRp1j2[indexHT_process]->Fill(dRp1j2_);
								}

								if(cut_check == 13) { //|dRp2j1|
									dRp2j1[indexHT_process]->Fill(dRp2j1_);
								}
				*/


			}

			if(cut_check ==0) {
				photonpt1_nocut[indexHT_process]->Fill(abs(QCD_phopt1_flt));
				photonpt2_nocut[indexHT_process]->Fill(abs(QCD_phopt2_flt));
				deltaetaj1j2_nocut[indexHT_process]->Fill(deltaetajj);
				invmassj1j2_nocut[indexHT_process]->Fill(Mjj);
			}

		} // jet >= 2 && pho >= 2


	}//end of event loop
	//string color = "kCyan" ;
	array_eff[indexHT_process] =  (float) pass_entry/total_events ;
	array_eff_photonpt[indexHT_process] = (float) pass_entry_photonpt/total_events ;

	cout << "-------------------------"<< endl;
	cout <<"pass_entry= " << pass_entry_photonpt   << "  || Total_entry=  "  << total_events <<endl;
	cout << "  eff_photonpt = " << array_eff_photonpt[indexHT_process]<<endl;
	cout << "------------------------" << endl;

	//not cut
	plot_setting(photonpt1_nocut[indexHT_process] , weighted_process ,color);
	plot_setting(photonpt2_nocut[indexHT_process] , weighted_process , color);
	plot_setting(deltaetaj1j2_nocut[indexHT_process] , weighted_process ,color);
	plot_setting(invmassj1j2_nocut[indexHT_process], weighted_process , color);

	plot_setting(photonpt1_cut[indexHT_process], weighted_process ,color);
	plot_setting(photonpt2_cut[indexHT_process], weighted_process ,color);
	plot_setting(jetpt1_cut[indexHT_process], weighted_process ,color);
	plot_setting(jetpt2_cut[indexHT_process], weighted_process ,color);

	plot_setting(jetdeltaeta[indexHT_process], weighted_process ,color);
	plot_setting(mjj_cut[indexHT_process], weighted_process ,color);
	plot_setting(mpp_cut[indexHT_process], weighted_process ,color);


	plot_setting(photonIDMVA1_origin[indexHT_process], weighted_process ,color);
	plot_setting(photonIDMVA2_origin[indexHT_process], weighted_process ,color);
	plot_setting(Zi_origin[indexHT_process], weighted_process ,color);
	plot_setting(DPHI_origin[indexHT_process], weighted_process ,color);

	plot_setting(photoneta1_cut[indexHT_process], weighted_process ,color);
	plot_setting(photoneta2_cut[indexHT_process], weighted_process ,color);
	plot_setting(jeteta1_cut[indexHT_process], weighted_process ,color);
	plot_setting(jeteta2_cut[indexHT_process], weighted_process ,color);

	plot_setting(dRp1p2[indexHT_process], weighted_process ,color);
	plot_setting(dRj1j2[indexHT_process], weighted_process ,color);
	plot_setting(dRp1j2[indexHT_process], weighted_process ,color);
	plot_setting(dRp2j1[indexHT_process], weighted_process ,color);

	//photonpt2_nocut_stack_qcd->Add(photonpt2_nocut[indexHT_process]);

}



void process_VBS(string string_inputroot_process , string string_inputtree_process , double weighted_process, int indexHT_process  , int cut_check , string color ) {
	// Open NanoAOD file


	TFile *file = new TFile(string_inputroot_process.c_str());
	//TFile *file = TFile::Open(string_inputroot_process);
	TTree *tree;
	file->GetObject(string_inputtree_process.c_str(),tree);
	//tree->Print();
	cout << string_inputroot_process << "  ||  " << string_inputtree_process <<endl;

	// Set up branches for GenPart and Photon collections

	std::vector<float> *GenPart_pt = 0, *GenPart_eta = 0, *GenPart_phi = 0;
	std::vector<int> *GenPart_pdgId = 0, *GenPart_status = 0 , *Photon_cutBased = 0 , *Jet_Id=0 , *Jet_puId=0 ;
	std::vector<float> *Photon_eta = 0, *Photon_phi = 0 , *Photon_pt = 0 , *Photon_mass = 0;
	std::vector<float> *r9 = 0, *hoe = 0 , *phiso = 0 , *chiso = 0 , *sieie = 0 , *Photon_mvaID=0 ;
	std::vector<float> *Jet_pt = 0, *Jet_eta = 0, *Jet_phi = 0 , *Jet_m = 0  ;

	tree->SetBranchAddress("GenPart_pt", &GenPart_pt);
	tree->SetBranchAddress("GenPart_eta", &GenPart_eta);
	tree->SetBranchAddress("GenPart_phi", &GenPart_phi);
	tree->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
	tree->SetBranchAddress("GenPart_status", &GenPart_status);

	tree->SetBranchAddress("Photon_eta", &Photon_eta);
	tree->SetBranchAddress("Photon_phi", &Photon_phi);
	tree->SetBranchAddress("Photon_mass",   &Photon_mass  );
	tree->SetBranchAddress("Photon_pt",   &Photon_pt  );
	tree->SetBranchAddress("Photon_r9",   &r9  );
	tree->SetBranchAddress("Photon_hoe",   &hoe  );
	tree->SetBranchAddress("Photon_pfRelIso03_all",   &phiso  );
	tree->SetBranchAddress("Photon_pfRelIso03_chg",   &chiso  );
	tree->SetBranchAddress("Photon_sieie",   &sieie  );
	tree->SetBranchAddress("Photon_mvaID",   &Photon_mvaID  );
	tree->SetBranchAddress("Photon_cutBasedBitmap",   &Photon_cutBased  );

	tree->SetBranchAddress("Jet_pt", &Jet_pt);
	tree->SetBranchAddress("Jet_eta", &Jet_eta);
	tree->SetBranchAddress("Jet_phi", &Jet_phi);
	tree->SetBranchAddress("Jet_mass", &Jet_m );
	tree->SetBranchAddress("Jet_jetId", &Jet_Id );
	tree->SetBranchAddress("Jet_puId", &Jet_puId );

	// Loop over events
	Long64_t nEntries = tree->GetEntries();

	int total_events = tree->GetEntries() ;
	cout << " total entry=  "   <<  total_events  << endl;
	int pass_entry  = 0 ;
	int total_events_ = 1000;
	int pass_entry_photonpt = 0 ;

	TLorentzVector photon1_t , photon2_t , jet1_t , jet2_t ;
	int match_number = 0 ;



	for (int evt = 0 ;  evt <  total_events ; ++evt) {

		photon1_t.SetPtEtaPhiM(0, 0, 0, 0) ,photon2_t.SetPtEtaPhiM(0, 0, 0, 0) , jet1_t.SetPtEtaPhiM(0, 0, 0, 0) , jet2_t.SetPtEtaPhiM(0, 0, 0, 0)   ;

		Photon_pt->clear();
		Photon_eta->clear();
		Photon_phi->clear();
		Photon_mass->clear();
		r9->clear();
		hoe->clear();
		phiso->clear();
		chiso->clear();
		sieie->clear();
		Photon_mvaID ->clear() ;
		Photon_cutBased->clear();

		GenPart_pt->clear() ;
		GenPart_eta->clear() ;
		GenPart_phi->clear();
		GenPart_pdgId->clear() ;
		GenPart_status->clear() ;

		Jet_pt->clear();
		Jet_eta->clear();
		Jet_phi->clear();
		Jet_m->clear();
		Jet_Id->clear();
		Jet_puId->clear();

		tree->GetEntry(evt);
		int nGenPart = GenPart_pt->size()  ;
		int nPhoton  = Photon_eta->size() ;
		int nJet = Jet_eta->size();


		//cout << "NJET= " << nJet <<" | nGenPart" << nGenPart <<  endl;

		// Collect generator-level photons
		std::vector<int> genPhotonIndices;
		for (int j = 0; j < nGenPart; ++j) {
			if (std::abs(GenPart_pdgId->at(j)) == 22 && GenPart_status->at(j) == 1) {  // Photon and final state
				genPhotonIndices.push_back(j);
			}
		}

		// Vector to store matched reconstructed photons
		std::vector<MatchedPhoton> matchedPhotons;

		// Match reconstructed photons to generator-level photons
		for (int k = 0; k < nPhoton; ++k) {
			float reco_eta = Photon_eta->at(k);
			float reco_phi = Photon_phi->at(k);

			float min_deltaR = 0.2 ;  // Define matching threshold
			//float min_deltaR = 1000000;
			int bestMatchIndex = -1;

			for (int genIdx : genPhotonIndices) {
				float gen_eta = GenPart_eta->at(genIdx);
				float gen_phi = GenPart_phi->at(genIdx);
				float dR = deltaR(reco_eta, reco_phi, gen_eta, gen_phi);

				if (dR < min_deltaR) {
					min_deltaR = dR;
					bestMatchIndex = genIdx;
				}
			}

			// If a match was found, save the matched photon
			if (bestMatchIndex != -1) {
				matchedPhotons.push_back({
					k,                       // Index in reconstructed Photon collection
					Photon_pt->at(k),        // pt of the reco photon
					Photon_eta->at(k),       // eta of the reco photon
					Photon_phi->at(k) ,       // mass of the reco photon
					Photon_mass->at(k)   ,     // r9 of the reco photon
					r9->at(k),
					hoe->at(k),
					phiso->at(k),
					chiso->at(k),
					sieie->at(k) ,
					Photon_mvaID->at(k)  ,
					Photon_cutBased->at(k)
				});
			}

			//h_minDeltaR_photon->Fill(min_deltaR);

		}// end of gen match

		std::sort( matchedPhotons.begin(), matchedPhotons.end(),
		[](const MatchedPhoton &a, const MatchedPhoton &b) {
			return a.pt > b.pt;
		} );

		nGenPart = GenPart_pt->size()  ;
		// Step 1: Collect generator-level quarks
		std::vector<int> genQuarksIndices;
		//cout << "begin q-jet " << endl;
		for (int j = 0; j < nGenPart; j++) {
			int absPdgId = abs(GenPart_pdgId->at(j));
			if ((absPdgId >= 1 && absPdgId <= 6) && GenPart_status->at(j) == 23) {  // Quarks in the hard process
				genQuarksIndices.push_back(j);
				//cout << "begin q-jet step2" << endl;
			}
		}



		// Step 2: Match reconstructed jets to generator-level quarks
		std::vector<MatchedObject> matchedJets;


		for (int k = 0; k < nJet; k++) {
			float reco_eta = Jet_eta->at(k);
			float reco_phi = Jet_phi->at(k);

			float minDeltaR = 0.4;  // Matching threshold for jets
			//float minDeltaR = 1000000 ;

			int bestMatchIndex = -1;

			for (int genIdx : genQuarksIndices ) {
				float gen_eta = GenPart_eta->at(genIdx);
				float gen_phi = GenPart_phi->at(genIdx);
				float dR = deltaR(reco_eta, reco_phi, gen_eta, gen_phi);
				if (dR < minDeltaR) {
					minDeltaR = dR ;
					bestMatchIndex = genIdx ;
				}
			}

			// If a match was found, store the matched jet
			if (bestMatchIndex != -1) {
				matchedJets.push_back({
					k,
					Jet_pt->at(k),
					Jet_eta->at(k),
					Jet_phi->at(k),
					Jet_m->at(k),
					Jet_Id->at(k),
					Jet_puId->at(k)
				});
			}


			//h_minDeltaR->Fill(minDeltaR);
			//cout <<"======= " << minDeltaR <<endl;

		} // end of gen matching

		// Step 3: Sort gen quarks and matched jets by pt

		std::sort(matchedJets.begin(), matchedJets.end(), [](const MatchedObject &a, const MatchedObject &b) {
			return a.pt > b.pt;
		});

		// Step 4: Select and print the two best gen quarks and reco jets
		//std::cout << "Event " << i << ":\n";

		int check_condition = 0;
		int ID = 0 ;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (matchedPhotons.size() >= 2  && matchedJets.size() >= 2 ) {
			match_number = match_number + 1 ;
			const MatchedPhoton &photon_first = matchedPhotons[0];
			const MatchedPhoton &photon_second = matchedPhotons[1];
			const MatchedObject &jet_first = matchedJets[0];
			const MatchedObject &jet_second = matchedJets[1];
			photon1_t.SetPtEtaPhiM(photon_first.pt, photon_first.eta, photon_first.phi, photon_first.m) ;
			photon2_t.SetPtEtaPhiM(photon_second.pt, photon_second.eta, photon_second.phi , photon_second.m) ;
			jet1_t.SetPtEtaPhiM(jet_first.pt, jet_first.eta, jet_first.phi, jet_first.m) ;
			jet2_t.SetPtEtaPhiM(jet_second.pt, jet_second.eta, jet_second.phi , jet_second.m) ;

			TLorentzVector dijet = jet1_t + jet2_t ;
			TLorentzVector diphoton = photon1_t + photon2_t ;
			double Mjj = dijet.M();
			double Mpp= diphoton.M() ;
			double deltaetajj = jet_first.eta - jet_second.eta  ;
			double deltaetapp = photon_first.eta - photon_second.eta ;
			double deltaphijj = jet1_t.DeltaPhi(jet2_t);
			double deltaphipp = photon1_t.DeltaPhi(photon2_t);
			double Z = Zi( jet_first.eta , jet_second.eta,  photon_first.eta,  photon_second.eta) ;
			double dphi = DPHI(jet_first.phi , jet_second.phi,  photon_first.phi,  photon_second.phi ) ;
			double dRp1p2_ = photon1_t.DeltaR(photon2_t);
			double dRj1j2_ = jet1_t.DeltaR(jet2_t);
			double dRp1j1_ = photon1_t.DeltaR(jet1_t);
			double dRp1j2_ = photon1_t.DeltaR(jet2_t);
			double dRp2j1_ = photon2_t.DeltaR(jet1_t);
			double dRp2j2_ = photon2_t.DeltaR(jet2_t);
			
			if(cut_check == -100  ) { //photonpt
				Zi_origin[indexHT_process]->Fill(Z);
				DPHI_origin[indexHT_process]->Fill(dphi);
				photonIDMVA1_origin[indexHT_process]->Fill(photon_first.mvaID);
				photonIDMVA2_origin[indexHT_process]->Fill(photon_second.mvaID);
			}


			// pt_cut
			ID_check(ID, jet_first.jetid , jet_second.jetid , jet_first.jetpuid , jet_second.jetpuid  , 0 , 0 , photon_first.phocutBased , photon_second.phocutBased , pass_entry );
			Filter_first(cut_check , check_condition  , photon_first.pt , photon_second.pt , jet_first.pt , jet_second.pt  );
			Filter_second(cut_check ,check_condition , jet_first.eta , jet_second.eta , Mjj ,  Mpp , deltaetajj , photon_first.eta ,  photon_second.eta );
			Filter_third(cut_check , check_condition ,  dRj1j2_ , dRp1j1_ , dRp1j2_ ,  dRp2j1_ , dRp2j2_);

			if(check_condition  ==3 && ID ==1 ) {

				if(cut_check == 1) { //photonpt
					//cout << "true=  " <<  abs(QCD_phopt1_flt) <<endl;
					photonpt1_cut[indexHT_process]->Fill(abs(photon_first.pt));
					photonpt2_cut[indexHT_process]->Fill(abs(photon_second.pt));
					cout << "+++++++++++++++"<<endl;
				}

				if(cut_check == 2) { //jetpt
					jetpt1_cut[indexHT_process]->Fill(abs(jet_first.pt));
					jetpt2_cut[indexHT_process]->Fill(abs(jet_second.pt));
				}

				if(cut_check == 3) { //etaj1
					jeteta1_cut[indexHT_process]->Fill(jet_first.eta);
				}

				if(cut_check == 4) { //etaj2
					jeteta2_cut[indexHT_process]->Fill(jet_second.eta);
				}

				if(cut_check == 5) { //Mjj
					mjj_cut[indexHT_process]->Fill(Mjj);
				}

				if(cut_check == 6) { //Mpp
					mpp_cut[indexHT_process]->Fill(Mpp);
				}

				if(cut_check == 7) { //deltaetajj
					jetdeltaeta[indexHT_process]->Fill(deltaetajj);
				}

				if(cut_check == 8) { //|etap1|
					photoneta1_cut[indexHT_process]->Fill(photon_first.eta);
				}

				if(cut_check == 9) { //|etap2|
					photoneta2_cut[indexHT_process]->Fill(photon_second.eta);
				}

				if(cut_check == 10) { //|dRp1p2|
					dRp1p2[indexHT_process]->Fill(dRp1p2_);
				}

				if(cut_check == 11) { //|dRj1j2|
					dRj1j2[indexHT_process]->Fill(dRj1j2_);
				}

				if(cut_check == 12) { //|dRp1j2|
					dRp1j2[indexHT_process]->Fill(dRp1j2_);
				}

				if(cut_check == 13) { //|dRp2j1|
					dRp2j1[indexHT_process]->Fill(dRp2j1_);
				}

			}

			if(cut_check ==0) {
				photonpt1_nocut[indexHT_process]->Fill(abs(photon_first.pt));
				photonpt2_nocut[indexHT_process]->Fill(abs(photon_second.pt));
				deltaetaj1j2_nocut[indexHT_process]->Fill(deltaetajj);
				invmassj1j2_nocut[indexHT_process]->Fill(Mjj);
			}
		}


	}// end of event loop

	//int color = 6 ;
	array_eff[indexHT_process] =  (float) pass_entry/match_number ;
	array_eff_photonpt[indexHT_process] = (float) pass_entry_photonpt/match_number ;



	cout << "-------------------------"<< endl;
	cout <<"pass_entry= " << pass_entry   << "  || Total_entry=  "  <<  match_number <<endl;
	cout << "  eff_photonpt = " << array_eff_photonpt[indexHT_process]<<endl;
	cout << "------------------------" << endl;
	//not cut
	plot_setting(photonpt1_nocut[indexHT_process] , weighted_process ,color);
	plot_setting(photonpt2_nocut[indexHT_process] , weighted_process , color);
	plot_setting(deltaetaj1j2_nocut[indexHT_process] , weighted_process ,color);
	plot_setting(invmassj1j2_nocut[indexHT_process], weighted_process , color);

	plot_setting(photonpt1_cut[indexHT_process], weighted_process ,color);
	plot_setting(photonpt2_cut[indexHT_process], weighted_process ,color);
	plot_setting(jetpt1_cut[indexHT_process], weighted_process ,color);
	plot_setting(jetpt2_cut[indexHT_process], weighted_process ,color);
	plot_setting(jetdeltaeta[indexHT_process], weighted_process ,color);
	plot_setting(mjj_cut[indexHT_process], weighted_process ,color);
	plot_setting(mpp_cut[indexHT_process], weighted_process ,color);

	plot_setting(photonIDMVA1_origin[indexHT_process], weighted_process ,color);
	plot_setting(photonIDMVA2_origin[indexHT_process], weighted_process ,color);
	plot_setting(Zi_origin[indexHT_process], weighted_process ,color);
	plot_setting(DPHI_origin[indexHT_process], weighted_process ,color);

	plot_setting(photoneta1_cut[indexHT_process], weighted_process ,color);
	plot_setting(photoneta2_cut[indexHT_process], weighted_process ,color);
	plot_setting(jeteta1_cut[indexHT_process], weighted_process ,color);
	plot_setting(jeteta2_cut[indexHT_process], weighted_process ,color);

	plot_setting(dRp1p2[indexHT_process], weighted_process ,color);
	plot_setting(dRj1j2[indexHT_process], weighted_process ,color);
	plot_setting(dRp1j2[indexHT_process], weighted_process ,color);
	plot_setting(dRp2j1[indexHT_process], weighted_process ,color);




	file->Close();

}




void process(int indexHT , int cut_check) {


	if(indexHT == 16 ) { //50to100
		const double weighted = 285450.3883 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd50to100_filter.root" ;
		string string_inputtree = "TreeB50to100" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color) ;
	}

	if(indexHT == 15 ) { //100to200
		const double weighted = 14986.56339 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd100to200_filter.root" ;
		string string_inputtree = "TreeB100to200" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 14 ) { //200to300
		const double weighted = 1710 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd200to300_filter.root" ;
		string string_inputtree = "TreeB200to300" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 13 ) { //300to500
		const double weighted = 353.2093356 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd300to500_filter.root" ;
		string string_inputtree = "TreeB300to500" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 12 ) { //500to700
		const double weighted = 32.62486728 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd500to700_filter.root" ;
		string string_inputtree = "TreeB500to700" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 11) { //700to1000
		const double weighted = 7.864342292 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd700to1000_filter.root" ;
		string string_inputtree = "TreeB700to1000" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 10) { //1000to1500
		const double weighted = 4.21514623 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd1000to1500_filter.root" ;
		string string_inputtree = "TreeB1000to1500" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check  , color ) ;
	}

	if(indexHT == 9) { //1500to2000
		const double weighted = 0.5436444275 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd1500to2000_filter.root" ;
		string string_inputtree = "TreeB1500to2000" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check  , color ) ;
	}

	if(indexHT == 8) { //2000toInf
		const double weighted = 0.2218711951 ;
		string color = "kCyan" ;
		string string_inputroot = "2qcd2000toInf_filter.root" ;
		string string_inputtree = "TreeB2000toInf" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 7 ) { //GJ 40to100 7
		const double weighted = 132.6333278 ;
		string color = "kYellow" ;
		string string_inputroot = "GJET_40to100_.root" ;
		string string_inputtree = "TreeB40to100" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color) ;
	}

	if(indexHT == 6 ) { //GJ 100to200 6
		const double weighted = 56.35388668 ;
		string color = "kYellow" ;
		string string_inputroot = "GJET_100to200_.root" ;
		string string_inputtree = "TreeB100to200" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color) ;
	}

	if(indexHT == 5 ) { //GJ 200to400 5
		const double weighted = 7.292141497 ;
		string color = "kYellow" ;
		string string_inputroot = "GJET_200to400_.root" ;
		string string_inputtree = "TreeB200to400" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color) ;
	}

	if(indexHT == 4 ) { //GJ 400to600 4
		const double weighted = 3.529953383 ;
		string color = "kYellow" ;
		string string_inputroot = "GJET_400to600_.root" ;
		string string_inputtree = "TreeB400to600" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color) ;
	}

	if(indexHT == 3 ) { //GJ 600toInf 3
		const double weighted = 1.119542569;
		string color = "kYellow" ;
		string string_inputroot = "GJET_600toInf_.root" ;
		string string_inputtree = "TreeB600toInf" ;
		process_HLT(string_inputroot , string_inputtree , weighted, indexHT  ,  cut_check , color) ;
	}

	if(indexHT == 2) { //qcd-diphotonjets 2
		//TFile *file = TFile::Open("vbs_filter.root"); TTree *tree = (TTree*)file->Get("Event");
		const double weighted = 2.327106 ;
		string color = "kGreen+3" ;
		string string_inputroot = "qcd_filter_cutHLT.root" ;
		string string_inputtree = "Event" ;
		process_VBS( string_inputroot  , string_inputtree  , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 1) { //interference 1
		//TFile *file = TFile::Open("vbs_filter.root"); TTree *tree = (TTree*)file->Get("Event");
		const double weighted = 0.000868038 ;
		string color = "kAzure" ;
		string string_inputroot = "inter_filter_cutHLT.root" ;
		string string_inputtree = "Event" ;
		process_VBS( string_inputroot  , string_inputtree  , weighted, indexHT  ,  cut_check , color ) ;
	}


	if(indexHT == 0) { //VBS 0
		//TFile *file = TFile::Open("vbs_filter.root"); TTree *tree = (TTree*)file->Get("Event");
		const double weighted =  0.005320464;
		string color = "kRed+1" ;
		string string_inputroot = "vbs_filter_cutHLT.root" ;
		string string_inputtree = "Event" ;
		process_VBS( string_inputroot  , string_inputtree  , weighted, indexHT  ,  cut_check , color ) ;
	}

}






int plot_f() {


	for(int k1 = 0 ; k1<20 ; k1++ ) {

		stack_plot[k1] = new THStack("stack_plot", "stack_plot");
	}

	TCanvas *canvas_photonpt1_nocut = new TCanvas("canvas", "Divided Canvas Example", 800, 600);
	canvas_photonpt1_nocut->Divide(1, 1);
	TCanvas *canvas_photonpt2_nocut = new TCanvas("canvas1", "Divided Canvas Example1", 800, 600);
	canvas_photonpt2_nocut->Divide(1, 1);

	TCanvas *canvas_photonpt1_stack_nocut = new TCanvas("canvas_photonpt1_stack_nocut", "canvas_photonpt1_stack_nocut", 800, 600);
	canvas_photonpt1_stack_nocut->Divide(1, 1);
	TCanvas *canvas_photonpt2_stack_nocut = new TCanvas("canvas_photonpt2_stack_nocut", "canvas_photonpt2_stack_nocut", 800, 600);
	canvas_photonpt2_stack_nocut->Divide(1, 1);

	TCanvas *canvas_deltaetajj_stack_nocut = new TCanvas("canvas_deltaetajj_stack_nocut", "canvas_deltaetajj_stack_nocut", 800, 600);
	canvas_deltaetajj_stack_nocut->Divide(1, 1);
	TCanvas *canvas_invmassjj_stack_nocut = new TCanvas("canvas_invmassjj_stack_nocut", "canvas_invmassjj_stack_nocut", 800, 600);
	canvas_invmassjj_stack_nocut->Divide(1, 1);

	TCanvas *canvas_photonpt1_stack_cut = new TCanvas("canvas_photonpt1_stack_cut", "canvas_photonpt1_stack_cut", 800, 600);
	canvas_photonpt1_stack_cut->Divide(1, 1);
	TCanvas *canvas_photonpt2_stack_cut = new TCanvas("canvas_photonpt2_stack_cut", "canvas_photonpt2_stack_cut", 800, 600);
	canvas_photonpt2_stack_cut->Divide(1, 1);
	TCanvas *canvas_jetpt1_stack_cut = new TCanvas("canvas_jetpt1_stack_cut", "canvas_jetpt1_stack_cut", 800, 600);
	canvas_jetpt1_stack_cut->Divide(1, 1);
	TCanvas *canvas_jetpt2_stack_cut = new TCanvas("canvas_jetpt2_stack_cut", "canvas_jetpt2_stack_cut", 800, 600);
	canvas_jetpt2_stack_cut->Divide(1, 1);

	TCanvas *canvas_photoneta1_stack_cut = new TCanvas("canvas_photoneta1_stack_cut", "canvas_photoneta1_stack_cut", 800, 600);
	canvas_photoneta1_stack_cut->Divide(1, 1);
	TCanvas *canvas_photoneta2_stack_cut = new TCanvas("canvas_photoneta2_stack_cut", "canvas_photoneta2_stack_cut", 800, 600);
	canvas_photoneta2_stack_cut->Divide(1, 1);
	TCanvas *canvas_jeteta1_stack_cut = new TCanvas("canvas_jeteta1_stack_cut", "canvas_jeteta1_stack_cut", 800, 600);
	canvas_jeteta1_stack_cut->Divide(1, 1);
	TCanvas *canvas_jeteta2_stack_cut = new TCanvas("canvas_jeteta2_stack_cut", "canvas_jeteta2_stack_cut", 800, 600);
	canvas_jeteta2_stack_cut->Divide(1, 1);

	TCanvas *canvas_deltaphijj_stack_cut = new TCanvas("canvas_deltaphijj_stack_cut", "canvas_deltaphijj_stack_cut", 800, 600);
	canvas_deltaphijj_stack_cut->Divide(1, 1);
	TCanvas *canvas_deltaetajj_stack_cut = new TCanvas("canvas_deltaetajj_stack_cut", "canvas_deltaetajj_stack_cut", 800, 600);
	canvas_deltaetajj_stack_cut->Divide(1, 1);
	TCanvas *canvas_deltaphipp_stack_cut = new TCanvas("canvas_deltaphipp_stack_cut", "canvas_deltaphipp_stack_cut", 800, 600);
	canvas_deltaphipp_stack_cut->Divide(1, 1);
	TCanvas *canvas_deltaetapp_stack_cut = new TCanvas("canvas_deltaetapp_stack_cut", "canvas_deltaetapp_stack_cut", 800, 600);
	canvas_deltaetapp_stack_cut->Divide(1, 1);

	TCanvas *canvas_invmassjj_stack_cut = new TCanvas("canvas_invmassjj_stack_cut", "canvas_invmassjj_stack_cut", 800, 600);
	canvas_invmassjj_stack_cut->Divide(1, 1);
	TCanvas *canvas_invmasspp_stack_cut = new TCanvas("canvas_invmasspp_stack_cut", "canvas_invmasspp_stack_cut", 800, 600);
	canvas_invmasspp_stack_cut->Divide(1, 1);

	TCanvas *canvas_DPHI_stack = new TCanvas("canvas_DPHI_stack", "canvas_DPHI_stack", 800, 600);
	canvas_DPHI_stack->Divide(1, 1);
	TCanvas *canvas_Zi_stack = new TCanvas("canvas_Zi_stack", "canvas_Zi_stack", 800, 600);
	canvas_Zi_stack->Divide(1, 1);
	TCanvas *canvas_photonIDMVA1_stack = new TCanvas("canvas_photonIDMVA1_stack", "canvas_photonIDMVA1_stack", 800, 600);
	canvas_photonIDMVA1_stack->Divide(1, 1);
	TCanvas *canvas_photonIDMVA2_stack = new TCanvas("canvas_photonIDMVA2_stack", "canvas_photonIDMVA2_stack", 800, 600);
	canvas_photonIDMVA2_stack->Divide(1, 1);

	//deltaR canvas
	TCanvas *canvas_dRp1p2_stack_cut = new TCanvas("canvas_dRp1p2_stack_cut", "canvas_dRp1p2_stack_cut", 800, 600);
	canvas_dRp1p2_stack_cut->Divide(1, 1);
	TCanvas *canvas_dRj1j2_stack_cut = new TCanvas("canvas_dRj1j2_stack_cut", "canvas_dRj1j2_stack_cut", 800, 600);
	canvas_dRj1j2_stack_cut->Divide(1, 1);
	TCanvas *canvas_dRp1j1_stack_cut = new TCanvas("canvas_dRp1j1_stack_cut", "canvas_dRp1j1_stack_cut", 800, 600);
	canvas_dRp1j1_stack_cut->Divide(1, 1);
	TCanvas *canvas_dRp1j2_stack_cut = new TCanvas("canvas_dRp1j2_stack_cut", "canvas_dRp1j2_stack_cut", 800, 600);
	canvas_dRp1j2_stack_cut->Divide(1, 1);
	TCanvas *canvas_dRp2j1_stack_cut = new TCanvas("canvas_dRp2j1_stack_cut", "canvas_dRp2j1_stack_cut", 800, 600);
	canvas_dRp2j1_stack_cut->Divide(1, 1);
	TCanvas *canvas_dRp2j2_stack_cut = new TCanvas("canvas_dRp2j2_stack_cut", "canvas_dRp2j2_stack_cut", 800, 600);
	canvas_dRp2j2_stack_cut->Divide(1, 1);

	//cut_all--------------------------------------------------------------------------
	TCanvas *canvas_photonpt1_stack_cut_all = new TCanvas("canvas_photonpt1_stack_cut_all", "canvas_photonpt1_stack_cut_all", 800, 600);
	canvas_photonpt1_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_photonpt2_stack_cut_all = new TCanvas("canvas_photonpt2_stack_cut_all", "canvas_photonpt2_stack_cut_all", 800, 600);
	canvas_photonpt2_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_jetpt1_stack_cut_all = new TCanvas("canvas_jetpt1_stack_cut_all", "canvas_jetpt1_stack_cut_all", 800, 600);
	canvas_jetpt1_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_jetpt2_stack_cut_all = new TCanvas("canvas_jetpt2_stack_cut_all", "canvas_jetpt2_stack_cut_all", 800, 600);
	canvas_jetpt2_stack_cut_all->Divide(1, 1);

	TCanvas *canvas_deltaphijj_stack_cut_all = new TCanvas("canvas_deltaphijj_stack_cut_all", "canvas_deltaphijj_stack_cut_all", 800, 600);
	canvas_deltaphijj_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_deltaetajj_stack_cut_all = new TCanvas("canvas_deltaetajj_stack_cut_all", "canvas_deltaetajj_stack_cut_all", 800, 600);
	canvas_deltaetajj_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_deltaphipp_stack_cut_all = new TCanvas("canvas_deltaphipp_stack_cut_all", "canvas_deltaphipp_stack_cut_all", 800, 600);
	canvas_deltaphipp_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_deltaetapp_stack_cut_all = new TCanvas("canvas_deltaetapp_stack_cut_all", "canvas_deltaetapp_stack_cut_all", 800, 600);
	canvas_deltaetapp_stack_cut_all->Divide(1, 1);

	TCanvas *canvas_invmassjj_stack_cut_all = new TCanvas("canvas_invmassjj_stack_cut_all", "canvas_invmassjj_stack_cut_all", 800, 600);
	canvas_invmassjj_stack_cut_all->Divide(1, 1);
	TCanvas *canvas_invmasspp_stack_cut_all = new TCanvas("canvas_invmasspp_stack_cut_all", "canvas_invmasspp_stack_cut_all", 800, 600);
	canvas_invmasspp_stack_cut_all->Divide(1, 1);
	int cut_check=-100;

	//for(cut_check = 1 ; cut_check < 15 ; cut_check++) {


	for(int i = 0 ; i < 17 ; i++) {
		topPad[i] = new TPad("topPad", "Top Pad", 0, 0.3, 1, 1);
		bottomPad[i] = new TPad("botPad", "Bot Pad", 0, 0, 1, 0.3);;
		
		photonpt1_nocut[i] = new TH1F(  "nocutphotonpt1" , "nocutphotonpt1", 20, 0, 700);
		photonpt2_nocut[i] = new TH1F(  "nocutphotonpt2" , "nocutphotonpt2", 20, 0, 700);
		deltaetaj1j2_nocut[i] = new TH1F(  "deltaetaj1j2_nocut" , "deltaetaj1j2_nocut", 100, -8, 8);
		invmassj1j2_nocut[i] = new TH1F(  "invmassj1j2_nocut" , "invmassj1j2_nocut", 100, 0, 5000);

		//cut pt
		photonpt1_cut[i] =  new TH1F(  "photonpt1" , "photonpt1", 20, 0, 500)  ;
		photonpt2_cut[i] = new TH1F(  "photonpt2" , "photonpt2", 20, 0, 200)  ;
		jetpt1_cut[i] = new TH1F(  "jetpt1" , "jetpt1", 20, 0, 1000) ;
		jetpt2_cut[i] = new TH1F(  "jetpt2" , "jetpt2", 20, 0, 800) ;

		jetdeltaeta[i] = new TH1F(  "jetdeltaeta" , "jetdeltaeta", 20, -6, 6) ;
		photoneta1_cut[i] = new TH1F(  "photoneta1" , "photoneta1", 20, -10, 10) ;
		photoneta2_cut[i] = new TH1F(  "photoneta2" , "photoneta2", 20, -10, 10) ;
		jeteta1_cut[i] = new TH1F(  "jeteta1" , "jeteta1", 20, -6, 6) ;
		jeteta2_cut[i] = new TH1F(  "jeteta2" , "jeteta2", 20, -6, 6) ;

		dRp1p2[i] = new TH1F(  "dRp1p2" , "dRp1p2", 20, 0, 10)  ;
		dRj1j2[i] = new TH1F(  "dRj1j2" , "dRj1j2", 20, 0, 10)  ;
		dRp1j1[i] = new TH1F(  "dRp1j1" , "dRp1j1", 20, 0, 10)  ;
		dRp1j2[i] = new TH1F(  "dRp1j2" , "dRp1j2", 20, 0, 10)  ;
		dRp2j1[i] = new TH1F(  "dRp2j1" , "dRp2j1", 20, 0, 10)  ;
		dRp2j2[i] = new TH1F(  "dRp2j2" , "dRp2j2", 20, 0, 10)  ;

		mpp_cut[i] = new TH1F(  "mpp" , "mpp", 20, 0, 800)  ;
		mjj_cut[i] = new TH1F(  "mjj" , "mjj", 20, 0, 5000)  ;

		Zi_origin[i] = new TH1F(  "Zi_origin" , "Zi_origin", 20, 0, 10) ;
		DPHI_origin[i] = new TH1F(  "DPHI_origin" , "DPHI_origin", 20, 0, 4) ;
		photonIDMVA1_origin[i] = new TH1F(  "photonIDMVA1_origin" , "photonIDMVA1_origin", 20, -1, 1) ;
		photonIDMVA2_origin[i] = new TH1F(  "photonIDMVA2_origin" , "photonIDMVA2_origin", 20, -1, 1) ;

		//cut_all
		photonpt1_cut_all[i] =  new TH1F(  "photonpt1_cut_all" , "photonpt1_cut_all", 100, 0, 1000)  ;
		photonpt2_cut_all[i] = new TH1F(  "photonpt2_cut_all" , "photonpt2_cut_all", 100, 0, 1000)  ;
		jetpt1_cut_all[i] = new TH1F(  "jetpt1_cut_all" , "jetpt1_cut_all", 130, 0, 1300) ;
		jetpt2_cut_all[i] = new TH1F(  "jetpt2_cut_all" , "jetpt2_cut_all", 100, 0, 1000) ;
		photondeltaeta_cut_all[i] = new TH1F(  "photondeltaeta_cut_all" , "photondeltaeta_cut_all", 20, -6, 6) ;
		photoneta2_cut_all[i] = new TH1F(  "photoneta2_cut_all" , "photoneta2_cut_all", 100, -10, 10) ;
		jetdeltaeta_cut_all[i] = new TH1F(  "jetdeltaeta_cut_all" , "jetdeltaeta_cut_all", 20, -10, 10)  ;
		jeteta2_cut_all[i] = new TH1F(  "jeteta2_cut_all" , "jeteta2_cut_all", 100, -10, 10)  ;
		photondeltaphi_cut_all[i] = new TH1F(  "photondeltaphi_cut_all" , "photondeltaphi_cut_all", 20, -6, 6)  ;
		photonphi2_cut_all[i] = new TH1F(  "photonphi2_cut_all" , "photonphi2_cut_all", 20, -6, 6)  ;
		jetdeltaphi[i] =  new TH1F(  "jetdeltaphi_cut_all" , "jetdeltaphi_cut_all", 50, -4, 4) ;
		jetphi2_cut_all[i] =  new TH1F(  "jetphi2_cut_all" , "jetphi2_cut_all", 20, -6, 6) ;
		mpp_cut_all[i] = new TH1F(  "mpp_cut_all" , "mpp_cut_all", 100, 0, 1500)  ;
		mjj_cut_all[i] = new TH1F(  "mjj_cut_all" , "mjj_cut_all", 20, 0, 5000)  ;
		
		 
//////////////////////////////////////////////////////////////////////////////////////////
		cout << "******************************************************" <<endl;
		cout << "                                                  " <<endl;
		cout << "**           " <<i << " "<< i << " "<< i << " " << i << "  "<< i << "         **" <<endl;


		/*
		if(cut_check == 1) { //photonpt
		if(cut_check == 2) { //jetpt
		if(cut_check == 3) { //etaj1
		if(cut_check == 4) { //etaj2
		if(cut_check == 5) { //Mjj
		if(cut_check == 6) { //Mpp
		if(cut_check == 7) { //dEta(j,j)
		if(cut_check == 8) { //|eta_photon1|
		if(cut_check == 9) { //|eta_photon2|
		if(cut_check == 10)		dR(pho, pho)>0.7
		if(cut_check == 11)  dR(j,j)>0.7 > 0.7
		if(cut_check == 12)  dR(pho1,jet2) > 0.7
		if(cut_check == 13)   dR(pho2,jet1) // in addition
		*/
		process(i , cut_check );

//////////////////////////////////////////////////////////////////////////////////////////
		//no cut

		//cut---------------------------
		Add_to_stack(photonpt1_cut[i], stack_plot[0], canvas_photonpt1_stack_cut , i );
		Add_to_stack(photonpt2_cut[i], stack_plot[1] , canvas_photonpt2_stack_cut , i );
		Add_to_stack(jetpt1_cut[i], stack_plot[2] , canvas_jetpt1_stack_cut , i );
		Add_to_stack(jetpt2_cut[i], stack_plot[3] , canvas_jetpt2_stack_cut , i );
		Add_to_stack(jeteta1_cut[i], stack_plot[4] , canvas_jeteta1_stack_cut , i );
		Add_to_stack(jeteta2_cut[i], stack_plot[5] , canvas_jeteta2_stack_cut , i );
		Add_to_stack(mjj_cut[i], stack_plot[6] , canvas_invmassjj_stack_cut , i );
		Add_to_stack(mpp_cut[i], stack_plot[7]  , canvas_invmasspp_stack_cut , i );
		Add_to_stack(jetdeltaeta[i], stack_plot[8]  , canvas_deltaetajj_stack_cut , i );
		Add_to_stack(photoneta1_cut[i], stack_plot[9] , canvas_photoneta1_stack_cut , i );
		Add_to_stack(photoneta2_cut[i], stack_plot[10] , canvas_photoneta2_stack_cut , i );

		//Add_to_stack(dRp1p2[i], stack_plot[11]  , canvas_dRp1p2_stack_cut , i );
		Add_to_stack(dRj1j2[i], stack_plot[11]  , canvas_dRj1j2_stack_cut , i );
		Add_to_stack(dRp1j1[i], stack_plot[12]  , canvas_dRp1j1_stack_cut , i );
		Add_to_stack(dRp1j2[i], stack_plot[13]  , canvas_dRp1j2_stack_cut , i );
		Add_to_stack(dRp2j1[i], stack_plot[14]  , canvas_dRp2j1_stack_cut , i );
		Add_to_stack(dRp2j2[i], stack_plot[15]  , canvas_dRp2j2_stack_cut , i );

		Add_to_stack(DPHI_origin[i], DPHI_origin_stack , canvas_DPHI_stack , i );
		Add_to_stack(Zi_origin[i], Zi_origin_stack , canvas_Zi_stack , i );
		Add_to_stack(photonIDMVA1_origin[i], photonIDMVA1_origin_stack , canvas_photonIDMVA1_stack , i );
		Add_to_stack(photonIDMVA2_origin[i], photonIDMVA2_origin_stack , canvas_photonIDMVA2_stack , i );
		
		
		if(i==0){
			hSignificance[0] = (TH1F*)Zi_origin[0]->Clone("hSignificance[0]");
		}
		else if(i==1){
		 	hback_total[0] = (TH1F*)Zi_origin[1]->Clone("hback_total");
		}
		else{
			hback_total[0]->Add(Zi_origin[i]);
		}
		

	}//end of forloop with stack


	float Max_h = 1E6 ;
	float Min_h = 8E-5 ;




	TH1F *frame = new TH1F("frame", "Frame", 100, 0, 700);
	frame->SetMinimum(10);
	frame->SetMaximum(3E7);  // Set desired Y-axis range
	///////////////////////////////////////////////////////////////
	
	/*
	hSignificance[0] = (TH1F*)Zi_origin[0]->Clone("hSignificance[0]");
		hback_total[0] = (TH1F*)Zi_origin[1] ->Clone("hback_total");
	*/
	if(cut_check == -100) {
		Max_h = 1E9 ;
		Min_h = 1E-6 ;

		TH1F *frame_Zi_origin_stack = new TH1F("frame_Zi_origin_stack", "Z", 20, 0, 10);

		Cal_Sig(Zi_origin[0] , hback_total[0] , hSignificance[0])	;
		Set_TPad( topPad[0] , bottomPad[0] , canvas_Zi_stack  );
		Frame_with_pad( canvas_Zi_stack , frame_Zi_origin_stack , Zi_origin_stack , Min_h , Max_h ,topPad[0]);
		Set_Ratio(canvas_Zi_stack   ,  hSignificance[0] ,  bottomPad[0] );	
		TLegend *legend_Zi_origin ;  
		Add_legend(legend_Zi_origin  , canvas_Zi_stack , Zi_origin[15] , Zi_origin[3] , Zi_origin[2] , Zi_origin[0] , Zi_origin[1]   );
		canvas_Zi_stack->SaveAs("Zi_origin_stack.pdf");
		
		TCanvas *tmp = new TCanvas("tmp", "tmp", 800, 600);
		tmp->Divide(1, 1);
		tmp->cd(1);
		hSignificance[0]->Draw();
		tmp->SaveAs("Zi_origin_stack_sig.pdf");
		/*
		
		
		TH1F *frame_DPHI_origin_stack = new TH1F("frame_photon1_cut_stack", "DPHI", 20, 0, 4);
		Frame( canvas_DPHI_stack , frame_DPHI_origin_stack , DPHI_origin_stack , Min_h , Max_h );
		TLegend *legend_DPHI_origin ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_DPHI_origin  , canvas_DPHI_stack , DPHI_origin[15] , DPHI_origin[3] , DPHI_origin[2] , DPHI_origin[0] , DPHI_origin[1]   );
		canvas_DPHI_stack->SaveAs("DPHI_origin_stack.pdf");

		Max_h = 1E8 ;
		Min_h = 1 ;

		TH1F *frame_photonIDMVA1_origin_stack = new TH1F("frame_photonIDMVA1_origin_stack", "PhotonIDMVA1", 20, -1, 1);
		Frame( canvas_photonIDMVA1_stack , frame_photonIDMVA1_origin_stack , photonIDMVA1_origin_stack , Min_h , Max_h );
		TLegend *legend_photonIDMVA1_origin ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_photonIDMVA1_origin  , canvas_photonIDMVA1_stack , photonIDMVA1_origin[15] , photonIDMVA1_origin[3] , photonIDMVA1_origin[2] , photonIDMVA1_origin[0] , photonIDMVA1_origin[1]   );
		canvas_photonIDMVA1_stack->SaveAs("photonIDMVA1_origin_stack.pdf");

		TH1F *frame_photonIDMVA2_origin_stack = new TH1F("frame_photonIDMVA2_origin_stack", "PhotonIDMVA2", 20, -1, 1);
		Frame( canvas_photonIDMVA2_stack , frame_photonIDMVA2_origin_stack , photonIDMVA2_origin_stack , Min_h , Max_h );
		TLegend *legend_photonIDMVA2_origin ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_photonIDMVA2_origin  , canvas_photonIDMVA2_stack , photonIDMVA2_origin[15] , photonIDMVA2_origin[3] , photonIDMVA2_origin[2] , photonIDMVA2_origin[0] , photonIDMVA2_origin[1]   );
		canvas_photonIDMVA2_stack->SaveAs("photonIDMVA2_origin_stack.pdf");
*/

	}

TLegend *legend_[100];
	/////////////////////////////////////////////////////////////////////////////////
	if( cut_check == 1 ) {
		Max_h = 3E4 ;
		Min_h = 1E-6 ;
		TH1F *frame_photonpt1_cut_stack = new TH1F("frame_photon1_cut_stack", "photonpt1", 20, 0, 500);
		Frame( canvas_photonpt1_stack_cut , frame_photonpt1_cut_stack , stack_plot[0] , Min_h , Max_h );
		TLegend *legend_photonpt1 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_photonpt1  , canvas_photonpt1_stack_cut , photonpt1_cut[15] , photonpt1_cut[3] , photonpt1_cut[2] , photonpt1_cut[0] , photonpt1_cut[1]   );
		canvas_photonpt1_stack_cut->SaveAs("photonpt1_cut_stack.pdf");

		TH1F *frame_photonpt2_cut_stack = new TH1F("frame_photon2_cut_stack", "photonpt2", 20, 0, 200);
		Frame( canvas_photonpt2_stack_cut , frame_photonpt2_cut_stack , stack_plot[1] , Min_h , Max_h );
		TLegend *legend_photonpt2 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_photonpt2  , canvas_photonpt2_stack_cut , photonpt2_cut[15] , photonpt2_cut[3] , photonpt2_cut[2] , photonpt2_cut[0] , photonpt2_cut[1]   );
		canvas_photonpt2_stack_cut->SaveAs("photonpt2_cut_stack.pdf");
	}

	if( cut_check == 2 ) {
		Max_h = 3E4 ;
		Min_h = 1E-6 ;
		TH1F *frame_jetpt1_cut_stack = new TH1F("frame_jet1_cut_stack", "jetpt1", 20, 0, 1000);
		Frame( canvas_jetpt1_stack_cut , frame_jetpt1_cut_stack , stack_plot[2] , Min_h , Max_h );
		TLegend *legend_jetpt1 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_jetpt1  , canvas_jetpt1_stack_cut , jetpt1_cut[15] , jetpt1_cut[3] , jetpt1_cut[2] , jetpt1_cut[0] , jetpt1_cut[1]   );
		canvas_jetpt1_stack_cut->SaveAs("jetpt1_cut_stack.pdf");

		TH1F *frame_jetpt2_cut_stack = new TH1F("frame_jet2_cut_stack", "jetpt2", 20, 0, 800);
		Frame( canvas_jetpt2_stack_cut , frame_jetpt2_cut_stack , stack_plot[3] , Min_h , Max_h );
		TLegend *legend_jetpt2 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_jetpt2  , canvas_jetpt2_stack_cut , jetpt2_cut[15] , jetpt2_cut[3] , jetpt2_cut[2] , jetpt2_cut[0] , jetpt2_cut[1]   );
		canvas_jetpt2_stack_cut->SaveAs("jetpt2_cut_stack.pdf");
	}

	if( cut_check == 3 ) {
		Max_h = 3E4 ;
		Min_h = 1E-6 ;
		TH1F *frame_jeteta1_cut = new TH1F("frame_jet1_cut_stack", "jeteta1", 20, -6, 6);
		Frame( canvas_jeteta1_stack_cut , frame_jeteta1_cut , stack_plot[4] , Min_h , Max_h );
		TLegend *legend_jeteta1 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_jeteta1  , canvas_jeteta1_stack_cut , jeteta1_cut[15] , jeteta1_cut[3] , jeteta1_cut[2] , jeteta1_cut[0] , jeteta1_cut[1]   );
		canvas_jeteta1_stack_cut->SaveAs("jeteta1_cut.pdf");
	}

	if( cut_check == 4 ) {
		Max_h = 3E4 ;
		Min_h = 1E-6 ;
		TH1F *frame_jeteta2_cut = new TH1F("frame_jet2_cut_stack", "jeteta2", 20, -6, 6);
		Frame( canvas_jeteta2_stack_cut , frame_jeteta2_cut , stack_plot[5] , Min_h , Max_h );
		TLegend *legend_jeteta2 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_jeteta2  , canvas_jeteta2_stack_cut , jeteta2_cut[15] , jeteta2_cut[3] , jeteta2_cut[2] , jeteta2_cut[0] , jeteta2_cut[1]   );
		canvas_jeteta2_stack_cut->SaveAs("jeteta2_cut.pdf");
	}

	if( cut_check == 5 ) {
		TH1F *frame_mjj_cut_stack = new TH1F("frame_mjj_nocut_stack", "mjj", 20, 0, 5000);
		Frame( canvas_invmassjj_stack_cut , frame_mjj_cut_stack , stack_plot[6] , Min_h , Max_h );
		TLegend *legend_invmassjj ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_invmassjj  , canvas_invmassjj_stack_cut , mjj_cut[15] , mjj_cut[3] , mjj_cut[2] , mjj_cut[0] , mjj_cut[1] );
		canvas_invmassjj_stack_cut->SaveAs("mjj_cut_stack.pdf");
	}

	if( cut_check == 6 ) {
		TH1F *frame_mpp_cut_stack = new TH1F("frame_mpp_nocut_stack", "mpp", 20, 0, 800);
		Frame( canvas_invmasspp_stack_cut , frame_mpp_cut_stack , stack_plot[7] , Min_h , Max_h );
		TLegend *legend_invmasspp ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_invmasspp  , canvas_invmasspp_stack_cut , mpp_cut[15] , mpp_cut[3] , mpp_cut[2] , mpp_cut[0] , mpp_cut[1] );
		canvas_invmasspp_stack_cut->SaveAs("mpp_cut_stack.pdf");
	}

//dEta(j,j)
	if( cut_check == 7 ) {
		TH1F *frame_jetdeltaeta_cut_stack = new TH1F("frame_jetdeltaeta_nocut_stack", "dEta(j,j)", 20, -6, 6);
		Frame( canvas_deltaetajj_stack_cut , frame_jetdeltaeta_cut_stack , stack_plot[8] , Min_h , Max_h );
		TLegend *legend_deltaetajj ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_deltaetajj  , canvas_deltaetajj_stack_cut , jetdeltaeta[15] , jetdeltaeta[3] , jetdeltaeta[2] , jetdeltaeta[0] , jetdeltaeta[1] );
		canvas_deltaetajj_stack_cut->SaveAs("jetdeltaeta_stack.pdf");
	}


	//|eta_photon1|
	if( cut_check == 8 ) {
		TH1F *frame_etap1_cut_stack = new TH1F("frame_etap1_cut_stack", "etap1", 20, -5, 5);
		Frame( canvas_photoneta1_stack_cut , frame_etap1_cut_stack , stack_plot[9] , Min_h , Max_h );
		TLegend *legend_etap1 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_etap1  , canvas_photoneta1_stack_cut , photoneta1_cut[15] , photoneta1_cut[3] , photoneta1_cut[2] , photoneta1_cut[0] , photoneta1_cut[1] );
		canvas_photoneta1_stack_cut->SaveAs("photoneta1_cut_stack.pdf");
	}

	//|eta_photon2|
	if( cut_check == 9 ) {
		TH1F *frame_etap2_cut_stack = new TH1F("frame_etap2_cut_stack", "etap2", 20, -5, 5);
		Frame( canvas_photoneta2_stack_cut , frame_etap2_cut_stack , stack_plot[10] , Min_h , Max_h );
		TLegend *legend_etap2 ;  // Position: (x1, y1, x2, y2)
		Add_legend(legend_etap2  ,canvas_photoneta2_stack_cut , photoneta2_cut[15] , photoneta2_cut[3] , photoneta2_cut[2] , photoneta2_cut[0] , photoneta2_cut[1] );
		canvas_photoneta2_stack_cut->SaveAs("photoneta2_cut_stack.pdf");
	}
	/*
	//|dRp1p2|
			if( cut_check == 10 ) {
				TH1F *frame_dRp1p2_cut_stack = new TH1F("frame_dRp1p2_cut_stack", "dRp1p2", 20, 0, 10);
				Frame( canvas_dRp1p2_stack_cut , frame_dRp1p2_cut_stack , stack_plot[11] , Min_h , Max_h );
				TLegend *legend_dRp1p2 ;  // Position: (x1, y1, x2, y2)
				Add_legend(legend_dRp1p2  , canvas_dRp1p2_stack_cut , dRp1p2[15] , dRp1p2[3] , dRp1p2[2] , dRp1p2[0] , dRp1p2[1] );
				canvas_dRp1p2_stack_cut->SaveAs("dRp1p2_stack.pdf");
			}

			//|dRj1j2|
			if( cut_check == 11 ) {
				TH1F *frame_dRj1j2_cut_stack = new TH1F("frame_dRj1j2_cut_stack", "dRj1j2",20, 0, 10);
				Frame( canvas_dRj1j2_stack_cut , frame_dRj1j2_cut_stack , stack_plot[12] , Min_h , Max_h );
				TLegend *legend_dRj1j2 ;  // Position: (x1, y1, x2, y2)
				Add_legend(legend_dRj1j2  ,canvas_dRj1j2_stack_cut , dRj1j2[15] , dRj1j2[3] , dRj1j2[2] , dRj1j2[0] , dRj1j2[1] );
				canvas_dRj1j2_stack_cut->SaveAs("dRj1j2_stack.pdf");
			}

			//|dRp1j2|
			if( cut_check == 12 ) {
				TH1F *frame_dRp1j2_cut_stack = new TH1F("frame_dRp1j2_cut_stack", "dRp1j2", 20, 0, 10);
				Frame( canvas_dRp1j2_stack_cut , frame_dRp1j2_cut_stack , stack_plot[13] , Min_h , Max_h );
				TLegend *legend_dRp1j2 ;  // Position: (x1, y1, x2, y2)
				Add_legend(legend_dRp1j2  , canvas_dRp1j2_stack_cut , dRp1j2[15] , dRp1j2[3] , dRp1j2[2] , dRp1j2[0] , dRp1j2[1] );
				canvas_dRj1j2_stack_cut->SaveAs("dRp1j2_stack.pdf");
			}

			//|dRp2j1|
			if( cut_check == 13 ) {
				TH1F *frame_dRp2j1_cut_stack = new TH1F("frame_dRp2j1_cut_stack", "dRp2j1", 20, 0, 10);
				Frame( canvas_dRp2j1_stack_cut , frame_dRp2j1_cut_stack , stack_plot[14] , Min_h , Max_h );
				TLegend *legend_dRp2j1 ;  // Position: (x1, y1, x2, y2)
				Add_legend(legend_dRp2j1  ,canvas_dRp2j1_stack_cut , dRp2j1[15] , dRp2j1[3] , dRp2j1[2] , dRp2j1[0] , dRp2j1[1] );
				canvas_dRp2j1_stack_cut->SaveAs("dRp2j1_stack.pdf");
			}
	*/

	//}
//	canvas_photonpt1_nocut->SaveAs("photonpt1_nocut.pdf");
//	canvas_photonpt2_nocut->SaveAs("photonpt2_nocut.pdf");


	return 0 ;
}
