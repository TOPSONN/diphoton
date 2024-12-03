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

const int nHists = 25;
// Create an array of TH1F histograms
TH1F *photonpt1_nocut[nHists] , *deltaetaj1j2_nocut[nHists];
TH1F *photonpt2_nocut[nHists] , *invmassj1j2_nocut[nHists]  ;

TH1F *photoneta1_cut[nHists] , *photoneta2_cut[nHists] ,  *jeteta1_cut[nHists] , *jeteta2_cut[nHists] ;
TH1F *photonpt1_cut[nHists] , *photonpt2_cut[nHists] ,  *jetpt1_cut[nHists] , *jetpt2_cut[nHists] ;
TH1F *photondeltaeta[nHists] ,   *jetdeltaeta[nHists]  ;
TH1F *photondeltaphi[nHists] , *photonphi2_cut[nHists] ,  *jetdeltaphi[nHists]  ;
TH1F *dRp1p2[nHists] , *dRj1j2[nHists] , *dRp1j2[nHists] , *dRp2j1[nHists] ;
TH1F *mpp_cut[nHists] , *mjj_cut[nHists] ;
TH1F *Zi_origin[nHists] , *DPHI_origin[nHists] ;

TH1F *photonIDMVA1_origin[nHists] , *photonIDMVA2_origin[nHists];

TH1F *photonpt1_cut_all[nHists] , *photonpt2_cut_all[nHists] ,  *jetpt1_cut_all[nHists] , *jetpt2_cut_all[nHists] ;
TH1F *photondeltaeta_cut_all[nHists] , *photoneta2_cut_all[nHists] ,  *jetdeltaeta_cut_all[nHists] , *jeteta2_cut_all[nHists] ;
TH1F *photondeltaphi_cut_all[nHists] , *photonphi2_cut_all[nHists] ,  *jetdeltaphi_cut_all[nHists] , *jetphi2_cut_all[nHists] ;
TH1F *mpp_cut_all[nHists] , *mjj_cut_all[nHists] ;

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

float array_HLT[10][30]= {0};
float total_entry_array[10] = {0};

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
	if (photon_phi > M_PI) photon_phi -= 2 * M_PI;
	if (photon_phi < -M_PI) photon_phi += 2 * M_PI;

	float jet_phi= jet_phi1 + jet_phi2 ;
	if (jet_phi > M_PI) jet_phi -= 2 * M_PI;
	if (jet_phi < -M_PI) jet_phi += 2 * M_PI;

	float dphi = (photon_phi-jet_phi) ;

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


TCanvas* Add_legend( TLegend* legend,  TCanvas* canvas , TH1F* QCD , TH1F* GJETS , TH1F* QCD_Diphoton , TH1F* INTER , TH1F* VBS ) {
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
		if( Z_ < 2.4 && DPHI_ > 1.9  ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 15 ) { //Z_
		if(  DPHI_ > 1.9  ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check == 16 ) { //DPHI_
		if( Z_ < 2.4 ) {
			check_condition = check_condition +1;
		}
	}

	if( cut_check > 16 ) { //cut_all
		if( Z_ < 3 && DPHI_ > 1.9  ) {
			check_condition = check_condition +1;
		}
	}


	return check_condition;
}


int ID_check(int &ID,int QCD_jetid1_flt ,int QCD_jetid2_flt ,int QCD_jetpuid1_flt ,int QCD_jetpuid2_flt , int QCD_phocutbased1_flt , int QCD_phocutbased2_flt , int &pass_entry) {

	if(  QCD_jetid1_flt == 6 && QCD_jetid2_flt ==6 && QCD_jetpuid1_flt == 7 && QCD_jetpuid2_flt == 7  && QCD_phocutbased1_flt ==3 && QCD_phocutbased2_flt ==3   ) {
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
	std::vector<bool> *HLT14 = 0 ,*HLT15 = 0, *HLT16 = 0, *HLT17 = 0 , *HLT18 = 0 ,*HLT19=0 ;


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

	tree->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT14);
	tree->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT15);
	tree->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", &HLT16);
	tree->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55", &HLT17 );
	tree->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55", &HLT18 );
	tree->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto", &HLT19 );

	// Loop over events
	Long64_t nEntries = tree->GetEntries();

	int total_events = tree->GetEntries() ;
	cout << " total entry=  "   <<  total_events  << endl;
	int pass_entry  = 0 ;
	int total_events_ = 100000;
	int pass_entry_photonpt = 0 ;

	TLorentzVector photon1_t , photon2_t , jet1_t , jet2_t ;
	int match_number = 0 ;



	for (int evt = 0 ;  evt <  total_events ; evt++) {

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

		HLT14->clear();
		HLT15->clear();
		HLT16->clear();
		HLT17->clear();
		HLT18->clear();
		HLT19->clear();

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

			// pt_cut
			//cout << " jet_first.jetid= " << jet_first.jetid <<endl;
			if(jet_first.jetid==2 && jet_second.jetpuid==2) cout << " jet_first.jetpuid= " << jet_first.jetpuid <<endl;
			ID_check(ID, jet_first.jetid , jet_second.jetid , jet_first.jetpuid , jet_second.jetpuid  , photon_first.phocutBased , photon_second.phocutBased , pass_entry );
			Filter_first(cut_check , check_condition  , photon_first.pt , photon_second.pt , jet_first.pt , jet_second.pt  );
			Filter_second(cut_check ,check_condition , jet_first.eta , jet_second.eta , Mjj ,  Mpp , deltaetajj , photon_first.eta ,  photon_second.eta );
			Filter_third(cut_check , check_condition ,  dRj1j2_ ,dRp1j1_ ,  dRp1j2_ ,   dRp2j1_ , dRp2j2_ );
			//Filter_fourth(cut_check , check_condition  ,  Z ,  dphi );

			if(check_condition == 3 && ID ==1 ) cout << " ID= " << ID  << "  |check_condition=  " << check_condition << endl;

			if(check_condition  ==3 && ID ==1 ) {

				total_entry_array[indexHT_process]  = total_entry_array[indexHT_process]  + 1;
				//cout << "++++++++++++++"   << "    "  << (*HLT14)[0] << endl;
				if( (*HLT14)[0] == true) array_HLT[indexHT_process][0] = array_HLT[indexHT_process][0] + 1 ;
				if((*HLT15)[0] == true) array_HLT[indexHT_process][1] = array_HLT[indexHT_process][1] + 1 ;
				if((*HLT16)[0] == true) array_HLT[indexHT_process][2] = array_HLT[indexHT_process][2] + 1 ;
				if((*HLT17)[0]== true) array_HLT[indexHT_process][3] = array_HLT[indexHT_process][3] + 1 ;
				if((*HLT18)[0] == true) array_HLT[indexHT_process][4] = array_HLT[indexHT_process][4] + 1 ;
				if((*HLT19)[0] == true) array_HLT[indexHT_process][5] = array_HLT[indexHT_process][5] + 1 ;

				if((*HLT14)[0] == true || (*HLT16)[0] == true) array_HLT[indexHT_process][6] = array_HLT[indexHT_process][6] + 1 ;
				if((*HLT14)[0] == true || (*HLT17)[0] == true) array_HLT[indexHT_process][7] = array_HLT[indexHT_process][7] + 1 ;
				if((*HLT14)[0] == true || (*HLT18)[0] == true) array_HLT[indexHT_process][8] = array_HLT[indexHT_process][8] + 1 ;
				if((*HLT14)[0] == true || (*HLT19)[0] == true) array_HLT[indexHT_process][9] = array_HLT[indexHT_process][9] + 1 ;

				if((*HLT15)[0] == true || (*HLT16)[0] == true) array_HLT[indexHT_process][10] = array_HLT[indexHT_process][10] + 1 ;
				if((*HLT15)[0] == true || (*HLT17)[0] == true ) array_HLT[indexHT_process][11] = array_HLT[indexHT_process][11] + 1 ;
				if((*HLT15)[0] == true || (*HLT18)[0] == true ) array_HLT[indexHT_process][12] = array_HLT[indexHT_process][12] + 1 ;
				if((*HLT15)[0] == true || (*HLT19)[0] == true) array_HLT[indexHT_process][13] = array_HLT[indexHT_process][13] + 1 ;

				if((*HLT16)[0] == true || (*HLT18)[0] == true) array_HLT[indexHT_process][14] = array_HLT[indexHT_process][14] + 1 ;
				if((*HLT16)[0] == true || (*HLT19)[0] == true) array_HLT[indexHT_process][15] = array_HLT[indexHT_process][15] + 1 ;

				if((*HLT17)[0] == true || (*HLT18)[0] == true) array_HLT[indexHT_process][16] = array_HLT[indexHT_process][16] + 1 ;
				if((*HLT17)[0] == true || (*HLT19)[0] == true) array_HLT[indexHT_process][17] = array_HLT[indexHT_process][17] + 1 ;

			}

		}


	}// end of event loop

	//int color = 6 ;
	array_eff[indexHT_process] =  (float) pass_entry/match_number ;
	array_eff_photonpt[indexHT_process] = (float) pass_entry_photonpt/match_number ;

	file->Close();

}




void process(int indexHT , int cut_check) {

	if(indexHT == 2) { //qcd-diphotonjets 2
		//TFile *file = TFile::Open("vbs_filter.root"); TTree *tree = (TTree*)file->Get("Event");
		const double weighted = 2.327106 ;
		string color = "kGreen+3" ;
		string string_inputroot = "qcd_filter_nocutHLT.root" ;
		string string_inputtree = "Event" ;
		process_VBS( string_inputroot  , string_inputtree  , weighted, indexHT  ,  cut_check , color ) ;
	}

	if(indexHT == 0) { //interference 0
		//TFile *file = TFile::Open("vbs_filter.root"); TTree *tree = (TTree*)file->Get("Event");
		const double weighted =0.000868038 ;
		string color = "kAzure" ;
		string string_inputroot = "inter_filter_nocutHLT.root" ;
		string string_inputtree = "Event" ;
		process_VBS( string_inputroot  , string_inputtree  , weighted, indexHT  ,  cut_check , color ) ;
	}


	if(indexHT == 1) { //VBS 1
		//TFile *file = TFile::Open("vbs_filter.root"); TTree *tree = (TTree*)file->Get("Event");
		const double weighted =  0.005320464 ;
		string color = "kRed+1" ;
		string string_inputroot = "vbs_filter_nocutHLT.root" ;
		string string_inputtree = "Event" ;
		process_VBS( string_inputroot  , string_inputtree  , weighted, indexHT  ,  cut_check , color ) ;
	}

}






int plot_HLT() {


	for(int k1 = 0 ; k1<20 ; k1++ ) {

		stack_plot[k1] = new THStack("stack_plot", "stack_plot");
	}


	int cut_check=30;


	for(int i = 0 ; i < 3 ; i++) {

		process(i , cut_check );



	}


	for(int i = 0 ; i < 3 ; i++) {

		if(i==0) {
			cout << "-------------------------------" <<endl;
			cout <<"            interference             " <<endl;
			cout <<"                                 " <<endl;
			cout <<"  total_entry_array[i]=   " <<   total_entry_array[i]   << endl;
			cout <<  array_HLT[i][0]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][1]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][2]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][3]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][4]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][5]/total_entry_array[i] <<endl;
			cout <<"              14                   " <<endl;
			cout <<  array_HLT[i][6]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][7]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][8]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][9]/total_entry_array[i] <<endl;
			cout <<"              15                   " <<endl;
			cout <<  array_HLT[i][10]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][11]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][12]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][13]/total_entry_array[i] <<endl;
			cout <<"              16        17            " <<endl;
			cout <<  array_HLT[i][14]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][15]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][16]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][17]/total_entry_array[i] <<endl;
			cout << "-------------------------------" <<endl;

		}

		if(i==1) {
			cout << "-------------------------------" <<endl;
			cout <<"           VBS            " <<endl;
			cout <<"                                 " <<endl;
			cout <<"  total_entry_array[i]=   " <<   total_entry_array[i]   << endl;
			cout <<  array_HLT[i][0]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][1]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][2]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][3]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][4]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][5]/total_entry_array[i] <<endl;
			cout <<"              14                   " <<endl;
			cout <<  array_HLT[i][6]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][7]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][8]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][9]/total_entry_array[i] <<endl;
			cout <<"              15                   " <<endl;
			cout <<  array_HLT[i][10]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][11]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][12]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][13]/total_entry_array[i] <<endl;
			cout <<"              16        17            " <<endl;
			cout <<  array_HLT[i][14]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][15]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][16]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][17]/total_entry_array[i] <<endl;
			cout << "-------------------------------" <<endl;

		}

		if(i==2) {
			cout << "-------------------------------" <<endl;
			cout <<"    QCD_Diphtonjets      " <<endl;
			cout <<"                                 " <<endl;
			cout <<"  total_entry_array[i]=   " <<   total_entry_array[i]   << endl;
			cout <<  array_HLT[i][0]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][1]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][2]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][3]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][4]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][5]/total_entry_array[i] <<endl;
			cout <<"              14                   " <<endl;
			cout <<  array_HLT[i][6]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][7]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][8]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][9]/total_entry_array[i] <<endl;
			cout <<"              15                   " <<endl;
			cout <<  array_HLT[i][10]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][11]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][12]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][13]/total_entry_array[i] <<endl;
			cout <<"              16        17            " <<endl;
			cout <<  array_HLT[i][14]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][15]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][16]/total_entry_array[i] <<endl;
			cout <<  array_HLT[i][17]/total_entry_array[i] <<endl;
			cout << "-------------------------------" <<endl;

		}

	}


	float Max_h = 1E6 ;
	float Min_h = 8E-5 ;





//	canvas_photonpt1_nocut->SaveAs("photonpt1_nocut.pdf");
//	canvas_photonpt2_nocut->SaveAs("photonpt2_nocut.pdf");


	return 0 ;
}
