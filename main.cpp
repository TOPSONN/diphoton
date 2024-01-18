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

#include <cstdio>
#include <stdio.h>
#include <string>
#include <cctype>
#include <algorithm>
#include <iterator>

void sort_(Float_t *a, int n) {

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (a[i] <= a[j]) {
				Float_t temp;
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}


}




int main() {

	int n1;
	cout << "max";
	cin >> n1  ;




	//////////////////////////// define plot
	TH1F *deltapt_p1_j1 = new TH1F("p1j1_pt","p1j1_pt",50,-150,150);
	TH1F *deltapt_p1_j2 = new TH1F("p1j2_pt","p1j2_pt",50,-150,150);
	TH1F *deltapt_p2_j1 = new TH1F("p2j1_pt","p2j1_pt",50,-150,150);
	TH1F *deltapt_p2_j2 = new TH1F("p2j2_pt","p2j2_pt",50,-150,150);

	TH1F *deltaeta_p1_j1 = new TH1F("p1j1_eta","p1j1_eta",40,-10,10);
	TH1F *deltaeta_p1_j2 = new TH1F("p1j2_eta","p1j2_eta",40,-10,10);
	TH1F *deltaeta_p2_j1 = new TH1F("p2j1_eta","p2j1_eta",40,-10,10);
	TH1F *deltaeta_p2_j2 = new TH1F("p2j2_eta","p2j2_eta",40,-10,10);

	TH1F *deltaphi_p1_j1 = new TH1F("p1j1_phi","p1j1_phi",50,-4,4);
	TH1F *deltaphi_p1_j2 = new TH1F("p1j2_phi","p1j2_phi",50,-4,4);
	TH1F *deltaphi_p2_j1 = new TH1F("p2j1_phi","p2j1_phi",50,-4,4);
	TH1F *deltaphi_p2_j2 = new TH1F("p2j2_phi","p2j2_phi",50,-4,4);


	TH1F *deltaphi_p1_p2 = new TH1F("p1p2_phi","p1p2_phi",30,-4,4);
	TH1F *deltaphi_j1_j2 = new TH1F("j1j2_phi","j1j2_phi",30,-4,4);
	TH1F *deltaphi_test= new TH1F("test","test",100,-5,5);


	TH1F *deltaeta_p1_p2 = new TH1F("p1p2_eta","p1p2_eta",40,-10,10);
	TH1F *deltaeta_j1_j2 = new TH1F("j1j2_eta","j1j2_eta",40,-10,10);

	TH1F *deltapt_p1_p2 = new TH1F("p1p2_pt","p1p2_pt",50,-150,150);
	TH1F *deltapt_j1_j2 = new TH1F("j1j2_pt","j1j2_pt",50,-150,150);

	TH1F *deltaR_p1_p2 = new TH1F("p1p2_R","p1p2_R",50,0,10);
	TH1F *deltaR_j1_j2 = new TH1F("j1j2_R","j1j2_R",50,0,10);

	TH1F *deltaR_photon = new TH1F("photon(gen)_R","photon(gen)_R",100,0,10);
	TH1F *deltaR_jet = new TH1F("jet(genjet)_R","jet(genjet)_R",100,0,10);

	TH1F *deltaR_p1_j1 = new TH1F("p1j1_R","p1j1_R",50,0,10);
	TH1F *deltaR_p1_j2 = new TH1F("p1j2_R","p1j2_R",50,0,10);
	TH1F *deltaR_p2_j1 = new TH1F("p2j1_R","p2j1_R",50,0,10);
	TH1F *deltaR_p2_j2 = new TH1F("p2j2_R","p2j2_R",50,0,10);




	for(int i=1 ; i<= n1 ; i++ ) {

		string file0 = "inter_nanoAOD"+ std::to_string(i)  +".root";
		cout << file0 <<endl;
		TFile *f = new TFile(file0.c_str());
		//	TFile *f = new TFile(file0,"read");

		TTree *tree;
		f->GetObject("Events;1",tree);
		cout << "--------------------" <<endl;
		cout << tree->GetEntries() <<endl;


//		vector<Float_t> *INTER_phopt=0;
//		vector<Float_t> *INTER_phoeta=0 ;
//		vector<Float_t> *INTER_phophi=0 ;
//		vector<Float_t> *INTER_phom=0 ;

		Float_t INTER_phopt[100]= { 0 }, INTER_phoeta[100]= { 0 }, INTER_phophi[100]= { 0 }, INTER_phom[100]= { 0 };
		Float_t INTER_genpt[100]= { 0 }, INTER_geneta[100]= { 0 }, INTER_genphi[100]= { 0 }, INTER_genm[100]= { 0 };
		Float_t INTER_jetpt[100]= { 0 }, INTER_jeteta[100]= { 0 }, INTER_jetphi[100]= { 0 }, INTER_jetm[100]= { 0 };
		Float_t INTER_gjetpt[100]= { 0 }, INTER_gjeteta[100]= { 0 }, INTER_gjetphi[100]= { 0 }, INTER_gjetm[100]= { 0 };
		Int_t INTER_jetid[100]= { 0 }, INTER_genid[100]= { 0 }, INTER_genstatus[100]= { 0 }, INTER_genmomid[100]= { 0 }, INTER_phoid[100]= { 0 };


		tree->SetBranchAddress("Photon_pt"  ,   &INTER_phopt  );
		tree->SetBranchAddress("Photon_eta"  ,   &INTER_phoeta  );
		tree->SetBranchAddress("Photon_phi"  ,   &INTER_phophi  );
		tree->SetBranchAddress("Photon_mass"  ,   &INTER_phom  );

		tree->SetBranchAddress("GenPart_pt"  ,   &INTER_genpt  );
		tree->SetBranchAddress("GenPart_eta"  ,   &INTER_geneta  );
		tree->SetBranchAddress("GenPart_phi"  ,   &INTER_genphi  );
		tree->SetBranchAddress("GenPart_mass"  ,   &INTER_genm  );

		tree->SetBranchAddress("Jet_pt"  ,   &INTER_jetpt  );
		tree->SetBranchAddress("Jet_eta"  ,   &INTER_jeteta  );
		tree->SetBranchAddress("Jet_phi"  ,   &INTER_jetphi  );
		tree->SetBranchAddress("Jet_mass"  ,   &INTER_jetm  );

		tree->SetBranchAddress("GenJet_pt"  ,   &INTER_gjetpt  );
		tree->SetBranchAddress("GenJet_eta"  ,   &INTER_gjeteta  );
		tree->SetBranchAddress("GenJet_phi"  ,   &INTER_gjetphi  );
		tree->SetBranchAddress("GenJet_mass"  ,   &INTER_gjetm  );


		tree->SetBranchAddress("Jet_jetId"  ,   &INTER_jetid  );
		tree->SetBranchAddress("GenPart_pdgId"  ,   &INTER_genid  );
		tree->SetBranchAddress("GenPart_status"  ,   &INTER_genstatus  );
		tree->SetBranchAddress("GenPart_genPartIdxMother"  ,   &INTER_genmomid  );
		tree->SetBranchAddress("Photon_pdgId"  ,   &INTER_phoid  );

		TLorentzVector INTER_photon1,INTER_photon2 , INTER_jet1 , INTER_jet2;
		TLorentzVector INTER_genpart1 , INTER_genpart2 ;

		TLorentzVector INTER_photon , INTER_genpart ;
		TLorentzVector INTER_jet , INTER_genjet ;

		TLorentzVector INTER_true_photon1 , INTER_true_photon2;
		TLorentzVector INTER_true_jet1 , INTER_true_jet2;


		for (int evt = 0; evt <tree->GetEntries() ; evt++) {

			double delta_eta_photon1_jet1=0;
			double delta_eta_photon1_jet2=0;
			double delta_eta_photon2_jet1=0;
			double delta_eta_photon2_jet2=0;
			double delta_phi_photon1_jet1=0;
			double delta_phi_photon1_jet2=0;
			double delta_phi_photon2_jet1=0;
			double delta_phi_photon2_jet2=0;
			double delta_pt_photon1_jet1=0;
			double delta_pt_photon1_jet2=0;
			double delta_pt_photon2_jet1=0;
			double delta_pt_photon2_jet2=0;

			double delta_pt_photon1_photon2=0;
			double delta_pt_jet1_jet2=0;
			double delta_phi_photon1_photon2=0;
			double delta_phi_jet1_jet2=0;
			double delta_eta_photon1_photon2=0;
			double delta_eta_jet1_jet2=0;


			double delta_R_photon1_photon2 = 0;
			double delta_R_jet1_jet2 = 0;
			double delta_R_photon1_jet1 =0;
			double delta_R_photon1_jet2 = 0;
			double delta_R_photon2_jet1 = 0;
			double delta_R_photon2_jet2 = 0;

			INTER_photon1.SetPtEtaPhiM(0,0,0,0) ;
			INTER_photon2.SetPtEtaPhiM(0,0,0,0) ;
			INTER_genpart1.SetPtEtaPhiM(0,0,0,0) ;
			INTER_genpart2.SetPtEtaPhiM(0,0,0,0) ;

			INTER_jet1.SetPtEtaPhiM(0,0,0,0) ;
			INTER_jet2.SetPtEtaPhiM(0,0,0,0) ;








			int pho_size =0;
			int jet_size =0;
			int genpart_size =0;
			int genjet_size =0;

			tree->GetEntry(evt);

			bool check =true;
			bool check_ =true;

			//-----------calcualte size of array photon-------------------//
			for(int a1=0 ; a1<sizeof(INTER_phopt) / sizeof(INTER_phopt[0]) ; a1++) {
				if(INTER_phopt[a1]!=0)  pho_size++ ;

			}


			//---------------calcualte size of array genpart---------------//
			for(int a2=0 ; a2<sizeof(INTER_genpt) / sizeof(INTER_genpt[0]) ; a2++) {
				if(INTER_genpt[a2]!=0)  genpart_size++ ;
			}


			//---------------calcualte size of array jet---------------//
			for(int a3=0 ; a3<sizeof(INTER_jetpt) / sizeof(INTER_jetpt[0]) ; a3++) {
				if(INTER_jetpt[a3]!=0)  jet_size++ ;
			}

			//---------------calcualte size of array genjet---------------//
			for(int a4=0 ; a4<sizeof(INTER_gjetpt) / sizeof(INTER_gjetpt[0]) ; a4++) {
				if(INTER_gjetpt[a4]!=0)  genjet_size++ ;
			}


			//----------------------------------------------------------------------//

			//------------------------------large to small sort photon----------------------------------------//

			for (int p1 = 0; p1 < pho_size; p1++) {
				for (int p2 = p1 + 1; p2 < pho_size; p2++) {
					if (INTER_phopt[p1] <= INTER_phopt[p2]) {
						Float_t temp = 0 ;
						temp = INTER_phopt[p1];
						INTER_phopt[p1] = INTER_phopt[p2];
						INTER_phopt[p2] = temp;

						temp = 0 ;
						temp = INTER_phoeta[p1];
						INTER_phoeta[p1] = INTER_phoeta[p2];
						INTER_phoeta[p2] = temp;

						temp = 0 ;
						temp = INTER_phophi[p1];
						INTER_phophi[p1] = INTER_phophi[p2];
						INTER_phophi[p2] = temp;

						temp = 0 ;
						temp = INTER_phom[p1];
						INTER_phom[p1] = INTER_phom[p2];
						INTER_phom[p2] = temp;



					}
				}
			}


			//----------------------------------------------------------------------//





			// for of pho ----------------------------------------------------------//
			for(int l1=0 ; l1<jet_size ; l1++) {
				//cout << l1 << " ; "<< INTER_phopt[l1] << endl;
				//			if() {
				//				INTER_photon1.SetPtEtaPhiM(0,0,0,0) ;
				//				INTER_photon2.SetPtEtaPhiM(0,0,0,0) ;

				//			}
				INTER_photon.SetPtEtaPhiM(INTER_phopt[l1],INTER_phoeta[l1],INTER_phophi[l1],INTER_phom[l1]) ;



				for(int l2=0 ; l2<genjet_size; l2++) { // for of genpart

					if(INTER_genstatus[l2] ==1 && abs(INTER_genid[l2]) ==22) {
						INTER_genpart.SetPtEtaPhiM(INTER_genpt[l2],INTER_geneta[l2],INTER_genphi[l2],INTER_genm[l2]) ;

						double deltaR = INTER_photon.DeltaR(INTER_genpart);
						deltaR_photon->Fill(deltaR);

						if(deltaR < 0.3) {
							if(l1 == 0 )  INTER_true_photon1.SetPtEtaPhiM(INTER_phopt[l1],INTER_phoeta[l1],INTER_phophi[l1],INTER_phom[l1]) ;
							if(l2 == 1 )  INTER_true_photon2.SetPtEtaPhiM(INTER_phopt[l1],INTER_phoeta[l1],INTER_phophi[l1],INTER_phom[l1]) ;

							break;
						}
					}

				}
			}
			
			

			for(int l1=0 ; l1<pho_size ; l1++) {
				//cout << l1 << " ; "<< INTER_phopt[l1] << endl;
				//			if() {
				//				INTER_photon1.SetPtEtaPhiM(0,0,0,0) ;
				//				INTER_photon2.SetPtEtaPhiM(0,0,0,0) ;

				//			}
				INTER_photon.SetPtEtaPhiM(INTER_phopt[l1],INTER_phoeta[l1],INTER_phophi[l1],INTER_phom[l1]) ;



				for(int l2=0 ; l2<genpart_size ; l2++) { // for of genpart

					if(INTER_genstatus[l2] ==1 && abs(INTER_genid[l2]) ==22) {
						INTER_genpart.SetPtEtaPhiM(INTER_genpt[l2],INTER_geneta[l2],INTER_genphi[l2],INTER_genm[l2]) ;

						double deltaR = INTER_photon.DeltaR(INTER_genpart);
						deltaR_photon->Fill(deltaR);

						if(deltaR < 0.3) {
							if(l1 == 0 )  INTER_true_photon1.SetPtEtaPhiM(INTER_phopt[l1],INTER_phoeta[l1],INTER_phophi[l1],INTER_phom[l1]) ;
							if(l2 == 1 )  INTER_true_photon2.SetPtEtaPhiM(INTER_phopt[l1],INTER_phoeta[l1],INTER_phophi[l1],INTER_phom[l1]) ;

							break;
						}
					}

				}
			}






		}




	}


	TCanvas *cl01 = new TCanvas("cl01","photon_deltaR");

	deltaR_photon->SetStats(0);
	deltaR_photon->SetLineColor(3);
	deltaR_photon->SetLineWidth(3);
//	deltaR_photon->Scale(1/deltaR_photon->Integral());
//	deltaR_photon->GetYaxis()->SetRangeUser(0.0001,100000);
	deltaR_photon->GetXaxis()->SetRangeUser(-1,20);
	deltaR_photon->Draw("hist");


	cl01->SaveAs("30_deltaR_photon.pdf");
	return 0 ;


	////////////////////////////////////////////////////////////













}
