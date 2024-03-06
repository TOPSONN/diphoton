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



void nano2() {


	/////////define plot INTER/////////////////// define plot INTER /////////////////////////////define plot INTER
	TH1F *INTER_deltaeta_j1_j2 = new TH1F("INTER_j1j2_eta","INTER_j1j2_eta",100,-10,10);
	TH1F *INTER_deltaeta_j1_j2_lhe = new TH1F("INTER_j1j2_eta_ori","INTER_j1j2_eta_lhe",100,-10,10);

	TH1F *INTER_deltaphi_j1_j2 = new TH1F("INTER_j1j2_phi","INTER_j1j2_phi",100,-4,4);
	TH1F *INTER_deltaphi_j1_j2_lhe = new TH1F("INTER_j1j2_phi_ori","INTER_j1j2_phi_lhe",100,-4,4);

	TH1F *INTER_deltapt_j1_j2 = new TH1F("INTER_p1j2_pt","INTER_p1j2_pt",100,-1,100);
	TH1F *INTER_deltapt_j1_j2_lhe = new TH1F("INTER_p1j2_pt_ori","INTER_p1j2_pt_lhe",100,-1,100);

	TH1F *INTER_deltaR_j1_j2 = new TH1F("INTER_jnanolhe","INTER_jnanolhe",50,0,10);
	TH1F *INTER_deltaR_j1_j2_lhe = new TH1F("INTER_j1j2_R_ori","INTER_j1j2_R_lhe",100,0,10);
	TH1F *INTER_deltaR_p = new TH1F("INTER_pnanolhe","INTER_pnanolhe",50,0,10);

	TH1F *INTER_deltaetagen= new TH1F("INTER_etagen","INTER_eta",100,-4,4);
	TH1F *INTER_deltaphigen= new TH1F("INTER_phigen","INTER_phi",100,-4,4);

	TH1F *INTER_deltaR_p_genreco = new TH1F("INTER_pgenreco","INTER_pgenreco",50,0,1);
	TH1F *INTER_deltaR_j_genreco = new TH1F("INTER_jgenreco","INTER_jgenreco",50,0,1);

	TH1F *INTER_genphotrue_eta_plot=new TH1F("INTER_genphotrue_eta_plot","INTER_genphotrue_eta_plot",20,-4,4);
	TH1F *INTER_genphotrue_phi_plot=new TH1F("INTER_genphotrue_phi_plot","INTER_genphotrue_phi_plot",20,-3.5,3.5);
	TH1F *INTER_genphotrue_pt_plot=new TH1F("INTER_genphotrue_pt_plot","INTER_genphotrue_pt_plot",50,0,100);
	TH1F *INTER_genphotrue_m_plot=new TH1F("INTER_genphotrue_m_plot","INTER_genphotrue_m_plot",50,0,50);

	TH1F *INTER_genjettrue_eta_plot=new TH1F("INTER_genjettrue_eta_plot","INTER_genjettrue_eta_plot",100,-10,10);
	TH1F *INTER_genjettrue_phi_plot=new TH1F("INTER_genjettrue_phi_plot","INTER_genjettrue_phi_plot",20,-3.5,3.5);
	TH1F *INTER_genjettrue_pt_plot=new TH1F("INTER_genjettrue_pt_plot","INTER_genjettrue_pt_plot",50,0,100);
	TH1F *INTER_genjettrue_m_plot=new TH1F("INTER_genjettrue_m_plot","INTER_genjettrue_m_plot",100,0,1000);

	TH1F *INTER_genjettrue_eta_j1j2_plot=new TH1F("INTER_genjettruej_eta_plot","INTER_genjettruej_eta_plot",20,-5,5);
	TH1F *INTER_genjettrue_phi_j1j2_plot=new TH1F("INTER_genjettruej_phi_plot","INTER_genjettruej_phi_plot",20,-3.5,3.5);
	TH1F *INTER_genjettrue_pt_j1j2_plot=new TH1F("INTER_genjettruej_pt_plot","INTER_genjettruej_pt_plot",50,-10,100);
	TH1F *INTER_genjettrue_deltaR_j1j2_plot=new TH1F("INTER_genjettruejR_plot","INTER_genjettruejR_plot",50,0,10);
	TH1F *INTER_genjettrue_mjj_plot=new TH1F("INTER_genjettruej_m_plot","INTER_genjettruej_m_plot",100,0,1000);

	TH1F *INTER_genphotrue_eta_p1p2_plot=new TH1F("INTER_genphotruep_eta_plot","INTER_genphotruep_eta_plot",20,-5,5);
	TH1F *INTER_genphotrue_phi_p1p2_plot=new TH1F("INTER_genphotruep_phi_plot","INTER_genphotruep_phi_plot",20,-3.5,3.5);
	TH1F *INTER_genphotrue_pt_p1p2_plot=new TH1F("INTER_genphotruep_pt_plot","INTER_genphotruep_pt_plot",50,-10,100);
	TH1F *INTER_genphotrue_deltaR_p1p2_plot=new TH1F("INTER_genphotruepR_plot","INTER_genphotruepR_plot",50,0,10);
	TH1F *INTER_genphotrue_mpp_plot=new TH1F("INTER_genphotruep_m_plot","INTER_genphotruep_m_plot",100,0,1000);



/////////define plot VBS/////////////////// define plot VBS/////////////////////////////define plot VBS
	TH1F *VBS_deltaeta_j1_j2 = new TH1F("VBS_j1j2_eta","VBS_j1j2_eta",100,-10,10);
	TH1F *VBS_deltaeta_j1_j2_lhe = new TH1F("VBS_j1j2_eta_ori","VBS_j1j2_eta_lhe",100,-10,10);

	TH1F *VBS_deltaphi_j1_j2 = new TH1F("VBS_j1j2_phi","VBS_j1j2_phi",100,-4,4);
	TH1F *VBS_deltaphi_j1_j2_lhe = new TH1F("VBS_j1j2_phi_ori","VBS_j1j2_phi_lhe",100,-4,4);

	TH1F *VBS_deltapt_j1_j2 = new TH1F("VBS_p1j2_pt","VBS_p1j2_pt",100,-1,100);
	TH1F *VBS_deltapt_j1_j2_lhe = new TH1F("VBS_p1j2_pt_ori","VBS_p1j2_pt_lhe",100,-1,100);

	TH1F *VBS_deltaR_j1_j2 = new TH1F("VBS_j1j2_R","VBS_j1j2_R",50,0,10);
	TH1F *VBS_deltaR_j1_j2_lhe = new TH1F("VBS_j1j2_R_ori","VBS_j1j2_R_lhe",50,0,10);
	TH1F *VBS_deltaR_p = new TH1F("VBS_p","VBS_p",50,0,10);

	TH1F *VBS_deltaR_p_genreco = new TH1F("VBS_pgenreco","VBS_pgenreco",50,0,1);
	TH1F *VBS_deltaR_j_genreco = new TH1F("VBS_jgenreco","VBS_jgenreco",50,0,1);

	TH1F *VBS_genphotrue_eta_plot=new TH1F("VBS_genphotrue_eta_plot","VBS_genphotrue_eta_plot",20,-4,4);
	TH1F *VBS_genphotrue_phi_plot=new TH1F("VBS_genphotrue_phi_plot","VBS_genphotrue_phi_plot",20,-3.5,3.5);
	TH1F *VBS_genphotrue_pt_plot=new TH1F("VBS_genphotrue_pt_plot","VBS_genphotrue_pt_plot",50,0,100);
	TH1F *VBS_genphotrue_m_plot=new TH1F("VBS_genphotrue_m_plot","VBS_genphotrue_m_plot",50,0,50);

	TH1F *VBS_genjettrue_eta_plot=new TH1F("VBS_genjettrue_eta_plot","VBS_genjettrue_eta_plot",50,-10,10);
	TH1F *VBS_genjettrue_phi_plot=new TH1F("VBS_genjettrue_phi_plot","VBS_genjettrue_phi_plot",20,-3.5,3.5);
	TH1F *VBS_genjettrue_pt_plot=new TH1F("VBS_genjettrue_pt_plot","VBS_genjettrue_pt_plot",50,0,100);
	TH1F *VBS_genjettrue_m_plot=new TH1F("VBS_genjettrue_m_plot","VBS_genjettrue_m_plot",100,0,1000);

	TH1F *VBS_genjettrue_eta_j1j2_plot=new TH1F("VBS_genjettruej_eta_plot","VBS_genjettruej_eta_plot",20,-5,5);
	TH1F *VBS_genjettrue_phi_j1j2_plot=new TH1F("VBS_genjettruej_phi_plot","VBS_genjettruej_phi_plot",20,-3.5,3.5);
	TH1F *VBS_genjettrue_pt_j1j2_plot=new TH1F("VBS_genjettruej_pt_plot","VBS_genjettruej_pt_plot",50,-10,100);
	TH1F *VBS_genjettrue_deltaR_j1j2_plot=new TH1F("VBS_genjettruej_pt_plot","VBS_genjettruej_pt_plot",50,0,10);
	TH1F *VBS_genjettrue_mjj_plot=new TH1F("VBS_genjettruej_m_plot","VBS_genjettruej_m_plot",100,0,1000);

	TH1F *VBS_genphotrue_eta_p1p2_plot=new TH1F("VBS_genphotruep_eta_plot","VBS_genphotruep_eta_plot",20,-5,5);
	TH1F *VBS_genphotrue_phi_p1p2_plot=new TH1F("VBS_genphotruep_phi_plot","VBS_genphotruep_phi_plot",20,-3.5,3.5);
	TH1F *VBS_genphotrue_pt_p1p2_plot=new TH1F("VBS_genphotruep_pt_plot","VBS_genphotruep_pt_plot",50,-10,100);
	TH1F *VBS_genphotrue_deltaR_p1p2_plot=new TH1F("VBS_genphotruepR_plot","VBS_genphotruepR_plot",50,0,10);
	TH1F *VBS_genphotrue_mpp_plot=new TH1F("VBS_genphotruep_m_plot","VBS_genphotruep_m_plot",100,0,1000);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(int i=1 ; i<= 9 ; i++ ) {

		string file0 = "inter_nanoAOD" + std::to_string(i)  +".root";
		cout << file0 <<endl;
		TFile *f = new TFile(file0.c_str());
		TTree *tree;
		f->GetObject("Events;1",tree);
		Float_t INTER_phopt[10]= { 0 }, INTER_phoeta[10]= { 0 }, INTER_phophi[10]= { 0 }, INTER_phom[10]= { 0 };
		Float_t INTER_genpt[30]= { 0 }, INTER_geneta[30]= { 0 }, INTER_genphi[30]= { 0 }, INTER_genm[30]= { 0 };
		Float_t INTER_jetpt[10]= { 0 }, INTER_jeteta[10]= { 0 }, INTER_jetphi[10]= { 0 }, INTER_jetm[10]= { 0 };
		Float_t INTER_gjetpt[10]= { 0 }, INTER_gjeteta[10]= { 0 }, INTER_gjetphi[10]= { 0 }, INTER_gjetm[10]= { 0 };
		Int_t INTER_jetid[10]= { 0 }, INTER_genid[10]= { 0 }, INTER_genstatus[10]= { 0 }, INTER_genmomid[10]= { 0 }, INTER_phoid[10]= { 0 };
		Float_t INTER_r9[10]= { 0 }, INTER_hoe[10]= { 0 }, INTER_phiso[100]= { 0 }, INTER_chiso[10]= { 0 }, INTER_sieie[10]= { 0 } ;
		Float_t INTER_gphopt[10]= { 0 }, INTER_gphoeta[10]= { 0 }, INTER_gphophi[100]= { 0 }, INTER_gphom[10]= { 0 };
		Float_t INTER_deltaR_j1[10]= { 0 }, INTER_deltaR_j2[10]= { 0 } ,INTER_deltaR_p1[10]= { 0 } , INTER_deltaR_p2[100]= { 0 };
		Float_t INTER_genphotrue_pt[10]= { 0 } , INTER_genphotrue_eta[10]= { 0 } , INTER_genphotrue_phi[10]= { 0 } , INTER_genphotrue_m[10]= { 0 };
		Float_t INTER_genjettrue_pt[10]= { 0 } , INTER_genjettrue_eta[10]= { 0 } , INTER_genjettrue_phi[10]= { 0 } , INTER_genjettrue_m[10]= { 0 };

		tree->SetBranchAddress("Photon_pt",   &INTER_phopt  );
		tree->SetBranchAddress("Photon_eta",   &INTER_phoeta  );
		tree->SetBranchAddress("Photon_phi",   &INTER_phophi  );
		tree->SetBranchAddress("Photon_mass",   &INTER_phom  );

		tree->SetBranchAddress("Photon_r9",   &INTER_r9  );
		tree->SetBranchAddress("Photon_hoe",   &INTER_hoe  );
		tree->SetBranchAddress("Photon_pfRelIso03_all",   &INTER_phiso  );
		tree->SetBranchAddress("Photon_pfRelIso03_chg",   &INTER_chiso  );
		tree->SetBranchAddress("Photon_sieie",   &INTER_sieie  );

		tree->SetBranchAddress("GenPart_pt",   &INTER_genpt  );
		tree->SetBranchAddress("GenPart_eta",   &INTER_geneta  );
		tree->SetBranchAddress("GenPart_phi",   &INTER_genphi  );
		tree->SetBranchAddress("GenPart_mass",   &INTER_genm  );

		tree->SetBranchAddress("Jet_pt",   &INTER_jetpt  );
		tree->SetBranchAddress("Jet_eta",   &INTER_jeteta  );
		tree->SetBranchAddress("Jet_phi",   &INTER_jetphi  );
		tree->SetBranchAddress("Jet_mass",   &INTER_jetm  );

		tree->SetBranchAddress("GenJet_pt",   &INTER_gjetpt  );
		tree->SetBranchAddress("GenJet_eta",   &INTER_gjeteta  );
		tree->SetBranchAddress("GenJet_phi",   &INTER_gjetphi  );
		tree->SetBranchAddress("GenJet_mass",   &INTER_gjetm  );


		tree->SetBranchAddress("Jet_jetId",   &INTER_jetid  );
		tree->SetBranchAddress("GenPart_pdgId",   &INTER_genid  );
		tree->SetBranchAddress("GenPart_status",   &INTER_genstatus  );
		tree->SetBranchAddress("GenPart_genPartIdxMother",   &INTER_genmomid  );
		tree->SetBranchAddress("Photon_pdgId",   &INTER_phoid  );

		TLorentzVector INTER_photon1,INTER_photon2, INTER_jet1, INTER_jet2;
		TLorentzVector INTER_genpart1, INTER_genpart2 ;

		TLorentzVector INTER_photon, INTER_genpart ;
		TLorentzVector INTER_jet, INTER_genjet ;

		TLorentzVector INTER_true_photon1, INTER_true_photon2;
		TLorentzVector INTER_true_jet1, INTER_true_jet2;

		TLorentzVector INTER_gjet1, INTER_gjet2;
		TLorentzVector INTER_gphoton1,INTER_gphoton2;

		TLorentzVector INTER_q1_gen_sph , INTER_q2_gen_sph ;

		////////////////////////////////////////////////////////
		string file01 = "vbs_nanoAOD"+ std::to_string(i)  +".root";
		cout << file01 <<endl;
		TFile *f1 = new TFile(file01.c_str());

		TTree *tree1;
		f1->GetObject("Events;1",tree1);


		Float_t VBS_phopt[10]= { 0 }, VBS_phoeta[10]= { 0 }, VBS_phophi[10]= { 0 }, VBS_phom[10]= { 0 };
		Float_t VBS_genpt[30]= { 0 }, VBS_geneta[30]= { 0 }, VBS_genphi[30]= { 0 }, VBS_genm[30]= { 0 };
		Float_t VBS_jetpt[10]= { 0 }, VBS_jeteta[10]= { 0 }, VBS_jetphi[10]= { 0 }, VBS_jetm[10]= { 0 };
		Float_t VBS_gjetpt[10]= { 0 }, VBS_gjeteta[10]= { 0 }, VBS_gjetphi[10]= { 0 }, VBS_gjetm[10]= { 0 };
		Int_t VBS_jetid[10]= { 0 }, VBS_genid[10]= { 0 }, VBS_genstatus[10]= { 0 }, VBS_genmomid[10]= { 0 }, VBS_phoid[10]= { 0 };
		Float_t VBS_r9[10]= { 0 }, VBS_hoe[10]= { 0 }, VBS_phiso[10]= { 0 }, VBS_chiso[10]= { 0 }, VBS_sieie[10]= { 0 } ;
		Float_t VBS_gphopt[10]= { 0 }, VBS_gphoeta[10]= { 0 }, VBS_gphophi[10]= { 0 }, VBS_gphom[10]= { 0 };
		Float_t VBS_deltaR_j1[10]= { 0 }, VBS_deltaR_j2[10]= { 0 }, VBS_deltaR_p1[10]= { 0 } , VBS_deltaR_p2[10]= { 0 };
		Float_t VBS_genphotrue_pt[10]= { 0 } , VBS_genphotrue_eta[10]= { 0 } , VBS_genphotrue_phi[10]= { 0 } , VBS_genphotrue_m[10]= { 0 };
		Float_t VBS_genjettrue_pt[10]= { 0 } , VBS_genjettrue_eta[10]= { 0 } , VBS_genjettrue_phi[10]= { 0 } , VBS_genjettrue_m[10]= { 0 };

		tree1->SetBranchAddress("Photon_pt",   &VBS_phopt  );
		tree1->SetBranchAddress("Photon_eta",   &VBS_phoeta  );
		tree1->SetBranchAddress("Photon_phi",   &VBS_phophi  );
		tree1->SetBranchAddress("Photon_mass",   &VBS_phom  );

		tree1->SetBranchAddress("Photon_r9",   &VBS_r9  );
		tree1->SetBranchAddress("Photon_hoe",   &VBS_hoe  );
		tree1->SetBranchAddress("Photon_pfRelIso03_all",   &VBS_phiso  );
		tree1->SetBranchAddress("Photon_pfRelIso03_chg",   &VBS_chiso  );
		tree1->SetBranchAddress("Photon_sieie",   &VBS_sieie  );

		tree1->SetBranchAddress("GenPart_pt",   &VBS_genpt  );
		tree1->SetBranchAddress("GenPart_eta",   &VBS_geneta  );
		tree1->SetBranchAddress("GenPart_phi",   &VBS_genphi  );
		tree1->SetBranchAddress("GenPart_mass",   &VBS_genm  );

		tree1->SetBranchAddress("Jet_pt",   &VBS_jetpt  );
		tree1->SetBranchAddress("Jet_eta",   &VBS_jeteta  );
		tree1->SetBranchAddress("Jet_phi",   &VBS_jetphi  );
		tree1->SetBranchAddress("Jet_mass",   &VBS_jetm  );

		tree1->SetBranchAddress("GenJet_pt",   &VBS_gjetpt  );
		tree1->SetBranchAddress("GenJet_eta",   &VBS_gjeteta  );
		tree1->SetBranchAddress("GenJet_phi",   &VBS_gjetphi  );
		tree1->SetBranchAddress("GenJet_mass",   &VBS_gjetm  );


		tree1->SetBranchAddress("Jet_jetId",   &VBS_jetid  );
		tree1->SetBranchAddress("GenPart_pdgId",   &VBS_genid  );
		tree1->SetBranchAddress("GenPart_status",   &VBS_genstatus  );
		tree1->SetBranchAddress("GenPart_genPartIdxMother",   &VBS_genmomid  );
		tree1->SetBranchAddress("Photon_pdgId",   &VBS_phoid  );

		TLorentzVector VBS_photon1,VBS_photon2, VBS_jet1, VBS_jet2;
		TLorentzVector VBS_genpart1, VBS_genpart2 ;

		TLorentzVector VBS_photon, VBS_genpart ;
		TLorentzVector VBS_jet, VBS_genjet ;

		TLorentzVector VBS_true_photon1, VBS_true_photon2;
		TLorentzVector VBS_true_jet1, VBS_true_jet2;

		TLorentzVector VBS_gjet1, VBS_gjet2;
		TLorentzVector VBS_gphoton1,VBS_gphoton2;

		TLorentzVector VBS_q1_gen_sph , VBS_q2_gen_sph ;

		///////////////////////////////////////////loop INTER/////////////////////////////////////////////////////////////////////////

		for( int evt=0 ; evt < tree->GetEntries() ; evt++ ) {


			TLorentzVector tmp;
			tmp.SetPxPyPzE(0,0,0,0) ;
			INTER_q1_gen_sph.SetPtEtaPhiM(0,0,0,0) , INTER_q2_gen_sph.SetPtEtaPhiM(0,0,0,0);
			TLorentzVector  INTER_p1_gen_sph , INTER_p2_gen_sph ;

			int pho_size =0;
			int jet_size =0;
			int genpart_size =0;
			int genjet_size =0;

			tree->GetEntry(evt);

			//-----------calcualte size of array photon-------------------//
			for(int a1=0 ; a1<sizeof(INTER_phopt) / sizeof(INTER_phopt[0]) ; a1++) {
				//cout << " a1= " << a1 << " || INTER_phopt" << INTER_phopt[a1]<<endl;
				if(INTER_phopt[a1]!=0)  pho_size++ ;

			}


			//---------------calcualte size of array genpart---------------//
			for(int a2=0 ; a2<sizeof(INTER_genpt) / sizeof(INTER_genpt[0]) ; a2++) {
				//cout << " a2= " << a2 << " || INTER_genpt" << INTER_genpt[a2]<<endl;
				if(INTER_genpt[a2]!=0)  genpart_size++ ;
			}


			//---------------calcualte size of array jet---------------//
			for(int a3=0 ; a3<sizeof(INTER_jetpt) / sizeof(INTER_jetpt[0]) ; a3++) {
				//cout << " a3= " << a3 << " || INTER_jetpt" << INTER_jetpt[a3]<<endl;
				if(INTER_jetpt[a3]!=0)  jet_size++ ;
			}

			//---------------calcualte size of array genjet---------------//
			for(int a4=0 ; a4<sizeof(INTER_gjetpt) / sizeof(INTER_gjetpt[0]) ; a4++) {
				//cout << " a4= " << a4 << " || INTER_gjetpt" << INTER_gjetpt[a4]<<endl;
				if(INTER_gjetpt[a4]!=0)  genjet_size++ ;
			}
			//cout << " --------- "<<genjet_size <<endl;

/////////////////////////// calculate gen-reco photon deltaR////////////////////////////////////////////
			TLorentzVector INTER_photon_gen_sph ;
			TLorentzVector INTER_photon_reco_sph ;
			for(int l2=0 ; l2<pho_size ; l2++) {
				INTER_photon_reco_sph.SetPtEtaPhiM(INTER_phopt[l2],INTER_phoeta[l2],INTER_phophi[l2],INTER_phom[l2]);
				double tmpdeltaR = 1000;
				for(int l3=0 ; l3<genpart_size ; l3++) {

					int status0 = INTER_genstatus[l3];
					int PID0 = INTER_genid[l3] ;
					int gmother =INTER_genmomid[l3];

					if(status0 >0 	 ) {

						if(abs(PID0) ==22 && gmother<=1 ) {
							INTER_photon_gen_sph.SetPtEtaPhiM(INTER_genpt[l3],INTER_geneta[l3],INTER_genphi[l3],INTER_genm[l3]);
							INTER_deltaR_p1[l3] = INTER_photon_gen_sph.DeltaR(INTER_photon_reco_sph);


							if(INTER_deltaR_p1[l3] < tmpdeltaR  ) {
								tmpdeltaR = INTER_deltaR_p1[l3];

							}
						}

					}
				}

				//cout << tmpdeltaR <<endl;
				INTER_deltaR_p_genreco->Fill(tmpdeltaR);

			}


/////////////////////////// determine gen-reco photon deltaR////////////////////////////////////////////
TLorentzVector INTER_genpho_1st;
			TLorentzVector INTER_genpho_2nd;
			int acc_3 = 0 ;


			for(int l2=0 ; l2<pho_size ; l2++) {
				
				INTER_photon_reco_sph.SetPtEtaPhiM(INTER_phopt[l2],INTER_phoeta[l2],INTER_phophi[l2],INTER_phom[l2]);
				double tmpdeltaR = 1000;
				for(int l3=0 ; l3<genpart_size ; l3++) {

					int status0 = INTER_genstatus[l3];
					int PID0 = INTER_genid[l3] ;
					int gmother =INTER_genmomid[l3];

					if(status0 >0 	 ) {

						if(abs(PID0) ==22 && gmother<=1 ) {
							INTER_photon_gen_sph.SetPtEtaPhiM(INTER_genpt[l3],INTER_geneta[l3],INTER_genphi[l3],INTER_genm[l3]);
							INTER_deltaR_p1[l3] = INTER_photon_gen_sph.DeltaR(INTER_photon_reco_sph);


							if(INTER_deltaR_p1[l3] <= 0.2  ) {
								acc_3 = acc_3 +1;
								INTER_genphotrue_pt[l3]= INTER_genpt[l3];
								INTER_genphotrue_pt_plot->Fill(INTER_genphotrue_pt[l3]);
								INTER_genphotrue_phi[l3]= INTER_genphi[l3];
								INTER_genphotrue_phi_plot->Fill(INTER_genphotrue_phi[l3]);
								INTER_genphotrue_eta[l3]= INTER_geneta[l3];
								INTER_genphotrue_eta_plot->Fill(INTER_genphotrue_eta[l3]);
								INTER_genphotrue_m[l3]= INTER_genm[l3];
								INTER_genphotrue_m_plot->Fill(INTER_genphotrue_m[l3]);
								if(acc_3 == 1)  INTER_genpho_1st.SetPtEtaPhiM(INTER_genpt[l3],INTER_geneta[l3],INTER_genphi[l3],INTER_genm[l3]);
						if(acc_3 == 2)   INTER_genpho_2nd.SetPtEtaPhiM(INTER_genpt[l3],INTER_geneta[l3],INTER_genphi[l3],INTER_genm[l3]);
								break;
							}
						}

					}
				}

			}

///////////////////////////true genpho////////////////////////////////////////////

double INTER_deltaR_p1_p2,INTER_deltaeta_p1_p2,INTER_deltaphi_p1_p2,INTER_deltapt_p1_p2,INTER_invmass_p1_p2;
			if(acc_3 >=2) {
				double tmp_eta1,tmp_eta2,tmp_phi1,tmp_phi2,tmp_pt1,tmp_pt2;
				double tmp_deltaR1,tmp_deltaR2, tmp_invmass1,tmp_invmass2;

				tmp_eta1 =INTER_genpho_1st.Eta();
				tmp_eta2 =INTER_genpho_2nd.Eta();
				tmp_pt1 =INTER_genpho_1st.Pt();
				tmp_pt2 =INTER_genpho_2nd.Pt();

				INTER_deltaeta_p1_p2 = tmp_eta1 - tmp_eta2;
				INTER_deltaphi_p1_p2 = INTER_genpho_1st.DeltaPhi(INTER_genpho_2nd);
				INTER_deltaR_p1_p2 = INTER_genpho_1st.DeltaR(INTER_genpho_2nd);
				INTER_deltapt_p1_p2 = tmp_pt1 - tmp_pt2;
				INTER_invmass_p1_p2 = (INTER_genpho_1st+INTER_genpho_2nd).M();

				INTER_genphotrue_eta_p1p2_plot->Fill(INTER_deltaeta_p1_p2);
				INTER_genphotrue_phi_p1p2_plot->Fill(INTER_deltaphi_p1_p2);
				INTER_genphotrue_pt_p1p2_plot->Fill(INTER_deltapt_p1_p2);
				INTER_genphotrue_deltaR_p1p2_plot->Fill(INTER_deltaR_p1_p2);
				INTER_genphotrue_mpp_plot->Fill(INTER_invmass_p1_p2);





				//cout << "eta ||" << " phi || " << " R || " << " pt || " << " invmass " <<endl;
				//cout << INTER_deltaeta_j1_j2  << " || " << INTER_deltaphi_j1_j2  <<" || " << INTER_deltaR_j1_j2  <<" || " << INTER_deltapt_j1_j2  <<" ||  " << INTER_invmass_j1_j2  <<endl;
			}















///////////////////////////calculate gen-reco jet deltaR////genjet_size////////////////////////////////////////

			TLorentzVector INTER_jet_gen_sph ;
			TLorentzVector INTER_jet_reco_sph ;
			cout << genjet_size <<endl;
			for(int l2=0 ; l2<genjet_size ; l2++) {
				//cout << "loop in reco" << "-------------------"<<endl;
				INTER_jet_gen_sph.SetPtEtaPhiM(INTER_gjetpt[l2],INTER_gjeteta[l2],INTER_gjetphi[l2],INTER_gjetm[l2]);
				double tmpdeltaR = 100000;

				for(int l3=0 ; l3<jet_size  ; l3++) {
					//cout << "loop in genjet" << "-------------------"<<endl;
					INTER_jet_reco_sph.SetPtEtaPhiM(INTER_jetpt[l3],INTER_jeteta[l3],INTER_jetphi[l3],INTER_jetm[l3]);


					INTER_deltaR_j1[l3]  = INTER_jet_gen_sph.DeltaR(INTER_jet_reco_sph);
					if(INTER_deltaR_j1[l3] < tmpdeltaR  ) {

						tmpdeltaR = INTER_deltaR_j1[l3];
						//cout <<" tmp in genloop " << tmpdeltaR <<" || min_index = " << l3  <<endl;
					}


				}
				//cout << "tmp" << tmpdeltaR <<endl;
				INTER_deltaR_j_genreco->Fill(tmpdeltaR);

			}

///////////////////////////determine gen-reco jet deltaR////////////////////////////////////////////
			TLorentzVector INTER_genjet_1st;
			TLorentzVector INTER_genjet_2nd;
			int acc_1 = 0 ;


			for(int l2=0 ; l2<jet_size ; l2++) {
				INTER_jet_reco_sph.SetPtEtaPhiM(INTER_jetpt[l2],INTER_jeteta[l2],INTER_jetphi[l2],INTER_jetm[l2]);
				double tmpdeltaR = 1000;
				for(int l3=0 ; l3<genjet_size ; l3++) {


					INTER_jet_gen_sph.SetPtEtaPhiM(INTER_gjetpt[l3],INTER_gjeteta[l3],INTER_gjetphi[l3],INTER_gjetm[l3]);
					INTER_deltaR_j1[l3] = INTER_jet_gen_sph.DeltaR(INTER_jet_reco_sph);


					if(INTER_deltaR_j1[l3] < 0.2 && INTER_gjetpt[l3] >20 ) {
						acc_1 = acc_1 +1;
						INTER_genjettrue_pt[l3]= INTER_gjetpt[l3];
						INTER_genjettrue_pt_plot->Fill(INTER_genjettrue_pt[l3]);
						INTER_genjettrue_phi[l3]= INTER_gjetphi[l3];
						INTER_genjettrue_phi_plot->Fill(INTER_genjettrue_phi[l3]);
						INTER_genjettrue_eta[l3]= INTER_gjeteta[l3];
						INTER_genjettrue_eta_plot->Fill(INTER_genjettrue_eta[l3]);
						INTER_genjettrue_m[l3]= INTER_gjetm[l3];
						INTER_genjettrue_m_plot->Fill(INTER_genjettrue_m[l3]);
						if(acc_1 == 1)  INTER_genjet_1st.SetPtEtaPhiM(INTER_gjetpt[l3],INTER_gjeteta[l3],INTER_gjetphi[l3],INTER_gjetm[l3]);
						if(acc_1 == 2)   INTER_genjet_2nd.SetPtEtaPhiM(INTER_gjetpt[l3],INTER_gjeteta[l3],INTER_gjetphi[l3],INTER_gjetm[l3]);

						break;
					}



				}



			}


///////////////////////////true genjet////////////////////////////////////////////
			double INTER_deltaR_j1_j2,INTER_deltaeta_j1_j2,INTER_deltaphi_j1_j2,INTER_deltapt_j1_j2,INTER_invmass_j1_j2;
			if(acc_1 >=2) {
				double tmp_eta1,tmp_eta2,tmp_phi1,tmp_phi2,tmp_pt1,tmp_pt2;
				double tmp_deltaR1,tmp_deltaR2, tmp_invmass1,tmp_invmass2;

				tmp_eta1 =INTER_genjet_1st.Eta();
				tmp_eta2 =INTER_genjet_2nd.Eta();
				tmp_pt1 =INTER_genjet_1st.Pt();
				tmp_pt2 =INTER_genjet_2nd.Pt();

				INTER_deltaeta_j1_j2 = tmp_eta1 - tmp_eta2;
				INTER_deltaphi_j1_j2 = INTER_genjet_1st.DeltaPhi(INTER_genjet_2nd);
				INTER_deltaR_j1_j2 = INTER_genjet_1st.DeltaR(INTER_genjet_2nd);
				INTER_deltapt_j1_j2 = tmp_pt1 - tmp_pt2;
				INTER_invmass_j1_j2 = (INTER_genjet_1st+INTER_genjet_2nd).M();

				INTER_genjettrue_eta_j1j2_plot->Fill(INTER_deltaeta_j1_j2);
				INTER_genjettrue_phi_j1j2_plot->Fill(INTER_deltaphi_j1_j2);
				INTER_genjettrue_pt_j1j2_plot->Fill(INTER_deltapt_j1_j2);
				INTER_genjettrue_deltaR_j1j2_plot->Fill(INTER_deltaR_j1_j2);
				INTER_genjettrue_mjj_plot->Fill(INTER_invmass_j1_j2);





				//cout << "eta ||" << " phi || " << " R || " << " pt || " << " invmass " <<endl;
				//cout << INTER_deltaeta_j1_j2  << " || " << INTER_deltaphi_j1_j2  <<" || " << INTER_deltaR_j1_j2  <<" || " << INTER_deltapt_j1_j2  <<" ||  " << INTER_invmass_j1_j2  <<endl;
			}





///////////////////////////initial////////////////////////////////////////////
			INTER_genjet_1st.SetPtEtaPhiM(0,0,0,0);
			INTER_genjet_2nd.SetPtEtaPhiM(0,0,0,0);
			INTER_deltaR_j1_j2=0 ,INTER_deltaeta_j1_j2=0 ,INTER_deltaphi_j1_j2=0 ,INTER_deltapt_j1_j2=0 ,INTER_invmass_j1_j2=0 ;
			memset(INTER_gjetpt, 0, sizeof(INTER_gjetpt));
			memset(INTER_jetpt, 0, sizeof(INTER_jetpt));
			memset(INTER_deltaR_j1, 0, sizeof(INTER_deltaR_j1));

			memset(INTER_genjettrue_pt, 0, sizeof(INTER_genjettrue_pt));
			memset(INTER_genjettrue_phi, 0, sizeof(INTER_genjettrue_phi));
			memset(INTER_genjettrue_eta, 0, sizeof(INTER_genjettrue_eta));
			memset(INTER_genjettrue_m, 0, sizeof(INTER_genjettrue_m));

			memset(INTER_genphotrue_pt, 0, sizeof(INTER_genphotrue_pt));
			memset(INTER_genphotrue_phi, 0, sizeof(INTER_genphotrue_phi));
			memset(INTER_genphotrue_eta, 0, sizeof(INTER_genphotrue_eta));
			memset(INTER_genphotrue_m, 0, sizeof(INTER_genphotrue_m));

			memset(INTER_gjetpt, 0, sizeof(INTER_gjetpt));
			memset(INTER_gjetphi, 0, sizeof(INTER_gjetphi));
			memset(INTER_gjeteta, 0, sizeof(INTER_gjeteta));
			memset(INTER_gjetm, 0, sizeof(INTER_gjetm));
			
			memset(INTER_phopt, 0, sizeof(INTER_phopt));
			memset(INTER_phophi, 0, sizeof(INTER_phophi));
			memset(INTER_phoeta, 0, sizeof(INTER_phoeta));
			memset(INTER_phom, 0, sizeof(INTER_phom));
			
			memset(INTER_genpt, 0, sizeof(INTER_genpt));
			memset(INTER_genphi, 0, sizeof(INTER_genphi));
			memset(INTER_geneta, 0, sizeof(INTER_geneta));
			memset(INTER_genm, 0, sizeof(INTER_genm));
			



		}

		///////////////////////////////////////////loop VBS/////////////////////////////////////////////////////////////////////////

		for( int evt=0 ; evt < tree1->GetEntries() ; evt++ ) {


			TLorentzVector tmp;
			tmp.SetPxPyPzE(0,0,0,0) ;
			VBS_q1_gen_sph.SetPtEtaPhiM(0,0,0,0) , VBS_q2_gen_sph.SetPtEtaPhiM(0,0,0,0);
			TLorentzVector  VBS_p1_gen_sph , VBS_p2_gen_sph ;

			int pho_size =0;
			int jet_size =0;
			int genpart_size =0;
			int genjet_size =0;

			tree1->GetEntry(evt);
			//cout << evt <<endl;

			//-----------calcualte size of array photon-------------------//
			for(int a1=0 ; a1<sizeof(VBS_phopt) / sizeof(VBS_phopt[0]) ; a1++) {
				if(VBS_phopt[a1]!=0)  pho_size++ ;

			}


			//---------------calcualte size of array genpart---------------//
			for(int a2=0 ; a2<sizeof(VBS_genpt) / sizeof(VBS_genpt[0]) ; a2++) {
				if(VBS_genpt[a2]!=0)  genpart_size++ ;
			}


			//---------------calcualte size of array jet---------------//
			for(int a3=0 ; a3<sizeof(VBS_jetpt) / sizeof(VBS_jetpt[0]) ; a3++) {
				if(VBS_jetpt[a3]!=0)  jet_size++ ;
			}

			//---------------calcualte size of array genjet---------------//
			for(int a4=0 ; a4<sizeof(VBS_gjetpt) / sizeof(VBS_gjetpt[0]) ; a4++) {
				if(VBS_gjetpt[a4]!=0)  genjet_size++ ;
			}


			/////////////////////////// calculate gen-reco photon deltaR////////////////////////////////////////////
			TLorentzVector VBS_photon_gen_sph ;
			TLorentzVector VBS_photon_reco_sph ;
			for(int l2=0 ; l2<pho_size ; l2++) {
				VBS_photon_reco_sph.SetPtEtaPhiM(VBS_phopt[l2],VBS_phoeta[l2],VBS_phophi[l2],VBS_phom[l2]);
				double tmpdeltaR = 1000;
				for(int l3=0 ; l3<genpart_size ; l3++) {

					int status0 = VBS_genstatus[l3];
					int PID0 = VBS_genid[l3] ;
					int gmother =VBS_genmomid[l3];

					if(status0 >0 	 ) {

						if(abs(PID0) ==22 && gmother<=1 ) {
							VBS_photon_gen_sph.SetPtEtaPhiM(VBS_genpt[l3],VBS_geneta[l3],VBS_genphi[l3],VBS_genm[l3]);
							VBS_deltaR_p1[l3] = VBS_photon_gen_sph.DeltaR(VBS_photon_reco_sph);


							if(VBS_deltaR_p1[l3] < tmpdeltaR  ) {
								tmpdeltaR = VBS_deltaR_p1[l3];

							}
						}

					}
				}

				VBS_deltaR_p_genreco->Fill(tmpdeltaR);

			}


/////////////////////////// determine gen-reco photon deltaR////////////////////////////////////////////

			for(int l2=0 ; l2<pho_size ; l2++) {
				VBS_photon_reco_sph.SetPtEtaPhiM(VBS_phopt[l2],VBS_phoeta[l2],VBS_phophi[l2],VBS_phom[l2]);
				double tmpdeltaR = 1000;
				for(int l3=0 ; l3<genpart_size ; l3++) {

					int status0 = VBS_genstatus[l3];
					int PID0 = VBS_genid[l3] ;
					int gmother =VBS_genmomid[l3];

					if(status0 >0 	 ) {

						if(abs(PID0) ==22 && gmother<=1 ) {
							VBS_photon_gen_sph.SetPtEtaPhiM(VBS_genpt[l3],VBS_geneta[l3],VBS_genphi[l3],VBS_genm[l3]);
							VBS_deltaR_p1[l3] = VBS_photon_gen_sph.DeltaR(VBS_photon_reco_sph);


							if(VBS_deltaR_p1[l3] <= 0.2  ) {
								VBS_genphotrue_pt[l3]= VBS_genpt[l3];
								VBS_genphotrue_pt_plot->Fill(VBS_genphotrue_pt[l3]);
								VBS_genphotrue_phi[l3]= VBS_genphi[l3];
								VBS_genphotrue_phi_plot->Fill(VBS_genphotrue_phi[l3]);
								VBS_genphotrue_eta[l3]= VBS_geneta[l3];
								VBS_genphotrue_eta_plot->Fill(VBS_genphotrue_eta[l3]);
								VBS_genphotrue_m[l3]= VBS_genm[l3];
								VBS_genphotrue_m_plot->Fill(VBS_genphotrue_m[l3]);
								break;
							}
						}

					}
				}



			}

///////////////////////////calculate gen-reco jet deltaR////genjet_size////////////////////////////////////////

			TLorentzVector VBS_jet_gen_sph ;
			TLorentzVector VBS_jet_reco_sph ;
			for(int l2=0 ; l2<genjet_size ; l2++) {
				//cout << "loop in reco" << "-------------------"<<endl;
				VBS_jet_gen_sph.SetPtEtaPhiM(VBS_gjetpt[l2],VBS_gjeteta[l2],VBS_gjetphi[l2],VBS_gjetm[l2]);
				//cout << VBS_gjetpt[l2] <<endl;
				double tmpdeltaR = 1000;

				for(int l3=0 ; l3<jet_size  ; l3++) {
					//cout << "loop in genjet" << "-------------------"<<endl;
					VBS_jet_reco_sph.SetPtEtaPhiM(VBS_jetpt[l3],VBS_jeteta[l3],VBS_jetphi[l3],VBS_jetm[l3]);


					VBS_deltaR_j1[l3]  = VBS_jet_gen_sph.DeltaR(VBS_jet_reco_sph);

					if(VBS_deltaR_j1[l3] < tmpdeltaR ) {

						tmpdeltaR = VBS_deltaR_j1[l3];
						//	cout <<" tmp in genloop " << tmpdeltaR <<" || pt = " << VBS_gjetpt[l2] <<endl;
					}


				}
				if(VBS_gjetpt[l2] > 30  ) //cout << "tmp" << VBS_gjetpt[l2] <<endl;
					if(VBS_gjetpt[l2] > 30  ) VBS_deltaR_j_genreco->Fill(tmpdeltaR);

			}


///////////////////////////determine gen-reco jet deltaR////////////////////////////////////////////
			TLorentzVector VBS_genjet_1st;
			TLorentzVector VBS_genjet_2nd;
			int acc_ = 0 ;
			for(int l2=0 ; l2<jet_size ; l2++) {

				VBS_jet_reco_sph.SetPtEtaPhiM(VBS_jetpt[l2],VBS_jeteta[l2],VBS_jetphi[l2],VBS_jetm[l2]);
				double tmpdeltaR = 1000;
				for(int l3=0 ; l3<genjet_size ; l3++) {


					VBS_jet_gen_sph.SetPtEtaPhiM(VBS_gjetpt[l3],VBS_gjeteta[l3],VBS_gjetphi[l3],VBS_gjetm[l3]);
					VBS_deltaR_j1[l3] = VBS_jet_gen_sph.DeltaR(VBS_jet_reco_sph);


					if(VBS_deltaR_j1[l3] < 0.2 && VBS_gjetpt[l3] >20   ) {
						acc_ = acc_ +1;
						VBS_genjettrue_pt[l3]= VBS_gjetpt[l3];
						VBS_genjettrue_pt_plot->Fill(VBS_genjettrue_pt[l3]);
						VBS_genjettrue_phi[l3]= VBS_gjetphi[l3];
						VBS_genjettrue_phi_plot->Fill(VBS_genjettrue_phi[l3]);
						VBS_genjettrue_eta[l3]= VBS_gjeteta[l3];
						VBS_genjettrue_eta_plot->Fill(VBS_genjettrue_eta[l3]);
						VBS_genjettrue_m[l3]= VBS_gjetm[l3];
						VBS_genjettrue_m_plot->Fill(VBS_genjettrue_m[l3]);
						if(acc_ == 1)  VBS_genjet_1st.SetPtEtaPhiM(VBS_gjetpt[l3],VBS_gjeteta[l3],VBS_gjetphi[l3],VBS_gjetm[l3]);
						if(acc_ == 2)   VBS_genjet_2nd.SetPtEtaPhiM(VBS_gjetpt[l3],VBS_gjeteta[l3],VBS_gjetphi[l3],VBS_gjetm[l3]);

						break;
					}
				}
			}




///////////////////////////true genjet////////////////////////////////////////////
			double VBS_deltaR_j1_j2,VBS_deltaeta_j1_j2,VBS_deltaphi_j1_j2,VBS_deltapt_j1_j2,VBS_invmass_j1_j2;
			if(acc_ >=2) {
				double tmp_eta1,tmp_eta2,tmp_phi1,tmp_phi2,tmp_pt1,tmp_pt2;
				double tmp_deltaR1,tmp_deltaR2, tmp_invmass1,tmp_invmass2;

				tmp_eta1 =VBS_genjet_1st.Eta();
				tmp_eta2 =VBS_genjet_2nd.Eta();
				tmp_pt1 =VBS_genjet_1st.Pt();
				tmp_pt2 =VBS_genjet_2nd.Pt();

				VBS_deltaeta_j1_j2 = tmp_eta1 - tmp_eta2;
				VBS_deltaphi_j1_j2 = VBS_genjet_1st.DeltaPhi(VBS_genjet_2nd);
				VBS_deltaR_j1_j2 = VBS_genjet_1st.DeltaR(VBS_genjet_2nd);
				VBS_deltapt_j1_j2 = tmp_pt1 - tmp_pt2;
				VBS_invmass_j1_j2 = (VBS_genjet_1st+VBS_genjet_2nd).M();

				VBS_genjettrue_eta_j1j2_plot->Fill(VBS_deltaeta_j1_j2);
				VBS_genjettrue_phi_j1j2_plot->Fill(VBS_deltaphi_j1_j2);
				VBS_genjettrue_pt_j1j2_plot->Fill(VBS_deltapt_j1_j2);
				VBS_genjettrue_deltaR_j1j2_plot->Fill(VBS_deltaR_j1_j2);
				VBS_genjettrue_mjj_plot->Fill(VBS_invmass_j1_j2);


				//cout << "eta ||" << " phi || " << " R || " << " pt || " << " invmass " <<endl;
				//cout << VBS_deltaeta_j1_j2  << " || " << VBS_deltaphi_j1_j2  <<" || " << VBS_deltaR_j1_j2  <<" || " << VBS_deltapt_j1_j2  <<" ||  " << VBS_invmass_j1_j2  <<endl;
			}





///////////////////////////initial////////////////////////////////////////////
			VBS_genjet_1st.SetPtEtaPhiM(0,0,0,0);
			VBS_genjet_2nd.SetPtEtaPhiM(0,0,0,0);
			VBS_deltaR_j1_j2=0 ,VBS_deltaeta_j1_j2=0 ,VBS_deltaphi_j1_j2=0 ,VBS_deltapt_j1_j2=0 ,VBS_invmass_j1_j2=0 ;
			memset(VBS_gjetpt, 0, sizeof(VBS_gjetpt));
			memset(VBS_jetpt, 0, sizeof(VBS_jetpt));
			memset(VBS_deltaR_j1, 0, sizeof(VBS_deltaR_j1));

			memset(VBS_genjettrue_pt, 0, sizeof(VBS_genjettrue_pt));
			memset(VBS_genjettrue_phi, 0, sizeof(VBS_genjettrue_phi));
			memset(VBS_genjettrue_eta, 0, sizeof(VBS_genjettrue_eta));
			memset(VBS_genjettrue_m, 0, sizeof(VBS_genjettrue_m));

			memset(VBS_genphotrue_pt, 0, sizeof(VBS_genphotrue_pt));
			memset(VBS_genphotrue_phi, 0, sizeof(VBS_genphotrue_phi));
			memset(VBS_genphotrue_eta, 0, sizeof(VBS_genphotrue_eta));
			memset(VBS_genphotrue_m, 0, sizeof(VBS_genphotrue_m));

			memset(VBS_gjetpt, 0, sizeof(VBS_gjetpt));
			memset(VBS_gjetphi, 0, sizeof(VBS_gjetphi));
			memset(VBS_gjeteta, 0, sizeof(VBS_gjeteta));
			memset(VBS_gjetm, 0, sizeof(VBS_gjetm));









		}

		///////////////////////////end of VBS loop////

	}

	//	f->Close();
	//	f1->Close();
	//	delete f;
	//	delete f1;

	///////////////////////////////


	TCanvas *cl27 = new TCanvas("cl27","recopho_gen[part]");

	INTER_deltaR_p_genreco->SetStats(0);
	INTER_deltaR_p_genreco->SetLineColor(4);
	INTER_deltaR_p_genreco->SetLineWidth(3);

	INTER_deltaR_p_genreco->Draw("hist");

	VBS_deltaR_p_genreco->SetStats(0);
	VBS_deltaR_p_genreco->SetLineColor(3);
	VBS_deltaR_p_genreco->SetLineWidth(3);

	VBS_deltaR_p_genreco->Draw("hist,same");
	cl27->SetTitle("hghg");
	cl27->SaveAs("!!2gen_2_deltaRgen_pho.pdf");


///////////////////////////////////////////////////////

	TCanvas *cl28 = new TCanvas("cl28","recojjet_genjet");

	INTER_deltaR_j_genreco->SetStats(0);
	INTER_deltaR_j_genreco->SetLineColor(4);
	INTER_deltaR_j_genreco->SetLineWidth(3);
	INTER_deltaR_j_genreco->Scale(1/INTER_deltaR_j_genreco->Integral());
	INTER_deltaR_j_genreco->GetYaxis()->SetRangeUser(0, 0.8);
	INTER_deltaR_j_genreco->Draw("hist");

	VBS_deltaR_j_genreco->SetStats(0);
	VBS_deltaR_j_genreco->SetLineColor(3);
	VBS_deltaR_j_genreco->SetLineWidth(3);
	VBS_deltaR_j_genreco->Scale(1/VBS_deltaR_j_genreco->Integral());
	VBS_deltaR_j_genreco->Draw("hist,same");

	cl28->SaveAs("!!2gen_2_deltaRgenjet_recojet.pdf");
///////////////////////////////////Gentrue_parameter//////////////////////////////////////////////////////////////////////

//------------------pho
/////pt////

	TCanvas *cl29 = new TCanvas("cl29","pt_t_genpho");

	INTER_genphotrue_pt_plot->SetStats(0);
	INTER_genphotrue_pt_plot->SetLineColor(4);
	INTER_genphotrue_pt_plot->SetLineWidth(3);

	INTER_genphotrue_pt_plot->Draw("hist");

	VBS_genphotrue_pt_plot->SetStats(0);
	VBS_genphotrue_pt_plot->SetLineColor(3);
	VBS_genphotrue_pt_plot->SetLineWidth(3);

	VBS_genphotrue_pt_plot->Draw("hist,same");

	cl29->SaveAs("!!3truegenpho_pt.pdf");



///////eta/////


	TCanvas *cl30 = new TCanvas("cl30","eta_t_genpho");

	INTER_genphotrue_eta_plot->SetStats(0);
	INTER_genphotrue_eta_plot->SetLineColor(4);
	INTER_genphotrue_eta_plot->SetLineWidth(3);

	INTER_genphotrue_eta_plot->Draw("hist");

	VBS_genphotrue_eta_plot->SetStats(0);
	VBS_genphotrue_eta_plot->SetLineColor(3);
	VBS_genphotrue_eta_plot->SetLineWidth(3);

	VBS_genphotrue_eta_plot->Draw("hist,same");

	cl30->SaveAs("!!3truegenpho_eta.pdf");

/////phi////


	TCanvas *cl31= new TCanvas("cl31","phi_t_genpho");

	INTER_genphotrue_phi_plot->SetStats(0);
	INTER_genphotrue_phi_plot->SetLineColor(4);
	INTER_genphotrue_phi_plot->SetLineWidth(3);

	INTER_genphotrue_phi_plot->Draw("hist");

	VBS_genphotrue_phi_plot->SetStats(0);
	VBS_genphotrue_phi_plot->SetLineColor(3);
	VBS_genphotrue_phi_plot->SetLineWidth(3);

	VBS_genphotrue_phi_plot->Draw("hist,same");

	cl31->SaveAs("!!3truegenpho_phi.pdf");


/////m/////

	TCanvas *cl32 = new TCanvas("cl32","m_t_genpho");

	INTER_genphotrue_m_plot->SetStats(0);
	INTER_genphotrue_m_plot->SetLineColor(4);
	INTER_genphotrue_m_plot->SetLineWidth(3);

	INTER_genphotrue_m_plot->Draw("hist");

	VBS_genphotrue_m_plot->SetStats(0);
	VBS_genphotrue_m_plot->SetLineColor(3);
	VBS_genphotrue_m_plot->SetLineWidth(3);

	VBS_genphotrue_m_plot->Draw("hist,same");

	cl32->SaveAs("!!3truegenpho_m.pdf");

///-------------------------------truejet

	/////pt////

	TCanvas *cl33 = new TCanvas("cl33","pt_t_genjet");

	INTER_genjettrue_pt_plot->SetStats(0);
	INTER_genjettrue_pt_plot->SetLineColor(4);
	INTER_genjettrue_pt_plot->SetLineWidth(3);

	INTER_genjettrue_pt_plot->Draw("hist");

	VBS_genjettrue_pt_plot->SetStats(0);
	VBS_genjettrue_pt_plot->SetLineColor(3);
	VBS_genjettrue_pt_plot->SetLineWidth(3);

	VBS_genjettrue_pt_plot->Draw("hist,same");

	cl33->SaveAs("!!3truegenjet_pt.pdf");



///////eta/////


	TCanvas *cl34 = new TCanvas("cl34","eta_t_genjet");

	INTER_genjettrue_eta_plot->SetStats(0);
	INTER_genjettrue_eta_plot->SetLineColor(4);
	INTER_genjettrue_eta_plot->SetLineWidth(3);

	INTER_genjettrue_eta_plot->Draw("hist");

	VBS_genjettrue_eta_plot->SetStats(0);
	VBS_genjettrue_eta_plot->SetLineColor(3);
	VBS_genjettrue_eta_plot->SetLineWidth(3);

	VBS_genjettrue_eta_plot->Draw("hist,same");

	cl34->SaveAs("!!3truegenjet_eta.pdf");

/////phi////


	TCanvas *cl35= new TCanvas("cl35","phi_t_genjet");

	INTER_genjettrue_phi_plot->SetStats(0);
	INTER_genjettrue_phi_plot->SetLineColor(4);
	INTER_genjettrue_phi_plot->SetLineWidth(3);

	INTER_genjettrue_phi_plot->Draw("hist");

	VBS_genjettrue_phi_plot->SetStats(0);
	VBS_genjettrue_phi_plot->SetLineColor(3);
	VBS_genjettrue_phi_plot->SetLineWidth(3);

	VBS_genjettrue_phi_plot->Draw("hist,same");

	cl35->SaveAs("!!3truegenjet_phi.pdf");


/////m/////

	TCanvas *cl36 = new TCanvas("cl36","m_t_genjet");

	INTER_genjettrue_m_plot->SetStats(0);
	INTER_genjettrue_m_plot->SetLineColor(4);
	INTER_genjettrue_m_plot->SetLineWidth(3);

	INTER_genjettrue_m_plot->Draw("hist");

	VBS_genjettrue_m_plot->SetStats(0);
	VBS_genjettrue_m_plot->SetLineColor(3);
	VBS_genjettrue_m_plot->SetLineWidth(3);

	VBS_genjettrue_m_plot->Draw("hist,same");

	cl36->SaveAs("!!3truegenjet_m.pdf");


///-------------------------------truejet delta

//////deltapt/////
	TCanvas *cl37 = new TCanvas("cl37","genjettrue_pt_j1j2_plot");

	INTER_genjettrue_pt_j1j2_plot->SetStats(0);
	INTER_genjettrue_pt_j1j2_plot->SetLineColor(4);
	INTER_genjettrue_pt_j1j2_plot->SetLineWidth(3);
	int b_max_etaj1j2 = INTER_genjettrue_pt_j1j2_plot->GetMaximumBin();
	double y_max_etaj1j2 = INTER_genjettrue_pt_j1j2_plot->GetBinContent(b_max_etaj1j2);

	INTER_genjettrue_pt_j1j2_plot->GetYaxis()->SetRangeUser(0, y_max_etaj1j2*1.2);
	INTER_genjettrue_pt_j1j2_plot->Draw("hist");

	VBS_genjettrue_pt_j1j2_plot->SetStats(0);
	VBS_genjettrue_pt_j1j2_plot->SetLineColor(3);
	VBS_genjettrue_pt_j1j2_plot->SetLineWidth(3);

	VBS_genjettrue_pt_j1j2_plot->Draw("hist,same");

	cl37->SaveAs("!!4genjettrue_pt_j1j2.pdf");

/////deltaeta/////
	TCanvas *cl38 = new TCanvas("cl38","m_t_genjet");

	INTER_genjettrue_eta_j1j2_plot->SetStats(0);
	INTER_genjettrue_eta_j1j2_plot->SetLineColor(4);
	INTER_genjettrue_eta_j1j2_plot->SetLineWidth(3);
	b_max_etaj1j2 = INTER_genjettrue_eta_j1j2_plot->GetMaximumBin();
	y_max_etaj1j2 = INTER_genjettrue_eta_j1j2_plot->GetBinContent(b_max_etaj1j2);

	INTER_genjettrue_eta_j1j2_plot->GetYaxis()->SetRangeUser(0, y_max_etaj1j2*1.2);
	INTER_genjettrue_eta_j1j2_plot->Draw("hist");

	VBS_genjettrue_eta_j1j2_plot->SetStats(0);
	VBS_genjettrue_eta_j1j2_plot->SetLineColor(3);
	VBS_genjettrue_eta_j1j2_plot->SetLineWidth(3);

	VBS_genjettrue_eta_j1j2_plot->Draw("hist,same");

	cl38->SaveAs("!!4genjettrue_eta_j1j2.pdf");

/////deltaphi//////
	TCanvas *cl39 = new TCanvas("cl39","genjettrue_phi_j1j2");

	INTER_genjettrue_phi_j1j2_plot->SetStats(0);
	INTER_genjettrue_phi_j1j2_plot->SetLineColor(4);
	INTER_genjettrue_phi_j1j2_plot->SetLineWidth(3);
	b_max_etaj1j2 = INTER_genjettrue_phi_j1j2_plot->GetMaximumBin();
	y_max_etaj1j2 = INTER_genjettrue_phi_j1j2_plot->GetBinContent(b_max_etaj1j2);

	INTER_genjettrue_phi_j1j2_plot->GetYaxis()->SetRangeUser(0, y_max_etaj1j2*1.2);
	INTER_genjettrue_phi_j1j2_plot->Draw("hist");

	VBS_genjettrue_phi_j1j2_plot->SetStats(0);
	VBS_genjettrue_phi_j1j2_plot->SetLineColor(3);
	VBS_genjettrue_phi_j1j2_plot->SetLineWidth(3);

	VBS_genjettrue_phi_j1j2_plot->Draw("hist,same");

	cl39->SaveAs("!!4genjettrue_phi_j1j2.pdf");

/////deltaR//////
	TCanvas *cl40 = new TCanvas("cl40","genjettrue_deltaR_j1j2");

	INTER_genjettrue_deltaR_j1j2_plot->SetStats(0);
	INTER_genjettrue_deltaR_j1j2_plot->SetLineColor(4);
	INTER_genjettrue_deltaR_j1j2_plot->SetLineWidth(3);
	b_max_etaj1j2 = INTER_genjettrue_deltaR_j1j2_plot->GetMaximumBin();
	y_max_etaj1j2 = INTER_genjettrue_deltaR_j1j2_plot->GetBinContent(b_max_etaj1j2);

	INTER_genjettrue_deltaR_j1j2_plot->GetYaxis()->SetRangeUser(0 , y_max_etaj1j2*1.2);
	INTER_genjettrue_deltaR_j1j2_plot->Draw("hist");

	VBS_genjettrue_deltaR_j1j2_plot->SetStats(0);
	VBS_genjettrue_deltaR_j1j2_plot->SetLineColor(3);
	VBS_genjettrue_deltaR_j1j2_plot->SetLineWidth(3);

	VBS_genjettrue_deltaR_j1j2_plot->Draw("hist,same");

	cl40->SaveAs("!!4genjettrue_deltaR_j1j2.pdf");


///////invmass/////
	TCanvas *cl41 = new TCanvas("cl41","genjettrue_mjj_j1j2");

	INTER_genjettrue_mjj_plot->SetStats(0);
	INTER_genjettrue_mjj_plot->SetLineColor(4);
	INTER_genjettrue_mjj_plot->SetLineWidth(3);
	b_max_etaj1j2 = INTER_genjettrue_mjj_plot->GetMaximumBin();
	y_max_etaj1j2 = INTER_genjettrue_mjj_plot->GetBinContent(b_max_etaj1j2);

	INTER_genjettrue_mjj_plot->GetYaxis()->SetRangeUser(0, y_max_etaj1j2*1.2);
	INTER_genjettrue_mjj_plot->Draw("hist");

	VBS_genjettrue_mjj_plot->SetStats(0);
	VBS_genjettrue_mjj_plot->SetLineColor(3);
	VBS_genjettrue_mjj_plot->SetLineWidth(3);

	VBS_genjettrue_mjj_plot->Draw("hist,same");

	cl41->SaveAs("!!4genjettrue_mjj_j1j2.pdf");















}


























