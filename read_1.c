#include<cstdio>
#include<cstdlib>
#include<iostream>
#include <vector>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include <cassert>
#include <future>
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
#include "TF1.h"
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
#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <string>
#include <cctype>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <pthread.h>
#include <thread>
#include "TStopwatch.h"


#include <time.h>


/*
[nanoAOD]       ------>  [calculate BDT]                    --------->    [  ]
write           ------>    tmva --> compare(root to text)   ---------> read(text to root)  --> read_1( read root [data and 3MC] )

data path "/data1/ggNtuples/V10_06_30_03/job_EGamma_Run2018A(B,C,D)_UL"

-----------------------------
HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 -> 14
HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95 -> 15

HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 -> 16
HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55 -> 17

HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto -> 38
HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55
-----------------------------

						        2015  2016   2017   2018   2015¡V2018  2016¡V2018
Recorded luminosity [fb-1]      2.3 | 35.9 | 41.5 | 59.7   | 139       | 137
Total uncertainty [%]           2.3 | 2.5  | 2.3  | 2.5    | 1.8       | 1.8
Calibration uncertainty [%]     1.8 | 1.5  | 1.6  | 2.1
Normalization uncertainty [%]   1.5 | 2.0  | 1.7  | 1.3
-----------------------------
Global selection


*/
std::mutex vbs_mutex, inter_mutex , qcd_mutex , data_mutex , _mutex; // vbs_mutex.lock();

//useful API

	/*
	 std::pair<type1, type2> myPair(100, "Tom");
	 std::cout << myPair.first << "\n";
    std::cout << myPair.second << "\n";
    vector<pair<type1,type2>> vec1
    vec1.push_back(make_pair(type1,type2))
	*/


template<typename T>
class Vector {
	public:
		~Vector() {
			delete[] elems_;
		}
		T & operator[](size_t n) {
			return elems_[n];
		}
		const T & operator[](size_t n) const {
			return elems_[n];
		}
		void PushBack(const  T & val) {
			size_++;
			elems_[size_-1] = val;
		}
		size_t Size() const {
			return size_;
		}
		size_t Capacity() const {
			return 1000;
		}
	private:
		size_t size_ = 0;
		T * elems_ = new T[1000];
};

template<typename T>
std::ostream & operator<<(std::ostream & os, const Vector<T> & v) {
	if (v.Size() == 0) return os << "[]";
	os << '[' << v[0];
	for (size_t i = 1; i < v.Size(); i++) {
		os << ", " << v[i];
	}
	return os << ']';
}



template<typename T>
class Set {
		using ValueType = T;
		using ConstReference = const T&;
		using SizeType = size_t;
	public:
		bool Contains(ConstReference val) {
			return count(begin(data_), end(data_), val) > 0;
		}
		void Insert(ConstReference val) {
			if (Contains(val)) return;
			data_.push_back(val);
		}
		void Erase(ConstReference val) {
			data_.erase(find(begin(data_), end(data_), val));
		}
		SizeType Size() const {
			return data_.size();
		}
	private:
		vector<ValueType> data_;
};

TH1F *significant_special_plt = new TH1F("_special_plt","_special_plt",100,-1,1);


struct signal_event {
	Long64_t signal_true;
	Long64_t signal_fake;
	Long64_t signal_total;
};

struct background_event {
	Long64_t background_true;
	Long64_t background_fake;
	Long64_t background_total;
};

////////////////////Global variabless
float alpha = 59.7*1000 ; //pb
float qcd_total_cross = 38.98 , inter_total_cross = 0.014475636607476456 , vbs_total_cross = 0.08880813263142302 ; // pb
float qcd_over_vbs = qcd_total_cross/vbs_total_cross ;
float vbs_weight=0 , qcd_weight=0 , inter_weight=0;
double global_va10[20000] = {} ;
float barral = 1.44 , endcap = 1.57 ;


double fitFunction( double &a0, double &a1, double &a2, double &a3 , float &mjj ) {
	double fitval = exp(a0 + a1*(mjj) + a2*(mjj*mjj) + a3*(mjj*mjj*mjj));
	return fitval;



}

////////////////////////////
TH1F *significant_plt = new TH1F("significant_plt","significant_plt",2002,-1,1);
TH1F *trueratio_plt = new TH1F("trueratio","trueratio",2000,-1,1);
TH1F *fakeratio_plt = new TH1F("fakeratio","fakeratio",2000,-1,1);

/////////////////////////////
//after cut//
TH1F *INTER_reco_invmass_j1_j2_plot_aftercut = new TH1F("INTER_reco_invmass_j1_j2_plot_aftercut","INTER_reco_invmass_j1_j2_plot_aftercut",1000,0,5000);
TH1F *INTER_reco_deltaeta_j1_j2_plot_aftercut = new TH1F("INTER_reco_j1j2_eta_plot_aftercut","INTER_reco_j1j2_eta_plot_aftercut",40,-10,10);

//cut1
TH1F *INTER_reco_deltaphi_j1_j2_plot_cut1 = new TH1F("INTER_reco_j1j2_phi_plot_cut1","INTER_reco_j1j2_phi_plot_cut1",40,-10,10);
TH1F *INTER_reco_deltaphi_p1_p2_plot_cut1 = new TH1F("INTER_reco_p1p2_phi_plot_cut1","INTER_reco_p1p2_phi_plot_cut1",40,-10,10);
TH1F *INTER_reco_deltaphi_p1_j1_plot_cut1 = new TH1F("INTER_reco_p1j1_phi_plot_cut1","INTER_reco_p1j1_phi_plot_cut1",40,-10,10);
TH1F *INTER_reco_deltaphi_p1_j2_plot_cut1 = new TH1F("INTER_reco_p1j2_phi_plot_cut1","INTER_reco_p1j2_phi_plot_cut1",40,-10,10);
TH1F *INTER_reco_deltaphi_p2_j1_plot_cut1 = new TH1F("INTER_reco_p2j1_phi_plot_cut1","INTER_reco_p2j1_phi_plot_cut1",40,-10,10);
TH1F *INTER_reco_deltaphi_p2_j2_plot_cut1 = new TH1F("INTER_reco_p2j2_phi_plot_cut1","INTER_reco_p2j2_phi_plot_cut1",40,-10,10);

//cut2
TH1F *INTER_reco_invmass_j1_j2_plot_cut2 = new TH1F("INTER_reco_invmass_j1_j2_plot_cut2","INTER_reco_invmass_j1_j2_plot_cut2",1000,0,5000);

//cut3
TH1F *INTER_reco_deltaeta_j1_j2_plot_cut3 = new TH1F("INTER_reco_deltaeta_j1_j2_plot_cut3","INTER_reco_deltaeta_j1_j2_plot_cut3",40,-10,10);


//BDT///
TH1F *INTER_BDT = new TH1F("INTER_BDT","INTER_BDT",100,-1,1);
TH1F *INTER_BDTG = new TH1F("INTER_BDTG","INTER_BDTG",40,-1,1);
//after BDTcut////////////////////////////////////
TH1F *INTER_reco_deltaeta_p1_j1_plot_true = new TH1F("INTER_reco_p1j1_eta_plot_true","INTER_reco_p1j1_eta_plot_true",50,0,10);
TH1F *INTER_reco_deltaeta_p1_j1_plot_fake = new TH1F("INTER_reco_p1j1_eta_plot_fake","INTER_reco_p1j1_eta_plot_fake",20,0,8);
TH1F *INTER_reco_invmass_j1_j2_plot_true = new TH1F("INTER_reco_invmass_j1_j2_plot_true","INTER_reco_invmass_j1_j2_plot_true",50,0,5000);
TH1F *INTER_reco_invmass_j1_j2_plot_fake = new TH1F("INTER_reco_invmass_j1_j2_plot_fake","INTER_reco_invmass_j1_j2_plot_fake",50,0,5000);
TH1F *INTER_reco_deltaphi_j1_j2_plot_true = new TH1F("INTER_reco_j1j2_phi_plot_true","INTER_reco_j1j2_phi_plot_true",50,0,5);
TH1F *INTER_reco_deltaphi_p1_p2_plot_true = new TH1F("INTER_reco_p1p2_phi_plot_true","INTER_reco_p1p2_phi_plot_true",50,0,5);
TH1F *INTER_reco_deltaeta_j1_j2_plot_true = new TH1F("INTER_reco_j1j2_eta_plot_true","INTER_reco_j1j2_eta_plot_true",50,0,10);
TH1F *INTER_reco_invmass_p1_p2_plot_true = new TH1F("INTER_reco_invmass_p1_p2_plot_true","INTER_reco_invmass_p1_p2_plot_true",100,0,2500);

/////shower-shape ////////////
TH1F *INTER_reco_r9_p1_bar = new TH1F("INTER_reco_r9_p1_true_bar","INTER_reco_r9_p1_true_bar",20,0,1);
TH1F *INTER_reco_r9_p2_bar = new TH1F("INTER_reco_r9_p2_true_bar","INTER_reco_r9_p2_true_bar",20,0,1);
TH1F *INTER_reco_sieie_p1_bar = new TH1F("INTER_reco_sieie_p1_true_bar","INTER_reco_sieie_p1_true_bar",50,0,0.08);
TH1F *INTER_reco_sieie_p2_bar = new TH1F("INTER_reco_sieie_p2_true_bar","INTER_reco_sieie_p2_true_bar",50,0,0.08);
TH1F *INTER_reco_hoe_p1_bar = new TH1F("INTER_reco_hoe_p1_true_bar","INTER_reco_hoe_p1_true_bar",30,0,0.15);
TH1F *INTER_reco_hoe_p2_bar = new TH1F("INTER_reco_hoe_p2_true_bar","INTER_reco_hoe_p2_true_bar",30,0,0.15);
TH1F *INTER_reco_chiso_p1_bar = new TH1F("INTER_reco_chiso_p1_true_bar","INTER_reco_chiso_p1_true_bar",20,0,5);
TH1F *INTER_reco_chiso_p2_bar = new TH1F("INTER_reco_chiso_p2_true_bar","INTER_reco_chiso_p2_true_bar",20,0,5);

TH1F *INTER_reco_r9_p1_end = new TH1F("INTER_reco_r9_p1_true_end","INTER_reco_r9_p1_true_end",20,0,1);
TH1F *INTER_reco_r9_p2_end = new TH1F("INTER_reco_r9_p2_true_end","INTER_reco_r9_p2_true_end",20,0,1);
TH1F *INTER_reco_sieie_p1_end = new TH1F("INTER_reco_sieie_p1_true_end","INTER_reco_sieie_p1_true_end",50,0,0.08);
TH1F *INTER_reco_sieie_p2_end = new TH1F("INTER_reco_sieie_p2_true_end","INTER_reco_sieie_p2_true_end",50,0,0.08);
TH1F *INTER_reco_hoe_p1_end = new TH1F("INTER_reco_hoe_p1_true_end","INTER_reco_hoe_p1_true_end",30,0,0.15);
TH1F *INTER_reco_hoe_p2_end = new TH1F("INTER_reco_hoe_p2_true_end","INTER_reco_hoe_p2_true_end",30,0,0.15);
TH1F *INTER_reco_chiso_p1_end = new TH1F("INTER_reco_chiso_p1_true_end","INTER_reco_chiso_p1_true_end",20,0,5);
TH1F *INTER_reco_chiso_p2_end = new TH1F("INTER_reco_chiso_p2_true_end","INTER_reco_chiso_p2_true_end",20,0,5);

TH1F *INTER_reco_mva_p1 = new TH1F("INTER_reco_mva_p1_true","INTER_reco_mva_p1_true",100,-1.1,1.1);
TH1F *INTER_reco_mva_p2 = new TH1F("INTER_reco_mva_p2_true","INTER_reco_mva_p2_true",100,-1.1,1.1);

/////////define plot VBS/////////////////// define plot VBS/////////////////////////////define plot VBS////////////////////////////
//after cut//
TH1F *VBS_reco_invmass_j1_j2_plot_aftercut = new TH1F("VBS_reco_invmass_j1_j2_plot_aftercut","VBS_reco_invmass_j1_j2_plot_aftercut",1000,0,5000);
TH1F *VBS_reco_deltaeta_j1_j2_plot_aftercut = new TH1F("VBS_reco_j1j2_eta_plot_aftercut","VBS_reco_j1j2_eta_plot_aftercut",40,-10,10);

//cut1
TH1F *VBS_reco_deltaphi_j1_j2_plot_cut1 = new TH1F("VBS_reco_j1j2_phi_plot_cut1","VBS_reco_j1j2_phi_plot_cut1",40,-10,10);
TH1F *VBS_reco_deltaphi_p1_p2_plot_cut1 = new TH1F("VBS_reco_p1p2_phi_plot_cut1","VBS_reco_p1p2_phi_plot_cut1",40,-10,10);
TH1F *VBS_reco_deltaphi_p1_j1_plot_cut1 = new TH1F("VBS_reco_p1j1_phi_plot_cut1","VBS_reco_p1j1_phi_plot_cut1",40,-10,10);
TH1F *VBS_reco_deltaphi_p1_j2_plot_cut1 = new TH1F("VBS_reco_p1j2_phi_plot_cut1","VBS_reco_p1j2_phi_plot_cut1",40,-10,10);
TH1F *VBS_reco_deltaphi_p2_j1_plot_cut1 = new TH1F("VBS_reco_p2j1_phi_plot_cut1","VBS_reco_p2j1_phi_plot_cut1",40,-10,10);
TH1F *VBS_reco_deltaphi_p2_j2_plot_cut1 = new TH1F("VBS_reco_p2j2_phi_plot_cut1","VBS_reco_p2j2_phi_plot_cut1",40,-10,10);

//cut2
TH1F *VBS_reco_invmass_j1_j2_plot_cut2 = new TH1F("VBS_reco_invmass_j1_j2_plot_cut2","VBS_reco_invmass_j1_j2_plot_cut2",1000,0,5000);

//cut3
TH1F *VBS_reco_deltaeta_j1_j2_plot_cut3 = new TH1F("VBS_reco_deltaeta_j1_j2_plot_cut3","VBS_reco_deltaeta_j1_j2_plot_cut3",40,-10,10);


//BDT///
TH1F *VBS_BDT = new TH1F("VBS_BDT","VBS_BDT",100,-1,1);
TH1F *VBS_BDTG = new TH1F("VBS_BDTG","VBS_BDTG",40,-1,1);
//after BDTcut
TH1F *VBS_reco_deltaeta_p1_j1_plot_true = new TH1F("VBS_reco_p1j1_eta_plot_true","VBS_reco_p1j1_eta_plot_true",50,0,10);
TH1F *VBS_reco_deltaeta_p1_j1_plot_fake = new TH1F("VBS_reco_p1j1_eta_plot_fake","VBS_reco_p1j1_eta_plot_fake",20,0,8);
TH1F *VBS_reco_invmass_j1_j2_plot_true = new TH1F("VBS_reco_invmass_j1_j2_plot_true","VBS_reco_invmass_j1_j2_plot_true",50,0,5000);
TH1F *VBS_reco_invmass_j1_j2_plot_fake = new TH1F("VBS_reco_invmass_j1_j2_plot_fake","VBS_reco_invmass_j1_j2_plot_fake",50,0,5000);
TH1F *VBS_reco_deltaphi_j1_j2_plot_true = new TH1F("VBS_reco_j1j2_phi_plot_true","VBS_reco_j1j2_phi_plot_true",50,0,5);
TH1F *VBS_reco_deltaphi_p1_p2_plot_true = new TH1F("VBS_reco_p1p2_phi_plot_true","VBS_reco_p1p2_phi_plot_true",50,0,5);
TH1F *VBS_reco_deltaeta_j1_j2_plot_true = new TH1F("VBS_reco_j1j2_eta_plot_true","VBS_reco_j1j2_eta_plot_true",50,0,10);
TH1F *VBS_reco_invmass_p1_p2_plot_true = new TH1F("VBS_reco_invmass_p1_p2_plot_true","VBS_reco_invmass_p1_p2_plot_true",100,0,2500);

/////shower-shape ////////////
TH1F *VBS_reco_r9_p1_bar = new TH1F("VBS_reco_r9_p1_true_bar","VBS_reco_r9_p1_true_bar",20,0,1);
TH1F *VBS_reco_r9_p2_bar = new TH1F("VBS_reco_r9_p2_true_bar","VBS_reco_r9_p2_true_bar",20,0,1);
TH1F *VBS_reco_sieie_p1_bar = new TH1F("VBS_reco_sieie_p1_true_bar","VBS_reco_sieie_p1_true_bar",50,0,0.08);
TH1F *VBS_reco_sieie_p2_bar = new TH1F("VBS_reco_sieie_p2_true_bar","VBS_reco_sieie_p2_true_bar",50,0,0.08);
TH1F *VBS_reco_hoe_p1_bar = new TH1F("VBS_reco_hoe_p1_true_bar","VBS_reco_hoe_p1_true_bar",30,0,0.15);
TH1F *VBS_reco_hoe_p2_bar = new TH1F("VBS_reco_hoe_p2_true_bar","VBS_reco_hoe_p2_true_bar",30,0,0.15);
TH1F *VBS_reco_chiso_p1_bar = new TH1F("VBS_reco_chiso_p1_true_bar","VBS_reco_chiso_p1_true_bar",20,0,5);
TH1F *VBS_reco_chiso_p2_bar = new TH1F("VBS_reco_chiso_p2_true_bar","VBS_reco_chiso_p2_true_bar",20,0,5);

TH1F *VBS_reco_r9_p1_end = new TH1F("VBS_reco_r9_p1_true_end","VBS_reco_r9_p1_true_end",20,0,1);
TH1F *VBS_reco_r9_p2_end = new TH1F("VBS_reco_r9_p2_true_end","VBS_reco_r9_p2_true_end",20,0,1);
TH1F *VBS_reco_sieie_p1_end = new TH1F("VBS_reco_sieie_p1_true_end","VBS_reco_sieie_p1_true_end",50,0,0.08);
TH1F *VBS_reco_sieie_p2_end = new TH1F("VBS_reco_sieie_p2_true_end","VBS_reco_sieie_p2_true_end",50,0,0.08);
TH1F *VBS_reco_hoe_p1_end = new TH1F("VBS_reco_hoe_p1_true_end","VBS_reco_hoe_p1_true_end",30,0,0.15);
TH1F *VBS_reco_hoe_p2_end = new TH1F("VBS_reco_hoe_p2_true_end","VBS_reco_hoe_p2_true_end",30,0,0.15);
TH1F *VBS_reco_chiso_p1_end = new TH1F("VBS_reco_chiso_p1_true_end","VBS_reco_chiso_p1_true_end",20,0,5);
TH1F *VBS_reco_chiso_p2_end = new TH1F("VBS_reco_chiso_p2_true_end","VBS_reco_chiso_p2_true_end",20,0,5);

TH1F *VBS_reco_mva_p1 = new TH1F("VBS_reco_mva_p1_true","VBS_reco_mva_p1_true",100,-1.1,1.1);
TH1F *VBS_reco_mva_p2 = new TH1F("VBS_reco_mva_p2_true","VBS_reco_mva_p2_true",100,-1.1,1.1);

///////////////////////////////////////////


/////////define plot QCD/////////////////// define plot QCD/////////////////////////////define plot QCD//////////////////////////////
//cut1
TH1F *QCD_reco_deltaphi_j1_j2_plot_cut1 = new TH1F("QCD_reco_j1j2_phi_plot_cut1","QCD_reco_j1j2_phi_plot_cut1",40,-10,10);
TH1F *QCD_reco_deltaphi_p1_p2_plot_cut1 = new TH1F("QCD_reco_p1p2_phi_plot_cut1","QCD_reco_p1p2_phi_plot_cut1",40,-10,10);
TH1F *QCD_reco_deltaphi_p1_j1_plot_cut1 = new TH1F("QCD_reco_p1j1_phi_plot_cut1","QCD_reco_p1j1_phi_plot_cut1",40,-10,10);
TH1F *QCD_reco_deltaphi_p1_j2_plot_cut1 = new TH1F("QCD_reco_p1j2_phi_plot_cut1","QCD_reco_p1j2_phi_plot_cut1",40,-10,10);
TH1F *QCD_reco_deltaphi_p2_j1_plot_cut1 = new TH1F("QCD_reco_p2j1_phi_plot_cut1","QCD_reco_p2j1_phi_plot_cut1",40,-10,10);
TH1F *QCD_reco_deltaphi_p2_j2_plot_cut1 = new TH1F("QCD_reco_p2j2_phi_plot_cut1","QCD_reco_p2j2_phi_plot_cut1",40,-10,10);

//cut2
TH1F *QCD_reco_invmass_j1_j2_plot_cut2 = new TH1F("QCD_reco_invmass_j1_j2_plot_cut2","QCD_reco_invmass_j1_j2_plot_cut2",100,0,5000);

//cut3
TH1F *QCD_reco_deltaeta_j1_j2_plot_cut3 = new TH1F("QCD_reco_deltaeta_j1_j2_plot_cut3","QCD_reco_deltaeta_j1_j2_plot_cut3",40,-10,10);

//BDT//
TH1F *QCD_BDT = new TH1F("QCD_BDT","QCD_BDT",100,-1,1);
TH1F *QCD_BDTG = new TH1F("QCD_BDTG","QCD_BDTG",40,-1,1);
//after BDTcut
TH1F *QCD_reco_deltaeta_p1_j1_plot_true = new TH1F("QCD_reco_p1j1_eta_plot_true","QCD_reco_p1j1_eta_plot_true",50,0,10);
TH1F *QCD_reco_deltaeta_p1_j1_plot_fake = new TH1F("QCD_reco_p1j1_eta_plot_fake","QCD_reco_p1j1_eta_plot_fake",20,0,8);
TH1F *QCD_reco_invmass_j1_j2_plot_true = new TH1F("QCD_reco_invmass_j1_j2_plot_true","QCD_reco_invmass_j1_j2_plot_true",50,0,5000);
TH1F *QCD_reco_invmass_j1_j2_plot_fake = new TH1F("QCD_reco_invmass_j1_j2_plot_fake","QCD_reco_invmass_j1_j2_plot_fake",50,0,5000);
TH1F *QCD_reco_deltaphi_j1_j2_plot_true = new TH1F("QCD_reco_j1j2_phi_plot_true","QCD_reco_j1j2_phi_plot_true",50,0,5);
TH1F *QCD_reco_deltaphi_p1_p2_plot_true = new TH1F("QCD_reco_p1p2_phi_plot_true","QCD_reco_p1p2_phi_plot_true",50,0,5);
TH1F *QCD_reco_deltaeta_j1_j2_plot_true = new TH1F("QCD_reco_j1j2_eta_plot_true","QCD_reco_j1j2_eta_plot_true",50,0,10);
TH1F *QCD_reco_invmass_p1_p2_plot_true = new TH1F("QCD_reco_invmass_p1_p2_plot_true","QCD_reco_invmass_p1_p2_plot_true",100,0,2500);

/////shower-shape ////////////
TH1F *QCD_reco_r9_p1_bar = new TH1F("QCD_reco_r9_p1_true_bar","QCD_reco_r9_p1_true_bar",20,0,1);
TH1F *QCD_reco_r9_p2_bar = new TH1F("QCD_reco_r9_p2_true_bar","QCD_reco_r9_p2_true_bar",20,0,1);
TH1F *QCD_reco_sieie_p1_bar = new TH1F("QCD_reco_sieie_p1_true_bar","QCD_reco_sieie_p1_true_bar",50,0,0.08);
TH1F *QCD_reco_sieie_p2_bar = new TH1F("QCD_reco_sieie_p2_true_bar","QCD_reco_sieie_p2_true_bar",50,0,0.08);
TH1F *QCD_reco_hoe_p1_bar = new TH1F("QCD_reco_hoe_p1_true_bar","QCD_reco_hoe_p1_true_bar",30,0,0.15);
TH1F *QCD_reco_hoe_p2_bar = new TH1F("QCD_reco_hoe_p2_true_bar","QCD_reco_hoe_p2_true_bar",30,0,0.15);
TH1F *QCD_reco_chiso_p1_bar = new TH1F("QCD_reco_chiso_p1_true_bar","QCD_reco_chiso_p1_true_bar",20,0,5);
TH1F *QCD_reco_chiso_p2_bar = new TH1F("QCD_reco_chiso_p2_true_bar","QCD_reco_chiso_p2_true_bar",20,0,5);

TH1F *QCD_reco_r9_p1_end = new TH1F("QCD_reco_r9_p1_true_end","QCD_reco_r9_p1_true_end",20,0,1);
TH1F *QCD_reco_r9_p2_end = new TH1F("QCD_reco_r9_p2_true_end","QCD_reco_r9_p2_true_end",20,0,1);
TH1F *QCD_reco_sieie_p1_end = new TH1F("QCD_reco_sieie_p1_true_end","QCD_reco_sieie_p1_true_end",50,0,0.08);
TH1F *QCD_reco_sieie_p2_end = new TH1F("QCD_reco_sieie_p2_true_end","QCD_reco_sieie_p2_true_end",50,0,0.08);
TH1F *QCD_reco_hoe_p1_end = new TH1F("QCD_reco_hoe_p1_true_end","QCD_reco_hoe_p1_true_end",30,0,0.15);
TH1F *QCD_reco_hoe_p2_end = new TH1F("QCD_reco_hoe_p2_true_end","QCD_reco_hoe_p2_true_end",30,0,0.15);
TH1F *QCD_reco_chiso_p1_end = new TH1F("QCD_reco_chiso_p1_true_end","QCD_reco_chiso_p1_true_end",20,0,5);
TH1F *QCD_reco_chiso_p2_end = new TH1F("QCD_reco_chiso_p2_true_end","QCD_reco_chiso_p2_true_end",20,0,5);

TH1F *QCD_reco_mva_p1 = new TH1F("QCD_reco_mva_p1_true","QCD_reco_mva_p1_true",100,-1.1,1.1);
TH1F *QCD_reco_mva_p2 = new TH1F("QCD_reco_mva_p2_true","QCD_reco_mva_p2_true",100,-1.1,1.1);

///////////Data///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//after cut//
TH1F *Data_reco_invmass_j1_j2_plot_aftercut = new TH1F("Data_reco_invmass_j1_j2_plot_aftercut","Data_reco_invmass_j1_j2_plot_aftercut",1000,0,5000);
TH1F *Data_reco_deltaeta_j1_j2_plot_aftercut = new TH1F("Data_reco_j1j2_eta_plot_aftercut","Data_reco_j1j2_eta_plot_aftercut",40,-10,10);

//cut1
TH1F *Data_reco_deltaphi_j1_j2_plot_cut1 = new TH1F("Data_reco_j1j2_phi_plot_cut1","Data_reco_j1j2_phi_plot_cut1",40,-10,10);
TH1F *Data_reco_deltaphi_p1_p2_plot_cut1 = new TH1F("Data_reco_p1p2_phi_plot_cut1","Data_reco_p1p2_phi_plot_cut1",40,-10,10);
TH1F *Data_reco_deltaphi_p1_j1_plot_cut1 = new TH1F("Data_reco_p1j1_phi_plot_cut1","Data_reco_p1j1_phi_plot_cut1",40,-10,10);
TH1F *Data_reco_deltaphi_p1_j2_plot_cut1 = new TH1F("Data_reco_p1j2_phi_plot_cut1","Data_reco_p1j2_phi_plot_cut1",40,-10,10);
TH1F *Data_reco_deltaphi_p2_j1_plot_cut1 = new TH1F("Data_reco_p2j1_phi_plot_cut1","Data_reco_p2j1_phi_plot_cut1",40,-10,10);
TH1F *Data_reco_deltaphi_p2_j2_plot_cut1 = new TH1F("Data_reco_p2j2_phi_plot_cut1","Data_reco_p2j2_phi_plot_cut1",40,-10,10);

//cut2
TH1F *Data_reco_invmass_j1_j2_plot_cut2 = new TH1F("Data_reco_invmass_j1_j2_plot_cut2","Data_reco_invmass_j1_j2_plot_cut2",1000,0,5000);

//cut3
TH1F *Data_reco_deltaeta_j1_j2_plot_cut3 = new TH1F("Data_reco_deltaeta_j1_j2_plot_cut3","Data_reco_deltaeta_j1_j2_plot_cut3",40,-10,10);


//BDT///
TH1F *Data_BDT = new TH1F("Data_BDT","Data_BDT",100,-1,1);
TH1F *Data_BDTG = new TH1F("Data_BDTG","Data_BDTG",100,-1,1);
//after BDTcut
TH1F *Data_reco_deltaeta_p1_j1_plot_true = new TH1F("Data_reco_p1j1_eta_plot_true","Data_reco_p1j1_eta_plot_true",50,0,10);
TH1F *Data_reco_deltaeta_p1_j1_plot_fake = new TH1F("Data_reco_p1j1_eta_plot_fake","Data_reco_p1j1_eta_plot_fake",20,0,8);
TH1F *Data_reco_invmass_j1_j2_plot_true = new TH1F("Data_reco_invmass_j1_j2_plot_true","Data_reco_invmass_j1_j2_plot_true",50,0,5000);
TH1F *Data_reco_invmass_j1_j2_plot_fake = new TH1F("Data_reco_invmass_j1_j2_plot_fake","Data_reco_invmass_j1_j2_plot_fake",50,0,5000);
TH1F *Data_reco_deltaphi_j1_j2_plot_true = new TH1F("Data_reco_j1j2_phi_plot_true","Data_reco_j1j2_phi_plot_true",50,0,5);
TH1F *Data_reco_deltaphi_p1_p2_plot_true = new TH1F("Data_reco_p1p2_phi_plot_true","Data_reco_p1p2_phi_plot_true",50,0,5);
TH1F *Data_reco_deltaeta_j1_j2_plot_true = new TH1F("Data_reco_j1j2_eta_plot_true","Data_reco_j1j2_eta_plot_true",50,0,8);
TH1F *Data_reco_invmass_p1_p2_plot_true = new TH1F("Data_reco_invmass_p1_p2_plot_true","Data_reco_invmass_p1_p2_plot_true",100,0,2500);

/////shower-shape ////////////
TH1F *Data_reco_r9_p1_bar = new TH1F("Data_reco_r9_p1_true_bar","Data_reco_r9_p1_true_bar",20,0,1);
TH1F *Data_reco_r9_p2_bar = new TH1F("Data_reco_r9_p2_true_bar","Data_reco_r9_p2_true_bar",20,0,1);
TH1F *Data_reco_sieie_p1_bar = new TH1F("Data_reco_sieie_p1_true_bar","Data_reco_sieie_p1_true_bar",50,0,0.08);
TH1F *Data_reco_sieie_p2_bar = new TH1F("Data_reco_sieie_p2_true_bar","Data_reco_sieie_p2_true_bar",50,0,0.08);
TH1F *Data_reco_hoe_p1_bar = new TH1F("Data_reco_hoe_p1_true_bar","Data_reco_hoe_p1_true_bar",30,0,0.15);
TH1F *Data_reco_hoe_p2_bar = new TH1F("Data_reco_hoe_p2_true_bar","Data_reco_hoe_p2_true_bar",30,0,0.15);
TH1F *Data_reco_chiso_p1_bar = new TH1F("Data_reco_chiso_p1_true_bar","Data_reco_chiso_p1_true_bar",20,0,5);
TH1F *Data_reco_chiso_p2_bar = new TH1F("Data_reco_chiso_p2_true_bar","Data_reco_chiso_p2_true_bar",20,0,5);

TH1F *Data_reco_r9_p1_end = new TH1F("Data_reco_r9_p1_true_end","Data_reco_r9_p1_true_end",20,0,1);
TH1F *Data_reco_r9_p2_end = new TH1F("Data_reco_r9_p2_true_end","Data_reco_r9_p2_true_end",20,0,1);
TH1F *Data_reco_sieie_p1_end = new TH1F("Data_reco_sieie_p1_true_end","Data_reco_sieie_p1_true_end",50,0,0.08);
TH1F *Data_reco_sieie_p2_end = new TH1F("Data_reco_sieie_p2_true_end","Data_reco_sieie_p2_true_end",50,0,0.08);
TH1F *Data_reco_hoe_p1_end = new TH1F("Data_reco_hoe_p1_true_end","Data_reco_hoe_p1_true_end",30,0,0.15);
TH1F *Data_reco_hoe_p2_end = new TH1F("Data_reco_hoe_p2_true_end","Data_reco_hoe_p2_true_end",30,0,0.15);
TH1F *Data_reco_chiso_p1_end = new TH1F("Data_reco_chiso_p1_true_end","Data_reco_chiso_p1_true_end",20,0,5);
TH1F *Data_reco_chiso_p2_end = new TH1F("Data_reco_chiso_p2_true_end","Data_reco_chiso_p2_true_end",20,0,5);

TH1F *Data_reco_mva_p1 = new TH1F("Data_reco_mva_p1_true","Data_reco_mva_p1_true",100,-1.1,1.1);
TH1F *Data_reco_mva_p2 = new TH1F("Data_reco_mva_p2_true","Data_reco_mva_p2_true",100,-1.1,1.1);

////////////Stack//////////
auto hs_deltap1j1 = new THStack("hs","etap1j1");
auto hs_BDTG = new THStack("hsBDTG ","BDTG ");
auto hs_invmass_j1_j2 = new THStack("hs_invmass_j1_j2","hs_invmass_j1_j2");
auto hs_invmass_p1_p2 = new THStack("hs_invmass_p1_p2","hs_invmass_p1_p2");




////////////////////FIT///////////////////////

//TF1 *fitFcn = new TF1("fitFcn","fitFcn",0,3,5000);

////////////////////cut selection///////////////////////



//////////////////////////////////////////

int vbs( int &events , int &true_signal , int &fake_signal , double &BDT_cut  ) {
	_mutex.lock();
	string file02 = "new_BDT_vbs.root" ;
	cout << file02 <<endl;
	TFile *f2 = new TFile(file02.c_str());
	TTree *TreeSNew;
	f2->GetObject("TreeSNew",TreeSNew);

	events = TreeSNew->GetEntries() ;



	float reco_invmass_j1j2_SNew;
	float reco_deltaeta_j1j2_SNew;
	float reco_deltaphi_j1j2_SNew;
	float reco_deltaR_j1j2_SNew;
	float reco_deltapt_j1j2_SNew;

	float reco_invmass_p1p2_SNew;
	float reco_deltaeta_p1p2_SNew;
	float reco_deltaphi_p1p2_SNew;
	float reco_deltaR_p1p2_SNew;
	float reco_deltapt_p1p2_SNew;

	float reco_deltaeta_p1j1_SNew;
	float reco_deltaeta_p1j2_SNew;
	float reco_deltaeta_p2j1_SNew;
	float reco_deltaeta_p2j2_SNew;

	float reco_deltaphi_p1j1_SNew;
	float reco_deltaphi_p1j2_SNew;
	float reco_deltaphi_p2j1_SNew;
	float reco_deltaphi_p2j2_SNew;

	float reco_deltaR_p1j1_SNew;
	float reco_deltaR_p1j2_SNew;
	float reco_deltaR_p2j1_SNew;
	float reco_deltaR_p2j2_SNew;

	float reco_deltapt_p1j1_SNew;
	float reco_deltapt_p1j2_SNew;
	float reco_deltapt_p2j1_SNew;
	float reco_deltapt_p2j2_SNew;
	float BDT_SNew;
	float BDTG_SNew;

	float reco_pho1_r9_SNew , reco_pho1_hoe_SNew , reco_pho1_chiso_SNew , reco_pho1_sieie_SNew , reco_pho1_mva_SNew ;
	float reco_pho2_r9_SNew , reco_pho2_hoe_SNew , reco_pho2_chiso_SNew , reco_pho2_sieie_SNew , reco_pho2_mva_SNew ;
	float reco_eta_p1_SNew , reco_eta_p2_SNew ;

	TreeSNew->SetBranchAddress("reco_invmass_j1j2", &reco_invmass_j1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaeta_j1j2", &reco_deltaeta_j1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaphi_j1j2", &reco_deltaphi_j1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaR_j1j2", &reco_deltaR_j1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltapt_j1j2", &reco_deltapt_j1j2_SNew);

	TreeSNew->SetBranchAddress("reco_invmass_p1p2", &reco_invmass_p1p2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaeta_p1p2", &reco_deltaeta_p1p2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaphi_p1p2", &reco_deltaphi_p1p2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaR_p1p2", &reco_deltaR_p1p2_SNew);
	TreeSNew->SetBranchAddress("reco_deltapt_p1p2", &reco_deltapt_p1p2_SNew);

	TreeSNew->SetBranchAddress("reco_deltaeta_p1j1", &reco_deltaeta_p1j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltaeta_p1j2", &reco_deltaeta_p1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaeta_p2j1", &reco_deltaeta_p2j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltaeta_p2j2", &reco_deltaeta_p2j2_SNew);

	TreeSNew->SetBranchAddress("reco_deltaphi_p1j1", &reco_deltaphi_p1j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltaphi_p1j2", &reco_deltaphi_p1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaphi_p2j1", &reco_deltaphi_p2j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltaphi_p2j2", &reco_deltaphi_p2j2_SNew);

	TreeSNew->SetBranchAddress("reco_deltaR_p1j1", &reco_deltaR_p1j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltaR_p1j2", &reco_deltaR_p1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltaR_p2j1", &reco_deltaR_p2j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltaR_p2j2", &reco_deltaR_p2j2_SNew);

	TreeSNew->SetBranchAddress("reco_deltapt_p1j1", &reco_deltapt_p1j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltapt_p1j2", &reco_deltapt_p1j2_SNew);
	TreeSNew->SetBranchAddress("reco_deltapt_p2j1", &reco_deltapt_p2j1_SNew);
	TreeSNew->SetBranchAddress("reco_deltapt_p2j2", &reco_deltapt_p2j2_SNew);
	TreeSNew->SetBranchAddress("BDT", &BDT_SNew);
	TreeSNew->SetBranchAddress("BDTG", &BDTG_SNew);

	TreeSNew->SetBranchAddress("reco_pho1_r9", &reco_pho1_r9_SNew);
	TreeSNew->SetBranchAddress("reco_pho1_hoe", &reco_pho1_hoe_SNew);
	TreeSNew->SetBranchAddress("reco_pho1_chiso", &reco_pho1_chiso_SNew);
	TreeSNew->SetBranchAddress("reco_pho1_sieie", &reco_pho1_sieie_SNew);
	TreeSNew->SetBranchAddress("reco_pho1_mva", &reco_pho1_mva_SNew);
	TreeSNew->SetBranchAddress("reco_eta_p1", &reco_eta_p1_SNew);

	TreeSNew->SetBranchAddress("reco_pho2_r9", &reco_pho2_r9_SNew);
	TreeSNew->SetBranchAddress("reco_pho2_hoe", &reco_pho2_hoe_SNew);
	TreeSNew->SetBranchAddress("reco_pho2_chiso", &reco_pho2_chiso_SNew);
	TreeSNew->SetBranchAddress("reco_pho2_sieie", &reco_pho2_sieie_SNew);
	TreeSNew->SetBranchAddress("reco_pho2_mva", &reco_pho2_mva_SNew);
	TreeSNew->SetBranchAddress("reco_eta_p2", &reco_eta_p2_SNew);

	for(int evt=0 ; evt < events ; evt++) {

		TreeSNew->GetEntry(evt);

		VBS_BDT->Fill(BDT_SNew,vbs_weight);
		VBS_BDTG->Fill(BDTG_SNew,vbs_weight);





		//true signal
		if(BDTG_SNew > BDT_cut && reco_invmass_p1p2_SNew > 10 && reco_invmass_j1j2_SNew > 400) {
			VBS_reco_deltaeta_p1_j1_plot_true->Fill(reco_deltaeta_p1j1_SNew,vbs_weight);
			VBS_reco_invmass_j1_j2_plot_true->Fill(reco_invmass_j1j2_SNew,vbs_weight);
			VBS_reco_deltaphi_j1_j2_plot_true->Fill(reco_deltaphi_j1j2_SNew,vbs_weight);
			VBS_reco_deltaphi_p1_p2_plot_true->Fill(reco_deltaphi_p1p2_SNew,vbs_weight);
			VBS_reco_deltaeta_j1_j2_plot_true->Fill(reco_deltaeta_j1j2_SNew,vbs_weight);
			VBS_reco_invmass_p1_p2_plot_true->Fill(reco_invmass_p1p2_SNew,vbs_weight);
			VBS_reco_mva_p1->Fill(reco_pho1_mva_SNew,vbs_weight);
			VBS_reco_mva_p2->Fill(reco_pho2_mva_SNew,vbs_weight);

			//barral
			if(abs(reco_eta_p1_SNew) > 0 && abs(reco_eta_p1_SNew) < 1.44) {
				VBS_reco_r9_p1_bar->Fill(reco_pho1_r9_SNew,vbs_weight);
				VBS_reco_hoe_p1_bar->Fill(reco_pho1_hoe_SNew,vbs_weight);
				VBS_reco_chiso_p1_bar->Fill(reco_pho1_chiso_SNew,vbs_weight);
				VBS_reco_sieie_p1_bar->Fill(reco_pho1_sieie_SNew,vbs_weight);
			}

			if(abs(reco_eta_p2_SNew) > 0 && abs(reco_eta_p2_SNew) < 1.44) {
				VBS_reco_r9_p2_bar->Fill(reco_pho2_r9_SNew,vbs_weight);
				VBS_reco_hoe_p2_bar->Fill(reco_pho2_hoe_SNew,vbs_weight);
				VBS_reco_chiso_p2_bar->Fill(reco_pho2_chiso_SNew,vbs_weight);
				VBS_reco_sieie_p2_bar->Fill(reco_pho2_sieie_SNew,vbs_weight);
			}

			//endcap
			if(abs(reco_eta_p1_SNew) > 1.57 && abs(reco_eta_p1_SNew) < 2.5) {
				VBS_reco_r9_p1_end->Fill(reco_pho1_r9_SNew,vbs_weight);
				VBS_reco_hoe_p1_end->Fill(reco_pho1_hoe_SNew,vbs_weight);
				VBS_reco_chiso_p1_end->Fill(reco_pho1_chiso_SNew,vbs_weight);
				VBS_reco_sieie_p1_end->Fill(reco_pho1_sieie_SNew,vbs_weight);
			}

			if(abs(reco_eta_p2_SNew) > 1.57 && abs(reco_eta_p2_SNew) < 2.5) {
				VBS_reco_r9_p2_end->Fill(reco_pho2_r9_SNew,vbs_weight);
				VBS_reco_hoe_p2_end->Fill(reco_pho2_hoe_SNew,vbs_weight);
				VBS_reco_chiso_p2_end->Fill(reco_pho2_chiso_SNew,vbs_weight);
				VBS_reco_sieie_p2_end->Fill(reco_pho2_sieie_SNew,vbs_weight);
			}

			//
			double a0=-12.5 , a1=11*10^(-4) , a2=-6*10^(-7) , a3=4*10^(-11) ;

			true_signal++;
		}


		//fake signal
		if(BDTG_SNew <= BDT_cut) {
			VBS_reco_deltaeta_p1_j1_plot_fake->Fill(reco_deltaeta_p1j1_SNew,vbs_weight);
			VBS_reco_invmass_j1_j2_plot_fake->Fill(reco_invmass_j1j2_SNew,vbs_weight);
			//		if( evt < 5) cout << "V_F " << reco_deltaeta_p1j1_SNew <<endl;
			fake_signal++;
		}


	}




	f2->Close();
	_mutex.unlock();
	//signal_event result =   { 0, 0, TreeSNew->GetEntries() };
	return 0 ;
}


int qcd( int &events , int &true_background , int &fake_background , double &BDT_cut ) {
	_mutex.lock();
	string file02 = "new_BDT_qcd.root" ;
	cout << file02 <<endl;
	TFile *f2 = new TFile(file02.c_str());
	TTree *TreeBNew;
	f2->GetObject("TreeBNew",TreeBNew);
	events = TreeBNew->GetEntries() ;


	float reco_invmass_j1j2_BNew;
	float reco_deltaeta_j1j2_BNew;
	float reco_deltaphi_j1j2_BNew;
	float reco_deltaR_j1j2_BNew;
	float reco_deltapt_j1j2_BNew;

	float reco_invmass_p1p2_BNew;
	float reco_deltaeta_p1p2_BNew;
	float reco_deltaphi_p1p2_BNew;
	float reco_deltaR_p1p2_BNew;
	float reco_deltapt_p1p2_BNew;

	float reco_deltaeta_p1j1_BNew;
	float reco_deltaeta_p1j2_BNew;
	float reco_deltaeta_p2j1_BNew;
	float reco_deltaeta_p2j2_BNew;

	float reco_deltaphi_p1j1_BNew;
	float reco_deltaphi_p1j2_BNew;
	float reco_deltaphi_p2j1_BNew;
	float reco_deltaphi_p2j2_BNew;

	float reco_deltaR_p1j1_BNew;
	float reco_deltaR_p1j2_BNew;
	float reco_deltaR_p2j1_BNew;
	float reco_deltaR_p2j2_BNew;

	float reco_deltapt_p1j1_BNew;
	float reco_deltapt_p1j2_BNew;
	float reco_deltapt_p2j1_BNew;
	float reco_deltapt_p2j2_BNew;
	float BDT_BNew;
	float BDTG_BNew;

	float reco_pho1_r9_BNew , reco_pho1_hoe_BNew , reco_pho1_chiso_BNew , reco_pho1_sieie_BNew , reco_pho1_mva_BNew ;
	float reco_pho2_r9_BNew , reco_pho2_hoe_BNew , reco_pho2_chiso_BNew , reco_pho2_sieie_BNew , reco_pho2_mva_BNew ;
	float reco_eta_p1_BNew , reco_eta_p2_BNew ;


	TreeBNew->SetBranchAddress("reco_invmass_j1j2", &reco_invmass_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_j1j2", &reco_deltaeta_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_j1j2", &reco_deltaphi_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_j1j2", &reco_deltaR_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_j1j2", &reco_deltapt_j1j2_BNew);

	TreeBNew->SetBranchAddress("reco_invmass_p1p2", &reco_invmass_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p1p2", &reco_deltaeta_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p1p2", &reco_deltaphi_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p1p2", &reco_deltaR_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p1p2", &reco_deltapt_p1p2_BNew);

	TreeBNew->SetBranchAddress("reco_deltaeta_p1j1", &reco_deltaeta_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p1j2", &reco_deltaeta_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p2j1", &reco_deltaeta_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p2j2", &reco_deltaeta_p2j2_BNew);

	TreeBNew->SetBranchAddress("reco_deltaphi_p1j1", &reco_deltaphi_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p1j2", &reco_deltaphi_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p2j1", &reco_deltaphi_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p2j2", &reco_deltaphi_p2j2_BNew);

	TreeBNew->SetBranchAddress("reco_deltaR_p1j1", &reco_deltaR_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p1j2", &reco_deltaR_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p2j1", &reco_deltaR_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p2j2", &reco_deltaR_p2j2_BNew);

	TreeBNew->SetBranchAddress("reco_deltapt_p1j1", &reco_deltapt_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p1j2", &reco_deltapt_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p2j1", &reco_deltapt_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p2j2", &reco_deltapt_p2j2_BNew);
	TreeBNew->SetBranchAddress("BDT", &BDT_BNew);
	TreeBNew->SetBranchAddress("BDTG", &BDTG_BNew);

	TreeBNew->SetBranchAddress("reco_pho1_r9", &reco_pho1_r9_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_hoe", &reco_pho1_hoe_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_chiso", &reco_pho1_chiso_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_sieie", &reco_pho1_sieie_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_mva", &reco_pho1_mva_BNew);
	TreeBNew->SetBranchAddress("reco_eta_p1", &reco_eta_p1_BNew);

	TreeBNew->SetBranchAddress("reco_pho2_r9", &reco_pho2_r9_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_hoe", &reco_pho2_hoe_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_chiso", &reco_pho2_chiso_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_sieie", &reco_pho2_sieie_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_mva", &reco_pho2_mva_BNew);
	TreeBNew->SetBranchAddress("reco_eta_p2", &reco_eta_p2_BNew);

	_mutex.unlock();

	for(int evt=0 ; evt < events ; evt++) {

		TreeBNew->GetEntry(evt);



		QCD_BDT->Fill(BDT_BNew,qcd_weight);
		QCD_BDTG->Fill(BDTG_BNew,qcd_weight);
		//true background
		if(BDTG_BNew >= BDT_cut && reco_invmass_p1p2_BNew > 10  && reco_invmass_j1j2_BNew > 400 ) {
			QCD_reco_deltaeta_p1_j1_plot_true->Fill(reco_deltaeta_p1j1_BNew,qcd_weight);
			QCD_reco_invmass_j1_j2_plot_true->Fill(reco_invmass_j1j2_BNew,qcd_weight);
			QCD_reco_deltaphi_j1_j2_plot_true->Fill(reco_deltaphi_j1j2_BNew,qcd_weight);
			QCD_reco_deltaphi_p1_p2_plot_true->Fill(reco_deltaphi_p1p2_BNew,qcd_weight);
			QCD_reco_deltaeta_j1_j2_plot_true->Fill(reco_deltaeta_j1j2_BNew,qcd_weight);
			QCD_reco_invmass_p1_p2_plot_true->Fill(reco_invmass_p1p2_BNew,qcd_weight);
			QCD_reco_mva_p1->Fill(reco_pho1_mva_BNew,qcd_weight);
			QCD_reco_mva_p2->Fill(reco_pho2_mva_BNew,qcd_weight);

			//barral
			if(abs(reco_eta_p1_BNew) > 0 && abs(reco_eta_p1_BNew) < 1.44) {
				QCD_reco_r9_p1_bar->Fill(reco_pho1_r9_BNew,qcd_weight);
				QCD_reco_hoe_p1_bar->Fill(reco_pho1_hoe_BNew,qcd_weight);
				QCD_reco_chiso_p1_bar->Fill(reco_pho1_chiso_BNew,qcd_weight);
				QCD_reco_sieie_p1_bar->Fill(reco_pho1_sieie_BNew,qcd_weight);
			}

			if(abs(reco_eta_p2_BNew) > 0 && abs(reco_eta_p2_BNew) < 1.44) {
				QCD_reco_r9_p2_bar->Fill(reco_pho2_r9_BNew,qcd_weight);
				QCD_reco_hoe_p2_bar->Fill(reco_pho2_hoe_BNew,qcd_weight);
				QCD_reco_chiso_p2_bar->Fill(reco_pho2_chiso_BNew,qcd_weight);
				QCD_reco_sieie_p2_bar->Fill(reco_pho2_sieie_BNew,qcd_weight);
			}

			//endcap
			if(abs(reco_eta_p1_BNew) > 1.57 && abs(reco_eta_p1_BNew) < 2.5) {
				QCD_reco_r9_p1_end->Fill(reco_pho1_r9_BNew,qcd_weight);
				QCD_reco_hoe_p1_end->Fill(reco_pho1_hoe_BNew,qcd_weight);
				QCD_reco_chiso_p1_end->Fill(reco_pho1_chiso_BNew,qcd_weight);
				QCD_reco_sieie_p1_end->Fill(reco_pho1_sieie_BNew,qcd_weight);
			}

			if(abs(reco_eta_p2_BNew) > 1.57 && abs(reco_eta_p2_BNew) < 2.5) {
				QCD_reco_r9_p2_end->Fill(reco_pho2_r9_BNew,qcd_weight);
				QCD_reco_hoe_p2_end->Fill(reco_pho2_hoe_BNew,qcd_weight);
				QCD_reco_chiso_p2_end->Fill(reco_pho2_chiso_BNew,qcd_weight);
				QCD_reco_sieie_p2_end->Fill(reco_pho2_sieie_BNew,qcd_weight);
			}

			//

			true_background++;
		}


		//fake background
		if(BDTG_BNew < BDT_cut) {
			QCD_reco_deltaeta_p1_j1_plot_fake->Fill(reco_deltaeta_p1j1_BNew , qcd_weight);
			QCD_reco_invmass_j1_j2_plot_fake->Fill(reco_invmass_j1j2_BNew , qcd_weight);
//			if( evt < 10) cout << "Q_T "<<reco_deltaeta_p1j1_BNew <<endl;
			fake_background++;
		}


	}







	/*
		background_event result =   { 0, 0, TreeBNew->GetEntries() };
		return result ;
		*/
	f2->Close();
	return 0 ;
}



int inter( int &events , int &true_background , int &fake_background , double &BDT_cut ) {
	_mutex.lock();
	string file02 = "new_BDT_inter.root" ;
	cout << file02 <<endl;
	TFile *f2 = new TFile(file02.c_str());
	TTree *TreeBNew;
	f2->GetObject("TreeBNew",TreeBNew);
	events = TreeBNew->GetEntries() ;


	float reco_invmass_j1j2_BNew;
	float reco_deltaeta_j1j2_BNew;
	float reco_deltaphi_j1j2_BNew;
	float reco_deltaR_j1j2_BNew;
	float reco_deltapt_j1j2_BNew;

	float reco_invmass_p1p2_BNew;
	float reco_deltaeta_p1p2_BNew;
	float reco_deltaphi_p1p2_BNew;
	float reco_deltaR_p1p2_BNew;
	float reco_deltapt_p1p2_BNew;

	float reco_deltaeta_p1j1_BNew;
	float reco_deltaeta_p1j2_BNew;
	float reco_deltaeta_p2j1_BNew;
	float reco_deltaeta_p2j2_BNew;

	float reco_deltaphi_p1j1_BNew;
	float reco_deltaphi_p1j2_BNew;
	float reco_deltaphi_p2j1_BNew;
	float reco_deltaphi_p2j2_BNew;

	float reco_deltaR_p1j1_BNew;
	float reco_deltaR_p1j2_BNew;
	float reco_deltaR_p2j1_BNew;
	float reco_deltaR_p2j2_BNew;

	float reco_deltapt_p1j1_BNew;
	float reco_deltapt_p1j2_BNew;
	float reco_deltapt_p2j1_BNew;
	float reco_deltapt_p2j2_BNew;
	float BDT_BNew;
	float BDTG_BNew;

	float reco_pho1_r9_BNew , reco_pho1_hoe_BNew , reco_pho1_chiso_BNew , reco_pho1_sieie_BNew , reco_pho1_mva_BNew ;
	float reco_pho2_r9_BNew , reco_pho2_hoe_BNew , reco_pho2_chiso_BNew , reco_pho2_sieie_BNew , reco_pho2_mva_BNew ;
	float reco_eta_p1_BNew , reco_eta_p2_BNew ;

	TreeBNew->SetBranchAddress("reco_invmass_j1j2", &reco_invmass_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_j1j2", &reco_deltaeta_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_j1j2", &reco_deltaphi_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_j1j2", &reco_deltaR_j1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_j1j2", &reco_deltapt_j1j2_BNew);

	TreeBNew->SetBranchAddress("reco_invmass_p1p2", &reco_invmass_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p1p2", &reco_deltaeta_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p1p2", &reco_deltaphi_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p1p2", &reco_deltaR_p1p2_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p1p2", &reco_deltapt_p1p2_BNew);

	TreeBNew->SetBranchAddress("reco_deltaeta_p1j1", &reco_deltaeta_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p1j2", &reco_deltaeta_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p2j1", &reco_deltaeta_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaeta_p2j2", &reco_deltaeta_p2j2_BNew);

	TreeBNew->SetBranchAddress("reco_deltaphi_p1j1", &reco_deltaphi_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p1j2", &reco_deltaphi_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p2j1", &reco_deltaphi_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaphi_p2j2", &reco_deltaphi_p2j2_BNew);

	TreeBNew->SetBranchAddress("reco_deltaR_p1j1", &reco_deltaR_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p1j2", &reco_deltaR_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p2j1", &reco_deltaR_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltaR_p2j2", &reco_deltaR_p2j2_BNew);

	TreeBNew->SetBranchAddress("reco_deltapt_p1j1", &reco_deltapt_p1j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p1j2", &reco_deltapt_p1j2_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p2j1", &reco_deltapt_p2j1_BNew);
	TreeBNew->SetBranchAddress("reco_deltapt_p2j2", &reco_deltapt_p2j2_BNew);
	TreeBNew->SetBranchAddress("BDT", &BDT_BNew);
	TreeBNew->SetBranchAddress("BDTG", &BDTG_BNew);

	TreeBNew->SetBranchAddress("reco_pho1_r9", &reco_pho1_r9_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_hoe", &reco_pho1_hoe_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_chiso", &reco_pho1_chiso_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_sieie", &reco_pho1_sieie_BNew);
	TreeBNew->SetBranchAddress("reco_pho1_mva", &reco_pho1_mva_BNew);
	TreeBNew->SetBranchAddress("reco_eta_p1", &reco_eta_p1_BNew);

	TreeBNew->SetBranchAddress("reco_pho2_r9", &reco_pho2_r9_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_hoe", &reco_pho2_hoe_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_chiso", &reco_pho2_chiso_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_sieie", &reco_pho2_sieie_BNew);
	TreeBNew->SetBranchAddress("reco_pho2_mva", &reco_pho2_mva_BNew);
	TreeBNew->SetBranchAddress("reco_eta_p2", &reco_eta_p2_BNew);

	_mutex.unlock();

	for(int evt=0 ; evt < events ; evt++) {

		TreeBNew->GetEntry(evt);



		INTER_BDT->Fill(BDT_BNew,inter_weight);
		INTER_BDTG->Fill(BDTG_BNew,inter_weight);
		//true background
		if(BDTG_BNew >= BDT_cut && reco_invmass_j1j2_BNew > 400 && reco_invmass_p1p2_BNew > 10 ) {
			INTER_reco_deltaeta_p1_j1_plot_true->Fill(reco_deltaeta_p1j1_BNew,inter_weight);
			INTER_reco_invmass_j1_j2_plot_true->Fill(reco_invmass_j1j2_BNew,inter_weight);
			INTER_reco_deltaphi_j1_j2_plot_true->Fill(reco_deltaphi_j1j2_BNew,inter_weight);
			INTER_reco_deltaphi_p1_p2_plot_true->Fill(reco_deltaphi_p1p2_BNew,inter_weight);
			INTER_reco_deltaeta_j1_j2_plot_true->Fill(reco_deltaeta_j1j2_BNew,inter_weight);
			INTER_reco_invmass_p1_p2_plot_true->Fill(reco_invmass_p1p2_BNew,inter_weight);
			INTER_reco_mva_p1->Fill(reco_pho1_mva_BNew,inter_weight);
			INTER_reco_mva_p2->Fill(reco_pho2_mva_BNew,inter_weight);

			//barral
			if(abs(reco_eta_p1_BNew) > 0 && abs(reco_eta_p1_BNew) < 1.44) {
				INTER_reco_r9_p1_bar->Fill(reco_pho1_r9_BNew,inter_weight);
				INTER_reco_hoe_p1_bar->Fill(reco_pho1_hoe_BNew,inter_weight);
				INTER_reco_chiso_p1_bar->Fill(reco_pho1_chiso_BNew,inter_weight);
				INTER_reco_sieie_p1_bar->Fill(reco_pho1_sieie_BNew,inter_weight);
			}

			if(abs(reco_eta_p2_BNew) > 0 && abs(reco_eta_p2_BNew) < 1.44) {
				INTER_reco_r9_p2_bar->Fill(reco_pho2_r9_BNew,inter_weight);
				INTER_reco_hoe_p2_bar->Fill(reco_pho2_hoe_BNew,inter_weight);
				INTER_reco_chiso_p2_bar->Fill(reco_pho2_chiso_BNew,inter_weight);
				INTER_reco_sieie_p2_bar->Fill(reco_pho2_sieie_BNew,inter_weight);
			}

			//endcap
			if(abs(reco_eta_p1_BNew) > 1.57 && abs(reco_eta_p1_BNew) < 2.5) {
				INTER_reco_r9_p1_end->Fill(reco_pho1_r9_BNew,inter_weight);
				INTER_reco_hoe_p1_end->Fill(reco_pho1_hoe_BNew,inter_weight);
				INTER_reco_chiso_p1_end->Fill(reco_pho1_chiso_BNew,inter_weight);
				INTER_reco_sieie_p1_end->Fill(reco_pho1_sieie_BNew,inter_weight);
			}

			if(abs(reco_eta_p2_BNew) > 1.57 && abs(reco_eta_p2_BNew) < 2.5) {
				INTER_reco_r9_p2_end->Fill(reco_pho2_r9_BNew,inter_weight);
				INTER_reco_hoe_p2_end->Fill(reco_pho2_hoe_BNew,inter_weight);
				INTER_reco_chiso_p2_end->Fill(reco_pho2_chiso_BNew,inter_weight);
				INTER_reco_sieie_p2_end->Fill(reco_pho2_sieie_BNew,inter_weight);
			}

			//





			true_background++;
		}


		//fake background
		if(BDTG_BNew < BDT_cut) {
			INTER_reco_deltaeta_p1_j1_plot_fake->Fill(reco_deltaeta_p1j1_BNew , inter_weight);
			INTER_reco_invmass_j1_j2_plot_fake->Fill(reco_invmass_j1j2_BNew , inter_weight);
//			if( evt < 10) cout << "Q_T "<<reco_deltaeta_p1j1_BNew <<endl;
			fake_background++;
		}


	}







	/*
		background_event result =   { 0, 0, TreeBNew->GetEntries() };
		return result ;
		*/
	f2->Close();

	return 0 ;
}



void data(int &events) {

	_mutex.lock();
	string file = "data_101.root" ;
	cout << file <<endl;
	TFile *f2 = new TFile(file.c_str());
	TTree *DataTree;
	f2->GetObject("ggNtuplizer/EventTree;41",DataTree);


	int nPho ;
	int nJet ;
	int lumis;
Float_t VBS_phopt[100]= { 0 }, VBS_phoeta[100]= { 0 }, VBS_phophi[100]= { 0 }, VBS_phom[100]= { 0 };
Float_t VBS_r9[100]= { 0 }, VBS_hoe[100]= { 0 }, VBS_phiso[100]= { 0 }, VBS_chiso[100]= { 0 }, VBS_sieie[100]= { 0 } , VBS_Photon_mvaID[100]= { 0 } ;
Float_t VBS_jetpt[100]= { 0 }, VBS_jeteta[100]= { 0 }, VBS_jetphi[100]= { 0 }, VBS_jetm[100]= { 0 };

Bool_t VBS_HLT_Diphoton_14 , VBS_HLT_Diphoton_15 , VBS_HLT_Diphoton_16 , VBS_HLT_Diphoton_17 ;
	DataTree->SetBranchAddress("Photon_pt",   &VBS_phopt  );
	DataTree->SetBranchAddress("Photon_eta",   &VBS_phoeta  );
	DataTree->SetBranchAddress("Photon_phi",   &VBS_phophi  );
	DataTree->SetBranchAddress("Photon_mass",   &VBS_phom  );

	DataTree->SetBranchAddress("Photon_r9",   &VBS_r9  );
	DataTree->SetBranchAddress("Photon_hoe",   &VBS_hoe  );
	DataTree->SetBranchAddress("Photon_pfRelIso03_all",   &VBS_phiso  );
	DataTree->SetBranchAddress("Photon_pfRelIso03_chg",   &VBS_chiso  );
	DataTree->SetBranchAddress("Photon_sieie",   &VBS_sieie  );
	DataTree->SetBranchAddress("Photon_mvaID",   &VBS_Photon_mvaID  );
	
	DataTree->SetBranchAddress("Jet_pt",   &VBS_jetpt  );
	DataTree->SetBranchAddress("Jet_eta",   &VBS_jeteta  );
	DataTree->SetBranchAddress("Jet_phi",   &VBS_jetphi  );
	DataTree->SetBranchAddress("Jet_mass",   &VBS_jetm  );
	
	DataTree->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",   &VBS_HLT_Diphoton_14  );
	DataTree->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55",   &VBS_HLT_Diphoton_17  );
	


	events = DataTree->GetEntries() ;
	cout << "      +++++++++++++++            " <<endl;
	cout << DataTree->GetEntries() <<endl;
	cout << "      +++++++++++++++            " <<endl;
	int events_loop = 20 ;
	TLorentzVector VBS_photon1,VBS_photon2, VBS_jet1, VBS_jet2;
	_mutex.unlock();

	for( int evt=0 ; evt < DataTree->GetEntries()  ; evt++ ) {
		VBS_photon1.SetPtEtaPhiE(0,0,0,0) ,  VBS_photon2.SetPtEtaPhiE(0,0,0,0)  , VBS_jet1.SetPtEtaPhiM(0,0,0,0)   , VBS_jet2.SetPtEtaPhiM(0,0,0,0)  ;
		int pho_size = 0 ;
		int jet_size = 0 ;
		DataTree->GetEntry(evt);

		float phoeta[10] = {} , phophi[10] = {} ,  phoe[10] = {} ,  phoet[10] = {} ;
		float jeteta[10] = {} , jetphi[10] = {} ,  jetpt[10] = {} ,  jetmt[10] = {} ;

		float deltapt_j1j2 = 0 , deltaphi_j1j2 = 0 , deltaeta_j1j2 = 0 , mjj = 0;
		float deltapt_p1p2 = 0 , deltaphi_p1p2 = 0 , deltaeta_p1p2 = 0 , mpp = 0;

		float deltapt_p1j1 = 0 , deltapt_p1j2 = 0 ,  deltapt_p2j1 = 0 ,  deltapt_p2j2 = 0 ;
		float deltaeta_p1j1 = 0 , deltaeta_p1j2 = 0 ,  deltaeta_p2j1 = 0 ,  deltaeta_p2j2 = 0 ;
		float deltaphi_p1j1 = 0 , deltaphi_p1j2 = 0 ,  deltaphi_p2j1 = 0 ,  deltaphi_p2j2 = 0 ;

		float phoR9_p1 = 0 , phoR9_p2 = 0 , phoHoverE_p1 = 0 , phoHoverE_p2 = 0 , phoIDMVA_p1 = 0 , phoIDMVA_p2 = 0 ;
		float phoPFChIso_p1 = 0 , phoPFChIso_p2 = 0 , phoSigmaIEtaIEta_p1 = 0 , phoSigmaIEtaIEta_p2 = 0 ;

		float eta_j1 = 0 , eta_j2 = 0 , eta_p1 = 0 , eta_p2 = 0 , pt_j1 = 0 , pt_j2 = 0 , pt_p1 = 0 , pt_p2 = 0 ;
		int chr = 0 ;

		if ( nPho >=2 && nJet >= 2 ) {

		

	


			if ( abs(pt_p1) > 40 && abs(pt_p2) > 20  &&  abs(pt_j1) > 100  && abs(pt_j2) > 30  &&
			        abs(eta_p1) < 3.5  && abs(eta_p2) <  3.5 &&  abs(eta_j1)  <  5  && abs(eta_j2) <  5 ) {
				chr = 1 ;
			}





			deltapt_p1j1 = ( VBS_photon1 - VBS_jet1 ).Pt()  , deltapt_p1j2 = ( VBS_photon1 - VBS_jet2 ).Pt() , deltapt_p2j1 = ( VBS_photon2 - VBS_jet1 ).Pt() , deltapt_p2j2 = ( VBS_photon2 - VBS_jet2 ).Pt() ;
			deltaeta_p1j1 = eta_p1-eta_p2  , deltaeta_p1j2 = ( VBS_photon1 - VBS_jet2 ).Eta() , deltaeta_p2j1 = ( VBS_photon2 - VBS_jet1 ).Eta() , deltaeta_p2j2 = ( VBS_photon2 - VBS_jet2 ).Eta() ;
			deltaphi_p1j1 = ( VBS_photon1 - VBS_jet1 ).Phi()  , deltaphi_p1j2 = ( VBS_photon1 - VBS_jet2 ).Phi() , deltaphi_p2j1 = ( VBS_photon2 - VBS_jet1 ).Phi() , deltaphi_p2j2 = ( VBS_photon2 - VBS_jet2 ).Phi() ;
			deltapt_p1p2 = ( VBS_photon1 - VBS_photon2 ).Pt() , deltaphi_p1p2 = VBS_photon1.DeltaPhi(VBS_photon2) , deltaeta_p1p2 = ( VBS_photon1 - VBS_photon2 ).Eta() , mpp = ( VBS_photon1 + VBS_photon2 ).M();
			deltapt_j1j2 = ( VBS_jet1 - VBS_jet2 ).Pt() , deltaphi_j1j2 = VBS_jet1.DeltaPhi(VBS_jet2) , deltaeta_j1j2 = eta_j1 - eta_j2 , mjj = ( VBS_jet1 + VBS_jet2 ).M();
}
		
			//select cut chr ==1



			//barral
			if(abs(eta_p1) > 0 && abs(eta_p1) < 1.44) {
				Data_reco_r9_p1_bar->Fill(phoR9_p1,1);
				Data_reco_hoe_p1_bar->Fill(phoHoverE_p1,1);
				Data_reco_chiso_p1_bar->Fill(phoPFChIso_p1,1);
				Data_reco_sieie_p1_bar->Fill(phoSigmaIEtaIEta_p1,1);
			}

			if(abs(eta_p2) > 0 && abs(eta_p2) < 1.44) {
				Data_reco_r9_p2_bar->Fill(phoR9_p2,1);
				Data_reco_hoe_p2_bar->Fill(phoHoverE_p2,1);
				Data_reco_chiso_p2_bar->Fill(phoPFChIso_p2,1);
				Data_reco_sieie_p2_bar->Fill(phoSigmaIEtaIEta_p2,1);
			}

			//endcap
			if(abs(eta_p1) > 1.57 && abs(eta_p1) < 2.5) {
				Data_reco_r9_p1_end->Fill(phoR9_p1,1);
				Data_reco_hoe_p1_end->Fill(phoHoverE_p1,1);
				Data_reco_chiso_p1_end->Fill(phoPFChIso_p1,1);
				Data_reco_sieie_p1_end->Fill(phoSigmaIEtaIEta_p1,1);
			}

			if(abs(eta_p2) > 1.57 && abs(eta_p2) < 2.5) {
				Data_reco_r9_p2_end->Fill(phoR9_p2,1);
				Data_reco_hoe_p2_end->Fill(phoHoverE_p2,1);
				Data_reco_chiso_p2_end->Fill(phoPFChIso_p2,1);
				Data_reco_sieie_p2_end->Fill(phoSigmaIEtaIEta_p2,1);
			}




			//	}





		}

	}


}





void read_1() {

////initial parameter//////////////



	TStopwatch t;
	t.Start();


	int x0 = 0 , x1 = 0 ,x2 = 0;
	int y0 = 0 , y1 = 0 , y2 = 0;
	int z0 = 0 , z1 = 0 , z2 = 0 ;
	int d0;
	double BDT__ ,BDT_cut = 0.12;
	signal_event signal_result ;
	background_event background_result;

	vbs_weight = 59.9/(1000000/(vbs_total_cross*1000));
	qcd_weight = 59.9/(1000000/(qcd_total_cross*1000));
	inter_weight = 59.9/(1000000/(inter_total_cross*1000));
	cout << 59.9/(1000000/(vbs_total_cross*1000)) << " " <<  59.9/(1000000/(qcd_total_cross*1000)) << " " << 59.9/(1000000/(inter_total_cross*1000)) <<endl;
	cout << "---------------"  << endl ;


	float vbs_true_numberofevent = (alpha*vbs_total_cross);
	float vbs_fake_numberofevent = (alpha*vbs_total_cross);
	float qcd_true_numberofevent = (alpha*qcd_total_cross);
	float qcd_fake_numberofevent = (alpha*qcd_total_cross);
	float inter_true_numberofevent = (alpha*inter_total_cross);
	float inter_fake_numberofevent = (alpha*inter_total_cross);
	double  BDT_ = 0  ;
	double  BDT_max = 0  ;
	float sig_max =0 ,  max_ = 0 , Z = 0  , qcd_true = 0 ,vbs_true = 0;
	clock_t tStart = clock();



	//for(int i = 0 ; i < 99; i++ ) {
	////////////////////////////////////
	x0 = 0 , x1 = 0 , x2 = 0;
	y0 = 0 , y1 = 0 , y2 = 0;
	z0 = 0 , z1 = 0 , z2 = 0 ;
	d0 = 0 ;
	int data_events = 0 ;

	/////////main function ///////   // { function , total , true , fake ,cut}



	int thre_detect = 1 ;


	if(thre_detect ==1) {
		std::thread t1 ([&]() {
			vbs( x0,y0 ,z0 , BDT_cut ) ;

		});

		std::thread t2 ([&]() {
			qcd( x1,y1 ,z1 , BDT_cut ) ;
		});

		std::thread t3 ([&]() {
			inter( x2 , y2 , z2 , BDT_cut ) ;
		});
		
	//	std::thread t4 ([&]() {
	//		data(d0);
	//	});

		t1.join();
		t2.join();
		t3.join();
		//t4.join();
	}


	if(thre_detect ==2) {

		vbs( x0,y0 ,z0 , BDT_cut ) ;
		qcd( x1,y1 ,z1 , BDT_cut ) ;
		inter( x2 , y2 , z2 , BDT_cut ) ;
		//	data(d0);
	}




	/////////////////////

	float  f_inter_total = x2 , f_vbs_total=x0 , f_qcd_total=x1 ;
	float  f_inter_true = y2 , f_vbs_true=y0 , f_qcd_true=y1 ;
	/*
		cout << "inter_total_entry= " << x2 << "  ||vbs_total_entry= " << x0 << "  ||qcd_total_entry= " << x1 <<endl;
		cout << "inter_total_events= " << x2*inter_weight << "  ||vbs_total_entry= " << x0*vbs_weight << "  ||qcd_total_entry= " << x1*qcd_weight <<endl;
		cout << "----------------  after BDT cut  --------------  "   << endl;
		cout << "inter_total_entry= " << y2 << "  ||vbs_total_entry= " << y0 << "  ||qcd_total_entry= " << y1 <<endl;
		cout << "inter_total_events= " << y2*inter_weight << "  ||vbs_total_entry= " << y0*vbs_weight << "  ||qcd_total_entry= " << y1*qcd_weight <<endl;
		cout << "inter_total_ratio= " << f_inter_true/f_inter_total << "  ||vbs_total_ratio= " << f_vbs_true/f_vbs_total << "  ||qcd_total_ratio= " << f_qcd_true/f_qcd_total <<endl;
		cout << "----------------  after weight  --------------  "   << endl;
		cout << " <><><><><><><><><><><><><><><><><><>"<<endl;

		cout << "inter= "<< inter_weight << " || vbs= " << vbs_weight << " || qcd= " << qcd_weight <<endl;
		*/
	float significant_weight = (vbs_weight*y0)/sqrt(qcd_weight*y1 + vbs_weight*y0 + inter_weight*y2 );

	//cout << "inter_total= "<< x2*inter_weight << " || vbs_total= " << x0*vbs_weight << " || qcd_total= " << x1*qcd_weight <<endl;
	//cout << " --------------   significant_weight = " <<  significant_weight << endl;

	if(significant_weight > max_ && qcd_weight > 100 ) {
		max_ = significant_weight;
		BDT_ = BDT_cut;

	}

	trueratio_plt->Fill(BDT_cut ,y0/x0 );
	fakeratio_plt->Fill(BDT_cut ,y1/x1 );
	significant_plt->Fill(BDT_cut ,significant_weight );


	BDT_cut = BDT_cut + 0.02;




	///////////SIG////////////
	cout << "-----------start cal SIG--------" <<endl;
	float vbs=vbs_weight*y0 , inter=inter_weight*y2  , qcd=qcd_weight*y1 ;
	Z = sqrt( 2*( ( vbs+inter+qcd )*log( (vbs+inter+qcd)/(inter+qcd) ) + (inter+qcd) - (vbs+inter+qcd) ) );
	cout << "qcd=  " << y1 << endl;
	significant_special_plt->Fill(BDT_cut,Z);
	if(Z > sig_max && y1 > 1000 ) {
		sig_max = Z ;
		BDT_max = BDT_cut ;
	}
	cout << "Z= " << Z <<"  sig_max "<<  sig_max << "  max_B= " << BDT_max <<endl;


	//////////////////////

	//}//end of calculate SIG



	cout << "---------------" << endl ;
	cout << "BDT_max= "<< BDT_max << " ||  sig_max= "<< sig_max << endl ;


///////////////fit///////////////
//"exp( [0] + [1]*(x) + [2]*(x*x) + [3]*(x*x*x) )";

///////invjj//////
	cout << "inter_total_entry= " << x2 << "  ||vbs_total_entry= " << x0 << "  ||qcd_total_entry= " << x1 << "  ||data_total_entry= " << d0 << endl;
	TF1 *func_vbs_inv = new TF1("func_vbs_inv", "exp( [0] + [1]*(x) + [2]*(x*x) + [3]*(x*x*x) )",500,5000);
	func_vbs_inv->SetParameters(-12.5 , 11e-4 , -6e-7 , 4e-11 );
	func_vbs_inv->SetLineColor(3);

	TF1 *func_qcd_inv = new TF1("func_qcd_inv", "exp( [0] + [1]*(x) + [2]*(x*x) + [3]*(x*x*x) )",500,5000);
	func_qcd_inv->SetParameters(-9.87 , -9.7e-4 , -9.8e-8 , 1e-11 );
	func_qcd_inv->SetLineColor(4);

	TF1 *func_inter_inv = new TF1("func_inter_inv", "exp( [0] + [1]*(x) + [2]*(x*x) + [3]*(x*x*x) )",500,5000);
	func_inter_inv->SetParameters(1 , 10e-4 , -7e-7 , 3e-11);
	func_inter_inv->SetLineColor(2);

///////deltaetajj//////
	TF1 *func_vbs_deltaetajj = new TF1("func_vbs_deltaetajj", "exp( [0] + [1]*(x) + [2]*(x*x)  )",0,8);
	func_vbs_deltaetajj->SetParameters(-3.17 , -0.7 , -19.4 );
	func_vbs_deltaetajj->SetLineColor(3);

	TF1 *func_qcd_deltaetajj = new TF1("func_qcd_deltaetajj", "exp( [0] + [1]*(x) + [2]*(x*x)  )",0,8);
	func_qcd_deltaetajj->SetParameters(-2.55 , 3.4 , -15.95  );
	func_qcd_deltaetajj->SetLineColor(4);

	TF1 *func_inter_deltaetajj = new TF1("func_inter_deltaetajj", "exp( [0] + [1]*(x) + [2]*(x*x)  )",0,8);
	func_inter_deltaetajj->SetParameters(-2.98 , -0.5 , -12.8 );
	func_inter_deltaetajj->SetLineColor(2);

//////////////////end of fit///////////////

//////////////////////////////e

	/////Tcanvas////////
	//BDT and BDTG/ c000~
	{
		TCanvas *c001= new TCanvas("c001","BDT");

		QCD_BDT->SetStats(0);
		QCD_BDT->SetLineColor(2);
		QCD_BDT->SetLineWidth(3);
		//	QCD_BDT->Scale(1/QCD_BDT->Integral());
		//	QCD_BDT->GetYaxis()->SetRangeUser(0,0.15);
		QCD_BDT->Draw("hist");

		VBS_BDT->SetStats(0);
		VBS_BDT->SetLineColor(8);
		VBS_BDT->SetLineWidth(3);
		//	VBS_BDT->Scale(1/VBS_BDT->Integral());
		VBS_BDT->Draw("hist,same");

		INTER_BDT->SetStats(0);
		INTER_BDT->SetLineColor(46);
		INTER_BDT->SetLineWidth(3);
		//	VBS_BDT->Scale(1/VBS_BDT->Integral());
		INTER_BDT->Draw("hist,same");




		gPad->SetLogy();
		TLegend *BDT_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
		BDT_0->AddEntry(VBS_BDT,"VBS_BDT")  ;
		BDT_0->Draw()  ;
		BDT_0->AddEntry(QCD_BDT,"QCD_BDT")  ;
		BDT_0->Draw()  ;
		BDT_0->AddEntry(INTER_BDT,"INTER_BDT")  ;
		BDT_0->Draw()  ;

		c001->SaveAs("!!20BDT.pdf");
////////////////////////////////////////////////////////////
		TCanvas *c002= new TCanvas("c002","BDTG");

		QCD_BDTG->SetStats(0);
		QCD_BDTG->SetLineColor(4);
		QCD_BDTG->SetLineWidth(3);
		//	QCD_BDTG->Scale(1/QCD_BDTG->Integral());
		QCD_BDTG->GetXaxis()->SetRangeUser(-1.1,1.1);
		QCD_BDTG->GetYaxis()->SetRangeUser(1,1000000000);
		QCD_BDTG->Draw("hist");

		VBS_BDTG->SetStats(0);
		VBS_BDTG->SetLineColor(8);
		VBS_BDTG->SetLineWidth(3);
		//	VBS_BDTG->Scale(1/VBS_BDTG->Integral());
		VBS_BDTG->Draw("hist,same");

		INTER_BDTG->SetStats(0);
		INTER_BDTG->SetLineColor(46);
		INTER_BDTG->SetLineWidth(3);
		//	INTER_BDTG->Scale(1/INTER_BDTG->Integral());
		INTER_BDTG->Draw("hist,same");

		gPad->SetLogy();



		TLegend *BDTG_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
		BDTG_0->AddEntry(VBS_BDTG,"VBS_BDTG")  ;
		BDTG_0->Draw()  ;
		BDTG_0->AddEntry(QCD_BDTG,"QCD_BDTG")  ;
		BDTG_0->Draw()  ;
		BDTG_0->AddEntry(INTER_BDTG,"INTER_BDTG")  ;
		BDTG_0->Draw()  ;

		c002->SaveAs("!!20BDTG.pdf");


		//BDTG stack

		TCanvas *cBTDG_stack= new TCanvas("cBTDG_stack","BTDG_stack");


		hs_BDTG->SetMinimum(0.00001);
		hs_BDTG->SetMaximum(2);
		//	gPad->SetLogy();

		QCD_BDTG->SetStats(0);
		QCD_BDTG->SetLineColor(46);
		QCD_BDTG->SetLineWidth(3);
		QCD_BDTG->SetFillColor(46);
		QCD_BDTG->SetFillStyle(4000);
		QCD_BDTG->GetYaxis()->SetRangeUser(0.001,5000000);
		QCD_BDTG->Scale(1/QCD_BDTG->Integral());
		QCD_BDTG->Scale(0.5);
		QCD_BDTG->GetXaxis()->SetRangeUser(-1.1,1.1);

		VBS_BDTG->SetStats(0);
		VBS_BDTG->SetLineColor(8);
		VBS_BDTG->SetFillColor(8);
		VBS_BDTG->SetFillStyle(4000);
		VBS_BDTG->SetLineWidth(3);
		VBS_BDTG->Scale(1/VBS_BDTG->Integral());

		INTER_BDTG->SetStats(0);
		INTER_BDTG->SetLineColor(46);
		INTER_BDTG->SetFillColor(46);
		INTER_BDTG->SetFillStyle(4000);
		INTER_BDTG->SetLineWidth(3);
		INTER_BDTG->Scale(1/INTER_BDTG->Integral());
		INTER_BDTG->Scale(0.5);


		VBS_BDTG->Draw("hist") ;
		hs_BDTG->Add(INTER_BDTG);
		hs_BDTG->Add(QCD_BDTG);
		hs_BDTG->Draw("hist,same");


		cBTDG_stack->SaveAs("!!20BDTG_stack.pdf");





//
		TCanvas *c003= new TCanvas("c003","SIGNIFICANT");

		significant_special_plt->SetStats(0);
		significant_special_plt->SetLineColor(1);
		significant_special_plt->SetLineWidth(3);
		significant_special_plt->Draw("hist");


		TLegend *significant_special_plt_0 = new TLegend(0.7, 0.6, 0.4, 0.4);
		significant_special_plt_0->AddEntry(significant_special_plt,"significant_special_plt")  ;
		significant_special_plt_0->Draw()  ;


		c003->SaveAs("!!20SIGNIFICANT.pdf");


//
		TCanvas *c004= new TCanvas("c004","ratio");

		trueratio_plt->SetStats(0);
		trueratio_plt->SetLineColor(4);
		trueratio_plt->SetLineWidth(3);
		trueratio_plt->Draw("hist");


		fakeratio_plt->SetStats(0);
		fakeratio_plt->SetLineColor(2);
		fakeratio_plt->SetLineWidth(3);
		fakeratio_plt->Draw("hist,same");


		TLegend *fakeratio_plt_0 = new TLegend(0.7, 0.6, 0.4, 0.4);
		fakeratio_plt_0->AddEntry(trueratio_plt,"S")  ;
		fakeratio_plt_0->Draw()  ;
		fakeratio_plt_0->AddEntry(fakeratio_plt,"B")  ;
		fakeratio_plt_0->Draw()  ;

		c004->SaveAs("!!20ratio.pdf");




	}
//delta eta p1j1 ,c101
	{

		TCanvas *c101= new TCanvas("c101","deltaeta");

		QCD_reco_deltaeta_p1_j1_plot_true->SetStats(0);
		QCD_reco_deltaeta_p1_j1_plot_true->SetLineColor(4);
		QCD_reco_deltaeta_p1_j1_plot_true->SetLineWidth(3);
//		QCD_reco_deltaeta_p1_j1_plot_true->Scale(qcd_true_numberofevent);
		QCD_reco_deltaeta_p1_j1_plot_true->GetYaxis()->SetRangeUser(0.0001, 10000000);
		QCD_reco_deltaeta_p1_j1_plot_true->Draw("");
		//	QCD_reco_deltaeta_p1_j1_plot_true->Fit("func_qcd_deltaetajj");

//		cout << 2*alpha*qcd_total_cross << endl;
		cout << QCD_reco_deltaeta_p1_j1_plot_true->Integral() << endl;
		/*
				QCD_reco_deltaeta_p1_j1_plot_fake->SetStats(0);
				QCD_reco_deltaeta_p1_j1_plot_fake->SetLineColor(9);
				QCD_reco_deltaeta_p1_j1_plot_fake->SetLineWidth(3);
				QCD_reco_deltaeta_p1_j1_plot_fake->Scale(qcd_fake_numberofevent);
				QCD_reco_deltaeta_p1_j1_plot_fake->Draw("hist,same");
		*/

		VBS_reco_deltaeta_p1_j1_plot_true->SetStats(0);
		VBS_reco_deltaeta_p1_j1_plot_true->SetLineColor(3);
		VBS_reco_deltaeta_p1_j1_plot_true->SetLineWidth(3);
		//	VBS_reco_deltaeta_p1_j1_plot_true->Scale(vbs_true_numberofevent);
		VBS_reco_deltaeta_p1_j1_plot_true->Draw("same");
		//	VBS_reco_deltaeta_p1_j1_plot_true->Fit("func_vbs_deltaetajj");

		/*
				VBS_reco_deltaeta_p1_j1_plot_fake->SetStats(0);
				VBS_reco_deltaeta_p1_j1_plot_fake->SetLineColor(8);
				VBS_reco_deltaeta_p1_j1_plot_fake->SetLineWidth(3);
				VBS_reco_deltaeta_p1_j1_plot_fake->Scale(vbs_fake_numberofevent);
				VBS_reco_deltaeta_p1_j1_plot_fake->Draw("hist,same");
				*/
		INTER_reco_deltaeta_p1_j1_plot_true->SetStats(0);
		INTER_reco_deltaeta_p1_j1_plot_true->SetLineColor(46);
		INTER_reco_deltaeta_p1_j1_plot_true->SetLineWidth(3);
		//	VBS_reco_deltaeta_p1_j1_plot_true->Scale(vbs_true_numberofevent);
		INTER_reco_deltaeta_p1_j1_plot_true->Draw("same");
		//	INTER_reco_deltaeta_p1_j1_plot_true->Fit("func_inter_deltaetajj");


		Data_reco_deltaeta_p1_j1_plot_true->SetStats(0);
		Data_reco_deltaeta_p1_j1_plot_true->SetLineColor(1);
		Data_reco_deltaeta_p1_j1_plot_true->SetLineWidth(3);
		//	VBS_reco_deltaeta_p1_j1_plot_true->Scale(vbs_true_numberofevent);
		Data_reco_deltaeta_p1_j1_plot_true->Draw("same");



		gPad->SetLogy();

		//QCD_reco_deltaeta_p1_j1_plot_true->Fit("func","R");

		TLegend *deltaeta_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
		deltaeta_0->AddEntry(VBS_reco_deltaeta_p1_j1_plot_true,"VBS_reco_deltaeta_p1_j1_plot_true")  ;
		deltaeta_0->Draw()  ;
//		deltaeta_0->AddEntry(VBS_reco_deltaeta_p1_j1_plot_fake,"VBS_reco_deltaeta_p1_j1_plot_fake")  ;
//		deltaeta_0->Draw()  ;
		deltaeta_0->AddEntry(QCD_reco_deltaeta_p1_j1_plot_true,"QCD_reco_deltaeta_p1_j1_plot_true")  ;
		deltaeta_0->Draw()  ;
//		deltaeta_0->AddEntry(QCD_reco_deltaeta_p1_j1_plot_fake,"QCD_reco_deltaeta_p1_j1_plot_fake")  ;
//		deltaeta_0->Draw()  ;
		deltaeta_0->AddEntry(INTER_reco_deltaeta_p1_j1_plot_true,"INTER_reco_deltaeta_p1_j1_plot_true")  ;
		deltaeta_0->Draw()  ;
		deltaeta_0->AddEntry(Data_reco_deltaeta_p1_j1_plot_true,"Data_reco_deltaeta_p1_j1_plot_true")  ;
		deltaeta_0->Draw()  ;

		c101->SaveAs("!!20deltaeta.pdf");


		////stack//////
		TCanvas *c101_stack= new TCanvas("c101_stack","deltaeta_stack");


		hs_deltap1j1->SetMinimum(0.01);
		hs_deltap1j1->SetMaximum(2000);
		gPad->SetLogy();


		QCD_reco_deltaeta_p1_j1_plot_true->SetFillColor(4);
		QCD_reco_deltaeta_p1_j1_plot_true->SetFillStyle(4000);
		hs_deltap1j1->Add(QCD_reco_deltaeta_p1_j1_plot_true);

		VBS_reco_deltaeta_p1_j1_plot_true->SetFillColor(3);
		VBS_reco_deltaeta_p1_j1_plot_true->SetFillStyle(4000);
		VBS_reco_deltaeta_p1_j1_plot_true->SetLineWidth(2);
		VBS_reco_deltaeta_p1_j1_plot_true->Scale(10);
		hs_deltap1j1->Add(VBS_reco_deltaeta_p1_j1_plot_true);


		INTER_reco_deltaeta_p1_j1_plot_true->SetFillColor(46);
		INTER_reco_deltaeta_p1_j1_plot_true->SetFillStyle(4000);
		INTER_reco_deltaeta_p1_j1_plot_true->SetLineWidth(4);
		INTER_reco_deltaeta_p1_j1_plot_true->Scale(100);
		hs_deltap1j1->Add(INTER_reco_deltaeta_p1_j1_plot_true);
		//gPad->SetFrameFillColor(1);
		//gPad->SetTheta(4);
		//gPad->SetPhi(3.14159);

		hs_deltap1j1->Draw("hist");

		Data_reco_deltaeta_p1_j1_plot_true->SetStats(0);
		Data_reco_deltaeta_p1_j1_plot_true->SetLineColor(1);
		Data_reco_deltaeta_p1_j1_plot_true->SetLineWidth(3);
		//	VBS_reco_deltaeta_p1_j1_plot_true->Scale(vbs_true_numberofevent);
		Data_reco_deltaeta_p1_j1_plot_true->Draw("same");


		TLegend *deltaeta_0_stack = new TLegend(0.7, 0.6, 0.9, 0.9);
		deltaeta_0_stack ->AddEntry(VBS_reco_deltaeta_p1_j1_plot_true,"VBS_reco_deltaeta_p1_j1_plot_true*10")  ;
		deltaeta_0_stack ->Draw()  ;

		deltaeta_0_stack ->AddEntry(QCD_reco_deltaeta_p1_j1_plot_true,"QCD_reco_deltaeta_p1_j1_plot_true")  ;
		deltaeta_0_stack ->Draw()  ;

		deltaeta_0_stack ->AddEntry(INTER_reco_deltaeta_p1_j1_plot_true,"INTER_reco_deltaeta_p1_j1_plot_true*100")  ;
		deltaeta_0_stack ->Draw()  ;
		deltaeta_0_stack->AddEntry(Data_reco_deltaeta_p1_j1_plot_true,"Data_reco_deltaeta_p1_j1_plot_true")  ;
		deltaeta_0_stack->Draw()  ;
		c101_stack->SaveAs("!!20deltaeta_stack.pdf");
	}


//delta eta j1j2 ,c104
	{
		TCanvas *c104= new TCanvas("c104","deltaetaj1j2");

		QCD_reco_deltaeta_j1_j2_plot_true->SetStats(0);
		QCD_reco_deltaeta_j1_j2_plot_true->SetLineColor(4);
		QCD_reco_deltaeta_j1_j2_plot_true->SetLineWidth(3);
		QCD_reco_deltaeta_j1_j2_plot_true->GetYaxis()->SetRangeUser(0.0001, 10000000);
		QCD_reco_deltaeta_j1_j2_plot_true->Draw("hist");

		VBS_reco_deltaeta_j1_j2_plot_true->SetStats(0);
		VBS_reco_deltaeta_j1_j2_plot_true->SetLineColor(3);
		VBS_reco_deltaeta_j1_j2_plot_true->SetLineWidth(3);
		VBS_reco_deltaeta_j1_j2_plot_true->Draw("same");

		INTER_reco_deltaeta_j1_j2_plot_true->SetStats(0);
		INTER_reco_deltaeta_j1_j2_plot_true->SetLineColor(46);
		INTER_reco_deltaeta_j1_j2_plot_true->SetLineWidth(3);
		INTER_reco_deltaeta_j1_j2_plot_true->Draw("same");

		Data_reco_deltaeta_j1_j2_plot_true->SetStats(0);
		Data_reco_deltaeta_j1_j2_plot_true->SetLineColor(1);
		Data_reco_deltaeta_j1_j2_plot_true->SetLineWidth(3);
		Data_reco_deltaeta_j1_j2_plot_true->Draw("same");

		gPad->SetLogy();

		TLegend *deltaeta_0_j1j2 = new TLegend(0.7, 0.6, 0.9, 0.9);
		deltaeta_0_j1j2->AddEntry(VBS_reco_deltaeta_j1_j2_plot_true,"VBS_reco_deltaeta_j1_j2_plot_true")  ;
		deltaeta_0_j1j2->Draw()  ;
		deltaeta_0_j1j2->AddEntry(QCD_reco_deltaeta_j1_j2_plot_true,"QCD_reco_deltaeta_j1_j2_plot_true")  ;
		deltaeta_0_j1j2->Draw()  ;
		deltaeta_0_j1j2->AddEntry(INTER_reco_deltaeta_j1_j2_plot_true,"INTER_reco_deltaeta_j1_j2_plot_true")  ;
		deltaeta_0_j1j2->Draw()  ;
		deltaeta_0_j1j2->AddEntry(Data_reco_deltaeta_j1_j2_plot_true,"Data_reco_deltaeta_j1_j2_plot_true")  ;
		deltaeta_0_j1j2->Draw()  ;

		c104->SaveAs("!!21deltaeta_j1j2.pdf");

	}


//delta phi p1p2 ,c106
	{
		TCanvas *c106= new TCanvas("c106","deltaphip1p2");

		QCD_reco_deltaphi_p1_p2_plot_true->SetStats(0);
		QCD_reco_deltaphi_p1_p2_plot_true->SetLineColor(4);
		QCD_reco_deltaphi_p1_p2_plot_true->SetLineWidth(3);
		QCD_reco_deltaphi_p1_p2_plot_true->GetYaxis()->SetRangeUser(0.0001, 10000000);
		QCD_reco_deltaphi_p1_p2_plot_true->Draw("hist");

		VBS_reco_deltaphi_p1_p2_plot_true->SetStats(0);
		VBS_reco_deltaphi_p1_p2_plot_true->SetLineColor(3);
		VBS_reco_deltaphi_p1_p2_plot_true->SetLineWidth(3);
		VBS_reco_deltaphi_p1_p2_plot_true->Draw("same");

		INTER_reco_deltaphi_p1_p2_plot_true->SetStats(0);
		INTER_reco_deltaphi_p1_p2_plot_true->SetLineColor(46);
		INTER_reco_deltaphi_p1_p2_plot_true->SetLineWidth(3);
		INTER_reco_deltaphi_p1_p2_plot_true->Draw("same");

		Data_reco_deltaphi_p1_p2_plot_true->SetStats(0);
		Data_reco_deltaphi_p1_p2_plot_true->SetLineColor(1);
		Data_reco_deltaphi_p1_p2_plot_true->SetLineWidth(3);
		Data_reco_deltaphi_p1_p2_plot_true->Draw("same");

		gPad->SetLogy();

		TLegend *deltaphi_0_p1p2 = new TLegend(0.7, 0.6, 0.9, 0.9);
		deltaphi_0_p1p2->AddEntry(VBS_reco_deltaphi_p1_p2_plot_true,"VBS_reco_deltaphi_p1_p2_plot_true")  ;
		deltaphi_0_p1p2->Draw()  ;
		deltaphi_0_p1p2->AddEntry(QCD_reco_deltaphi_p1_p2_plot_true,"QCD_reco_deltaphi_p1_p2_plot_true")  ;
		deltaphi_0_p1p2->Draw()  ;
		deltaphi_0_p1p2->AddEntry(INTER_reco_deltaphi_p1_p2_plot_true,"INTER_reco_deltaphi_p1_p2_plot_true")  ;
		deltaphi_0_p1p2->Draw()  ;
		deltaphi_0_p1p2->AddEntry(Data_reco_deltaphi_p1_p2_plot_true,"Data_reco_deltaphi_p1_p2_plot_true")  ;
		deltaphi_0_p1p2->Draw()  ;

		c106->SaveAs("!!21deltaphi_p1p2.pdf");
	}


//delta phi j1j2 ,c105
	{
		TCanvas *c105= new TCanvas("c105","deltaphij1j2");

		QCD_reco_deltaphi_j1_j2_plot_true->SetStats(0);
		QCD_reco_deltaphi_j1_j2_plot_true->SetLineColor(4);
		QCD_reco_deltaphi_j1_j2_plot_true->SetLineWidth(3);
		QCD_reco_deltaphi_j1_j2_plot_true->GetYaxis()->SetRangeUser(0.0001, 10000000);
		QCD_reco_deltaphi_j1_j2_plot_true->Draw("hist");

		VBS_reco_deltaphi_j1_j2_plot_true->SetStats(0);
		VBS_reco_deltaphi_j1_j2_plot_true->SetLineColor(3);
		VBS_reco_deltaphi_j1_j2_plot_true->SetLineWidth(3);
		VBS_reco_deltaphi_j1_j2_plot_true->Draw("same");

		INTER_reco_deltaphi_j1_j2_plot_true->SetStats(0);
		INTER_reco_deltaphi_j1_j2_plot_true->SetLineColor(46);
		INTER_reco_deltaphi_j1_j2_plot_true->SetLineWidth(3);
		INTER_reco_deltaphi_j1_j2_plot_true->Draw("same");

		Data_reco_deltaphi_j1_j2_plot_true->SetStats(0);
		Data_reco_deltaphi_j1_j2_plot_true->SetLineColor(1);
		Data_reco_deltaphi_j1_j2_plot_true->SetLineWidth(3);
		Data_reco_deltaphi_j1_j2_plot_true->Draw("same");

		gPad->SetLogy();

		TLegend *deltaphi_0_j1j2 = new TLegend(0.7, 0.6, 0.9, 0.9);
		deltaphi_0_j1j2->AddEntry(VBS_reco_deltaphi_j1_j2_plot_true,"VBS_reco_deltaphi_j1_j2_plot_true")  ;
		deltaphi_0_j1j2->Draw()  ;
		deltaphi_0_j1j2->AddEntry(QCD_reco_deltaphi_j1_j2_plot_true,"QCD_reco_deltaphi_j1_j2_plot_true")  ;
		deltaphi_0_j1j2->Draw()  ;
		deltaphi_0_j1j2->AddEntry(INTER_reco_deltaphi_j1_j2_plot_true,"INTER_reco_deltaphi_j1_j2_plot_true")  ;
		deltaphi_0_j1j2->Draw()  ;
		deltaphi_0_j1j2->AddEntry(Data_reco_deltaphi_j1_j2_plot_true,"Data_reco_deltaphi_j1_j2_plot_true")  ;
		deltaphi_0_j1j2->Draw()  ;

		c105->SaveAs("!!21deltaphi_j1j2.pdf");
	}

//invmass p1p2 , c107
	{
		TCanvas *c107= new TCanvas("c107","invmass mpp");

		QCD_reco_invmass_p1_p2_plot_true->SetStats(0);
		QCD_reco_invmass_p1_p2_plot_true->SetLineColor(4);
		QCD_reco_invmass_p1_p2_plot_true->SetLineWidth(3);
		QCD_reco_invmass_p1_p2_plot_true->GetYaxis()->SetRangeUser(0.0001, 100000000);
		QCD_reco_invmass_p1_p2_plot_true->Draw("hist");

		VBS_reco_invmass_p1_p2_plot_true->SetStats(0);
		VBS_reco_invmass_p1_p2_plot_true->SetLineColor(3);
		VBS_reco_invmass_p1_p2_plot_true->SetLineWidth(3);
		VBS_reco_invmass_p1_p2_plot_true->Draw("hist,same");

		INTER_reco_invmass_p1_p2_plot_true->SetStats(0);
		INTER_reco_invmass_p1_p2_plot_true->SetLineColor(46);
		INTER_reco_invmass_p1_p2_plot_true->SetLineWidth(3);
		INTER_reco_invmass_p1_p2_plot_true->Draw("hist,same");

		Data_reco_invmass_p1_p2_plot_true->SetStats(0);
		Data_reco_invmass_p1_p2_plot_true->SetLineColor(1);
		Data_reco_invmass_p1_p2_plot_true->SetLineWidth(3);
		Data_reco_invmass_p1_p2_plot_true->Draw("hist,same");
		gPad->SetLogy();

		TLegend *invmasspp_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
		invmasspp_0->AddEntry(VBS_reco_invmass_p1_p2_plot_true,"VBS_reco_invmass_p1_p2_plot_true")  ;
		invmasspp_0->Draw()  ;
		invmasspp_0->AddEntry(QCD_reco_invmass_p1_p2_plot_true,"QCD_reco_invmass_p1_p2_plot_true")  ;
		invmasspp_0->Draw()  ;
		invmasspp_0->AddEntry(INTER_reco_invmass_p1_p2_plot_true,"INTER_reco_invmass_p1_p2_plot_true")  ;
		invmasspp_0->Draw()  ;
		invmasspp_0->AddEntry(Data_reco_invmass_p1_p2_plot_true,"Data_reco_invmass_p1_p2_plot_true")  ;
		invmasspp_0->Draw()  ;

		c107->SaveAs("!!21invmasspp.pdf");
		
		////stack//////
		TCanvas *c107_stack= new TCanvas("c107_stack","invmassp1p2_stack");
		hs_invmass_p1_p2->SetMinimum(0.1);
		hs_invmass_p1_p2->SetMaximum(1000000);
		gPad->SetLogy();

		QCD_reco_invmass_p1_p2_plot_true->SetFillColor(4);
		QCD_reco_invmass_p1_p2_plot_true->SetFillStyle(4000);
		hs_invmass_p1_p2->Add(QCD_reco_invmass_p1_p2_plot_true);

		VBS_reco_invmass_p1_p2_plot_true->SetFillColor(3);
		VBS_reco_invmass_p1_p2_plot_true->SetFillStyle(4000);
		VBS_reco_invmass_p1_p2_plot_true->SetLineWidth(2);
		VBS_reco_invmass_p1_p2_plot_true->Scale(1);
		hs_invmass_p1_p2->Add(VBS_reco_invmass_p1_p2_plot_true);

		INTER_reco_invmass_p1_p2_plot_true->SetFillColor(46);
		INTER_reco_invmass_p1_p2_plot_true->SetFillStyle(4000);
		INTER_reco_invmass_p1_p2_plot_true->SetLineWidth(4);
		INTER_reco_invmass_p1_p2_plot_true->Scale(1);
		hs_invmass_p1_p2->Add(INTER_reco_invmass_p1_p2_plot_true);
		hs_invmass_p1_p2->Draw("hist");

		Data_reco_invmass_p1_p2_plot_true->SetStats(0);
		Data_reco_invmass_p1_p2_plot_true->SetLineColor(1);
		Data_reco_invmass_p1_p2_plot_true->SetLineWidth(3);
		Data_reco_invmass_p1_p2_plot_true->Draw("same");

		TLegend *invmass_p1_p2_0_stack = new TLegend(0.7, 0.6, 0.9, 0.9);
		invmass_p1_p2_0_stack->AddEntry(VBS_reco_invmass_p1_p2_plot_true,"VBS_reco_invmass_p1_p2*10")  ;
		invmass_p1_p2_0_stack->Draw()  ;
		invmass_p1_p2_0_stack->AddEntry(QCD_reco_invmass_p1_p2_plot_true,"QCD_reco_invmass_p1_p2")  ;
		invmass_p1_p2_0_stack->Draw()  ;
		invmass_p1_p2_0_stack->AddEntry(INTER_reco_invmass_p1_p2_plot_true,"INTER_reco_invmass_p1_p2*100")  ;
		invmass_p1_p2_0_stack->Draw()  ;
		invmass_p1_p2_0_stack->AddEntry(Data_reco_invmass_p1_p2_plot_true,"Data_reco_invmass_p1_p2")  ;
		invmass_p1_p2_0_stack->Draw()  ;
		c107_stack->SaveAs("!!20invmass_p1_p2_stack.pdf");
		
		
	}


//invmass j1j2 ,c102 c103

	{
		TCanvas *c102= new TCanvas("c102","invmass");


		QCD_reco_invmass_j1_j2_plot_true->SetLineColor(4);
		QCD_reco_invmass_j1_j2_plot_true->SetLineWidth(3);
		QCD_reco_invmass_j1_j2_plot_true->GetYaxis()->SetRangeUser(0.0001, 100000000);
		QCD_reco_invmass_j1_j2_plot_true->Draw("");

		VBS_reco_invmass_j1_j2_plot_true->SetStats(0);
		VBS_reco_invmass_j1_j2_plot_true->SetLineColor(3);
		VBS_reco_invmass_j1_j2_plot_true->SetLineWidth(3);
		VBS_reco_invmass_j1_j2_plot_true->Draw("same");

		INTER_reco_invmass_j1_j2_plot_true->SetStats(0);
		INTER_reco_invmass_j1_j2_plot_true->SetLineColor(46);
		INTER_reco_invmass_j1_j2_plot_true->SetLineWidth(3);
		INTER_reco_invmass_j1_j2_plot_true->Draw("same");

		Data_reco_invmass_j1_j2_plot_true->SetStats(0);
		Data_reco_invmass_j1_j2_plot_true->SetLineColor(1);
		Data_reco_invmass_j1_j2_plot_true->SetLineWidth(3);
		Data_reco_invmass_j1_j2_plot_true->Draw("same");

		gPad->SetLogy();

		TLegend *invmassjj_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
		invmassjj_0->AddEntry(VBS_reco_invmass_j1_j2_plot_true,"VBS_reco_invmass_j1_j2_plot_true")  ;
		invmassjj_0->Draw()  ;
		invmassjj_0->AddEntry(QCD_reco_invmass_j1_j2_plot_true,"QCD_reco_invmass_j1_j2_plot_true")  ;
		invmassjj_0->Draw()  ;
		invmassjj_0->AddEntry(INTER_reco_invmass_j1_j2_plot_true,"INTER_reco_invmass_j1_j2_plot_true")  ;
		invmassjj_0->Draw()  ;
		invmassjj_0->AddEntry(Data_reco_invmass_j1_j2_plot_true,"Data_reco_invmass_j1_j2_plot_true")  ;
		invmassjj_0->Draw()  ;

		c102->SaveAs("!!21invmassjj.pdf");
		
		
		////stack//////
		TCanvas *c102_stack= new TCanvas("c102_stack","invmassj1j2_stack");
		hs_invmass_j1_j2->SetMinimum(0.1);
		hs_invmass_j1_j2->SetMaximum(1000000);
		gPad->SetLogy();

		QCD_reco_invmass_j1_j2_plot_true->SetFillColor(4);
		QCD_reco_invmass_j1_j2_plot_true->SetFillStyle(4000);
		hs_invmass_j1_j2->Add(QCD_reco_invmass_j1_j2_plot_true);

		VBS_reco_invmass_j1_j2_plot_true->SetFillColor(3);
		VBS_reco_invmass_j1_j2_plot_true->SetFillStyle(4000);
		VBS_reco_invmass_j1_j2_plot_true->SetLineWidth(2);
		VBS_reco_invmass_j1_j2_plot_true->Scale(1);
		hs_invmass_j1_j2->Add(VBS_reco_invmass_j1_j2_plot_true);

		INTER_reco_invmass_j1_j2_plot_true->SetFillColor(46);
		INTER_reco_invmass_j1_j2_plot_true->SetFillStyle(4000);
		INTER_reco_invmass_j1_j2_plot_true->SetLineWidth(4);
		INTER_reco_invmass_j1_j2_plot_true->Scale(1);
		hs_invmass_j1_j2->Add(INTER_reco_invmass_j1_j2_plot_true);
		hs_invmass_j1_j2->Draw("hist");

		Data_reco_invmass_j1_j2_plot_true->SetStats(0);
		Data_reco_invmass_j1_j2_plot_true->SetLineColor(1);
		Data_reco_invmass_j1_j2_plot_true->SetLineWidth(3);
		Data_reco_invmass_j1_j2_plot_true->Draw("same");

		TLegend *invmass_j1_j2_0_stack = new TLegend(0.7, 0.6, 0.9, 0.9);
		invmass_j1_j2_0_stack->AddEntry(VBS_reco_invmass_j1_j2_plot_true,"VBS_reco_invmass_j1_j2*10")  ;
		invmass_j1_j2_0_stack->Draw()  ;
		invmass_j1_j2_0_stack->AddEntry(QCD_reco_invmass_j1_j2_plot_true,"QCD_reco_invmass_j1_j2")  ;
		invmass_j1_j2_0_stack->Draw()  ;
		invmass_j1_j2_0_stack->AddEntry(INTER_reco_invmass_j1_j2_plot_true,"INTER_reco_invmass_j1_j2*100")  ;
		invmass_j1_j2_0_stack->Draw()  ;
		invmass_j1_j2_0_stack->AddEntry(Data_reco_invmass_j1_j2_plot_true,"Data_reco_invmass_j1_j2")  ;
		invmass_j1_j2_0_stack->Draw()  ;
		c102_stack->SaveAs("!!20invmass_j1_j2_stack.pdf");
/*
		TCanvas *c103= new TCanvas("c103","invmass0");

		VBS_reco_invmass_j1_j2_plot_true->SetStats(0);
		VBS_reco_invmass_j1_j2_plot_true->SetLineColor(3);
		VBS_reco_invmass_j1_j2_plot_true->SetLineWidth(3);
		VBS_reco_invmass_j1_j2_plot_true->GetYaxis()->SetRangeUser(0.0001, 10000000);
		VBS_reco_invmass_j1_j2_plot_true->Draw();

		QCD_reco_invmass_j1_j2_plot_true->SetStats(0);
		QCD_reco_invmass_j1_j2_plot_true->SetLineColor(4);
		QCD_reco_invmass_j1_j2_plot_true->SetLineWidth(2);
		QCD_reco_invmass_j1_j2_plot_true->Draw("same");

		VBS_reco_invmass_j1_j2_plot_true->Fit("func_vbs_inv");
		QCD_reco_invmass_j1_j2_plot_true->Fit("func_qcd_inv");
		
*/		


	}



//shower shape
	{
		//r9 photon1 c300
		{
			//INTER_reco_r9_p1_bar
			TCanvas *c300= new TCanvas("c300","C_300");

			QCD_reco_r9_p1_bar->SetStats(0);
			QCD_reco_r9_p1_bar->SetLineColor(4);
			QCD_reco_r9_p1_bar->SetLineWidth(2);
			QCD_reco_r9_p1_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_r9_p1_bar->Draw("hist");

			VBS_reco_r9_p1_bar->SetStats(0);
			VBS_reco_r9_p1_bar->SetLineColor(3);
			VBS_reco_r9_p1_bar->SetLineWidth(2);
			VBS_reco_r9_p1_bar->Draw("hist,same");

			INTER_reco_r9_p1_bar->SetStats(0);
			INTER_reco_r9_p1_bar->SetLineColor(2);
			INTER_reco_r9_p1_bar->SetLineWidth(2);
			INTER_reco_r9_p1_bar->Draw("hist,same");

			Data_reco_r9_p1_bar->SetStats(0);
			Data_reco_r9_p1_bar->SetLineColor(1);
			Data_reco_r9_p1_bar->SetLineWidth(2);
			Data_reco_r9_p1_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_r9_p1_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_r9_p1_bar_0->AddEntry(VBS_reco_r9_p1_bar,"VBS_reco_r9_p1_bar")  ;
			reco_r9_p1_bar_0->Draw()  ;
			reco_r9_p1_bar_0->AddEntry(QCD_reco_r9_p1_bar,"QCD_reco_r9_p1_bar")  ;
			reco_r9_p1_bar_0->Draw()  ;
			reco_r9_p1_bar_0->AddEntry(INTER_reco_r9_p1_bar,"INTER_reco_r9_p1_bar")  ;
			reco_r9_p1_bar_0->Draw()  ;
			reco_r9_p1_bar_0->AddEntry(Data_reco_r9_p1_bar,"Data_reco_r9_p1_bar")  ;
			reco_r9_p1_bar_0->Draw()  ;


			c300->SaveAs("!!31reco_r9_p1_bar.pdf");


			TCanvas *c300_0= new TCanvas("c300_0","C_300_0");

			QCD_reco_r9_p1_end->SetStats(0);
			QCD_reco_r9_p1_end->SetLineColor(4);
			QCD_reco_r9_p1_end->SetLineWidth(2);
			QCD_reco_r9_p1_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_r9_p1_end->Draw("hist");

			VBS_reco_r9_p1_end->SetStats(0);
			VBS_reco_r9_p1_end->SetLineColor(3);
			VBS_reco_r9_p1_end->SetLineWidth(2);
			VBS_reco_r9_p1_end->Draw("hist,same");

			INTER_reco_r9_p1_end->SetStats(0);
			INTER_reco_r9_p1_end->SetLineColor(2);
			INTER_reco_r9_p1_end->SetLineWidth(2);
			INTER_reco_r9_p1_end->Draw("hist,same");

			Data_reco_r9_p1_end->SetStats(0);
			Data_reco_r9_p1_end->SetLineColor(1);
			Data_reco_r9_p1_end->SetLineWidth(2);
			Data_reco_r9_p1_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_r9_p1_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_r9_p1_end_0->AddEntry(VBS_reco_r9_p1_end,"VBS_reco_r9_p1_end")  ;
			reco_r9_p1_end_0->Draw()  ;
			reco_r9_p1_end_0->AddEntry(QCD_reco_r9_p1_end,"QCD_reco_r9_p1_end")  ;
			reco_r9_p1_end_0->Draw()  ;
			reco_r9_p1_end_0->AddEntry(INTER_reco_r9_p1_end,"INTER_reco_r9_p1_end")  ;
			reco_r9_p1_end_0->Draw()  ;
			reco_r9_p1_end_0->AddEntry(Data_reco_r9_p1_end,"Data_reco_r9_p1_end")  ;
			reco_r9_p1_end_0->Draw()  ;

			c300_0->SaveAs("!!31reco_r9_p1_end.pdf");
		}

		//r9 photon2 c301
		{
			//INTER_reco_r9_p1_bar
			TCanvas *c301= new TCanvas("c301","C_301");

			QCD_reco_r9_p2_bar->SetStats(0);
			QCD_reco_r9_p2_bar->SetLineColor(4);
			QCD_reco_r9_p2_bar->SetLineWidth(2);
			QCD_reco_r9_p2_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_r9_p2_bar->Draw("hist");

			VBS_reco_r9_p2_bar->SetStats(0);
			VBS_reco_r9_p2_bar->SetLineColor(3);
			VBS_reco_r9_p2_bar->SetLineWidth(2);
			VBS_reco_r9_p2_bar->Draw("hist,same");

			INTER_reco_r9_p2_bar->SetStats(0);
			INTER_reco_r9_p2_bar->SetLineColor(2);
			INTER_reco_r9_p2_bar->SetLineWidth(2);
			INTER_reco_r9_p2_bar->Draw("hist,same");

			Data_reco_r9_p2_bar->SetStats(0);
			Data_reco_r9_p2_bar->SetLineColor(1);
			Data_reco_r9_p2_bar->SetLineWidth(2);
			Data_reco_r9_p2_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_r9_p2_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_r9_p2_bar_0->AddEntry(VBS_reco_r9_p2_bar,"VBS_reco_r9_p2_bar")  ;
			reco_r9_p2_bar_0->Draw()  ;
			reco_r9_p2_bar_0->AddEntry(QCD_reco_r9_p2_bar,"QCD_reco_r9_p2_bar")  ;
			reco_r9_p2_bar_0->Draw()  ;
			reco_r9_p2_bar_0->AddEntry(INTER_reco_r9_p2_bar,"INTER_reco_r9_p2_bar")  ;
			reco_r9_p2_bar_0->Draw()  ;
			reco_r9_p2_bar_0->AddEntry(Data_reco_r9_p2_bar,"Data_reco_r9_p2_bar")  ;
			reco_r9_p2_bar_0->Draw()  ;


			c301->SaveAs("!!31reco_r9_p2_bar.pdf");


			TCanvas *c301_0= new TCanvas("c301_0","C_301_0");

			QCD_reco_r9_p2_end->SetStats(0);
			QCD_reco_r9_p2_end->SetLineColor(4);
			QCD_reco_r9_p2_end->SetLineWidth(2);
			QCD_reco_r9_p2_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_r9_p2_end->Draw("hist");

			VBS_reco_r9_p2_end->SetStats(0);
			VBS_reco_r9_p2_end->SetLineColor(3);
			VBS_reco_r9_p2_end->SetLineWidth(2);
			VBS_reco_r9_p2_end->Draw("hist,same");

			INTER_reco_r9_p2_end->SetStats(0);
			INTER_reco_r9_p2_end->SetLineColor(2);
			INTER_reco_r9_p2_end->SetLineWidth(2);
			INTER_reco_r9_p2_end->Draw("hist,same");

			Data_reco_r9_p2_end->SetStats(0);
			Data_reco_r9_p2_end->SetLineColor(1);
			Data_reco_r9_p2_end->SetLineWidth(2);
			Data_reco_r9_p2_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_r9_p2_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_r9_p2_end_0->AddEntry(VBS_reco_r9_p2_end,"VBS_reco_r9_p2_end")  ;
			reco_r9_p2_end_0->Draw()  ;
			reco_r9_p2_end_0->AddEntry(QCD_reco_r9_p2_end,"QCD_reco_r9_p2_end")  ;
			reco_r9_p2_end_0->Draw()  ;
			reco_r9_p2_end_0->AddEntry(INTER_reco_r9_p2_end,"INTER_reco_r9_p2_end")  ;
			reco_r9_p2_end_0->Draw()  ;
			reco_r9_p2_end_0->AddEntry(Data_reco_r9_p2_end,"Data_reco_r9_p2_end")  ;
			reco_r9_p2_end_0->Draw()  ;

			c301_0->SaveAs("!!31reco_r9_p2_end.pdf");
		}

		//sieie photon1 c302
		{

			//INTER_reco_r9_p1_bar
			TCanvas *c302= new TCanvas("c302","C_302");

			QCD_reco_sieie_p1_bar->SetStats(0);
			QCD_reco_sieie_p1_bar->SetLineColor(4);
			QCD_reco_sieie_p1_bar->SetLineWidth(2);
			QCD_reco_sieie_p1_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_sieie_p1_bar->Draw("hist");

			VBS_reco_sieie_p1_bar->SetStats(0);
			VBS_reco_sieie_p1_bar->SetLineColor(3);
			VBS_reco_sieie_p1_bar->SetLineWidth(2);
			VBS_reco_sieie_p1_bar->Draw("hist,same");

			INTER_reco_sieie_p1_bar->SetStats(0);
			INTER_reco_sieie_p1_bar->SetLineColor(2);
			INTER_reco_sieie_p1_bar->SetLineWidth(2);
			INTER_reco_sieie_p1_bar->Draw("hist,same");

			Data_reco_sieie_p1_bar->SetStats(0);
			Data_reco_sieie_p1_bar->SetLineColor(1);
			Data_reco_sieie_p1_bar->SetLineWidth(2);
			Data_reco_sieie_p1_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_sieie_p1_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_sieie_p1_bar_0->AddEntry(VBS_reco_sieie_p1_bar,"VBS_reco_sieie_p1_bar")  ;
			reco_sieie_p1_bar_0->Draw()  ;
			reco_sieie_p1_bar_0->AddEntry(QCD_reco_sieie_p1_bar,"QCD_reco_sieie_p1_bar")  ;
			reco_sieie_p1_bar_0->Draw()  ;
			reco_sieie_p1_bar_0->AddEntry(INTER_reco_sieie_p1_bar,"INTER_reco_sieie_p1_bar")  ;
			reco_sieie_p1_bar_0->Draw()  ;
			reco_sieie_p1_bar_0->AddEntry(Data_reco_sieie_p1_bar,"Data_reco_sieie_p1_bar")  ;
			reco_sieie_p1_bar_0->Draw()  ;


			c302->SaveAs("!!31reco_sieie_p1_bar.pdf");


			TCanvas *c302_0= new TCanvas("c302_0","C_302_0");

			QCD_reco_sieie_p1_end->SetStats(0);
			QCD_reco_sieie_p1_end->SetLineColor(4);
			QCD_reco_sieie_p1_end->SetLineWidth(2);
			QCD_reco_sieie_p1_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_sieie_p1_end->Draw("hist");

			VBS_reco_sieie_p1_end->SetStats(0);
			VBS_reco_sieie_p1_end->SetLineColor(3);
			VBS_reco_sieie_p1_end->SetLineWidth(2);
			VBS_reco_sieie_p1_end->Draw("hist,same");

			INTER_reco_sieie_p1_end->SetStats(0);
			INTER_reco_sieie_p1_end->SetLineColor(2);
			INTER_reco_sieie_p1_end->SetLineWidth(2);
			INTER_reco_sieie_p1_end->Draw("hist,same");

			Data_reco_sieie_p1_end->SetStats(0);
			Data_reco_sieie_p1_end->SetLineColor(1);
			Data_reco_sieie_p1_end->SetLineWidth(2);
			Data_reco_sieie_p1_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_sieie_p1_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_sieie_p1_end_0->AddEntry(VBS_reco_sieie_p1_end,"VBS_reco_sieie_p1_end")  ;
			reco_sieie_p1_end_0->Draw()  ;
			reco_sieie_p1_end_0->AddEntry(QCD_reco_sieie_p1_end,"QCD_reco_sieie_p1_end")  ;
			reco_sieie_p1_end_0->Draw()  ;
			reco_sieie_p1_end_0->AddEntry(INTER_reco_sieie_p1_end,"INTER_reco_sieie_p1_end")  ;
			reco_sieie_p1_end_0->Draw()  ;
			reco_sieie_p1_end_0->AddEntry(Data_reco_sieie_p1_end,"Data_reco_sieie_p1_end")  ;
			reco_sieie_p1_end_0->Draw()  ;

			c302_0->SaveAs("!!31reco_sieie_p1_end.pdf");
		}

		//sieie photon2 c303
		{
			//INTER_reco_r9_p2_bar
			TCanvas *c303= new TCanvas("c303","C_303");

			QCD_reco_sieie_p2_bar->SetStats(0);
			QCD_reco_sieie_p2_bar->SetLineColor(4);
			QCD_reco_sieie_p2_bar->SetLineWidth(2);
			QCD_reco_sieie_p2_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_sieie_p2_bar->Draw("hist");

			VBS_reco_sieie_p2_bar->SetStats(0);
			VBS_reco_sieie_p2_bar->SetLineColor(3);
			VBS_reco_sieie_p2_bar->SetLineWidth(2);
			VBS_reco_sieie_p2_bar->Draw("hist,same");

			INTER_reco_sieie_p2_bar->SetStats(0);
			INTER_reco_sieie_p2_bar->SetLineColor(2);
			INTER_reco_sieie_p2_bar->SetLineWidth(2);
			INTER_reco_sieie_p2_bar->Draw("hist,same");

			Data_reco_sieie_p2_bar->SetStats(0);
			Data_reco_sieie_p2_bar->SetLineColor(1);
			Data_reco_sieie_p2_bar->SetLineWidth(2);
			Data_reco_sieie_p2_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_sieie_p2_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_sieie_p2_bar_0->AddEntry(VBS_reco_sieie_p2_bar,"VBS_reco_sieie_p2_bar")  ;
			reco_sieie_p2_bar_0->Draw()  ;
			reco_sieie_p2_bar_0->AddEntry(QCD_reco_sieie_p2_bar,"QCD_reco_sieie_p2_bar")  ;
			reco_sieie_p2_bar_0->Draw()  ;
			reco_sieie_p2_bar_0->AddEntry(INTER_reco_sieie_p2_bar,"INTER_reco_sieie_p2_bar")  ;
			reco_sieie_p2_bar_0->Draw()  ;
			reco_sieie_p2_bar_0->AddEntry(Data_reco_sieie_p2_bar,"Data_reco_sieie_p2_bar")  ;
			reco_sieie_p2_bar_0->Draw()  ;


			c303->SaveAs("!!31reco_sieie_p2_bar.pdf");


			TCanvas *c303_0= new TCanvas("c303_0","C_303_0");

			QCD_reco_sieie_p2_end->SetStats(0);
			QCD_reco_sieie_p2_end->SetLineColor(4);
			QCD_reco_sieie_p2_end->SetLineWidth(2);
			QCD_reco_sieie_p2_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_sieie_p2_end->Draw("hist");

			VBS_reco_sieie_p2_end->SetStats(0);
			VBS_reco_sieie_p2_end->SetLineColor(3);
			VBS_reco_sieie_p2_end->SetLineWidth(2);
			VBS_reco_sieie_p2_end->Draw("hist,same");

			INTER_reco_sieie_p2_end->SetStats(0);
			INTER_reco_sieie_p2_end->SetLineColor(2);
			INTER_reco_sieie_p2_end->SetLineWidth(2);
			INTER_reco_sieie_p2_end->Draw("hist,same");

			Data_reco_sieie_p2_end->SetStats(0);
			Data_reco_sieie_p2_end->SetLineColor(1);
			Data_reco_sieie_p2_end->SetLineWidth(2);
			Data_reco_sieie_p2_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_sieie_p2_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_sieie_p2_end_0->AddEntry(VBS_reco_sieie_p2_end,"VBS_reco_sieie_p2_end")  ;
			reco_sieie_p2_end_0->Draw()  ;
			reco_sieie_p2_end_0->AddEntry(QCD_reco_sieie_p2_end,"QCD_reco_sieie_p2_end")  ;
			reco_sieie_p2_end_0->Draw()  ;
			reco_sieie_p2_end_0->AddEntry(INTER_reco_sieie_p2_end,"INTER_reco_sieie_p2_end")  ;
			reco_sieie_p2_end_0->Draw()  ;
			reco_sieie_p2_end_0->AddEntry(Data_reco_sieie_p2_end,"Data_reco_sieie_p2_end")  ;
			reco_sieie_p2_end_0->Draw()  ;

			c303_0->SaveAs("!!31reco_sieie_p2_end.pdf");
		}

		//hoe photon1 c304
		{
			//INTER_reco_r9_p1_bar
			TCanvas *c304= new TCanvas("c304","C_304");

			QCD_reco_hoe_p1_bar->SetStats(0);
			QCD_reco_hoe_p1_bar->SetLineColor(4);
			QCD_reco_hoe_p1_bar->SetLineWidth(2);
			QCD_reco_hoe_p1_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_hoe_p1_bar->Draw("hist");

			VBS_reco_hoe_p1_bar->SetStats(0);
			VBS_reco_hoe_p1_bar->SetLineColor(3);
			VBS_reco_hoe_p1_bar->SetLineWidth(2);
			VBS_reco_hoe_p1_bar->Draw("hist,same");

			INTER_reco_hoe_p1_bar->SetStats(0);
			INTER_reco_hoe_p1_bar->SetLineColor(2);
			INTER_reco_hoe_p1_bar->SetLineWidth(2);
			INTER_reco_hoe_p1_bar->Draw("hist,same");

			Data_reco_hoe_p1_bar->SetStats(0);
			Data_reco_hoe_p1_bar->SetLineColor(1);
			Data_reco_hoe_p1_bar->SetLineWidth(2);
			Data_reco_hoe_p1_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_hoe_p1_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_hoe_p1_bar_0->AddEntry(VBS_reco_hoe_p1_bar,"VBS_reco_hoe_p1_bar")  ;
			reco_hoe_p1_bar_0->Draw()  ;
			reco_hoe_p1_bar_0->AddEntry(QCD_reco_hoe_p1_bar,"QCD_reco_hoe_p1_bar")  ;
			reco_hoe_p1_bar_0->Draw()  ;
			reco_hoe_p1_bar_0->AddEntry(INTER_reco_hoe_p1_bar,"INTER_reco_hoe_p1_bar")  ;
			reco_hoe_p1_bar_0->Draw()  ;
			reco_hoe_p1_bar_0->AddEntry(Data_reco_hoe_p1_bar,"Data_reco_hoe_p1_bar")  ;
			reco_hoe_p1_bar_0->Draw()  ;


			c304->SaveAs("!!31reco_hoe_p1_bar.pdf");


			TCanvas *c304_0= new TCanvas("c304_0","C_304_0");

			QCD_reco_hoe_p1_end->SetStats(0);
			QCD_reco_hoe_p1_end->SetLineColor(4);
			QCD_reco_hoe_p1_end->SetLineWidth(2);
			QCD_reco_hoe_p1_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_hoe_p1_end->Draw("hist");

			VBS_reco_hoe_p1_end->SetStats(0);
			VBS_reco_hoe_p1_end->SetLineColor(3);
			VBS_reco_hoe_p1_end->SetLineWidth(2);
			VBS_reco_hoe_p1_end->Draw("hist,same");

			INTER_reco_hoe_p1_end->SetStats(0);
			INTER_reco_hoe_p1_end->SetLineColor(2);
			INTER_reco_hoe_p1_end->SetLineWidth(2);
			INTER_reco_hoe_p1_end->Draw("hist,same");

			Data_reco_hoe_p1_end->SetStats(0);
			Data_reco_hoe_p1_end->SetLineColor(1);
			Data_reco_hoe_p1_end->SetLineWidth(2);
			Data_reco_hoe_p1_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_hoe_p1_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_hoe_p1_end_0->AddEntry(VBS_reco_hoe_p1_end,"VBS_reco_hoe_p1_end")  ;
			reco_hoe_p1_end_0->Draw()  ;
			reco_hoe_p1_end_0->AddEntry(QCD_reco_hoe_p1_end,"QCD_reco_hoe_p1_end")  ;
			reco_hoe_p1_end_0->Draw()  ;
			reco_hoe_p1_end_0->AddEntry(INTER_reco_hoe_p1_end,"INTER_reco_hoe_p1_end")  ;
			reco_hoe_p1_end_0->Draw()  ;
			reco_hoe_p1_end_0->AddEntry(Data_reco_hoe_p1_end,"Data_reco_hoe_p1_end")  ;
			reco_hoe_p1_end_0->Draw()  ;

			c304_0->SaveAs("!!31reco_hoe_p1_end.pdf");
		}

		//hoe photon2 c305
		{
			//INTER_reco_r9_p2_bar
			TCanvas *c305= new TCanvas("c305","C_305");

			QCD_reco_hoe_p2_bar->SetStats(0);
			QCD_reco_hoe_p2_bar->SetLineColor(4);
			QCD_reco_hoe_p2_bar->SetLineWidth(2);
			QCD_reco_hoe_p2_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_hoe_p2_bar->Draw("hist");

			VBS_reco_hoe_p2_bar->SetStats(0);
			VBS_reco_hoe_p2_bar->SetLineColor(3);
			VBS_reco_hoe_p2_bar->SetLineWidth(2);
			VBS_reco_hoe_p2_bar->Draw("hist,same");

			INTER_reco_hoe_p2_bar->SetStats(0);
			INTER_reco_hoe_p2_bar->SetLineColor(2);
			INTER_reco_hoe_p2_bar->SetLineWidth(2);
			INTER_reco_hoe_p2_bar->Draw("hist,same");

			Data_reco_hoe_p2_bar->SetStats(0);
			Data_reco_hoe_p2_bar->SetLineColor(1);
			Data_reco_hoe_p2_bar->SetLineWidth(2);
			Data_reco_hoe_p2_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_hoe_p2_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_hoe_p2_bar_0->AddEntry(VBS_reco_hoe_p2_bar,"VBS_reco_hoe_p2_bar")  ;
			reco_hoe_p2_bar_0->Draw()  ;
			reco_hoe_p2_bar_0->AddEntry(QCD_reco_hoe_p2_bar,"QCD_reco_hoe_p2_bar")  ;
			reco_hoe_p2_bar_0->Draw()  ;
			reco_hoe_p2_bar_0->AddEntry(INTER_reco_hoe_p2_bar,"INTER_reco_hoe_p2_bar")  ;
			reco_hoe_p2_bar_0->Draw()  ;
			reco_hoe_p2_bar_0->AddEntry(Data_reco_hoe_p2_bar,"Data_reco_hoe_p2_bar")  ;
			reco_hoe_p2_bar_0->Draw()  ;


			c305->SaveAs("!!31reco_hoe_p2_bar.pdf");


			TCanvas *c305_0= new TCanvas("c305_0","C_305_0");

			QCD_reco_hoe_p2_end->SetStats(0);
			QCD_reco_hoe_p2_end->SetLineColor(4);
			QCD_reco_hoe_p2_end->SetLineWidth(2);
			QCD_reco_hoe_p2_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_hoe_p2_end->Draw("hist");

			VBS_reco_hoe_p2_end->SetStats(0);
			VBS_reco_hoe_p2_end->SetLineColor(3);
			VBS_reco_hoe_p2_end->SetLineWidth(2);
			VBS_reco_hoe_p2_end->Draw("hist,same");

			INTER_reco_hoe_p2_end->SetStats(0);
			INTER_reco_hoe_p2_end->SetLineColor(2);
			INTER_reco_hoe_p2_end->SetLineWidth(2);
			INTER_reco_hoe_p2_end->Draw("hist,same");

			Data_reco_hoe_p2_end->SetStats(0);
			Data_reco_hoe_p2_end->SetLineColor(1);
			Data_reco_hoe_p2_end->SetLineWidth(2);
			Data_reco_hoe_p2_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_hoe_p2_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_hoe_p2_end_0->AddEntry(VBS_reco_hoe_p2_end,"VBS_reco_hoe_p2_end")  ;
			reco_hoe_p2_end_0->Draw()  ;
			reco_hoe_p2_end_0->AddEntry(QCD_reco_hoe_p2_end,"QCD_reco_hoe_p2_end")  ;
			reco_hoe_p2_end_0->Draw()  ;
			reco_hoe_p2_end_0->AddEntry(INTER_reco_hoe_p2_end,"INTER_reco_hoe_p2_end")  ;
			reco_hoe_p2_end_0->Draw()  ;
			reco_hoe_p2_end_0->AddEntry(Data_reco_hoe_p2_end,"Data_reco_hoe_p2_end")  ;
			reco_hoe_p2_end_0->Draw()  ;

			c305_0->SaveAs("!!31reco_hoe_p2_end.pdf");
		}

		//chiso photon1 c306
		{
			
			//INTER_reco_r9_p1_bar
			TCanvas *c306= new TCanvas("c306","C_306");

			QCD_reco_chiso_p1_bar->SetStats(0);
			QCD_reco_chiso_p1_bar->SetLineColor(4);
			QCD_reco_chiso_p1_bar->SetLineWidth(2);
			QCD_reco_chiso_p1_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_chiso_p1_bar->Draw("hist");

			VBS_reco_chiso_p1_bar->SetStats(0);
			VBS_reco_chiso_p1_bar->SetLineColor(3);
			VBS_reco_chiso_p1_bar->SetLineWidth(2);
			VBS_reco_chiso_p1_bar->Draw("hist,same");

			INTER_reco_chiso_p1_bar->SetStats(0);
			INTER_reco_chiso_p1_bar->SetLineColor(2);
			INTER_reco_chiso_p1_bar->SetLineWidth(2);
			INTER_reco_chiso_p1_bar->Draw("hist,same");

			Data_reco_chiso_p1_bar->SetStats(0);
			Data_reco_chiso_p1_bar->SetLineColor(1);
			Data_reco_chiso_p1_bar->SetLineWidth(2);
			Data_reco_chiso_p1_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_chiso_p1_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_chiso_p1_bar_0->AddEntry(VBS_reco_chiso_p1_bar,"VBS_reco_chiso_p1_bar")  ;
			reco_chiso_p1_bar_0->Draw()  ;
			reco_chiso_p1_bar_0->AddEntry(QCD_reco_chiso_p1_bar,"QCD_reco_chiso_p1_bar")  ;
			reco_chiso_p1_bar_0->Draw()  ;
			reco_chiso_p1_bar_0->AddEntry(INTER_reco_chiso_p1_bar,"INTER_reco_chiso_p1_bar")  ;
			reco_chiso_p1_bar_0->Draw()  ;
			reco_chiso_p1_bar_0->AddEntry(Data_reco_chiso_p1_bar,"Data_reco_chiso_p1_bar")  ;
			reco_chiso_p1_bar_0->Draw()  ;


			c306->SaveAs("!!31reco_chiso_p1_bar.pdf");


			TCanvas *c306_0= new TCanvas("c306_0","C_306_0");

			QCD_reco_chiso_p1_end->SetStats(0);
			QCD_reco_chiso_p1_end->SetLineColor(4);
			QCD_reco_chiso_p1_end->SetLineWidth(2);
			QCD_reco_chiso_p1_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_chiso_p1_end->Draw("hist");

			VBS_reco_chiso_p1_end->SetStats(0);
			VBS_reco_chiso_p1_end->SetLineColor(3);
			VBS_reco_chiso_p1_end->SetLineWidth(2);
			VBS_reco_chiso_p1_end->Draw("hist,same");

			INTER_reco_chiso_p1_end->SetStats(0);
			INTER_reco_chiso_p1_end->SetLineColor(2);
			INTER_reco_chiso_p1_end->SetLineWidth(2);
			INTER_reco_chiso_p1_end->Draw("hist,same");

			Data_reco_chiso_p1_end->SetStats(0);
			Data_reco_chiso_p1_end->SetLineColor(1);
			Data_reco_chiso_p1_end->SetLineWidth(2);
			Data_reco_chiso_p1_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_chiso_p1_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_chiso_p1_end_0->AddEntry(VBS_reco_chiso_p1_end,"VBS_reco_chiso_p1_end")  ;
			reco_chiso_p1_end_0->Draw()  ;
			reco_chiso_p1_end_0->AddEntry(QCD_reco_chiso_p1_end,"QCD_reco_chiso_p1_end")  ;
			reco_chiso_p1_end_0->Draw()  ;
			reco_chiso_p1_end_0->AddEntry(INTER_reco_chiso_p1_end,"INTER_reco_chiso_p1_end")  ;
			reco_chiso_p1_end_0->Draw()  ;
			reco_chiso_p1_end_0->AddEntry(Data_reco_chiso_p1_end,"Data_reco_chiso_p1_end")  ;
			reco_chiso_p1_end_0->Draw()  ;

			c306_0->SaveAs("!!31reco_chiso_p1_end.pdf");
		}

		//chiso photon2 c307
	{
			
			//INTER_reco_r9_p1_bar
			TCanvas *c307= new TCanvas("c307","C_307");

			QCD_reco_chiso_p2_bar->SetStats(0);
			QCD_reco_chiso_p2_bar->SetLineColor(4);
			QCD_reco_chiso_p2_bar->SetLineWidth(2);
			QCD_reco_chiso_p2_bar->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_chiso_p2_bar->Draw("hist");

			VBS_reco_chiso_p2_bar->SetStats(0);
			VBS_reco_chiso_p2_bar->SetLineColor(3);
			VBS_reco_chiso_p2_bar->SetLineWidth(2);
			VBS_reco_chiso_p2_bar->Draw("hist,same");

			INTER_reco_chiso_p2_bar->SetStats(0);
			INTER_reco_chiso_p2_bar->SetLineColor(2);
			INTER_reco_chiso_p2_bar->SetLineWidth(2);
			INTER_reco_chiso_p2_bar->Draw("hist,same");

			Data_reco_chiso_p2_bar->SetStats(0);
			Data_reco_chiso_p2_bar->SetLineColor(1);
			Data_reco_chiso_p2_bar->SetLineWidth(2);
			Data_reco_chiso_p2_bar->Draw("E1,same");

			gPad->SetLogy();

			TLegend *reco_chiso_p2_bar_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_chiso_p2_bar_0->AddEntry(VBS_reco_chiso_p2_bar,"VBS_reco_chiso_p2_bar")  ;
			reco_chiso_p2_bar_0->Draw()  ;
			reco_chiso_p2_bar_0->AddEntry(QCD_reco_chiso_p2_bar,"QCD_reco_chiso_p2_bar")  ;
			reco_chiso_p2_bar_0->Draw()  ;
			reco_chiso_p2_bar_0->AddEntry(INTER_reco_chiso_p2_bar,"INTER_reco_chiso_p2_bar")  ;
			reco_chiso_p2_bar_0->Draw()  ;
			reco_chiso_p2_bar_0->AddEntry(Data_reco_chiso_p2_bar,"Data_reco_chiso_p2_bar")  ;
			reco_chiso_p2_bar_0->Draw()  ;


			c307->SaveAs("!!31reco_chiso_p2_bar.pdf");


			TCanvas *c307_0= new TCanvas("c307_0","C_307_0");

			QCD_reco_chiso_p2_end->SetStats(0);
			QCD_reco_chiso_p2_end->SetLineColor(4);
			QCD_reco_chiso_p2_end->SetLineWidth(2);
			QCD_reco_chiso_p2_end->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_chiso_p2_end->Draw("hist");

			VBS_reco_chiso_p2_end->SetStats(0);
			VBS_reco_chiso_p2_end->SetLineColor(3);
			VBS_reco_chiso_p2_end->SetLineWidth(2);
			VBS_reco_chiso_p2_end->Draw("hist,same");

			INTER_reco_chiso_p2_end->SetStats(0);
			INTER_reco_chiso_p2_end->SetLineColor(2);
			INTER_reco_chiso_p2_end->SetLineWidth(2);
			INTER_reco_chiso_p2_end->Draw("hist,same");

			Data_reco_chiso_p2_end->SetStats(0);
			Data_reco_chiso_p2_end->SetLineColor(1);
			Data_reco_chiso_p2_end->SetLineWidth(2);
			Data_reco_chiso_p2_end->Draw("E1,same");

			gPad->SetLogy();
			TLegend *reco_chiso_p2_end_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_chiso_p2_end_0->AddEntry(VBS_reco_chiso_p2_end,"VBS_reco_chiso_p2_end")  ;
			reco_chiso_p2_end_0->Draw()  ;
			reco_chiso_p2_end_0->AddEntry(QCD_reco_chiso_p2_end,"QCD_reco_chiso_p2_end")  ;
			reco_chiso_p2_end_0->Draw()  ;
			reco_chiso_p2_end_0->AddEntry(INTER_reco_chiso_p2_end,"INTER_reco_chiso_p2_end")  ;
			reco_chiso_p2_end_0->Draw()  ;
			reco_chiso_p2_end_0->AddEntry(Data_reco_chiso_p2_end,"Data_reco_chiso_p2_end")  ;
			reco_chiso_p2_end_0->Draw()  ;

			c307_0->SaveAs("!!31reco_chiso_p2_end.pdf");
		}

		//mva photon1 c308 INTER_reco_mva_p1
		{

			TCanvas *c308= new TCanvas("c308","c_308");

			QCD_reco_mva_p1->SetStats(0);
			QCD_reco_mva_p1->SetLineColor(4);
			QCD_reco_mva_p1->SetLineWidth(2);
			QCD_reco_mva_p1->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_mva_p1->Draw("hist");

			VBS_reco_mva_p1->SetStats(0);
			VBS_reco_mva_p1->SetLineColor(3);
			VBS_reco_mva_p1->SetLineWidth(2);
			VBS_reco_mva_p1->Draw("hist,same");

			INTER_reco_mva_p1->SetStats(0);
			INTER_reco_mva_p1->SetLineColor(46);
			INTER_reco_mva_p1->SetLineWidth(2);
			INTER_reco_mva_p1->Draw("hist,same");

			Data_reco_mva_p1->SetStats(0);
			Data_reco_mva_p1->SetLineColor(1);
			Data_reco_mva_p1->SetLineWidth(2);
			Data_reco_mva_p1->Draw("hist,same");
			gPad->SetLogy();

			TLegend *reco_mva_p1_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_mva_p1_0->AddEntry(VBS_reco_mva_p1,"VBS_reco_mva_p1")  ;
			reco_mva_p1_0->Draw()  ;
			reco_mva_p1_0->AddEntry(QCD_reco_mva_p1,"QCD_reco_mva_p1")  ;
			reco_mva_p1_0->Draw()  ;
			reco_mva_p1_0->AddEntry(INTER_reco_mva_p1,"INTER_reco_mva_p1")  ;
			reco_mva_p1_0->Draw()  ;
			reco_mva_p1_0->AddEntry(Data_reco_mva_p1,"Data_reco_mva_p1")  ;
			reco_mva_p1_0->Draw()  ;

			c308->SaveAs("!!21reco_mva_p1.pdf");
		}

		//mva photon2 c309
		{
			TCanvas *c309= new TCanvas("c309","c_309");

			QCD_reco_mva_p2->SetStats(0);
			QCD_reco_mva_p2->SetLineColor(4);
			QCD_reco_mva_p2->SetLineWidth(2);
			QCD_reco_mva_p2->GetYaxis()->SetRangeUser(0.0001, 100000000);
			QCD_reco_mva_p2->Draw("hist");

			VBS_reco_mva_p2->SetStats(0);
			VBS_reco_mva_p2->SetLineColor(3);
			VBS_reco_mva_p2->SetLineWidth(2);
			VBS_reco_mva_p2->Draw("hist,same");

			INTER_reco_mva_p2->SetStats(0);
			INTER_reco_mva_p2->SetLineColor(46);
			INTER_reco_mva_p2->SetLineWidth(2);
			INTER_reco_mva_p2->Draw("hist,same");

			Data_reco_mva_p2->SetStats(0);
			Data_reco_mva_p2->SetLineColor(1);
			Data_reco_mva_p2->SetLineWidth(2);
			Data_reco_mva_p2->Draw("hist,same");
			gPad->SetLogy();

			TLegend *reco_mva_p2_0 = new TLegend(0.7, 0.6, 0.9, 0.9);
			reco_mva_p2_0->AddEntry(VBS_reco_mva_p2,"VBS_reco_mva_p2")  ;
			reco_mva_p2_0->Draw()  ;
			reco_mva_p2_0->AddEntry(QCD_reco_mva_p2,"QCD_reco_mva_p2")  ;
			reco_mva_p2_0->Draw()  ;
			reco_mva_p2_0->AddEntry(INTER_reco_mva_p2,"INTER_reco_mva_p2")  ;
			reco_mva_p2_0->Draw()  ;
			reco_mva_p2_0->AddEntry(Data_reco_mva_p2,"Data_reco_mva_p2")  ;
			reco_mva_p2_0->Draw()  ;

			c309->SaveAs("!!21reco_mva_p2.pdf");
		}



		gPad->SetLogy();
		TTimeStamp ts;
		cout << ts.AsString() << endl;

		t.Stop();
		t.Print();

		cout<<"Time taken: "<<(double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
		cout << " .................. " <<endl;

	}

}





