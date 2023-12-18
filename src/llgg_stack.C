#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"
#include "classes/SortableObject.h"
#include "modules/Delphes.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTask.h"
#include <string.h>
#include <TVector.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <time.h>
#include <TVector3.h>
#include <vector>
#include <algorithm>
#include <cmath>
#endif

bool compareTVector3(const TVector3& a, const TVector3& b) {
    return a.Mag() > b.Mag();
}
double calculateThrust(const std::vector<TLorentzVector>& momenta) {
    double t;
    TVector3 totalMomentum(0, 0, 0);
    for (const auto& particle : momenta) {
	totalMomentum += particle.Vect();
    }

    double totalEnergy = 0;
    for (const auto& particle : momenta) {
	totalEnergy += particle.E();
    }
    TVector3 beta = -totalMomentum * (1.0 / totalEnergy);
    
    std::vector<TVector3> momentaCp;
    for (const auto& particle : momenta) {
	TLorentzVector tmpPar = particle;
	tmpPar.Boost(beta);
	momentaCp.push_back(tmpPar.Vect());
    }
    unsigned  int n = momentaCp.size();
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;

    assert(n >= 3);
    n = 3;
    if (momentaCp.size() == 3) n = 3;

    std::vector<TVector3> tvec;
    std::vector<double> tval;
    std::sort(momentaCp.begin(), momentaCp.end(), compareTVector3);

    /*
    t = 0;
    TVector3 tmpAxis = momentaCp[0].Unit();

    for (unsigned int k = 0; k < momentaCp.size(); k++) {
	t += fabs(tmpAxis.Dot(momentaCp[k]));
    }
    */
    
    for (unsigned int i = 0; i < pow(2, n - 1); i++) {
	TVector3 foo(0, 0, 0);
	int sign = i;
	for (unsigned int k = 0; k < n; k++) {
	    (sign % 2) == 1 ? foo += momentaCp[k] : foo -= momentaCp[k];
	    sign /= 2;
	}
	foo = foo.Unit();

	double diff = 999.;
	while (diff > 1e-5) {
	    TVector3 foobar(0, 0, 0);
	    for (unsigned int k = 0; k < momentaCp.size(); k++) {
		foo.Dot(momentaCp[k]) > 0 ? foobar += momentaCp[k] : foobar -= momentaCp[k];
	    }
	    diff = (foo - foobar.Unit()).Mag();
	    foo = foobar.Unit();
	}

	t = 0.;
	for (unsigned int k = 0; k < momentaCp.size(); k++) {
	    t += fabs(foo.Dot(momentaCp[k]));
	}

	tval.push_back(t);
	tvec.push_back(foo);
    }
    
    double norm = 0.0;
    for (const auto& p : momentaCp) {
        norm += p.Mag();
    }
    
    t = 0.;
    for (unsigned int i = 0; i < tvec.size(); i++) {
        if (tval[i] > t) {
     	    t = tval[i];
	    std::cout << "new biggest thrust = " << t/norm << std::endl;
	    //taxis = tvec[i];
	}
    }
    
    return t / norm;
}
void llgg(const char *inputFile, TH1D *Hmass) {

    std::cout << "loading Delphes library" << std::endl;
    gSystem -> Load("/eos/user/h/hjia/ATLAS_HH_bbgg/Delphes-3.5.0/libDelphes.so");
    //gSystem->Load("libDelphes");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;

    std::cout << "open up file" << std::endl;
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t nEntries = treeReader->GetEntries();

    TClonesArray *branchsmallJet = treeReader->UseBranch("UniqueJet");
    TClonesArray *branchJet = treeReader->UseBranch("PFJet10");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchGenPar = treeReader->UseBranch("Particle");
    TClonesArray *branchEvent = treeReader->UseBranch("Event");

    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    //TH1D *Zmass = new TH1D("Zmass", "lepton pair invariant mass", 100, 30, 220);
    //TH1D *Hmass = new TH1D("Hmass", "jet pair invariant mass", 30, 30, 220);
    //TH1D *tau21 = new TH1D("tau21", "higgs tau21", 30, 0, 1);
    //TH1D *nneutralhisto = new TH1D("nneutralhisto", "number of neutral constituents", 30, 0, 70);
    //TH1D *pthisto = new TH1D("pthisto", "higgs pt", 30, 0, 1000);
    //TH1D *ChargeFrachisto = new TH1D("ChargeFrachisto", "higgs Charged Energy Fraction", 30, 0, 1);
    //TH1D *btaghisto = new TH1D("btaghisto", "higgs jet btag", 4, 0, 2);
    //TH1D *zhdishisto = new TH1D("zhdishisto", "z and higgs jet distance", 30, 0, 4);
    //TH1D *truejethisto = new TH1D("truejethisto", "true higgs and higgs jet distance", 30, 0, 4);

    std::cout << "Start loop ever tree" << std::endl;

    Int_t total = 0;
    Int_t find = 0;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
	bool higgstag = false;
	treeReader->ReadEntry(entry);

	Int_t nJet = branchJet->GetEntriesFast();
	Int_t ne = branchElectron->GetEntries();
	Int_t nMu = branchMuon->GetEntries();
	Int_t nGenPar = branchGenPar->GetEntries();

	TObject *object;

	Int_t emuflag = 0; 

	if (ne == 2 and nMu == 0) {
	    emuflag = 1;
	}
	if (nMu == 2 and ne == 0) {
	    emuflag = 2;
	}
	if (emuflag == 0) {
	    continue;
	}
	if (nJet == 0) {
	    continue;
	}


	Electron *elec1, *elec2;
	Muon *mu1, *mu2;

	Double_t l1pt;
	Double_t l1eta;
	Double_t l1phi;
	Double_t l1charge;

	Double_t l2pt;
	Double_t l2eta;
	Double_t l2phi;
	Double_t l2charge;

	Double_t lmass;

	Double_t hpt;
	Double_t heta;
	Double_t hphi;
	Double_t hmass;
	Double_t htau21;
	Double_t htau42;
	Int_t hNCharged;
	Int_t hNNeutral;
	Double_t hChargeFrac;
	Int_t hbtag = 2;

	TLorentzVector l1;
	TLorentzVector l2;
	TLorentzVector h;
	l1.SetPtEtaPhiM(0,0,0,0);
	l2.SetPtEtaPhiM(0,0,0,0);
	h.SetPtEtaPhiM(0,0,0,0);

	switch (emuflag) {
	    case 1:
		elec1 = (Electron *) branchElectron->At(0);
		elec2 = (Electron *) branchElectron->At(1);

	        l1pt = elec1 -> PT;
	        l1phi = elec1 -> Phi;
	        l1eta = elec1 -> Eta;
	        l1charge = elec1 -> Charge;

	        l2pt = elec2 -> PT;
	        l2phi = elec2 -> Phi;
	        l2eta = elec2 -> Eta;
	        l2charge = elec2 -> Charge;

	        lmass = 0.000511;
	        break;
	    case 2:
	        mu1 = (Muon *) branchMuon->At(0);
		mu2 = (Muon *) branchMuon->At(1);

	        l1pt = mu1 -> PT;
	        l1phi = mu1 -> Phi;
	        l1eta = mu1 -> Eta;
	        l1charge = mu1 -> Charge;

	        l2pt = mu2 -> PT;
	        l2phi = mu2 -> Phi;
	        l2eta = mu2 -> Eta;
	        l2charge = mu2 -> Charge;

		lmass = 0.10566;
	        break;
	    default:
		std::cout << "error!" << std::endl;
	}

	if (l1charge + l2charge != 0) {
	    continue;
	}

	l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, lmass);
	l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, lmass);

	TLorentzVector z_ll;
	z_ll.SetPtEtaPhiM(0, 0, 0, 0);
	TLorentzVector h_gg;
	h_gg.SetPtEtaPhiM(0, 0, 0, 0);

	z_ll = l1 + l2;
	if (l1.DeltaR(l2) ==0) {
	    continue;
	}
	Double_t mass_min = 100000000;
	Jet *jet;
	Jet *hjet;
	Int_t hentry = 0;
	bool find_jet = false;
	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    jet = (Jet *) branchJet -> At(jentry);
	    
	    if (fabs(jet -> Eta) > 3.0) {
		continue;
	    }
	
	    if (TMath::Abs(jet->P4().DeltaR(z_ll) - TMath::Pi())>0.2) {
	    	continue;
	    }
	    if (TMath::Abs(jet -> Mass - 125) < mass_min) {
			mass_min = TMath::Abs(jet -> Mass - 125);
	        //hjet = (Jet *) branchJet -> At(jentry);
			hentry = jentry;
			find_jet = true; 
	    }
	}
	if (find_jet == false) {
	    continue;
	}

	//std::cout << hentry << std::endl;
	hjet = (Jet *) branchJet -> At(hentry);
	if (hjet == 0) {
	    std::cout << "unaccessible" << std::endl;
	    continue;
	}
	/*
	if (hjet -> BTag == 1){
	    continue;
	}
	*/

	hmass = hjet -> Mass;
	hpt = hjet -> PT;
	heta = hjet -> Eta;
	hphi = hjet -> Phi;
	htau21 = (hjet -> Tau[1])/(hjet -> Tau[0]);
	htau42 = (hjet -> Tau[3])/(hjet -> Tau[1]);
	hNNeutral = hjet -> NNeutrals;
	hNCharged = hjet -> NCharged;
	hChargeFrac = hjet -> ChargedEnergyFraction;
	hbtag = hjet -> BTag;
	
	/*
	if (htau42 < 0.5) {
	    continue;
	}
	if (htau21 > 0.5) {
	    continue;
	}
	*/
	/*
	if (hNCharged < 20) {
	    continue;
	}
	if (hNNeutral < 15) {
	    continue;
	}
	*/
	// Uncomment to get jet constituents info
	/*
	for(int i = 0; i < hjet->Constituents.GetEntriesFast(); ++i) {
	    object = hjet->Constituents.At(i);
	    if(object == 0) continue;
	        if(object->IsA() == Tower::Class()) {
	  	    Tower *tower;
		    tower = (Tower*) object;
		    std::cout << "is tower with Eem = " << tower->Eem << "and Ehad = " <<tower->Ehad;
	    	}
	    	if(object->IsA() == GenParticle::Class()) {
		    std::cout << "is genpar";
	    	} 
		if(object->IsA() == Track::Class()) {
		    Track *track;
		    track = (Track*) object;
		    std::cout << "is track";
		}
	    	if(object->IsA() == Candidate::Class()) {
		    std::cout << "is candidate";
	    	}
	    	std::cout << std::endl;
	}
	*/

	h_gg.SetPtEtaPhiM(hpt, heta, hphi, hmass);

	if (z_ll.M() < 76 or z_ll.M() > 106) {
	    continue;
	}
	if (h_gg.M() < 0 or h_gg.M() > 250) {
	    continue;
	}
	if (z_ll.Pt() < 250 or h_gg.Pt() < 250) {
	    continue;
	}
	
	Double_t lpairmass = z_ll.Mag();
	Double_t gpairmass = h_gg.Mag();

	Int_t nsmallJet = branchsmallJet->GetEntriesFast();
	bool smallBtag = false;
	Int_t nsmallB = 0;
	Jet *smalljet;
	for (int j=0; j < nsmallJet; ++j) {
	    smalljet = (Jet *) branchsmallJet -> At(j);
	    if (h_gg.DeltaR(smalljet -> P4()) < 1.0) {
		nsmallB++;
		if (smalljet -> BTag >= 1) {
		    smallBtag = true;
		}
	    }
	}
	
	if (smallBtag == true) {
	    continue;
	}
	/*
	if (TMath::Abs(h_gg.DeltaR(z_ll) - TMath::Pi())>0.5) {
	    continue;
	}
	*/
	//std::cout << "# of small jet inside" << nsmallB << std::endl;
        /*	
	if (hpt < 400) {
	    continue;
	}
	*/
	/*
	if (gpairmass > 85){
	    continue;
	}
	*/
	GenParticle *trueHiggs;
	GenParticle *genpar;
	/*
	for (int gentry=0; gentry < nGenPar; ++gentry){
	    genpar = (GenParticle *) branchGenPar -> At(gentry);
	    if (genpar -> PID == 25){
		GenParticle *D1;
	        GenParticle *D2;
		trueHiggs = (GenParticle *) branchGenPar -> At(gentry);
		D1 = (GenParticle *) branchGenPar -> At(trueHiggs->D1);
		D2 = (GenParticle *) branchGenPar -> At(trueHiggs->D2);
		if (D1 -> PID == 21 and D2 -> PID == 21) {
		    TLorentzVector HiggsP4 = trueHiggs -> P4();
		    if (h_gg.DeltaR(D1 -> P4()) > 1.0) {
			//std::cout << "Find D1 have distance " << h_gg.DeltaR(D1 -> P4()) << " from jet."<< std::endl;
		    }
		    if (h_gg.DeltaR(D2 -> P4()) > 1.0) {
			//std::cout << "Find D2 have distance " << h_gg.DeltaR(D2 -> P4()) << " from jet."<< std::endl;
		    }
		    if (h_gg.DeltaR(D1 -> P4()) <= 1.0 and h_gg.DeltaR(D2 -> P4()) <= 1.0) {
		        std::cout << "find higgs near jet with mass " << gpairmass << std::endl;
			find++;
			higgstag = true;
			//truejethisto -> Fill(h_gg.DeltaR(HiggsP4));
			break;
		    }
		}
	    }
	}
	*/
	//std::cout << "z_ll mass: " << lpairmass << "; h_gg mass: " << gpairmass<< std::endl;
	
	/*
	std::vector<TLorentzVector> particlesArr;	
	for(int i = 0; i < hjet->Constituents.GetEntriesFast(); ++i) {
	    object = hjet->Constituents.At(i);
	    if(object == 0) continue;
	    if(object->IsA() == Tower::Class()) {
	        Tower *tower = (Tower*) object;
		TLorentzVector p4 = tower->P4();
		particlesArr.push_back(p4); 
	    }
	    if(object->IsA() == Track::Class()) {
	        Track *track = (Track*) object;
	        TLorentzVector p4 = track->P4();
	        particlesArr.push_back(p4);
	    }
	}

	double thrust = calculateThrust(particlesArr);
	//std::cout << "Thrust: " << thrust << std::endl;
	if (thrust < 0.8) {
	    //continue;
	}
	total++;
	*/
	/*
	if (higgstag == false) {
	    continue;
	}
	*/
	//std::cout << "event number " << entry << std::endl;
	//Zmass -> Fill(lpairmass);
	Hmass -> Fill(h_gg.Mag());
	//tau21 -> Fill(htau21);
	//tau32 -> Fill(htau32);
	//nchargehisto -> Fill(hNCharged);
	//nneutralhisto -> Fill(hNNeutral);
	//pthisto -> Fill(hpt);
	//ChargeFrachisto -> Fill(hChargeFrac);
	//btaghisto -> Fill(hjet->BTag);
	//zhdishisto -> Fill(h_gg.DeltaR(z_ll));
	//std::cout << nJet << std::endl;
	//std::cout <<  z_ll.DeltaR(h_gg) << std::endl;
    }
    //std::cout << "find total of " << total << " event." << std::endl;
    //std::cout << "find " << find << " event match wich higgs to gluglu." << std::endl;
    Int_t totalnum = Hmass -> GetEntries();
    std::cout << "total number of event pass requirement is " << totalnum << std::endl;
    std::cout << "total is " << total << std::endl;

}
void llgg_stack(const char *sigFile, const char *bkgzjjFile, const char *bkgzhbbFile, const char *bkgzhccFile, const char *bkgzhwwFile) {
    TFile *output = new TFile("ZH_analysis.root", "RECREATE");
    TTree *tree_output = new TTree("tree_output", "Delphes");
	int histmin = 0;
	int histmax = 250;
    TH1D *sigHmass = new TH1D("Hmass", "jet pair invariant mass", 30, histmin, histmax);
    TH1D *bkgzjjHmass = new TH1D("Hmass", "jet pair invariant mass", 30, histmin, histmax);
    TH1D *bkgzhbbHmass = new TH1D("Hmass", "jet pair invariant mass", 30, histmin, histmax);
    TH1D *bkgzhccHmass = new TH1D("Hmass", "jet pair invariant mass", 30, histmin, histmax);
    TH1D *bkgzhwwHmass = new TH1D("Hmass", "jet pair invariant mass", 30, histmin, histmax);
    llgg(sigFile, sigHmass);
    llgg(bkgzjjFile, bkgzjjHmass);
    llgg(bkgzhbbFile, bkgzhbbHmass);
    llgg(bkgzhccFile, bkgzhccHmass);
    llgg(bkgzhwwFile, bkgzhwwHmass);

    THStack* Hmass = new THStack("Hmass", ";Reconstructed Higgs Invariant Mass [GeV]; Events");

    TH1* sigHmass_clone = new TH1D(*sigHmass);
    
    Double_t sig_w = 0.0093396; 
    Double_t sig_w_enh = 1000 * sig_w;
    Double_t bkgzjj_w = 556.6953864;
    Double_t bkgzhbb_w = 0.06324864;
    Double_t bkgzhcc_w = 0.003139626;
    Double_t bkgzhww_w = 0.01182654;

    Int_t sig = sigHmass -> GetEntries();
    Int_t bkgzjj = bkgzjjHmass -> GetEntries();
    Int_t bkgzhbb = bkgzhbbHmass -> GetEntries();
    Int_t bkgzhcc = bkgzhccHmass -> GetEntries();
    Int_t bkgzhww = bkgzhwwHmass -> GetEntries();
    Double_t significance = sig*sig_w/TMath::Sqrt(bkgzjj*bkgzjj_w + bkgzhbb*bkgzhbb_w+ bkgzhcc*bkgzhcc_w + bkgzhww*bkgzhww_w);
    std::cout << "signal expected" << sig*sig_w << std::endl; 
    std::cout << "zjj expected" << bkgzjj*bkgzjj_w << std::endl; 
    std::cout << "zhbb expected" << bkgzhbb*bkgzhbb_w << std::endl; 
    std::cout << "zhcc expected" << bkgzhcc*bkgzhcc_w << std::endl; 
    std::cout << "zhww expected" << bkgzhww*bkgzhww_w << std::endl; 
    std::cout << "Expected significance is" << significance << std::endl;
    sigHmass_clone -> SetLineColor(kPink-1);
    sigHmass_clone -> Scale(sig_w);
    sigHmass_clone -> SetLineWidth(4);
    
    bkgzhccHmass -> SetFillColor(kBlue+3);
    bkgzhccHmass -> Scale(bkgzhcc_w);
    Hmass -> Add(bkgzhccHmass);
    bkgzhbbHmass -> SetFillColor(kAzure);
    bkgzhbbHmass -> Scale(bkgzhbb_w);
    Hmass -> Add(bkgzhbbHmass);
    bkgzhwwHmass -> SetFillColor(kAzure+7);
    bkgzhwwHmass -> Scale(bkgzhww_w);
    Hmass -> Add(bkgzhwwHmass);
    bkgzjjHmass -> SetFillColor(kAzure+10);
    bkgzjjHmass -> Scale(bkgzjj_w);
    Hmass -> Add(bkgzjjHmass);
    sigHmass -> SetFillColor(kPink-1);
    sigHmass -> Scale(sig_w_enh);
    Hmass -> Add(sigHmass);

    Hmass -> SetMaximum(1e6);
    Hmass -> SetMinimum(10);
    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    totalcanvas -> SetLogy();
    TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend -> AddEntry(sigHmass, "signal #times 1000", "f");
    legend -> AddEntry(sigHmass_clone, "signal", "l");
    legend -> AddEntry(bkgzjjHmass, "bkg Z+jj", "f");
    legend -> AddEntry(bkgzhwwHmass, "bkg Z+ww", "f");
    legend -> AddEntry(bkgzhbbHmass, "bkg Z+bb", "f");
    legend -> AddEntry(bkgzhccHmass, "bkg Z+cc", "f");
    Hmass -> Draw("hist");
    sigHmass_clone -> Draw("same");
    legend -> Draw("same");
    totalcanvas -> SaveAs("ZH_stack.png");

    output -> cd();
    sigHmass -> Write();
    bkgzjjHmass -> Write();
    bkgzhbbHmass -> Write();
    bkgzhccHmass -> Write();
    bkgzhwwHmass -> Write();
    Hmass -> Write();
    tree_output -> Write();
    output -> Close();
}
