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
#include <utility>
#endif
	
bool compareTVector3(const TVector3& a, const TVector3& b) {
    return a.Mag() > b.Mag();
}

std::pair<double, double> calculateThrust(const std::vector<TLorentzVector>& momenta) {
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
    if (n == 0) return std::make_pair(0.0, 0.0);
    if (n == 1) return std::make_pair(1.0, 0.0);

    assert(n >= 3);
    n = 3;
    if (momentaCp.size() == 3) n = 3;

    std::vector<TVector3> tvec;
    std::vector<double> tval;
    TVector3 taxis(0, 0, 0);
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
	    taxis = tvec[i];
	}
    }
    double tm;
    tm = 0.;
    for (unsigned int k = 0; k < momentaCp.size(); k++) {
        tm += fabs(taxis.Cross(momentaCp[k]).Mag());
    }
    return std::make_pair(t / norm, tm / norm);
}

void llgg(const char *inputFile) {

    std::cout << "loading Delphes library" << std::endl;
    gSystem -> Load("/eos/user/h/hjia/ATLAS_HH_bbgg/Delphes-3.5.0/libDelphes.so");
    //gSystem->Load("libDelphes");
    //gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;

    std::cout << "open up file" << std::endl;
    TChain chain("Delphes");
    chain.Add(inputFile);

    std::string input(inputFile);
    size_t pos = input.find("_250GeV_14TeV.root");
    std::string channelName;
    if (pos != std::string::npos) {
        channelName = input.substr(0, pos);
    }
    std::string prefix = "../../../../d/hjia625/ATLAS_H_gg/data/Delphes/";
    size_t posprefix = channelName.find(prefix);
    if (pos != std::string::npos) {
        channelName.erase(posprefix, prefix.length());
    }

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t nEntries = treeReader->GetEntries();

    TClonesArray *branchJet = treeReader->UseBranch("PFJet10");
    TClonesArray *branchsmallJet = treeReader->UseBranch("UniqueJet");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchGenPar = treeReader->UseBranch("Particle");
    TClonesArray *branchEvent = treeReader->UseBranch("Event");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");
    TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");

    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");


    TH1D *Zmass = new TH1D("Zmass", "lepton pair invariant mass", 100, 30, 220);
    TH1D *Hmass = new TH1D("Hmass", "jet pair invariant mass", 50, 30, 220);
    TH1D *tau21 = new TH1D("tau21", "higgs tau21", 30, 0, 1);
    TH1D *tau42 = new TH1D("tau32", "higgs tau42", 30, 0, 1);
    TH1D *nchargehisto = new TH1D("nchargehisto", "number of charge constituents", 30, 0, 70);
    TH1D *nneutralhisto = new TH1D("nneutralhisto", "number of neutral constituents", 30, 0, 70);
    TH1D *pthisto = new TH1D("pthisto", "higgs pt", 30, 0, 1000);
    TH1D *ChargeFrachisto = new TH1D("ChargeFrachisto", "higgs Charged Energy Fraction", 30, 0, 1);
    TH1D *btaghisto = new TH1D("btaghisto", "higgs jet btag", 4, 0, 2);
    TH1D *zhdishisto = new TH1D("zhdishisto", "z and higgs jet distance", 30, 0, 4);
    TH1D *truejethisto = new TH1D("truejethisto", "true higgs and higgs jet distance", 30, 0, 4);
    TH1D *METhisto = new TH1D("METhisto", "MissingET", 30, 0, 250);
    TH1D *METhDishisto = new TH1D("METhDishisto", "MissingET higgs Dis", 30, 0, 6.28);
    TH1D *nsmallhisto = new TH1D("nsmallhisto", "number of small jet", 4, 0, 4);
    TH1D *nsmallmasshisto = new TH1D("nsmallmasshisto", "number of small jet", 50, 0, 200);
    TH1D *thrusthisto = new TH1D("thrusthisto", "jet thrust", 100, 0.5, 1);
    TH1D *thrustmhisto = new TH1D("thrustmhisto", "jet thrust minor", 100, 0, 0.8);

    std::cout << "Start loop ever tree" << std::endl;

    Int_t total = 0;
    Int_t find = 0;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
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
	if (nMu == 2) {
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
	if (l1.DeltaR(l2) < 0.05) {
	    continue;
	}
	Double_t mass_min = 100000000;
	Jet *jet;
	Jet *hjet;
	Int_t hentry;
	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    jet = (Jet *) branchJet -> At(jentry);
	    
	    /*
	    if (jet -> BTag == 1) {
		continue;
	    }
	    */
	    if (TMath::Abs(jet -> Mass - 125) < mass_min) {
		mass_min = TMath::Abs(jet -> Mass - 125);
	        //hjet = (Jet *) branchJet -> At(jentry);
		hentry = jentry; 
	    }
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
		
	if (hpt < 250){
	    continue;
	}

	
	TLorentzVector z_ll;
	z_ll.SetPtEtaPhiM(0, 0, 0, 0);
	TLorentzVector h_gg;
	h_gg.SetPtEtaPhiM(0, 0, 0, 0);
	
	z_ll = l1 + l2;
	h_gg.SetPtEtaPhiM(hpt, heta, hphi, hmass);
        if (z_ll.DeltaR(h_gg) < 1.5) {
	    continue;
	}
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
	Double_t Summass = 0;
	Jet *smalljet;
	for (int j=0; j < nsmallJet; ++j) {
	    smalljet = (Jet *) branchsmallJet -> At(j);
	    nsmallB++;
	    Summass = smalljet -> Mass + Summass;
	    if (h_gg.DeltaR(smalljet -> P4()) < 1.0) {
		if (smalljet -> BTag >= 1) {
		    smallBtag = true;
		}
	    }
	}
	
	if (smallBtag == true) {
	    continue;
	}
	
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

	//true tag *********************************************
	
	MissingET *met;
	met = (MissingET *) branchMET -> At(0);
	GenParticle *trueHiggs;
	GenParticle *genpar;
	GenParticle *muon;
	bool higgstag = false;
	//bool muontag = false;
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
	    /*
	    if (abs(abs(genpar -> PID) - 13) == 1) {
		muon = (GenParticle *) branchGenPar -> At(gentry);
		if (h_gg.DeltaR(muon -> P4()) <= 1) {
		    muontag = true;
		}
	    }
	    */
	}
	
	/*
	if (muontag == false) {
	    continue;
	}
	*/
	/*
	if (hpt > 250) {
	    continue;
	}
	*/
	//std::cout << "z_ll mass: " << lpairmass << "; h_gg mass: " << gpairmass<< std::endl;
	
	
	if (higgstag == false) {
	    //continue;
	}
	
	/*
	if (hjet->BTag >= 1){
	    //continue;
	    std::cout << "weird:" << hbtag << " at event " << entry <<std::endl;
	}
	*/

// Uncomment to get jet constituents info
        std::vector<TLorentzVector> particlesArr;	
	for(int i = 0; i < hjet->Constituents.GetEntriesFast(); ++i) {
	    object = hjet->Constituents.At(i);
	    if(object == 0) continue;
	        if(object->IsA() == Tower::Class()) {
	  	    Tower *tower = (Tower*) object;
		    TLorentzVector p4 = tower->P4();
		    particlesArr.push_back(p4); 
	    	}
		/*
	    	if(object->IsA() == GenParticle::Class()) {
		    std::cout << "is genpar";
	    	}
	        */	
		if(object->IsA() == Track::Class()) {
		    Track *track = (Track*) object;
		    TLorentzVector p4 = track->P4();
		    particlesArr.push_back(p4);
		}
		/*
	    	if(object->IsA() == Candidate::Class()) {
		    std::cout << "is candidate";
	    	}
		*/
	}

	std::pair<double, double> thrust_result = calculateThrust(particlesArr);
	double thrust = thrust_result.first;
	double thrust_minor = thrust_result.second;
	std::cout << "Thrust: " << thrust << std::endl;
	std::cout << "Thrust_minor: " << thrust_minor << std::endl;

	total++;
	//std::cout << "event number " << entry << std::endl;
	Zmass -> Fill(lpairmass);
	Hmass -> Fill(h_gg.Mag());
	tau21 -> Fill(htau21);
	tau42 -> Fill(htau42);
	nchargehisto -> Fill(hNCharged);
	nneutralhisto -> Fill(hNNeutral);
	pthisto -> Fill(hpt);
	ChargeFrachisto -> Fill(hChargeFrac);
	btaghisto -> Fill(hjet->BTag);
	zhdishisto -> Fill(h_gg.DeltaR(z_ll));
	//METhisto -> Fill(met->MET);
	//METhDishisto -> Fill(h_gg.DeltaR(met->P4()));
	nsmallhisto -> Fill(nsmallB);
	nsmallmasshisto -> Fill(Summass/nsmallB);
	thrusthisto -> Fill(thrust);
	thrustmhisto -> Fill(thrust_minor);
	std::cout << "event number" << entry << " hmass = " << h_gg.Mag() << std::endl;
	//std::cout << nJet << std::endl;
	//std::cout <<  z_ll.DeltaR(h_gg) << std::endl;
    }
    //std::cout << "find " << find << " event match wich higgs to gluglu." << std::endl;
    Int_t totalnum = Hmass -> GetEntries();
    std::cout << "find total of " << totalnum << " event." << std::endl;

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    //Hmass -> SetMaximum(9000);
    //Zmass -> SetMaximum(9000);
    Zmass -> SetLineColor(kPink);
    Hmass -> SetLineColor(kBlue);
    TLegend *legend = new TLegend(0.15, 0.65, 0.35, 0.85);
    legend -> AddEntry(Zmass, "lepton pair mass", "l");
    legend -> AddEntry(Hmass, "gluon pair mass", "l");
    Zmass -> Draw("HIST");
    Hmass -> Draw("same");
    legend -> Draw("same");
    totalcanvas -> SaveAs((channelName + "_anti.png").c_str());

    TCanvas *tau21canvas = new TCanvas("21canvas", "Canvas", 1400, 1400, 1400, 1400);
    tau21canvas -> SetWindowSize(1204,1228);
    tau21canvas -> SetCanvasSize(1200,1200);
    tau21 -> Draw("HIST");
    tau21canvas -> SaveAs((channelName + "_anti_tau21.png").c_str());

    TCanvas *tau42canvas = new TCanvas("42canvas", "Canvas", 1400, 1400, 1400, 1400);
    tau42canvas -> SetWindowSize(1204,1228);
    tau42canvas -> SetCanvasSize(1200,1200);
    tau42 -> Draw("HIST");
    tau42canvas -> SaveAs((channelName + "_anti_tau42.png").c_str());

    TCanvas *chargecanvas = new TCanvas("chargecanvas", "Canvas", 1400, 1400, 1400, 1400);
    chargecanvas -> SetWindowSize(1204,1228);
    chargecanvas -> SetCanvasSize(1200,1200);
    nchargehisto -> Draw("HIST");
    chargecanvas -> SaveAs((channelName + "_anti_ncharge.png").c_str());

    TCanvas *neutralcanvas = new TCanvas("neutralcanvas", "Canvas", 1400, 1400, 1400, 1400);
    neutralcanvas -> SetWindowSize(1204,1228);
    neutralcanvas -> SetCanvasSize(1200,1200);
    nneutralhisto -> Draw("HIST");
    neutralcanvas -> SaveAs((channelName + "_anti_nneutral.png").c_str());

    TCanvas *ptcanvas = new TCanvas("ptcanvas", "Canvas", 1400, 1400, 1400, 1400);
    ptcanvas -> SetWindowSize(1204,1228);
    ptcanvas -> SetCanvasSize(1200,1200);
    pthisto -> Draw("HIST");
    ptcanvas -> SaveAs((channelName + "_anti_pt.png").c_str());
/*
    TCanvas *metcanvas = new TCanvas("metcanvas", "Canvas", 1400, 1400, 1400, 1400);
    metcanvas -> SetWindowSize(1204,1228);
    metcanvas -> SetCanvasSize(1200,1200);
    METhisto -> Draw("HIST");
    metcanvas -> SaveAs((channelName + "_anti_met.png").c_str());
*/
    TCanvas *ChargeFraccanvas = new TCanvas("ChargeFraccanvas", "Canvas", 1400, 1400, 1400, 1400);
    ChargeFraccanvas -> SetWindowSize(1204,1228);
    ChargeFraccanvas -> SetCanvasSize(1200,1200);
    ChargeFrachisto -> Draw("HIST");
    ChargeFraccanvas -> SaveAs((channelName + "_anti_ChargeFrac.png").c_str());

    TCanvas *btagcanvas = new TCanvas("btagcanvas", "Canvas", 1400, 1400, 1400, 1400);
    btagcanvas -> SetWindowSize(1204,1228);
    btagcanvas -> SetCanvasSize(1200,1200);
    btaghisto -> Draw("HIST");
    btagcanvas -> SaveAs((channelName + "_anti_btag.png").c_str());

    TCanvas *zhdiscanvas = new TCanvas("zhdiscanvas", "Canvas", 1400, 1400, 1400, 1400);
    zhdiscanvas -> SetWindowSize(1204,1228);
    zhdiscanvas -> SetCanvasSize(1200,1200);
    zhdishisto -> Draw("HIST");
    zhdiscanvas -> SaveAs((channelName + "_anti_zhdis.png").c_str());

/*
    TCanvas *METhDiscanvas = new TCanvas("METhDiscanvas", "Canvas", 1400, 1400, 1400, 1400);
    METhDiscanvas -> SetWindowSize(1204,1228);
    METhDiscanvas -> SetCanvasSize(1200,1200);
    METhDishisto -> Draw("HIST");
    METhDiscanvas -> SaveAs((channelName + "_anti_METhDis.png").c_str());

    TCanvas *truejetdiscanvas = new TCanvas("truejetdiscanvas", "Canvas", 1400, 1400, 1400, 1400);
    truejetdiscanvas -> SetWindowSize(1204,1228);
    truejetdiscanvas -> SetCanvasSize(1200,1200);
    truejethisto -> Draw("HIST");
    truejetdiscanvas -> SaveAs((channelName + "_anti_truejet.png").c_str());
*/    
    TCanvas *nsmallcanvas = new TCanvas("nsmallcanvas", "Canvas", 1400, 1400, 1400, 1400);
    nsmallcanvas -> SetWindowSize(1204,1228);
    nsmallcanvas -> SetCanvasSize(1200,1200);
    nsmallhisto -> Draw("HIST");
    nsmallcanvas -> SaveAs((channelName + "_anti_nsmall.png").c_str());

    TCanvas *nsmallmasscanvas = new TCanvas("nsmallmasscanvas", "Canvas", 1400, 1400, 1400, 1400);
    nsmallmasscanvas -> SetWindowSize(1204,1228);
    nsmallmasscanvas -> SetCanvasSize(1200,1200);
    nsmallmasshisto -> Draw("HIST");
    nsmallmasscanvas -> SaveAs((channelName + "_anti_nsmallmass.png").c_str());

    TCanvas *thrustcanvas = new TCanvas("thrustcanvas", "Canvas", 1400, 1400, 1400, 1400);
    thrustcanvas -> SetWindowSize(1204,1228);
    thrustcanvas -> SetCanvasSize(1200,1200);
    thrusthisto -> Draw("HIST");
    thrustcanvas -> SaveAs((channelName + "_anti_thrust.png").c_str());

    TCanvas *thrustmcanvas = new TCanvas("thrustmcanvas", "Canvas", 1400, 1400, 1400, 1400);
    thrustmcanvas -> SetWindowSize(1204,1228);
    thrustmcanvas -> SetCanvasSize(1200,1200);
    thrustmhisto -> Draw("HIST");
    thrustmcanvas -> SaveAs((channelName + "_anti_thrust_minor.png").c_str());
}
