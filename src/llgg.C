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
#endif
	
void llgg(const char *inputFile) {

    std::cout << "loading Delphes library" << std::endl;
    gSystem -> Load("/eos/user/h/hjia/ATLAS_HH_bbgg/Delphes-3.5.0/libDelphes.so");
    //gSystem->Load("libDelphes");
    //gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;

    std::cout << "open up file" << std::endl;
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t nEntries = treeReader->GetEntries();

    TClonesArray *branchJet = treeReader->UseBranch("PFJet10");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchGenPar = treeReader->UseBranch("Particle");
    TClonesArray *branchEvent = treeReader->UseBranch("Event");

    TH1D *Zmass = new TH1D("Zmass", "lepton pair invariant mass", 100, 30, 220);
    TH1D *Hmass = new TH1D("Hmass", "jet pair invariant mass", 30, 30, 220);
    TH1D *tau21 = new TH1D("tau21", "higgs tau21", 30, 0, 1);
    //TH1D *tau32 = new TH1D("tau32", "higgs tau32", 30, 0, 1);
    TH1D *nchargehisto = new TH1D("nchargehisto", "number of charge constituents", 30, 0, 70);
    TH1D *nneutralhisto = new TH1D("nneutralhisto", "number of neutral constituents", 30, 0, 70);
    TH1D *pthisto = new TH1D("pthisto", "higgs pt", 30, 0, 1000);
    TH1D *ChargeFrachisto = new TH1D("ChargeFrachisto", "higgs Charged Energy Fraction", 30, 0, 1);
    TH1D *btaghisto = new TH1D("btaghisto", "higgs jet btag", 2, 0, 1);
    TH1D *zhdishisto = new TH1D("zhdishisto", "z and higgs jet distance", 30, 0, 4);
    TH1D *truejethisto = new TH1D("truejethisto", "true higgs and higgs jet distance", 30, 0, 4);

    std::cout << "Start loop ever tree" << std::endl;

    Int_t total = 0;
    Int_t find = 0;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
	bool higgstag = false;
	treeReader->ReadEntry(entry);

	Int_t nJet = branchJet->GetEntries();
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
	Double_t htau32;
	Int_t hNCharged;
	Int_t hNNeutral;
	Double_t hChargeFrac;
	Int_t hbtag;

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
	Jet *hjet;
	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    hjet = (Jet *) branchJet -> At(jentry);
	    /*
	    if (hjet -> BTag == 1) {
		continue;
	    }
	    */
	    if (TMath::Abs(hjet -> Mass - 125) < mass_min) {
		mass_min = TMath::Abs(hjet -> Mass - 125);
		hmass = hjet -> Mass;
		hpt = hjet -> PT;
		heta = hjet -> Eta;
		hphi = hjet -> Phi;
		htau21 = (hjet -> Tau[1])/(hjet -> Tau[0]);
		//htau32 = (hjet -> Tau[2])/(hjet -> Tau[1]);
		hNNeutral = hjet -> NNeutrals;
		hNCharged = hjet -> NCharged;
		hChargeFrac = hjet -> ChargedEnergyFraction;
		hbtag = hjet -> BTag;
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
	    }
	}
	TLorentzVector z_ll;
	z_ll.SetPtEtaPhiM(0, 0, 0, 0);
	TLorentzVector h_gg;
	h_gg.SetPtEtaPhiM(0, 0, 0, 0);
	
	z_ll = l1 + l2;
	h_gg.SetPtEtaPhiM(hpt, heta, hphi, hmass);
        if (z_ll.DeltaR(h_gg) < 1) {
	    continue;
	}
	if (z_ll.M() < 76 or z_ll.M() > 106) {
	    continue;
	}
	
	Double_t lpairmass = z_ll.Mag();
	Double_t gpairmass = h_gg.Mag();
	
	/*
	if (gpairmass >= 85){
	    continue;
	}
	*/
	GenParticle *trueHiggs;
	GenParticle *genpar;
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
			truejethisto -> Fill(h_gg.DeltaR(HiggsP4));
			break;
		    }
		}
	    }
	}
	//std::cout << "z_ll mass: " << lpairmass << "; h_gg mass: " << gpairmass<< std::endl;
	total++;
	/*
	if (higgstag == false) {
	    continue;
	}
	*/
	Zmass -> Fill(lpairmass);
	Hmass -> Fill(h_gg.Mag());
	tau21 -> Fill(htau21);
	//tau32 -> Fill(htau32);
	nchargehisto -> Fill(hNCharged);
	nneutralhisto -> Fill(hNNeutral);
	pthisto -> Fill(hpt);
	ChargeFrachisto -> Fill(hChargeFrac);
	btaghisto -> Fill(hbtag);
	zhdishisto -> Fill(h_gg.DeltaR(z_ll));
	//std::cout << nJet << std::endl;
	//std::cout <<  z_ll.DeltaR(h_gg) << std::endl;
    }
    std::cout << "find total of " << total << " event which have bad mass reconstruction." << std::endl;
    std::cout << "find " << find << " event which have bad mass reconstruction but do match wich higgs to gluglu." << std::endl;
    Int_t totalnum = Hmass -> GetEntries();

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    //Hmass -> SetMaximum(200);
    //Zmass -> SetMaximum(200)
    Zmass -> SetLineColor(kPink);
    Hmass -> SetLineColor(kBlue);
    TLegend *legend = new TLegend(0.15, 0.65, 0.35, 0.85);
    legend -> AddEntry(Zmass, "lepton pair mass", "l");
    legend -> AddEntry(Hmass, "gluon pair mass", "l");
    Zmass -> Draw("HIST");
    Hmass -> Draw("same");
    legend -> Draw("same");
    totalcanvas -> SaveAs("ZH_sig.png");

    TCanvas *tau21canvas = new TCanvas("21canvas", "Canvas", 1400, 1400, 1400, 1400);
    tau21canvas -> SetWindowSize(1204,1228);
    tau21canvas -> SetCanvasSize(1200,1200);
    tau21 -> Draw("HIST");
    tau21canvas -> SaveAs("ZH_sig_tau21.png");
/*
    TCanvas *tau32canvas = new TCanvas("32canvas", "Canvas", 1400, 1400, 1400, 1400);
    tau32canvas -> SetWindowSize(1204,1228);
    tau32canvas -> SetCanvasSize(1200,1200);
    tau32 -> Draw("HIST");
    tau32canvas -> SaveAs("ZH_gen_tau32.png");
*/
    TCanvas *chargecanvas = new TCanvas("chargecanvas", "Canvas", 1400, 1400, 1400, 1400);
    chargecanvas -> SetWindowSize(1204,1228);
    chargecanvas -> SetCanvasSize(1200,1200);
    nchargehisto -> Draw("HIST");
    chargecanvas -> SaveAs("ZH_sig_ncharge.png");

    TCanvas *neutralcanvas = new TCanvas("neutralcanvas", "Canvas", 1400, 1400, 1400, 1400);
    neutralcanvas -> SetWindowSize(1204,1228);
    neutralcanvas -> SetCanvasSize(1200,1200);
    nneutralhisto -> Draw("HIST");
    neutralcanvas -> SaveAs("ZH_sig_nneutral.png");

    TCanvas *ptcanvas = new TCanvas("ptcanvas", "Canvas", 1400, 1400, 1400, 1400);
    ptcanvas -> SetWindowSize(1204,1228);
    ptcanvas -> SetCanvasSize(1200,1200);
    pthisto -> Draw("HIST");
    ptcanvas -> SaveAs("ZH_sig_pt.png");

    TCanvas *ChargeFraccanvas = new TCanvas("ChargeFraccanvas", "Canvas", 1400, 1400, 1400, 1400);
    ChargeFraccanvas -> SetWindowSize(1204,1228);
    ChargeFraccanvas -> SetCanvasSize(1200,1200);
    ChargeFrachisto -> Draw("HIST");
    ChargeFraccanvas -> SaveAs("ZH_sig_ChargeFrac.png");

    TCanvas *btagcanvas = new TCanvas("btagcanvas", "Canvas", 1400, 1400, 1400, 1400);
    btagcanvas -> SetWindowSize(1204,1228);
    btagcanvas -> SetCanvasSize(1200,1200);
    btaghisto -> Draw("HIST");
    btagcanvas -> SaveAs("ZH_sig_btag.png");

    TCanvas *zhdiscanvas = new TCanvas("zhdiscanvas", "Canvas", 1400, 1400, 1400, 1400);
    zhdiscanvas -> SetWindowSize(1204,1228);
    zhdiscanvas -> SetCanvasSize(1200,1200);
    zhdishisto -> Draw("HIST");
    zhdiscanvas -> SaveAs("ZH_sig_zhdis.png");

    TCanvas *truejetdiscanvas = new TCanvas("truejetdiscanvas", "Canvas", 1400, 1400, 1400, 1400);
    truejetdiscanvas -> SetWindowSize(1204,1228);
    truejetdiscanvas -> SetCanvasSize(1200,1200);
    truejethisto -> Draw("HIST");
    truejetdiscanvas -> SaveAs("ZH_sig_truejet.png");
}