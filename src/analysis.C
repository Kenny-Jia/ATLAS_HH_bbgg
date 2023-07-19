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
#endif
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

bool HiggsBTag(Int_t nPar, TLeaf *Par_eta, TLeaf *Par_phi, TLeaf *Par_pid, TLeaf *Par_M1, Double_t eta, Double_t phi) {
    Int_t gluentry = 0;
    Double_t minDeltaR = 0.4;
    for (Int_t parentry = 0; parentry < nPar; parentry++) {
	if (TMath::Abs(Par_pid -> GetValue(parentry)) != 5) {
	    continue;
	}
	Double_t beta = Par_eta -> GetValue(parentry);
	Double_t bphi = Par_phi -> GetValue(parentry);
	Double_t deltaR = TMath::Sqrt((beta - eta) * (beta - eta) + (bphi - phi) * (bphi - phi));
	if (deltaR < minDeltaR) {
	    Int_t mother = Par_M1 -> GetValue(parentry);
	    if (Par_pid -> GetValue(mother) == 25) {
		return true;
	    }
	}
    }
    return false;
}
	
bool HiggsGluTag(Int_t nPar, TLeaf *Par_eta, TLeaf *Par_phi, TLeaf *Par_pid, TLeaf *Par_M1, Double_t eta, Double_t phi) {
    Int_t gluentry = 0;
    Double_t minDeltaR = 0.4;
    for (Int_t parentry = 0; parentry < nPar; parentry++) {
	if (Par_pid -> GetValue(parentry) != 21) {
	    continue;
	}
	Double_t glueta = Par_eta -> GetValue(parentry);
	Double_t gluphi = Par_phi -> GetValue(parentry);
	Double_t deltaR = TMath::Sqrt((glueta - eta) * (glueta - eta) + (gluphi - phi) * (gluphi - phi));
	if (deltaR < minDeltaR) {
	    Int_t mother = Par_M1 -> GetValue(parentry);
	    if (Par_pid -> GetValue(mother) == 25) {
		return true;
	    }
	}
    }
    return false;
}
	
void analysis(const char *inputFile) {
    gSystem -> Load("/eos/user/h/hjia/ATLAS_HH_bbgg/Delphes-3.5.0/libDelphes.so");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;
    TFile *file_sig = new TFile(inputFile);
    TTree *tree_sig = (TTree*) file_sig -> Get("Delphes");

    TLeaf *Jet_size = tree_sig -> GetLeaf("Jet_size");
    TLeaf *Jet_pt = tree_sig -> GetLeaf("Jet.PT");
    TLeaf *Jet_eta = tree_sig -> GetLeaf("Jet.Eta");
    TLeaf *Jet_phi = tree_sig -> GetLeaf("Jet.Phi");
    TLeaf *Jet_mass = tree_sig -> GetLeaf("Jet.Mass");
    TLeaf *Jet_btag = tree_sig -> GetLeaf("Jet.BTag");

    TLeaf *Par_size = tree_sig -> GetLeaf("Particle_size");
    TLeaf *Par_eta = tree_sig -> GetLeaf("Particle.Eta");
    TLeaf *Par_phi = tree_sig -> GetLeaf("Particle.Phi");
    TLeaf *Par_M1 = tree_sig -> GetLeaf("Particle.M1");
    TLeaf *Par_pid = tree_sig -> GetLeaf("Particle.PID");

    Int_t nEntries = tree_sig -> GetEntries();


    TH1D *Hbbmass = new TH1D("Hbbmass", "jet pair invariant mass", 30, 0, 250);
    TH1D *Hggmass = new TH1D("Hggmass", "jet pair invariant mass", 30, 0, 250);

    Int_t lesscnt = 0;
    Int_t lessbcnt = 0;
    Int_t lessgcnt = 0;
    Int_t morebcnt = 0;
    Int_t moregcnt = 0;
    for (Long64_t entry = 0; entry < nEntries; entry++) {

	Int_t glucnt = 0;
	Int_t bcnt = 0;
        tree_sig -> GetEntry(entry);
	Jet_size -> GetBranch() -> GetEntry(entry);
	Jet_eta -> GetBranch() -> GetEntry(entry);
	Jet_phi -> GetBranch() -> GetEntry(entry);
	Jet_pt -> GetBranch() -> GetEntry(entry);
	Jet_mass -> GetBranch() -> GetEntry(entry);
	Jet_btag -> GetBranch() -> GetEntry(entry);

	Par_size -> GetBranch() -> GetEntry(entry);
	Par_eta -> GetBranch() -> GetEntry(entry);
	Par_phi -> GetBranch() -> GetEntry(entry);
	Par_M1 -> GetBranch() -> GetEntry(entry);
	Par_pid -> GetBranch() -> GetEntry(entry);

	Int_t nJet = Jet_size -> GetValue();
	Int_t nPar = Par_size -> GetValue();

	if (nJet < 4) {
	    lesscnt++;
	    continue;
	}

	Jet_pt -> GetBranch() -> GetEntry(entry);
	Jet_eta -> GetBranch() -> GetEntry(entry);
	Jet_phi -> GetBranch() -> GetEntry(entry);
	Jet_btag -> GetBranch() -> GetEntry(entry);
	Jet_mass -> GetBranch() -> GetEntry(entry);
	
	Int_t nbJet = 0;

	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    if (Jet_btag -> GetValue(jentry) == 1) {
                nbJet++;
	    }
	}
	/*
	if (nbJet != 2) {
	    continue;
	}
        */

	bool b1find = false;
	Double_t b1pt;
	Double_t b1eta;
	Double_t b1phi;
	Double_t b1mass;

	Double_t b2pt;
	Double_t b2eta;
	Double_t b2phi;
	Double_t b2mass;

	TLorentzVector bjet1;
	TLorentzVector bjet2;
	/*
        for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    if (Jet_btag -> GetValue(jentry) == 1) {
		if (b1find == false){
		    b1pt = Jet_pt -> GetValue(jentry);
		    b1eta = Jet_eta -> GetValue(jentry);
		    b1phi = Jet_phi -> GetValue(jentry);
		    b1mass = Jet_mass -> GetValue(jentry);
		    bjet1.SetPtEtaPhiM(b1pt, b1eta, b1phi, b1mass);
		    b1find = true;
		} else {
		    b2pt = Jet_pt -> GetValue(jentry);
		    b2eta = Jet_eta -> GetValue(jentry);
		    b2phi = Jet_phi -> GetValue(jentry);
		    b2mass = Jet_mass -> GetValue(jentry);
		    bjet2.SetPtEtaPhiM(b2pt, b2eta, b2phi, b2mass);
		    break;
		}
	    }
	}
	*/

	//TLorentzVector h_bb = bjet1 + bjet2;
	TLorentzVector h_bb;

	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    Double_t eta = Jet_eta -> GetValue(jentry);
	    Double_t phi = Jet_phi -> GetValue(jentry);
	    if (HiggsBTag(nPar, Par_eta, Par_phi, Par_pid, Par_M1, eta, phi) == true) {
		TLorentzVector bjet;
		Double_t bpt = Jet_pt -> GetValue(jentry);
		if (bpt < 40) {
		    continue;
		}
		Double_t bmass = Jet_mass -> GetValue(jentry);
		bjet.SetPtEtaPhiM(bpt, eta, phi, bmass);
		h_bb = h_bb + bjet;
		bcnt++;
	    }
	}
	if (bcnt < 2) {
	    lessbcnt++;
	    continue;
	}
	if (bcnt > 2) {
	    morebcnt++;
	    continue;
	}
	Double_t bpairmass = h_bb.Mag();
	TLorentzVector h_gg;
	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    if (Jet_btag -> GetValue(jentry) == 1) {
		continue;		
	    } 
	    Double_t eta = Jet_eta -> GetValue(jentry);
	    Double_t phi = Jet_phi -> GetValue(jentry);
	    if (HiggsGluTag(nPar, Par_eta, Par_phi, Par_pid, Par_M1, eta, phi) == true) {
		TLorentzVector gjet;
		Double_t gpt = Jet_pt -> GetValue(jentry);
		if (gpt < 40) {
		    continue;
		}
		Double_t gmass = Jet_mass -> GetValue(jentry);
		gjet.SetPtEtaPhiM(gpt, eta, phi, gmass);
		h_gg = h_gg + gjet;
		glucnt++;
	    }
	}
	Double_t gpairmass = h_gg.Mag();

	if (glucnt < 2) {
	    lessgcnt++;
	    continue;
	}
	if (glucnt > 2) {
	    moregcnt++;
	}
	Hbbmass -> Fill(bpairmass);
	Hggmass -> Fill(gpairmass);
	std::cout << "event " << entry << " mass: " << bpairmass << ", " << gpairmass << std::endl;
    }

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    Hbbmass -> SetMaximum(2500);
    Hggmass -> SetMaximum(2500);
    Hbbmass -> SetLineColor(kPink);
    Hggmass -> SetLineColor(kBlue);
    TLegend *legend = new TLegend(0.6, 0.65, 0.85, 0.85);
    legend -> AddEntry(Hbbmass, "b pair mass", "l");
    legend -> AddEntry(Hggmass, "gluon pair mass", "l");
    Hbbmass -> Draw("HIST");
    Hggmass -> Draw("same");
    legend -> Draw("same");
    totalcanvas -> SaveAs("JetPairs.png");
    file_sig -> Close();
    std::cout << "less jet " << lesscnt << std::endl;
    std::cout << "less b " << lessbcnt << std::endl;
    std::cout << "more b " << morebcnt << std::endl;
    std::cout << "less g " << lessgcnt << std::endl;
    std::cout << "more g " << moregcnt << std::endl;
}
