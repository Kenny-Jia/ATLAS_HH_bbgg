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

TLorentzVector HiggsBTag(Int_t nPar, TLeaf *Par_eta, TLeaf *Par_phi, TLeaf *Par_pid, TLeaf *Par_M1, TLeaf *Par_pt, Double_t eta, Double_t phi, Double_t pt) {
    Double_t minDeltaR = 0.4;
    TLorentzVector truthbjet;
    truthbjet.SetPtEtaPhiM(0, 0, 0, 1000);
    for (Int_t parentry = 0; parentry < nPar; parentry++) {
	if (TMath::Abs(Par_pid -> GetValue(parentry)) != 5) {
	    continue;
	}
	Double_t beta = Par_eta -> GetValue(parentry);
	Double_t bphi = Par_phi -> GetValue(parentry);
	Double_t bpt = Par_pt -> GetValue(parentry);
	Double_t deltaR = TMath::Sqrt((beta - eta) * (beta - eta) + (bphi - phi) * (bphi - phi));
	if (deltaR < minDeltaR) {
	    Int_t mother = Par_M1 -> GetValue(parentry);
	    if (Par_pid -> GetValue(mother) == 25) {
		truthbjet.SetPtEtaPhiM(bpt, beta, bphi, 4.18);
		//std::cout << "find b:" << bpt/pt << std::endl;
		break;
	    }
	}
    }
    return truthbjet;
}

TLorentzVector HiggsGluTag(Int_t nPar, TLeaf *Par_eta, TLeaf *Par_phi, TLeaf *Par_pid, TLeaf *Par_M1, TLeaf *Par_pt, Double_t eta, Double_t phi, Double_t pt) {
    Double_t minDeltaR = 0.4;
    TLorentzVector truthgjet;
    truthgjet.SetPtEtaPhiM(0, 0, 0, 10000);
    for (Int_t parentry = 0; parentry < nPar; parentry++) {
	if (Par_pid -> GetValue(parentry) != 21) {
	    continue;
	}
	Double_t glueta = Par_eta -> GetValue(parentry);
	Double_t gluphi = Par_phi -> GetValue(parentry);
	Double_t glupt = Par_pt -> GetValue(parentry);
	Double_t deltaR = TMath::Sqrt((glueta - eta) * (glueta - eta) + (gluphi - phi) * (gluphi - phi));
	if (deltaR < minDeltaR) {
	    Int_t mother = Par_M1 -> GetValue(parentry);
	    if (Par_pid -> GetValue(mother) == 25) {
		truthgjet.SetPtEtaPhiM(glupt, glueta, gluphi, 0);
		//std::cout << "find g!" << std::endl;
		//std::cout << "find glu:" << glupt/pt << std::endl;
		break;
	    }
	}
    }
    return truthgjet;
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
    TLeaf *Par_pt = tree_sig -> GetLeaf("Particle.PT");
    TLeaf *Par_M1 = tree_sig -> GetLeaf("Particle.M1");
    TLeaf *Par_pid = tree_sig -> GetLeaf("Particle.PID");

    Int_t nEntries = tree_sig -> GetEntries();

    TH1D *Hbbmass = new TH1D("Hbbmass", "jet pair invariant mass", 30, 0, 250);
    TH1D *truthHbbmass = new TH1D("truthHbbmass", "particle pair invariant mass", 30, 124.7, 125.1);
    TH1D *Hggmass = new TH1D("Hggmass", "jet pair invariant mass", 30, 0, 250);
    TH1D *truthHggmass = new TH1D("truthHggmass", "particle pair invariant mass", 30, 124.7, 125.1);

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
	Par_pt -> GetBranch() -> GetEntry(entry);
	Par_M1 -> GetBranch() -> GetEntry(entry);
	Par_pid -> GetBranch() -> GetEntry(entry);

	Int_t nJet = Jet_size -> GetValue();
	Int_t nPar = Par_size -> GetValue();

	if (nJet < 4) {
	    lesscnt++;
	    continue;
	}

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
	h_bb.SetPtEtaPhiM(0, 0, 0, 0);
	TLorentzVector truthh_bb;
	truthh_bb.SetPtEtaPhiM(0, 0, 0, 0);

	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    Double_t eta = Jet_eta -> GetValue(jentry);
	    Double_t phi = Jet_phi -> GetValue(jentry);
	    Double_t bpt = Jet_pt -> GetValue(jentry);
	    if (HiggsBTag(nPar, Par_eta, Par_phi, Par_pid, Par_M1, Par_pt, eta, phi, bpt).M() != 1000) {
		TLorentzVector bjet;
	        TLorentzVector truthbjet = HiggsBTag(nPar, Par_eta, Par_phi, Par_pid, Par_M1, Par_pt, eta, phi, bpt);
		if (bpt < 20) {
		    continue;
		}
		Double_t bmass = Jet_mass -> GetValue(jentry);
		bjet.SetPtEtaPhiM(bpt, eta, phi, bmass);
		h_bb = h_bb + bjet;
		truthh_bb = truthh_bb + truthbjet;
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
	Double_t truthbpairmass = truthh_bb.Mag();

	TLorentzVector h_gg;
	h_gg.SetPtEtaPhiM(0, 0, 0, 0);
	TLorentzVector truthh_gg;
	truthh_gg.SetPtEtaPhiM(0, 0, 0, 0);
	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    if (Jet_btag -> GetValue(jentry) == 1) {
		continue;		
	    } 
	    Double_t eta = Jet_eta -> GetValue(jentry);
	    Double_t phi = Jet_phi -> GetValue(jentry);

	    Double_t gpt = Jet_pt -> GetValue(jentry);
	    if (HiggsGluTag(nPar, Par_eta, Par_phi, Par_pid, Par_M1, Par_pt, eta, phi, gpt).M() != 10000) {
		TLorentzVector gjet;
		TLorentzVector truthgjet = HiggsGluTag(nPar, Par_eta, Par_phi, Par_pid, Par_M1, Par_pt, eta, phi, gpt);
		if (gpt < 20) {
		    continue;
		}
		Double_t gmass = Jet_mass -> GetValue(jentry);
		gjet.SetPtEtaPhiM(gpt, eta, phi, gmass);
		h_gg = h_gg + gjet;
		truthh_gg = truthh_gg + truthgjet;
		glucnt++;
	    }
	}
	Double_t gpairmass = h_gg.Mag();
	Double_t truthgpairmass = truthh_gg.Mag();

	if (glucnt < 2) {
	    lessgcnt++;
	    continue;
	}
	if (glucnt > 2) {
	    moregcnt++;
	    continue;
	}
	if (abs(truthbpairmass - (4.18 + 4.18)) < 1e-2) { 
	    std::cout << "for event" << entry << " truth bpairmas wrong. Same parton matched" << std::endl;
	    continue;
	}
	if (truthgpairmass < 1e-2) { 
	    std::cout << "for event" << entry << " truth gpairmas wrong. Same parton matched" << std::endl;
	    continue;
	}

	Hbbmass -> Fill(bpairmass);
	truthHbbmass -> Fill(truthbpairmass);
	Hggmass -> Fill(gpairmass);
	truthHggmass -> Fill(truthgpairmass);
	std::cout << "event " << entry << " truth mass: " << truthbpairmass << ", " << truthgpairmass << std::endl;
    }

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    Hbbmass -> SetMaximum(4000);
    Hggmass -> SetMaximum(4000);
    Hbbmass -> SetLineColor(kPink);
    Hggmass -> SetLineColor(kBlue);
    TLegend *legend = new TLegend(0.15, 0.65, 0.35, 0.85);
    legend -> AddEntry(Hbbmass, "b pair mass", "l");
    legend -> AddEntry(Hggmass, "gluon pair mass", "l");
    Hbbmass -> Draw("HIST");
    Hggmass -> Draw("same");
    legend -> Draw("same");
    totalcanvas -> SaveAs("JetPairs.png");

    TCanvas *truthtotalcanvas = new TCanvas("truthtotalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    truthtotalcanvas -> cd();
    truthtotalcanvas -> SetWindowSize(1204,1228);
    truthtotalcanvas -> SetCanvasSize(1200,1200);
    truthHbbmass -> SetMaximum(24000);
    truthHggmass -> SetMaximum(24000);
    truthHbbmass -> SetLineColor(kPink);
    truthHggmass -> SetLineColor(kBlue);
    TLegend *truthlegend = new TLegend(0.2, 0.65, 0.4, 0.85);
    truthlegend -> AddEntry(truthHbbmass, "b pair mass", "l");
    truthlegend -> AddEntry(truthHggmass, "gluon pair mass", "l");
    truthHbbmass -> Draw("HIST");
    truthHggmass -> Draw("same");
    truthlegend -> Draw("same");
    truthtotalcanvas -> SaveAs("ParticlePairs.png");
    file_sig -> Close();
    std::cout << "less jet " << lesscnt << std::endl;
    std::cout << "less b " << lessbcnt << std::endl;
    std::cout << "more b " << morebcnt << std::endl;
    std::cout << "less g " << lessgcnt << std::endl;
    std::cout << "more g " << moregcnt << std::endl;
}
