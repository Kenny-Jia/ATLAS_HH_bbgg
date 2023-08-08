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
	
void truth_analysis(const char *inputFile) {
    gSystem -> Load("/eos/user/h/hjia/ATLAS_HH_bbgg/Delphes-3.5.0/libDelphes.so");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;
    TFile *file_sig = new TFile(inputFile);
    TTree *tree_sig = (TTree*) file_sig -> Get("Delphes");

    TLeaf *Par_size = tree_sig -> GetLeaf("Particle_size");
    TLeaf *Par_eta = tree_sig -> GetLeaf("Particle.Eta");
    TLeaf *Par_phi = tree_sig -> GetLeaf("Particle.Phi");
    TLeaf *Par_pt = tree_sig -> GetLeaf("Particle.PT");
    TLeaf *Par_M1 = tree_sig -> GetLeaf("Particle.M1");
    TLeaf *Par_pid = tree_sig -> GetLeaf("Particle.PID");

    Int_t nEntries = tree_sig -> GetEntries();

    TH1D *truthHbbmass = new TH1D("truthHbbmass", "particle pair invariant mass", 30, 124.7, 125.1);
    TH1D *truthHbbR = new TH1D("truthHbbR", "particle pair deltaR", 30, 0, 5);
    TH1D *truthHggmass = new TH1D("truthHggmass", "particle pair invariant mass", 30, 124.7, 125.1);
    TH1D *truthHggR = new TH1D("truthHggR", "particle pair deltaR", 30, 0, 5);
    TH1D *truthHHmass = new TH1D("truthHHmass", "particle pair invariant mass", 30, 0, 1500);

    for (Long64_t entry = 0; entry < nEntries; entry++) {
	Par_size -> GetBranch() -> GetEntry(entry);
	Par_eta -> GetBranch() -> GetEntry(entry);
	Par_phi -> GetBranch() -> GetEntry(entry);
	Par_pt -> GetBranch() -> GetEntry(entry);
	Par_M1 -> GetBranch() -> GetEntry(entry);
	Par_pid -> GetBranch() -> GetEntry(entry);

	Int_t nPar = Par_size -> GetValue();

	bool b1find = false;
	bool g1find = false;
	Double_t b1pt;
	Double_t b1eta;
	Double_t b1phi;

	Double_t b2pt;
	Double_t b2eta;
	Double_t b2phi;

	Double_t g1pt;
	Double_t g1eta;
	Double_t g1phi;

	Double_t g2pt;
	Double_t g2eta;
	Double_t g2phi;

	TLorentzVector b1;
	TLorentzVector b2;
	TLorentzVector g1;
	TLorentzVector g2;
	
	TLorentzVector h_bb;
	TLorentzVector h_gg;
	TLorentzVector hh;
	h_bb.SetPtEtaPhiM(0, 0, 0, 0);
	h_gg.SetPtEtaPhiM(0, 0, 0, 0);

	for (Int_t parentry = 0; parentry < nPar; parentry++) {
	    if ((TMath::Abs(Par_pid -> GetValue(parentry)) != 5) and (TMath::Abs(Par_pid -> GetValue(parentry)) != 21)) {
		continue;
	    }
	    Int_t mother = Par_M1 -> GetValue(parentry);
	    if (Par_pid -> GetValue(mother) == 25) {
		if (TMath::Abs(Par_pid -> GetValue(parentry)) == 5) {
		    if (b1find == false) {
		        b1eta = Par_eta -> GetValue(parentry);
		        b1phi = Par_phi -> GetValue(parentry);
		        b1pt = Par_pt -> GetValue(parentry);
			b1.SetPtEtaPhiM(b1pt, b1eta, b1phi, 4.18);
			b1find = true;
		    } else {
			b2eta = Par_eta -> GetValue(parentry);
		        b2phi = Par_phi -> GetValue(parentry);
		        b2pt = Par_pt -> GetValue(parentry);
			b2.SetPtEtaPhiM(b2pt, b2eta, b2phi, 4.18);
		    }
		} else {
		    if (g1find == false) {
		        g1eta = Par_eta -> GetValue(parentry);
		        g1phi = Par_phi -> GetValue(parentry);
		        g1pt = Par_pt -> GetValue(parentry);
			g1.SetPtEtaPhiM(g1pt, g1eta, g1phi, 0);
			g1find = true;
		    } else {
			g2eta = Par_eta -> GetValue(parentry);
		        g2phi = Par_phi -> GetValue(parentry);
		        g2pt = Par_pt -> GetValue(parentry);
			g2.SetPtEtaPhiM(g2pt, g2eta, g2phi, 0);
		    }
		}
	    }
    	}
	h_bb = b1 + b2;
	h_gg = g1 + g2;

	hh = h_bb + h_gg;
	truthHHmass -> Fill(hh.Mag());
	truthHbbmass -> Fill(h_bb.Mag());
	truthHggmass -> Fill(h_gg.Mag());
	truthHbbR -> Fill(b1.DeltaR(b2));
	truthHggR -> Fill(g1.DeltaR(g2));

	std::cout << "event " << entry << " truth mass: " << h_bb.Mag() << ", " << h_gg.Mag() << std::endl;
    }
    TCanvas *truthtotalcanvas = new TCanvas("truthtotalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    truthtotalcanvas -> cd();
    truthtotalcanvas -> SetWindowSize(1204,1228);
    truthtotalcanvas -> SetCanvasSize(1200,1200);
    truthHbbmass -> SetMaximum(100000);
    truthHggmass -> SetMaximum(100000);
    truthHbbmass -> SetLineColor(kPink);
    truthHggmass -> SetLineColor(kBlue);
    TLegend *truthlegend = new TLegend(0.2, 0.65, 0.4, 0.85);
    truthlegend -> AddEntry(truthHbbmass, "b pair mass", "l");
    truthlegend -> AddEntry(truthHggmass, "gluon pair mass", "l");
    truthHbbmass -> Draw("HIST");
    truthHggmass -> Draw("same");
    truthlegend -> Draw("same");
    truthtotalcanvas -> SaveAs("ParticlePairs_all.png");

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> cd();
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    truthHbbR -> SetMaximum(20000);
    truthHggR -> SetMaximum(20000);
    truthHbbR -> SetLineColor(kPink);
    truthHggR -> SetLineColor(kBlue);
    TLegend *legend = new TLegend(0.2, 0.65, 0.4, 0.85);
    legend -> AddEntry(truthHbbR, "b pair delta R", "l");
    legend -> AddEntry(truthHggR, "gluon pair delta R", "l");
    truthHbbR-> Draw("HIST");
    truthHggR -> Draw("same");
    legend -> Draw("same");
    totalcanvas -> SaveAs("ParticlePairs_deltaR_all.png");

    TCanvas *fourcanvas = new TCanvas("fourcanvas", "Canvas", 1400, 1400, 1400, 1400);
    fourcanvas -> cd();
    fourcanvas -> SetWindowSize(1204,1228);
    fourcanvas -> SetCanvasSize(1200,1200);
    truthHHmass -> Draw("same");
    fourcanvas -> SaveAs("ParticlePairs_HH_all.png");





    file_sig -> Close();
}
