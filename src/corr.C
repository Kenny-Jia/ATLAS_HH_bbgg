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

void corr(const char* sigFilename, const char* bkgFilename) {
    gStyle->SetOptStat(0);
    TFile* sigFile = TFile::Open(sigFilename);
    TFile* bkgFile = TFile::Open(bkgFilename);

    if (!sigFile || !bkgFile || sigFile->IsZombie() || bkgFile->IsZombie()) {
        std::cout << "Error opening files." << std::endl;
        return;
    }

    TTree* sigTree = dynamic_cast<TTree*>(sigFile->Get("HistData"));
    TTree* bkgTree = dynamic_cast<TTree*>(bkgFile->Get("HistData"));

    if (!sigTree || !bkgTree) {
        std::cout << "Error retrieving TTrees from the files." << std::endl;
        sigFile->Close();
        bkgFile->Close();
        return;
    }

    // Create 2D histograms
    TH2D* h_sig_PE_tau21 = new TH2D("h_sig_PE_tau21", "Signal: PE vs. tau21", 20, 0, 1, 20, 1.5, 5);
    TH2D* h_bkg_PE_tau21 = new TH2D("h_bkg_PE_tau21", "Background: PE vs. tau21", 20, 0, 1, 20, 1.5, 5);
    TH2D* h_sig_PE_asym = new TH2D("h_sig_PE_asym", "Signal: PE vs. asym", 20, 0, 1, 20, 1.5, 5);
    TH2D* h_sig_PE_mass = new TH2D("h_sig_PE_mass", "Signal: PE vs. mass", 30, 100, 150, 20, 1.5, 5);
    TH2D* h_bkg_PE_asym = new TH2D("h_bkg_PE_asym", "Background: PE vs. asym", 20, 0, 1, 20, 1.5, 5);
    TH2D* h_bkg_PE_mass = new TH2D("h_bkg_PE_mass", "Background: PE vs. mass", 30, 100, 150, 20, 1.5, 5);
    TH2D* h_sig_PE_D2 = new TH2D("h_sig_PE_D2", "Signal: PE vs. D2", 20, 0, 30, 20, 1.5, 5);
    TH2D* h_bkg_PE_D2 = new TH2D("h_bkg_PE_D2", "Background: PE vs. D2", 20, 0, 30, 20, 1.5, 5);
    TH2D* h_sig_PE_ntotal = new TH2D("h_sig_PE_ntotal", "Signal: PE vs. ntotal", 20, 0, 200, 20, 1.5, 5);
    TH2D* h_bkg_PE_ntotal = new TH2D("h_bkg_PE_ntotal", "Background: PE vs. ntotal", 20, 0, 200, 20, 1.5, 5);

    /*hmass > 100 && hmass < 150*/
    // Create 2D histograms
    sigTree->Draw("hPE:hasym>>h_sig_PE_asym", "hmass > 100 && hmass < 150", "goff");
    sigTree->Draw("hPE:hmass>>h_sig_PE_mass", "hmass > 100 && hmass < 150", "goff");
    bkgTree->Draw("hPE:hasym>>h_bkg_PE_asym", "hmass > 100 && hmass < 150", "goff");
    bkgTree->Draw("hPE:hmass>>h_bkg_PE_mass", "hmass > 100 && hmass < 150", "goff");

    sigTree->Draw("hPE:hntotal>>h_sig_PE_ntotal", "hmass > 100 && hmass < 150", "goff");
    sigTree->Draw("hPE:htau21>>h_sig_PE_tau21", "hmass > 100 && hmass < 150", "goff");
    bkgTree->Draw("hPE:hntotal>>h_bkg_PE_ntotal", "hmass > 100 && hmass < 150", "goff");
    bkgTree->Draw("hPE:htau21>>h_bkg_PE_tau21", "hmass > 100 && hmass < 150", "goff");
    sigTree->Draw("hPE:hD2>>h_sig_PE_D2", "hmass > 100 && hmass < 150", "goff");
    bkgTree->Draw("hPE:hD2>>h_bkg_PE_D2", "hmass > 100 && hmass < 150", "goff");

    // Create canvases and draw 2D histograms
    TCanvas* c_sig_PE_asym = new TCanvas("c_sig_PE_asym", "Signal: PE vs. asym", 800, 600);
    h_sig_PE_asym->Draw("CONT4Z");
    c_sig_PE_asym->SaveAs("sig_PE_asym.png");

    TCanvas* c_sig_PE_mass = new TCanvas("c_sig_PE_mass", "Signal: PE vs. mass", 800, 600);
    h_sig_PE_mass->Draw("CONT4Z");
    c_sig_PE_mass->SaveAs("sig_PE_mass.png");

    TCanvas* c_bkg_PE_asym = new TCanvas("c_bkg_PE_asym", "Background: PE vs. asym", 800, 600);
    h_bkg_PE_asym->Draw("CONT4Z");
    c_bkg_PE_asym->SaveAs("bkg_PE_asym.png");

    TCanvas* c_bkg_PE_mass = new TCanvas("c_bkg_PE_mass", "Background: PE vs. mass", 800, 600);
    h_bkg_PE_mass->Draw("CONT4Z");
    c_bkg_PE_mass->SaveAs("bkg_PE_mass.png");

    TCanvas* c_sig_PE_D2 = new TCanvas("c_sig_PE_D2", "Signal: PE vs. D2", 800, 600);
    h_sig_PE_D2->Draw("CONT4Z");
    c_sig_PE_D2->SaveAs("sig_PE_D2.png");

    TCanvas* c_bkg_PE_D2 = new TCanvas("c_bkg_PE_D2", "Background: PE vs. D2", 800, 600);
    h_bkg_PE_D2->Draw("CONT4Z");
    c_bkg_PE_D2->SaveAs("bkg_PE_D2.png");

    TCanvas* c_sig_PE_tau21 = new TCanvas("c_sig_PE_tau21", "Signal: PE vs. tau21", 800, 600);
    h_sig_PE_tau21->Draw("CONT4Z");
    c_sig_PE_tau21->SaveAs("sig_PE_tau21.png");

    TCanvas* c_bkg_PE_tau21 = new TCanvas("c_bkg_PE_tau21", "Background: PE vs. tau21", 800, 600);
    h_bkg_PE_tau21->Draw("CONT4Z");
    c_bkg_PE_tau21->SaveAs("bkg_PE_tau21.png");

    TCanvas* c_sig_PE_ntotal = new TCanvas("c_sig_PE_ntotal", "Signal: PE vs. ntotal", 800, 600);
    h_sig_PE_ntotal->Draw("CONT4Z");
    c_sig_PE_ntotal->SaveAs("sig_PE_ntotal.png");

    TCanvas* c_bkg_PE_ntotal = new TCanvas("c_bkg_PE_ntotal", "Background: PE vs. ntotal", 800, 600);
    h_bkg_PE_ntotal->Draw("CONT4Z");
    c_bkg_PE_ntotal->SaveAs("bkg_PE_ntotal.png");

    // Calculate correlation coefficients
    double corr_sig_PE_asym = h_sig_PE_asym->GetCorrelationFactor();
    double corr_sig_PE_mass = h_sig_PE_mass->GetCorrelationFactor();
    double corr_bkg_PE_asym = h_bkg_PE_asym->GetCorrelationFactor();
    double corr_bkg_PE_mass = h_bkg_PE_mass->GetCorrelationFactor();
    double corr_sig_PE_D2 = h_sig_PE_D2->GetCorrelationFactor();
    double corr_bkg_PE_D2 = h_bkg_PE_D2->GetCorrelationFactor();
    double corr_sig_PE_tau21 = h_sig_PE_tau21->GetCorrelationFactor();
    double corr_bkg_PE_tau21 = h_bkg_PE_tau21->GetCorrelationFactor();
    double corr_sig_PE_ntotal = h_sig_PE_ntotal->GetCorrelationFactor();
    double corr_bkg_PE_ntotal = h_bkg_PE_ntotal->GetCorrelationFactor();

    // Print correlation coefficients
    std::cout << "Correlation coefficients:" << std::endl;
    std::cout << "Signal: PE vs. asym: " << corr_sig_PE_asym << std::endl;
    std::cout << "Signal: PE vs. mass: " << corr_sig_PE_mass << std::endl;
    std::cout << "Background: PE vs. asym: " << corr_bkg_PE_asym << std::endl;
    std::cout << "Background: PE vs. mass: " << corr_bkg_PE_mass << std::endl;
    std::cout << "Signal: PE vs. tau21: " << corr_sig_PE_tau21 << std::endl;
    std::cout << "Background: PE vs. tau21: " << corr_bkg_PE_tau21 << std::endl;
    std::cout << "Signal: PE vs. D2: " << corr_sig_PE_D2 << std::endl;
    std::cout << "Background: PE vs. D2: " << corr_bkg_PE_D2 << std::endl;
    std::cout << "Signal: PE vs. ntotal: " << corr_sig_PE_ntotal << std::endl;
    std::cout << "Background: PE vs. ntotal: " << corr_bkg_PE_ntotal << std::endl;

    sigFile->Close();
    bkgFile->Close();
}