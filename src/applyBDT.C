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
TH1F* getBDTScoreHist(TTree* tree, const char* histName, const char* weightsFilename) {
    Double_t hpt, htau21, htau42, hthrust, hthrust_minor, hsubratio, hasym, hRjj, hsubjetwidth, hD2, hEntropy, hPE, hmass, heta, hphi;
    Double_t hsub1eta, hsub1phi, hsub1pt, hsub1mass, hsub2eta, hsub2phi, hsub2pt, hsub2mass;
    Double_t lep1eta, lep1phi, lep1pt, lep1mass, lep2eta, lep2phi, lep2pt, lep2mass;
    Double_t hsub1charge, hsub2charge;
    Int_t hntotal;
    
    tree->SetBranchAddress("hpt", &hpt);
    tree->SetBranchAddress("htau21", &htau21);
    tree->SetBranchAddress("htau42", &htau42);
    tree->SetBranchAddress("hthrust", &hthrust);
    tree->SetBranchAddress("hthrust_minor", &hthrust_minor);
    tree->SetBranchAddress("hsubratio", &hsubratio);
    tree->SetBranchAddress("hasym", &hasym);
    tree->SetBranchAddress("hRjj", &hRjj);
    tree->SetBranchAddress("hsubjetwidth", &hsubjetwidth);
    tree->SetBranchAddress("hD2", &hD2);
    tree->SetBranchAddress("hEntropy", &hEntropy);
    tree->SetBranchAddress("hPE", &hPE);
    tree->SetBranchAddress("hmass", &hmass);
    tree->SetBranchAddress("heta", &heta);
    tree->SetBranchAddress("hphi", &hphi);
    tree->SetBranchAddress("hsub1eta", &hsub1eta);
    tree->SetBranchAddress("hsub1phi", &hsub1phi);
    tree->SetBranchAddress("hsub1pt", &hsub1pt);
    tree->SetBranchAddress("hsub1mass", &hsub1mass);
    tree->SetBranchAddress("hsub2eta", &hsub2eta);
    tree->SetBranchAddress("hsub2phi", &hsub2phi);
    tree->SetBranchAddress("hsub2pt", &hsub2pt);
    tree->SetBranchAddress("hsub2mass", &hsub2mass);
    tree->SetBranchAddress("lep1eta", &lep1eta);
    tree->SetBranchAddress("lep1phi", &lep1phi);
    tree->SetBranchAddress("lep1pt", &lep1pt);
    tree->SetBranchAddress("lep1mass", &lep1mass);
    tree->SetBranchAddress("lep2eta", &lep2eta);
    tree->SetBranchAddress("lep2phi", &lep2phi);
    tree->SetBranchAddress("lep2pt", &lep2pt);
    tree->SetBranchAddress("lep2mass", &lep2mass);
    tree->SetBranchAddress("hsub1charge", &hsub1charge);
    tree->SetBranchAddress("hsub2charge", &hsub2charge);
    tree->SetBranchAddress("hntotal", &hntotal);
    
    Float_t f_hpt, f_htau21, f_htau42, f_hthrust, f_hthrust_minor, f_hsubratio, f_hasym, f_hRjj, f_hsubjetwidth, f_hD2, f_hEntropy, f_hPE, f_hmass, f_heta, f_hphi;
    Float_t f_hsub1eta, f_hsub1phi, f_hsub1pt, f_hsub1mass, f_hsub2eta, f_hsub2phi, f_hsub2pt, f_hsub2mass;
    Float_t f_lep1eta, f_lep1phi, f_lep1pt, f_lep1mass, f_lep2eta, f_lep2phi, f_lep2pt, f_lep2mass;
    Float_t f_hsub1charge, f_hsub2charge;
    Float_t f_hntotal;

    // Load the TMVA library
    TMVA::Tools::Instance();

    // Create a TMVA reader
    TMVA::Reader* reader = new TMVA::Reader("Color:!Silent");

    // Add the input variables to the reader
    reader->AddVariable("hpt", &f_hpt);
    reader->AddVariable("htau21", &f_htau21);
    reader->AddVariable("htau42", &f_htau42);
    reader->AddVariable("hthrust", &f_hthrust);
    reader->AddVariable("hthrust_minor", &f_hthrust_minor);
    reader->AddVariable("hsubratio", &f_hsubratio);
    reader->AddVariable("hasym", &f_hasym);
    reader->AddVariable("hRjj", &f_hRjj);
    reader->AddVariable("hsubjetwidth", &f_hsubjetwidth);
    reader->AddVariable("hD2", &f_hD2);
    reader->AddVariable("hEntropy", &f_hEntropy);
    reader->AddVariable("hPE", &f_hPE);
    reader->AddVariable("hmass", &f_hmass);
    reader->AddVariable("heta", &f_heta);
    reader->AddVariable("hphi", &f_hphi);
    reader->AddVariable("hsub1eta", &f_hsub1eta);
    reader->AddVariable("hsub1phi", &f_hsub1phi);
    reader->AddVariable("hsub1pt", &f_hsub1pt);
    reader->AddVariable("hsub1mass", &f_hsub1mass);
    reader->AddVariable("hsub2eta", &f_hsub2eta);
    reader->AddVariable("hsub2phi", &f_hsub2phi);
    reader->AddVariable("hsub2pt", &f_hsub2pt);
    reader->AddVariable("hsub2mass", &f_hsub2mass);
    reader->AddVariable("lep1eta", &f_lep1eta);
    reader->AddVariable("lep1phi", &f_lep1phi);
    reader->AddVariable("lep1pt", &f_lep1pt);
    reader->AddVariable("lep1mass", &f_lep1mass);
    reader->AddVariable("lep2eta", &f_lep2eta);
    reader->AddVariable("lep2phi", &f_lep2phi);
    reader->AddVariable("lep2pt", &f_lep2pt);
    reader->AddVariable("lep2mass", &f_lep2mass);
    reader->AddVariable("hsub1charge", &f_hsub1charge);
    reader->AddVariable("hsub2charge", &f_hsub2charge);
    reader->AddVariable("hntotal", &f_hntotal);

    // Book the trained BDT
    reader->BookMVA("BDT::BDTG", weightsFilename);
    TH1F* hist = new TH1F(histName, "", 100, -1, 1);
    Long64_t numEvents = tree->GetEntries();
    for (Long64_t i = 0; i < numEvents; i++) {
        tree->GetEntry(i);
        // Assign the values of the Double_t variables to the Float_t variables
        f_hpt = hpt;
        f_htau21 = htau21;
        f_htau42 = htau42;
        f_hthrust = hthrust;
        f_hthrust_minor = hthrust_minor;
        f_hsubratio = hsubratio;
        f_hasym = hasym;
        f_hRjj = hRjj;
        f_hsubjetwidth = hsubjetwidth;
        f_hD2 = hD2;
        f_hEntropy = hEntropy;
        f_hPE = hPE;
        f_hmass = hmass;
        f_heta = heta;
        f_hphi = hphi;
        f_hsub1eta = hsub1eta;
        f_hsub1phi = hsub1phi;
        f_hsub1pt = hsub1pt;
        f_hsub1mass = hsub1mass;
        f_hsub2eta = hsub2eta;
        f_hsub2phi = hsub2phi;
        f_hsub2pt = hsub2pt;
        f_hsub2mass = hsub2mass;
        f_lep1eta = lep1eta;
        f_lep1phi = lep1phi;
        f_lep1pt = lep1pt;
        f_lep1mass = lep1mass;
        f_lep2eta = lep2eta;
        f_lep2phi = lep2phi;
        f_lep2pt = lep2pt;
        f_lep2mass = lep2mass;
        f_hsub1charge = hsub1charge;
        f_hsub2charge = hsub2charge;
        f_hntotal = hntotal;

        Float_t bdtScore = reader->EvaluateMVA("BDT::BDTG");
        Double_t mvaErr = reader->GetMVAError();
        hist->Fill(bdtScore);
        //std::cout << bdtScore << std::endl;
        //std::cout << f_hntotal << std::endl;
    }
    
    // Clean up
    delete reader;

    return hist;
}
void applyBDT(const char* sigFilename, const char* bkgFilename, const char* weightsFilename){
    gStyle->SetOptStat(0);
    TFile* sigFile = TFile::Open(sigFilename);
    TFile* bkgFile = TFile::Open(bkgFilename);

    // Retrieve the "HistData" TTrees from the files
    TTree* sigTree = dynamic_cast<TTree*>(sigFile->Get("HistData"));
    TTree* bkgTree = dynamic_cast<TTree*>(bkgFile->Get("HistData"));

    // Declare the input variables
    TH1F* sigHist = getBDTScoreHist(sigTree, "sigHist", weightsFilename);
    TH1F* bkgHist = getBDTScoreHist(bkgTree, "bkgHist", weightsFilename);

    // Create a canvas and draw the histograms
    TCanvas* canvas = new TCanvas("canvas", "BDT Score", 800, 600);
    Double_t sigWeight = 0.0031132;
    Double_t bkgWeight = 92.79567345;
    // Normalize the histograms
    // sigHist->Scale(sigWeight);
    // bkgHist->Scale(bkgWeight);

    sigHist->SetLineColor(kBlue);
    sigHist->SetLineWidth(2);
    bkgHist->SetLineColor(kRed);
    bkgHist->SetLineWidth(2);

    // Draw the histograms
    bkgHist->Draw("hist");
    sigHist->Draw("hist same");
    

    // Find the best threshold
    Double_t bestThreshold = -1;
    Double_t maxSignificance = 0.00000;
    Int_t nBins = sigHist->GetNbinsX();
    for (Int_t i = 40; i <= nBins; i++) {
        Double_t threshold = sigHist->GetBinLowEdge(i);
        Double_t sig = sigHist->Integral(i, nBins);
        Double_t bkg = bkgHist->Integral(i, nBins);
        Double_t significance = sigWeight * sig / TMath::Sqrt(sigWeight * sig + bkgWeight * bkg);
        if (significance >= maxSignificance) {
            maxSignificance = significance;
            bestThreshold = threshold;
        }
    }

    // Draw a vertical line at the best threshold
    TLine* line = new TLine(bestThreshold, 0, bestThreshold, sigHist->GetMaximum());
    line->SetLineColor(kGreen);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    //line->Draw("same");

    // Add the best threshold and weights to the legend
    TLegend* legend = new TLegend(0.75, 0.75, 0.95, 0.95);
    legend->AddEntry(sigHist, Form("Signal (Weight: %.4f)", sigWeight), "l");
    legend->AddEntry(bkgHist, Form("Background (Weight: %.4f)", bkgWeight), "l");
    //legend->AddEntry(line, Form("Best Threshold: %.4f", bestThreshold), "l");
    legend->Draw("same");

    // Save the canvas as an image file
    canvas->SaveAs("bdt_score_histograms_with_threshold_and_weights.png");

    // Print the best threshold, maximum significance, and weights
    std::cout << "Best Threshold: " << bestThreshold << std::endl;
    std::cout << "Maximum Significance: " << maxSignificance << std::endl;
    std::cout << "Signal Weight: " << sigWeight << std::endl;
    std::cout << "Background Weight: " << bkgWeight << std::endl;


    // Clean up
    delete sigHist;
    delete bkgHist;
    delete canvas;
    delete legend;
    sigFile->Close();
    bkgFile->Close();
}