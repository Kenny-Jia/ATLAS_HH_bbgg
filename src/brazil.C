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

#include <TGraph.h>

#include <TGraphAsymmErrors.h>

#include <TCanvas.h>

#include <TLegend.h>

#include <TLine.h>

#include <TLatex.h>

#include <RooStats/NumberCountingUtils.h>

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
    reader->BookMVA("BDT::BDT", weightsFilename);
    TH1F* hist = new TH1F(histName, "", 60, -1, 1);
    Long64_t numEvents = tree->GetEntries();
    for (Long64_t i = 0; i < numEvents; i++) {
        tree->GetEntry(i);

        if  (hmass >= 100 && hmass <= 150) {
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

            Float_t bdtScore = reader->EvaluateMVA("BDT::BDT");
            Double_t mvaErr = reader->GetMVAError();
            hist->Fill(bdtScore);
        }
        //std::cout << bdtScore << std::endl;
        //std::cout << f_hntotal << std::endl;
    }
    
    // Clean up
    delete reader;

    return hist;
}

using namespace RooStats;

void CalculateCLs(const TH1F* h_sig, const TH1F* h_bkg, const std::vector<double>& mu_vals,
                  std::vector<double>& obs_cls, std::vector<double>& exp_cls,
                  std::vector<double>& exp_cls_err_up, std::vector<double>& exp_cls_err_down) {

    obs_cls.clear();
    exp_cls.clear();
    exp_cls_err_up.clear();
    exp_cls_err_down.clear();

    double bkg_total = h_bkg->Integral();
    double bkg_err = std::sqrt(bkg_total);

    for (double mu : mu_vals) {

        double sig_total = mu * h_sig->Integral();

        double obs_events = bkg_total + sig_total;
        double obs_cl = 1.0;
        for (int i = 0; i < obs_events; ++i) {
            obs_cl -= TMath::Poisson(i, bkg_total);
        }
        obs_cls.push_back(obs_cl);

        double exp_cl = 0.5;
        exp_cls.push_back(exp_cl);

        double exp_cl_err_up = 1.0;
        for (int i = 0; i < bkg_total + bkg_err; ++i) {
            exp_cl_err_up -= TMath::Poisson(i, bkg_total);
        }
        exp_cls_err_up.push_back(exp_cl_err_up - exp_cl);

        double exp_cl_err_down = 1.0;
        for (int i = 0; i < bkg_total - bkg_err; ++i) {
            exp_cl_err_down -= TMath::Poisson(i, bkg_total);
        }
        exp_cls_err_down.push_back(exp_cl - exp_cl_err_down);
    }
}


void PlotBrazil(const std::vector<double>& mu_vals, const std::vector<double>& obs_cls,
                const std::vector<double>& exp_cls, const std::vector<double>& exp_cls_err_up,
                const std::vector<double>& exp_cls_err_down) {

    TGraph* gr_obs = new TGraph(mu_vals.size(), &mu_vals[0], &obs_cls[0]);
    TGraph* gr_exp = new TGraph(mu_vals.size(), &mu_vals[0], &exp_cls[0]);
    TGraphAsymmErrors* gr_exp_err = new TGraphAsymmErrors(mu_vals.size(), &mu_vals[0], &exp_cls[0],
                                                           nullptr, nullptr, &exp_cls_err_down[0], &exp_cls_err_up[0]);

    gr_obs->SetLineColor(kBlack);
    gr_obs->SetLineWidth(2);
    gr_exp->SetLineColor(kRed);
    gr_exp->SetLineWidth(2);
    gr_exp->SetLineStyle(2);
    gr_exp_err->SetFillColor(kGreen);
    gr_exp_err->SetFillStyle(1001);

    TCanvas* c = new TCanvas("c", "Brazil Plot", 800, 600);
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    gr_exp_err->Draw("a3");
    gr_exp->Draw("same");
    gr_obs->Draw("same");

    leg->AddEntry(gr_obs, "Observed", "l");
    leg->AddEntry(gr_exp, "Expected", "l");
    leg->AddEntry(gr_exp_err, "Expected #pm 1#sigma", "f");
    leg->Draw();

    gr_exp_err->GetXaxis()->SetTitle("Signal Strength (#mu)");
    gr_exp_err->GetYaxis()->SetTitle("Confidence Level");
    gr_exp_err->SetTitle("Brazil Plot");

    TLine* line95 = new TLine(mu_vals.front(), 0.05, mu_vals.back(), 0.05);
    line95->SetLineColor(kBlack);
    line95->SetLineWidth(2);
    line95->SetLineStyle(7);
    line95->Draw();
    TLatex* text95 = new TLatex(mu_vals.back()*0.8, 0.96, "95% CL");
    text95->SetTextAlign(33);
    text95->Draw();

    c->SaveAs("brazil_plot.png");
}

void brazil(const char* weightsFilename,
                const char* EXP_sigFilename,
                    const char* EXP_zjetFilename,
                        const char* EXP_zhbbFilename,
                            const char* EXP_zhccFilename,
                                const char* EXP_zh4qFilename,
                                    const char* EXP_zzqqFilename,
                                        const char* EXP_zwqqFilename,
                                            const char* EXP_ttFilename) {
    gStyle->SetOptStat(0);

    TFile* EXP_sigFile = TFile::Open(EXP_sigFilename);
    TFile* EXP_zjetFile = TFile::Open(EXP_zjetFilename);
    TFile* EXP_zhbbFile = TFile::Open(EXP_zhbbFilename);
    TFile* EXP_zhccFile = TFile::Open(EXP_zhccFilename);
    TFile* EXP_zh4qFile = TFile::Open(EXP_zh4qFilename);
    TFile* EXP_zzqqFile = TFile::Open(EXP_zzqqFilename);
    TFile* EXP_zwqqFile = TFile::Open(EXP_zwqqFilename);
    TFile* EXP_ttFile = TFile::Open(EXP_ttFilename);

    TTree* EXP_sigTree = dynamic_cast<TTree*>(EXP_sigFile->Get("HistData"));
    TTree* EXP_zjetTree = dynamic_cast<TTree*>(EXP_zjetFile->Get("HistData"));
    TTree* EXP_zhbbTree = dynamic_cast<TTree*>(EXP_zhbbFile->Get("HistData"));
    TTree* EXP_zhccTree = dynamic_cast<TTree*>(EXP_zhccFile->Get("HistData"));
    TTree* EXP_zh4qTree = dynamic_cast<TTree*>(EXP_zh4qFile->Get("HistData"));
    TTree* EXP_zzqqTree = dynamic_cast<TTree*>(EXP_zzqqFile->Get("HistData"));
    TTree* EXP_zwqqTree = dynamic_cast<TTree*>(EXP_zwqqFile->Get("HistData"));
    TTree* EXP_ttTree = dynamic_cast<TTree*>(EXP_ttFile->Get("HistData"));

    TH1F* h_exp_sig = getBDTScoreHist(EXP_sigTree, "h_exp_sig", weightsFilename);
    TH1F* h_exp_zjet = getBDTScoreHist(EXP_zjetTree, "h_exp_zjet", weightsFilename);
    TH1F* h_exp_zhbb = getBDTScoreHist(EXP_zhbbTree, "h_exp_zhbb", weightsFilename);
    TH1F* h_exp_zhcc = getBDTScoreHist(EXP_zhccTree, "h_exp_zhcc", weightsFilename);
    TH1F* h_exp_zh4q = getBDTScoreHist(EXP_zh4qTree, "h_exp_zh4q", weightsFilename);
    TH1F* h_exp_zzqq = getBDTScoreHist(EXP_zzqqTree, "h_exp_zzqq", weightsFilename);
    TH1F* h_exp_zwqq = getBDTScoreHist(EXP_zwqqTree, "h_exp_zwqq", weightsFilename);
    TH1F* h_exp_tt = getBDTScoreHist(EXP_ttTree, "h_exp_tt", weightsFilename);

    TFile* outputFile = new TFile("BDT_score.root", "RECREATE");
    h_exp_sig->Write();
    h_exp_zjet->Write();
    h_exp_zhbb->Write();
    h_exp_zhcc->Write();
    h_exp_zh4q->Write();
    h_exp_zzqq->Write();
    h_exp_zwqq->Write();
    h_exp_tt->Write();
    outputFile->Close();

    delete h_exp_sig;
    delete h_exp_zjet;
    delete h_exp_zhbb;
    delete h_exp_zhcc;
    delete h_exp_zh4q;
    delete h_exp_zzqq;
    delete h_exp_zwqq;
    delete h_exp_tt;
}