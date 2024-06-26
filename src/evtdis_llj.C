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
	
void evtdis_llj(const char *inputFile, int evt_entry) {

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
    if (evt_entry >= nEntries) {
        std::cout << "Out of index!"<< std::endl;
        return;
    }

    TClonesArray *branchJet = treeReader->UseBranch("PFJet10");
    TClonesArray *branchGenPar = treeReader->UseBranch("Particle");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchEvent = treeReader->UseBranch("Event");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    treeReader -> ReadEntry(evt_entry);
    std::cout << "Reading the " << evt_entry << "th entry..." << std::endl;
    Int_t nJet = branchJet -> GetEntries();
    std::cout << nJet << " jets in this event." << std::endl;
    Int_t ne = branchElectron->GetEntries();
    std::cout << ne << " electrons in this event." << std::endl;
    Int_t nMu = branchMuon->GetEntries();
    std::cout << nMu << " muons in this event." << std::endl;
    Int_t nGenPar = branchGenPar->GetEntries();
    std::cout << nGenPar << " generator level particles in this event." << std::endl;

    GenParticle *trueZ;
    GenParticle *lep1;
    GenParticle *lep2;
    GenParticle *genpar;
    TLorentzVector ZP4;
    TLorentzVector lep1P4;
    TLorentzVector lep2P4;
    std::vector<TLorentzVector> final_products;
    bool emuflag = 0;
    for (int gentry=0; gentry < nGenPar; ++gentry){
        genpar = (GenParticle *) branchGenPar -> At(gentry);
        if (genpar -> PID == 23){
            trueZ = (GenParticle *) branchGenPar -> At(gentry);
            lep1 = (GenParticle *) branchGenPar -> At(trueZ->D1);
            lep2 = (GenParticle *) branchGenPar -> At(trueZ->D2);
            if (abs(abs(lep1 -> PID)-12) == 1 and abs(abs(lep2 -> PID)-12) == 1 and lep1->PID + lep2->PID == 0) {
                ZP4 = trueZ -> P4();
            lep1P4 = lep1 -> P4();
            lep2P4 = lep2 -> P4();
            }
        }
    }
    //std::cout << "Get Truth Info" << std::endl;
    const Double_t zeta[] = {ZP4.Eta()};
    const Double_t zphi[] = {ZP4.Phi()};
    const Double_t lepeta[] = {lep1P4.Eta(), lep2P4.Eta()};
    const Double_t lepphi[] = {lep1P4.Phi(), lep2P4.Phi()};
    TGraph *zscat = new TGraph(1, zeta, zphi);
    TGraph *lepscat = new TGraph(2, lepeta, lepphi);
    std::cout << "First lepton at " << lep1P4.Eta()<< ", " << lep1P4.Phi() << std::endl;
    std::cout << "Second lepton at " << lep2P4.Eta()<< ", " << lep2P4.Phi() << std::endl;
    Double_t mass_min = 100000000;
    Jet *jet;
    Jet *hjet;
    TLorentzVector jetP4;
    TLorentzVector tmpP4;
    for (Int_t jentry = 0; jentry < nJet; ++jentry) {
        jet = (Jet *) branchJet -> At(jentry);
        if (jet -> PT < 300) {
            continue;
        }
        tmpP4 = jet -> P4();
        if (tmpP4.DeltaR(ZP4) < 1.0) {
            continue;
        }
        if (TMath::Abs(jet -> Mass - 125) < mass_min) {
            mass_min = TMath::Abs(jet -> Mass - 125);
            hjet = (Jet *) branchJet -> At(jentry);
        }    
    }

    jetP4 = hjet -> P4();

    if (jetP4.DeltaR(ZP4) < 1.0) {
	    std::cout << "Z and H too close!" << std::endl;
	    return;
    }
    for (int gentry=0; gentry < nGenPar; ++gentry){
        genpar = (GenParticle *) branchGenPar -> At(gentry);
        if (genpar -> Status == 1) {
            if (jetP4.DeltaR(genpar -> P4()) <= 1.0) {
                final_products.push_back(genpar -> P4());
            }
        }
    }
    TEllipse* jetcirc = new TEllipse(jetP4.Eta(), jetP4.Phi() ,1.0, 1.0);
    jetcirc -> SetFillStyle(4000);
    jetcirc -> SetLineColor(2);
    jetcirc -> SetLineWidth(5);

    TObject *object;

    //std::cout << "Find jet" << std::endl;

    Int_t towercnt = 0;
    Int_t trackcnt = 0;
    std::cout <<  hjet -> Constituents.GetEntries() << " constituents in jet" << endl;
    for (int i = 0; i < hjet->Constituents.GetEntries(); ++i) {
        object = hjet->Constituents.At(i);
        if (object == 0) {
            std::cout << "Constituents is unaccessible! Check memory loading!" << std::endl;
        }
        if(object->IsA() == Tower::Class()) {
	        towercnt++;
    	}
        if(object->IsA() == Track::Class()) {
            trackcnt++;
        }
    }
    const Int_t towersize = towercnt;
    const Int_t tracksize = trackcnt;
    std::cout << towersize << " towers"<< " and " << tracksize << " tracks." << std::endl;
    if (towersize == 0 or tracksize == 0) {
	    std::cout << "Not enough tower or track!" << std::endl;
	    return;
    }
    Double_t towereta[towersize];
    Double_t tracketa[tracksize];
    Double_t towerphi[towersize];
    Double_t trackphi[tracksize];
    towercnt = 0;
    trackcnt = 0;
    //std::cout << "Store jet constituents" << std::endl;
    for(int i = 0; i < hjet->Constituents.GetEntries(); ++i) {
        object = hjet->Constituents.At(i);
	if(object == 0) {
	    std::cout << "Constituents is unaccessible! Check memory loading!" << std::endl;
	}
        if(object->IsA() == Tower::Class()) {
            Tower *tower;
            tower = (Tower*) object;
            TLorentzVector towerP4;
            towerP4 = tower -> P4();
            towereta[towercnt] = towerP4.Eta();
            towerphi[towercnt] = towerP4.Phi();
            towercnt++;
    	}
        if(object->IsA() == Track::Class()) {
            Track *track;
            track = (Track*) object;
            TLorentzVector trackP4;
            trackP4 = track -> P4();
            tracketa[trackcnt] = trackP4.Eta();
            trackphi[trackcnt] = trackP4.Phi();
            trackcnt++;
        }
    }
    std::cout << "Reconstruct higgs mass " << jetP4.Mag()<<std::endl;
    std::cout << "Reconstruct higgs pt " << jetP4.Pt()<<std::endl;

    const Double_t* consttowereta = towereta;
    const Double_t* consttowerphi = towerphi;
    const Double_t* consttracketa = tracketa;
    const Double_t* consttrackphi = trackphi;
    TGraph *towerscat = new TGraph(towersize, consttowereta, consttowerphi);
    TGraph *trackscat = new TGraph(tracksize, consttracketa, consttrackphi);
    //TH2D *towerhisto = new TH2D("towerhisto", "towerhisto", )
//Plot truth products
    const Int_t finalProductSize = final_products.size();

    // Arrays to store eta and phi values
    Double_t finalProductEta[finalProductSize];
    Double_t finalProductPhi[finalProductSize];

    // Fill the arrays
    for (int i = 0; i < finalProductSize; ++i) {
        finalProductEta[i] = final_products[i].Eta();
        finalProductPhi[i] = final_products[i].Phi();
    }
    TGraph *finalProductScat = new TGraph(finalProductSize, finalProductEta, finalProductPhi);


    auto mg = new TMultiGraph();

    TCanvas *totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    totalcanvas -> SetWindowSize(1204,1228);
    totalcanvas -> SetCanvasSize(1200,1200);
    zscat -> SetMarkerColor(kBlue);
    zscat -> SetMarkerStyle(kCircle);
    zscat -> SetMarkerSize(2);
    lepscat -> SetMarkerColor(kRed);
    lepscat -> SetMarkerStyle(kPlus);
    lepscat -> SetMarkerSize(2);
    towerscat -> SetMarkerColor(kGreen);
    towerscat -> SetMarkerStyle(kStar);
    towerscat -> SetMarkerSize(3);
    trackscat -> SetMarkerColor(kMagenta);
    trackscat -> SetMarkerStyle(kStar);
    trackscat -> SetMarkerSize(3);
    finalProductScat -> SetMarkerColor(kOrange);
    finalProductScat -> SetMarkerStyle(kCircle);
    finalProductScat -> SetMarkerSize(3);
    TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend -> AddEntry(zscat, "truth Z boson", "p");
    legend -> AddEntry(lepscat, "truth lepton pair", "p");
    legend -> AddEntry(towerscat, "tower in jet", "p");
    legend -> AddEntry(trackscat, "track in jet", "p");
    legend -> AddEntry(finalProductScat, "truth constituents inside jet", "p");
    legend -> AddEntry(jetcirc, "reco higgs jet", "l");
    //towerscat -> Draw("");
    //trackscat -> Draw("same");
    //higgsscat -> Draw("same");
    //zscat -> Draw("same");
    //gluonscat -> Draw("same");
    //lepscat -> Draw("same");
    //jetcirc -> Draw("same");
    //legend -> Draw("same");
    mg -> Add(finalProductScat);
    mg -> Add(towerscat);
    mg -> Add(trackscat);
    mg -> Add(lepscat);
    mg -> Add(zscat);
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("#phi");
    gPad->Modified();
    mg->GetXaxis()->SetLimits(-4,4);
    mg->SetMinimum(-3.14);
    mg->SetMaximum(3.14);
    jetcirc -> Draw("same");
    legend -> Draw("same");
    std::ostringstream filename;
    filename << "eventdisplay_bkg_" << evt_entry << ".png";
    totalcanvas -> SaveAs(filename.str().c_str());
}
