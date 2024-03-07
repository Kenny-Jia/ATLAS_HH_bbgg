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

bool compareTVector3(const TVector3 & a,
    const TVector3 & b) {
    return a.Mag() > b.Mag();
}

#include <vector>
#include <TLorentzVector.h>

double calculateSubJetWidth(const std::vector<TLorentzVector>& constituents, const std::vector<TLorentzVector>& subjets) {
    if (constituents.empty()) return 0; // Guard against division by zero

    double numerator = 0.0;
    double denominator = 0.0;

    for (const auto& constituent : constituents) {
        double minDeltaR = std::numeric_limits<double>::max();
        for (const auto& subjet : subjets) {
            double deltaR = constituent.DeltaR(subjet);
            if (deltaR < minDeltaR) {
                minDeltaR = deltaR;
            }
        }
        numerator += constituent.Pt() * minDeltaR;
        denominator += constituent.Pt();
    }

    return denominator != 0 ? numerator / denominator : 0;
}


std::pair < double, double > calculateThrust(const std::vector < TLorentzVector > & momenta) {
    double t;
    TVector3 totalMomentum(0, 0, 0);
    for (const auto & particle: momenta) {
        totalMomentum += particle.Vect();
    }

    double totalEnergy = 0;
    for (const auto & particle: momenta) {
        totalEnergy += particle.E();
    }
    TVector3 beta = -totalMomentum * (1.0 / totalEnergy);

    std::vector < TVector3 > momentaCp;
    for (const auto & particle: momenta) {
        TLorentzVector tmpPar = particle;
        tmpPar.Boost(beta);
        momentaCp.push_back(tmpPar.Vect());
    }
    unsigned int n = momentaCp.size();
    if (n == 0) return std::make_pair(0.0, 0.0);
    if (n == 1) return std::make_pair(1.0, 0.0);

    assert(n >= 3);
    n = 3;
    if (momentaCp.size() == 3) n = 3;

    std::vector < TVector3 > tvec;
    std::vector < double > tval;
    TVector3 taxis(0, 0, 0);
    std::sort(momentaCp.begin(), momentaCp.end(), compareTVector3);

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
    for (const auto & p: momentaCp) {
        norm += p.Mag();
    }

    t = 0.;
    for (unsigned int i = 0; i < tvec.size(); i++) {
        if (tval[i] > t) {
            t = tval[i];
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

void llgg(const char * inputFile,
    TH1D * Hmass,
        TH1D * tau21,
            TH1D * tau42,
                TH1D * thrust,
                    TH1D * thrust_minor,
                        TH1D * ntotal,
                            TH1D * pt,
                                TH1D * truth_total,
                                    TH1D * subratio,
                                        TH1D * asym,
                                            TH1D * Rjj,
                                                TH1D * subjetwidth) {

    Double_t ptcut = 350;
    std::cout << "loading Delphes library" << std::endl;
    gSystem -> Load("/eos/user/h/hjia/ATLAS_HH_bbgg/Delphes-3.5.0/libDelphes.so");
    //gSystem->Load("libDelphes");
    gROOT -> ProcessLine("gErrorIgnoreLevel = kWarning;");
    gErrorIgnoreLevel = kWarning;

    std::cout << "open up file" << std::endl;
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader * treeReader = new ExRootTreeReader( & chain);
    Long64_t nEntries = treeReader -> GetEntries();

    TClonesArray * branchsmallJet = treeReader -> UseBranch("UniqueJet");
    TClonesArray * branchJet = treeReader -> UseBranch("PFJet10");
    TClonesArray * branchElectron = treeReader -> UseBranch("Electron");
    TClonesArray * branchMuon = treeReader -> UseBranch("Muon");
    TClonesArray * branchGenPar = treeReader -> UseBranch("Particle");
    TClonesArray * branchMET = treeReader -> UseBranch("MissingET");
    TClonesArray * branchEvent = treeReader -> UseBranch("Event");

    TClonesArray * branchTrack = treeReader -> UseBranch("Track");
    TClonesArray * branchTower = treeReader -> UseBranch("Tower");
    TClonesArray * branchEFlowTrack = treeReader -> UseBranch("EFlowTrack");
    TClonesArray * branchEFlowPhoton = treeReader -> UseBranch("EFlowPhoton");
    TClonesArray * branchEFlowNeutralHadron = treeReader -> UseBranch("EFlowNeutralHadron");

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
        treeReader -> ReadEntry(entry);

        Int_t nJet = branchJet -> GetEntriesFast();
        Int_t ne = branchElectron -> GetEntries();
        Int_t nMu = branchMuon -> GetEntries();
        Int_t nGenPar = branchGenPar -> GetEntries();
        MissingET *met;
        met = (MissingET*) branchMET->At(0);


        TObject * object;

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
        if (met->MET >= 40) {
            continue;
        }

        Electron * elec1, * elec2;
        Muon * mu1, * mu2;

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
        l1.SetPtEtaPhiM(0, 0, 0, 0);
        l2.SetPtEtaPhiM(0, 0, 0, 0);
        h.SetPtEtaPhiM(0, 0, 0, 0);

        switch (emuflag) {
        case 1:
            elec1 = (Electron * ) branchElectron -> At(0);
            elec2 = (Electron * ) branchElectron -> At(1);

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
            mu1 = (Muon * ) branchMuon -> At(0);
            mu2 = (Muon * ) branchMuon -> At(1);

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
        if (z_ll.Pt() < ptcut) {
            //continue;
        }
        Double_t mass_min = 100000000;
        Jet * jet;
        Jet * hjet;
        Int_t hentry = 0;
        bool find_jet = false;
        for (Int_t jentry = 0; jentry < nJet; jentry++) {
            jet = (Jet * ) branchJet -> At(jentry);

            if (fabs(jet -> Eta) > 2.0) {
                continue;
            }

            if (jet -> PT < ptcut) {
                continue;
            }
            if (TMath::Abs(jet -> P4().DeltaR(z_ll) - TMath::Pi()) > 0.5) {
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
        hjet = (Jet * ) branchJet -> At(hentry);
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
        htau21 = (hjet -> Tau[1]) / (hjet -> Tau[0]);
        htau42 = (hjet -> Tau[3]) / (hjet -> Tau[1]);
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
        if (z_ll.Pt() < ptcut or h_gg.Pt() < ptcut) {
            //continue;
        }

        Double_t lpairmass = z_ll.Mag();
        Double_t gpairmass = h_gg.Mag();

        Int_t nsmallJet = branchsmallJet -> GetEntriesFast();
        bool smallBtag = false;
        Int_t nsmallB = 0;
        Jet * smalljet;
        Double_t ratio;
        Double_t ptasym;
        std::vector<Double_t> subjetPTArr;
        std::vector<TLorentzVector> subjetP4Arr;
        for (int j = 0; j < nsmallJet; ++j) {
            smalljet = (Jet * ) branchsmallJet -> At(j);
            if (h_gg.DeltaR(smalljet -> P4()) < 1.0) {
                nsmallB++;
                subjetPTArr.push_back(smalljet -> PT);
                subjetP4Arr.push_back(smalljet -> P4());
                if (smalljet -> BTag >= 1) {
                    smallBtag = true;
                }
            }
        }

        if (smallBtag == true) {
            continue;
        }

        if (subjetPTArr.size() < 2) {
            ratio = 1;
            ptasym = 1;
        } else {
            std::sort(subjetPTArr.begin(), subjetPTArr.end());
            std::reverse(subjetPTArr.begin(), subjetPTArr.end());
            ratio = subjetPTArr[1] * 1.0 / subjetPTArr[0];
            ptasym = (subjetPTArr[0] - subjetPTArr[1]) * 1.0 / (subjetPTArr[0] + subjetPTArr[1]);
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
        GenParticle * trueHiggs;
        GenParticle * genpar;
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
        Int_t ntotal_truth = 0;
        for (int gentry=0; gentry < nGenPar; ++gentry){
        	genpar = (GenParticle *) branchGenPar -> At(gentry);
        	if (genpar -> Status == 1) {
                if(genpar -> P4().DeltaR(h_gg) < 1) {
                    ntotal_truth++;
                }
            }
        }

        
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

        std::pair < double, double > thrust_result = calculateThrust(particlesArr);
        
        Double_t SubJetWidth = calculateSubJetWidth(particlesArr, subjetP4Arr);
        double hthrust = thrust_result.first;
        double hthrust_minor = thrust_result.second;

        if (hNNeutral+hNCharged < 40) { 
            //continue;
        }
        if (htau21 > 0.65) {
            //continue;
        }
        if (htau21 < 0.2) {
            //continue;
        }
        if (hthrust_minor < 0.15) {
            //continue;
        }
        if (hthrust > 0.97) {
            //continue;
        }
        //std::cout << "Thrust: " << thrust << std::endl;

        //total++;
        
        /*
        if (higgstag == false) {
        	continue;
        }
        */
        //std::cout << "event number " << entry << std::endl;
        //Zmass -> Fill(lpairmass);
        Hmass -> Fill(h_gg.Mag());
        tau21 -> Fill(htau21);
        tau42 -> Fill(htau42);
        thrust -> Fill(hthrust);
        thrust_minor -> Fill(hthrust_minor);
        ntotal -> Fill(hNNeutral+hNCharged);
        pt -> Fill(hpt);
        truth_total -> Fill(ntotal_truth);
        subratio -> Fill(ratio);
        asym -> Fill(ptasym);
        if (subjetPTArr.size() < 2) {
            Rjj -> Fill (0);
        } else {
            Rjj -> Fill (h_gg.Mag()/(TMath::Sqrt(subjetPTArr[0] * subjetPTArr[1])));
        }
        subjetwidth -> Fill(SubJetWidth);
        
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
void llgg_stack(const char * sigFile,
    const char * bkgzjetsFile,
        const char * bkgzhbbFile,
            const char * bkgzhccFile,
                const char * bkgzh4qFile,
                    const char * bkgzzqqFile,
                        const char * bkgzwqqFile,
                            const char * bkgttFile) {
    TFile * output = new TFile("ZH_analysis.root", "RECREATE");
    TTree * tree_output = new TTree("tree_output", "Delphes");

    int histmin = 0;
    int histmax = 250;
    int bins = 50;
    int taubins = 30;
    int thrustbins = 30;
    int ntotalbins = 30;
    int ptbins = 30;
    int ratiobins = 30;
    int asymbins = 30;
    int Rjjbins = 30;
    int subjetwidthbins = 30;

    TH1D * sigHmass = new TH1D("sigHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzjetsHmass = new TH1D("bkgzjetsHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzhbbHmass = new TH1D("bkgzhbbHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzhccHmass = new TH1D("bkgzhccHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzh4qHmass = new TH1D("bkgzh4qHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzzqqHmass = new TH1D("bkgzzqqHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzwqqHmass = new TH1D("bkgzwqqHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgttHmass = new TH1D("bkgttHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);

    TH1D * sigtau21 = new TH1D("sigtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzjetstau21 = new TH1D("bkgzjetstau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzhbbtau21 = new TH1D("bkgzhbbtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzhcctau21 = new TH1D("bkgzhcctau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzh4qtau21 = new TH1D("bkgzh4qtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzzqqtau21 = new TH1D("bkgzzqqtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzwqqtau21 = new TH1D("bkgzwqqtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgtttau21 = new TH1D("bkgtttau21", "Reconstucted Higgs tau21", taubins, 0, 1);

    TH1D * sigtau42 = new TH1D("sigtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzjetstau42 = new TH1D("bkgzjetstau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzhbbtau42 = new TH1D("bkgzhbbtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzhcctau42 = new TH1D("bkgzhcctau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzh4qtau42 = new TH1D("bkgzh4qtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzzqqtau42 = new TH1D("bkgzzqqtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzwqqtau42 = new TH1D("bkgzwqqtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgtttau42 = new TH1D("bkgtttau42", "Reconstucted Higgs tau42", taubins, 0, 1);

    TH1D * sigthrust = new TH1D("sigthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzjetsthrust = new TH1D("bkgzjetsthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzhbbthrust = new TH1D("bkgzhbbthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzhccthrust = new TH1D("bkgzhccthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzh4qthrust = new TH1D("bkgzh4qthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzzqqthrust = new TH1D("bkgzzqqthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzwqqthrust = new TH1D("bkgzwqqthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgttthrust = new TH1D("bkgttthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);

    TH1D * sigthrust_minor = new TH1D("sigthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzjetsthrust_minor = new TH1D("bkgzjetsthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzhbbthrust_minor = new TH1D("bkgzhbbthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzhccthrust_minor = new TH1D("bkgzhccthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzh4qthrust_minor = new TH1D("bkgzh4qthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzzqqthrust_minor = new TH1D("bkgzzqqthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzwqqthrust_minor = new TH1D("bkgzwqqthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgttthrust_minor = new TH1D("bkgttthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);

    TH1D * signtotal = new TH1D("signtotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgzjetsntotal = new TH1D("bkgzjetsntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgzhbbntotal = new TH1D("bkgzhbbntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgzhccntotal = new TH1D("bkgzhccntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgzh4qntotal = new TH1D("bkgzh4qntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgzzqqntotal = new TH1D("bkgzzqqntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgzwqqntotal = new TH1D("bkgzwqqntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);
    TH1D * bkgttntotal = new TH1D("bkgttntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 120);

    TH1D * sigtruth_total = new TH1D("sigtruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzjetstruth_total = new TH1D("bkgzjetstruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzhbbtruth_total = new TH1D("bkgzhbbtruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzhcctruth_total = new TH1D("bkgzhcctruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzh4qtruth_total = new TH1D("bkgzh4qtruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzzqqtruth_total = new TH1D("bkgzzqqtruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzwqqtruth_total = new TH1D("bkgzwqqtruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgtttruth_total = new TH1D("bkgtttruth_total", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);

    TH1D * sigpt = new TH1D("signpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzjetspt = new TH1D("bkgzjetsnpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzhbbpt = new TH1D("bkgzhbbpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzhccpt = new TH1D("bkgzhccpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzh4qpt = new TH1D("bkgzh4qpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzzqqpt = new TH1D("bkgzzqqpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzwqqpt = new TH1D("bkgzwqqpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgttpt = new TH1D("bkgttpt", "Reconstucted Higgs pt", ptbins, 200, 1000);

    TH1D * sigratio = new TH1D("signratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgzjetsratio = new TH1D("bkgzjetsnratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgzhbbratio = new TH1D("bkgzhbbratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgzhccratio = new TH1D("bkgzhccratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgzh4qratio = new TH1D("bkgzh4qratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgzzqqratio = new TH1D("bkgzzqqratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgzwqqratio = new TH1D("bkgzwqqratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);
    TH1D * bkgttratio = new TH1D("bkgttratio", "Reconstucted Higgs ratio", ratiobins, 0, 1);

    TH1D * sigasym = new TH1D("signasym", "Reconstucted Higgs asymmetry", asymbins, 0, 1);
    TH1D * bkgzjetsasym = new TH1D("bkgzjetsnasym", "Reconstucted Higgs asymmetry", asymbins, 0, 1);
    TH1D * bkgzhbbasym = new TH1D("bkgzhbbasym", "Reconstucted Higgs asymetry", asymbins, 0, 1);
    TH1D * bkgzhccasym = new TH1D("bkgzhccasym", "Reconstucted Higgs asymetry", asymbins, 0, 1);
    TH1D * bkgzh4qasym = new TH1D("bkgzh4qasym", "Reconstucted Higgs asymetry", asymbins, 0, 1);
    TH1D * bkgzzqqasym = new TH1D("bkgzzqqasym", "Reconstucted Higgs asymetry", asymbins, 0, 1);
    TH1D * bkgzwqqasym = new TH1D("bkgzwqqasym", "Reconstucted Higgs asymetry", asymbins, 0, 1);
    TH1D * bkgttasym = new TH1D("bkgttasym", "Reconstucted Higgs asymetry", asymbins, 0, 1);
    
    TH1D * sigRjj = new TH1D("signRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgzjetsRjj = new TH1D("bkgzjetsnRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgzhbbRjj = new TH1D("bkgzhbbRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgzhccRjj = new TH1D("bkgzhccRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgzh4qRjj = new TH1D("bkgzh4qRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgzzqqRjj = new TH1D("bkgzzqqRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgzwqqRjj = new TH1D("bkgzwqqRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);
    TH1D * bkgttRjj = new TH1D("bkgttRjj", "Reconstucted Higgs Rjj", Rjjbins, 0.4, 1.5);

    TH1D * sigsubjetwidth = new TH1D("signsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgzjetssubjetwidth = new TH1D("bkgzjetsnsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgzhbbsubjetwidth = new TH1D("bkgzhbbsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgzhccsubjetwidth = new TH1D("bkgzhccsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgzh4qsubjetwidth = new TH1D("bkgzh4qsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgzzqqsubjetwidth = new TH1D("bkgzzqqsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgzwqqsubjetwidth = new TH1D("bkgzwqqsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);
    TH1D * bkgttsubjetwidth = new TH1D("bkgttsubjetwidth", "Reconstucted Higgs subjetwidth", subjetwidthbins, 0, 0.3);

    llgg(sigFile, sigHmass, sigtau21, sigtau42, sigthrust, sigthrust_minor, signtotal, sigpt, sigtruth_total, sigratio, sigasym, sigRjj, sigsubjetwidth);
    llgg(bkgzjetsFile, bkgzjetsHmass, bkgzjetstau21, bkgzjetstau42, bkgzjetsthrust, bkgzjetsthrust_minor, bkgzjetsntotal, bkgzjetspt, bkgzjetstruth_total, bkgzjetsratio, bkgzjetsasym, bkgzjetsRjj, bkgzjetssubjetwidth);
    llgg(bkgzhbbFile, bkgzhbbHmass, bkgzhbbtau21, bkgzhbbtau42, bkgzhbbthrust, bkgzhbbthrust_minor, bkgzhbbntotal, bkgzhbbpt, bkgzhbbtruth_total, bkgzhbbratio, bkgzhbbasym, bkgzhbbRjj, bkgzhbbsubjetwidth);
    llgg(bkgzhccFile, bkgzhccHmass, bkgzhcctau21, bkgzhcctau42, bkgzhccthrust, bkgzhccthrust_minor, bkgzhccntotal, bkgzhccpt, bkgzhcctruth_total, bkgzhccratio, bkgzhccasym, bkgzhccRjj, bkgzhccsubjetwidth);
    llgg(bkgzh4qFile, bkgzh4qHmass, bkgzh4qtau21, bkgzh4qtau42, bkgzh4qthrust, bkgzh4qthrust_minor, bkgzh4qntotal, bkgzh4qpt, bkgzh4qtruth_total, bkgzh4qratio, bkgzh4qasym, bkgzh4qRjj, bkgzh4qsubjetwidth);
    llgg(bkgzzqqFile, bkgzzqqHmass, bkgzzqqtau21, bkgzzqqtau42, bkgzzqqthrust, bkgzzqqthrust_minor, bkgzzqqntotal, bkgzzqqpt, bkgzzqqtruth_total, bkgzzqqratio, bkgzzqqasym, bkgzzqqRjj, bkgzzqqsubjetwidth);
    llgg(bkgzwqqFile, bkgzwqqHmass, bkgzwqqtau21, bkgzwqqtau42, bkgzwqqthrust, bkgzwqqthrust_minor, bkgzwqqntotal, bkgzwqqpt, bkgzwqqtruth_total, bkgzwqqratio, bkgzwqqasym, bkgzwqqRjj, bkgzwqqsubjetwidth);
    llgg(bkgttFile, bkgttHmass, bkgtttau21, bkgtttau42, bkgttthrust, bkgttthrust_minor, bkgttntotal, bkgttpt, bkgtttruth_total, bkgttratio, bkgttasym, bkgttRjj, bkgttsubjetwidth);
    
    THStack * Hmass = new THStack("Hmass", ";Reconstructed Higgs Invariant Mass [GeV]; Events/5 GeV");

    TH1 * sigHmass_clone = new TH1D( * sigHmass);

    Double_t sig_w = 0.0031132;
    Double_t sig_w_enh = 1000 * sig_w;
    Double_t bkgzjets_w = 92.79567345; //556.6953864;
    Double_t bkgzhbb_w = 0.06324864;
    Double_t bkgzhcc_w = 0.003139626;
    Double_t bkgzh4q_w = 0.01182654;
    Double_t bkgzzqq_w = 0.867;
    Double_t bkgzwqq_w = 0.9042;
    Double_t bkgtt_w = 231.579;

    Int_t sig = sigHmass -> GetEntries();
    Int_t bkgzjets = bkgzjetsHmass -> GetEntries();
    Int_t bkgzhbb = bkgzhbbHmass -> GetEntries();
    Int_t bkgzhcc = bkgzhccHmass -> GetEntries();
    Int_t bkgzh4q = bkgzh4qHmass -> GetEntries();
    Int_t bkgzzqq = bkgzzqqHmass -> GetEntries();
    Int_t bkgzwqq = bkgzwqqHmass -> GetEntries();
    Int_t bkgtt = bkgttHmass -> GetEntries();
    Double_t significance = sig * sig_w / TMath::Sqrt(bkgzjets * bkgzjets_w + bkgzhbb * bkgzhbb_w + bkgzhcc * bkgzhcc_w + bkgzh4q * bkgzh4q_w + bkgzzqq * bkgzzqq_w + bkgzwqq * bkgzwqq_w+ bkgtt * bkgtt_w);
    std::cout << "signal expected " << sig * sig_w << std::endl;
    std::cout << "zjets expected " << bkgzjets * bkgzjets_w << std::endl;
    std::cout << "zhbb expected " << bkgzhbb * bkgzhbb_w << std::endl;
    std::cout << "zhcc expected " << bkgzhcc * bkgzhcc_w << std::endl;
    std::cout << "zh4q expected " << bkgzh4q * bkgzh4q_w << std::endl;
    std::cout << "zzqq expected " << bkgzzqq * bkgzzqq_w << std::endl;
    std::cout << "zwqq expected " << bkgzwqq * bkgzwqq_w << std::endl;
    std::cout << "tt expected " << bkgtt * bkgtt_w << std::endl;
    std::cout << "Expected significance is " << significance << std::endl;
    sigHmass_clone -> SetLineColor(kPink - 1);
    sigHmass_clone -> Scale(sig_w);
    sigHmass_clone -> SetLineWidth(4);
    bkgzhccHmass -> SetFillColor(kAzure + 6);
    bkgzhccHmass -> Scale(bkgzhcc_w, "width");
    Hmass -> Add(bkgzhccHmass);
    bkgzh4qHmass -> SetFillColor(kSpring -9);
    bkgzh4qHmass -> Scale(bkgzh4q_w, "width");
    Hmass -> Add(bkgzh4qHmass);
    bkgzhbbHmass -> SetFillColor(kMagenta - 10);
    bkgzhbbHmass -> Scale(bkgzhbb_w, "width");
    Hmass -> Add(bkgzhbbHmass);
    bkgzzqqHmass -> SetFillColor(kViolet - 6);
    bkgzzqqHmass -> Scale(bkgzzqq_w, "width");
    Hmass -> Add(bkgzzqqHmass);
    bkgzwqqHmass -> SetFillColor(kOrange + 7);
    bkgzwqqHmass -> Scale(bkgzwqq_w, "width");
    Hmass -> Add(bkgzwqqHmass);
    bkgttHmass -> SetFillColor(kGray + 3);
    bkgttHmass -> Scale(bkgtt_w, "width");
    Hmass -> Add(bkgttHmass);
    bkgzjetsHmass -> SetFillColor(kOrange - 2);
    bkgzjetsHmass -> Scale(bkgzjets_w, "width");
    Hmass -> Add(bkgzjetsHmass);
    sigHmass -> SetFillColor(kPink - 1);
    sigHmass -> Scale(sig_w_enh, "width");
    Hmass -> Add(sigHmass);
    Hmass -> SetMaximum(1e8);
    Hmass -> SetMinimum(1e-3);
    
    TCanvas * totalcanvas = new TCanvas("totalcanvas", "Canvas", 1400, 1400, 1600, 1600);
    totalcanvas -> SetWindowSize(1800, 1400);
    totalcanvas -> SetCanvasSize(1200, 1200);
    totalcanvas -> SetLogy();
    TLegend * legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend -> AddEntry(sigHmass, "signal #times 1000", "f");
    legend -> AddEntry(sigHmass_clone, "signal", "l");
    legend -> AddEntry(bkgzjetsHmass, "bkg Z+jets", "f");
    legend -> AddEntry(bkgzhbbHmass, "bkg Z+bb", "f");
    legend -> AddEntry(bkgzh4qHmass, "bkg Z+H->4q", "f");
    legend -> AddEntry(bkgzhccHmass, "bkg Z+cc", "f");
    legend -> AddEntry(bkgzzqqHmass, "bkg ZZ", "f");
    legend -> AddEntry(bkgzwqqHmass, "bkg ZW", "f");
    legend -> AddEntry(bkgttHmass, "bkg tt", "f");
    Hmass -> Draw("HIST");
    Hmass -> GetHistogram() -> GetYaxis() -> SetTitleOffset(1.5);
    gPad -> Update();
    sigHmass_clone -> Draw("same");
    legend -> Draw("same");

    significance = 0.07;
    
    TLatex latexUp(25, 1e7, "#font[52]{#sqrt{s}} #font[42]{=} #font[52]{14 TeV}");
    TLatex latexDown(25, 1e6, "#font[52]{L} #font[42]{=} #font[52]{3000 fb^{-1}}");
    string init("#font[52]{#frac{S}{#sqrt{B}}} #font[42]{=} #font[52]{");
    string add = to_string(significance);
    string end("} #font[52]{(100#leqm_{H}#leq150)}");
    init = init + add.substr(0, 5) + end;
    const char * latex = init.c_str();
    TLatex latexSigStrength(25, 1e5, latex);
    latexUp.SetTextSize(0.025);
    latexDown.SetTextSize(0.025);
    latexSigStrength.SetTextSize(0.025);
    latexUp.Draw("same");
    latexDown.Draw("same");
    latexSigStrength.Draw("same");

    totalcanvas -> SaveAs("ZH_stack.png");

// tau21
    bkgzhcctau21 -> SetLineColor(kMagenta + 2);
    bkgzhcctau21 -> SetLineWidth(4);
    bkgzhcctau21 -> Scale(1.0/bkgzhcctau21->Integral());
    bkgtttau21 -> SetLineColor(kGray + 3);
    bkgtttau21 -> SetLineWidth(4);
    bkgtttau21 -> Scale(1.0/bkgtttau21->Integral());
    bkgzzqqtau21 -> SetLineColor(kViolet - 6);
    bkgzzqqtau21 -> SetLineWidth(4);
    bkgzzqqtau21 -> Scale(1.0/bkgzzqqtau21->Integral());
    bkgzwqqtau21 -> SetLineColor(kOrange + 7);
    bkgzwqqtau21 -> SetLineWidth(4);
    bkgzwqqtau21 -> Scale(1.0/bkgzwqqtau21->Integral());
    bkgzh4qtau21 -> SetLineColor(kTeal + 3);
    bkgzh4qtau21 -> SetLineWidth(4);
    bkgzh4qtau21 -> Scale(1.0/bkgzh4qtau21->Integral());
    bkgzhbbtau21 -> SetLineColor(kAzure + 2);
    bkgzhbbtau21 -> SetLineWidth(4);
    bkgzhbbtau21 -> Scale(1.0/bkgzhbbtau21->Integral());
    bkgzjetstau21 -> SetLineColor(kOrange - 2);
    bkgzjetstau21 -> SetLineWidth(4);
    bkgzjetstau21 -> Scale(1.0/bkgzjetstau21->Integral());
    sigtau21 -> SetLineColor(kPink - 1);
    sigtau21 -> SetLineWidth(4);
    sigtau21 -> Scale(1.0/sigtau21->Integral());

    TCanvas * tau21canvas_zjet = new TCanvas("tau21canvas_zjet", "Canvas", 1400, 1400, 1400, 1400);
    tau21canvas_zjet -> cd();
    tau21canvas_zjet -> SetWindowSize(1248, 1228);
    tau21canvas_zjet -> SetCanvasSize(1200, 1200);
    tau21canvas_zjet-> SetLogy(0);
    TLegend * legendtau21_zjet = new TLegend(0.1, 0.65, 0.3, 0.85);
    legendtau21_zjet -> AddEntry(sigtau21, "signal", "l");
    legendtau21_zjet -> AddEntry(bkgzjetstau21, "bkg Z+jets", "l");
    //legendtau21_zjet -> AddEntry(bkgzh4qtau21, "bkg Z+qqqq", "l");
    //legendtau21_zjet -> AddEntry(bkgzhbbtau21, "bkg Z+bb", "l");
    //legendtau21_zjet -> AddEntry(bkgzhcctau21, "bkg Z+cc", "l");
    //legendtau21_zjet -> AddEntry(bkgzzqqtau21, "bkg ZZ", "l");
    //legendtau21_zjet -> AddEntry(bkgzwqqtau21, "bkg ZW", "l");
    //legendtau21_zjet -> AddEntry(bkgtttau21, "bkg tt", "l");
    sigtau21 -> Draw("hist");
    sigtau21 -> SetStats(0);
    sigtau21->GetXaxis()->SetTitle("#tau_{21}");
    // bkgzh4qtau21 -> Draw("hist");
    // bkgzh4qtau21 -> SetStats(0);
    // bkgzh4qtau21->GetXaxis()->SetTitle("#tau_{21}");
    bkgzjetstau21 -> Draw("same hist");
    //bkgzjetstau21 -> Draw("same hist");
    //bkgzhbbtau21 -> Draw("same hist");
    //bkgzhcctau21 -> Draw("same hist");
    //bkgzzqqtau21 -> Draw("same hist");
    //bkgzwqqtau21 -> Draw("same hist");
    //bkgtttau21 -> Draw("same hist");
    legendtau21_zjet -> Draw("same hist");

    tau21canvas_zjet -> SaveAs("tau21_zjet_stack.png");

    TCanvas * tau21canvas_2boson = new TCanvas("tau21canvas_2boson", "Canvas", 1400, 1400, 1400, 1400);
    tau21canvas_2boson -> cd();
    tau21canvas_2boson -> SetWindowSize(1248, 1228);
    tau21canvas_2boson -> SetCanvasSize(1200, 1200);
    tau21canvas_2boson-> SetLogy(0);
    
    TH1D *bkg2bosontau21 = (TH1D*)bkgzzqqtau21->Clone("bkg2bosontau21");
    bkg2bosontau21 -> Add(bkgzwqqtau21);
    bkg2bosontau21 -> SetLineColor(kGray + 3);
    //bkg2bosontau21 -> SetLineWidth(4);
    //bkg2bosontau21 -> SetLineStyle(kDotted);
    bkg2bosontau21 -> Scale(1.0/bkg2bosontau21->Integral());

    TLegend * legendtau21_2boson = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendtau21_2boson -> AddEntry(sigtau21, "signal", "l");
    //legendtau21_2boson -> AddEntry(bkgzjetstau21, "bkg Z+jets", "l");
    legendtau21_2boson -> AddEntry(bkgzh4qtau21, "bkg Z+H->4q", "l");
    legendtau21_2boson -> AddEntry(bkg2bosontau21, "bkg ZZ/ZW", "l");
    //legendtau21_2boson -> AddEntry(bkgzhbbtau21, "bkg Z+bb", "l");
    //legendtau21_2boson -> AddEntry(bkgzhcctau21, "bkg Z+cc", "l");
    // legendtau21_2boson -> AddEntry(bkgzzqqtau21, "bkg ZZ", "l");
    // legendtau21_2boson -> AddEntry(bkgzwqqtau21, "bkg ZW", "l");
    //legendtau21_2boson -> AddEntry(bkgtttau21, "bkg tt", "l");

    bkgzh4qtau21 -> Draw("hist");
    bkgzh4qtau21 -> SetStats(0);
    bkgzh4qtau21->GetXaxis()->SetTitle("#tau_{21}");
    sigtau21 -> Draw("same hist");
    //bkgzjetstau21 -> Draw("same hist");
    //bkgzhbbtau21 -> Draw("same hist");
    //bkgzhcctau21 -> Draw("same hist");
    bkg2bosontau21 -> Draw("same hist");
    // bkgzzqqtau21 -> Draw("same hist");
    // bkgzwqqtau21 -> Draw("same hist");
    //bkgtttau21 -> Draw("same hist");
    legendtau21_2boson -> Draw("same hist");

    tau21canvas_2boson -> SaveAs("tau21_2boson_stack.png");

// tau42
    bkgzhcctau42 -> SetLineColor(kMagenta + 2);
    bkgzhcctau42 -> SetLineWidth(4);
    bkgzhcctau42 -> Scale(1.0/bkgzhcctau42->Integral());
    bkgtttau42 -> SetLineColor(kGray + 3);
    bkgtttau42 -> SetLineWidth(4);
    bkgtttau42 -> Scale(1.0/bkgtttau42->Integral());
    bkgzzqqtau42 -> SetLineColor(kViolet - 6);
    bkgzzqqtau42 -> SetLineWidth(4);
    bkgzzqqtau42 -> Scale(1.0/bkgzzqqtau42->Integral());
    bkgzwqqtau42 -> SetLineColor(kOrange + 7);
    bkgzwqqtau42 -> SetLineWidth(4);
    bkgzwqqtau42 -> Scale(1.0/bkgzwqqtau42->Integral());
    bkgzh4qtau42 -> SetLineColor(kTeal + 3);
    bkgzh4qtau42 -> SetLineWidth(4);
    bkgzh4qtau42 -> Scale(1.0/bkgzh4qtau42->Integral());
    bkgzhbbtau42 -> SetLineColor(kAzure + 2);
    bkgzhbbtau42 -> SetLineWidth(4);
    bkgzhbbtau42 -> Scale(1.0/bkgzhbbtau42->Integral());
    bkgzjetstau42 -> SetLineColor(kOrange - 2);
    bkgzjetstau42 -> SetLineWidth(4);
    bkgzjetstau42 -> Scale(1.0/bkgzjetstau42->Integral());
    sigtau42 -> SetLineColor(kPink - 1);
    sigtau42 -> SetLineWidth(4);
    sigtau42 -> Scale(1.0/sigtau42->Integral());

    TCanvas * tau42canvas = new TCanvas("tau42canvas", "Canvas", 1400, 1400, 1400, 1400);
    tau42canvas -> SetWindowSize(1248, 1228);
    tau42canvas -> SetCanvasSize(1200, 1200);
    tau42canvas -> SetLogy(0);
    tau42canvas -> cd();
    TLegend * legendtau42 = new TLegend(0.15, 0.65, 0.35, 0.85);
    legendtau42 -> AddEntry(sigtau42, "signal", "l");
    //legendtau42 -> AddEntry(bkgzjetstau42, "bkg Z+jets", "l");
    //legendtau42 -> AddEntry(bkgzhbbtau42, "bkg Z+bb", "l");
    legendtau42 -> AddEntry(bkgzh4qtau42, "bkg Z+qqqq", "l");
    //legendtau42 -> AddEntry(bkgzhcctau42, "bkg Z+cc", "l");
    //legendtau42 -> AddEntry(bkgzzqqtau42, "bkg ZZ", "l");
    //legendtau42 -> AddEntry(bkgzwqqtau42, "bkg ZW", "l");
    //legendtau42 -> AddEntry(bkgtttau42, "bkg tt", "l");
    sigtau42 -> Draw("hist");
    sigtau42 -> SetStats(0);
    sigtau42 ->GetXaxis()->SetTitle("#tau_{42}");
    //bkgzjetstau42 -> Draw("same hist");
    //bkgzhbbtau42 -> Draw("same hist");
    //bkgzhcctau42 -> Draw("same hist");
    bkgzh4qtau42 -> Draw("same hist");
    //bkgzzqqtau42 -> Draw("same hist");
    //bkgzwqqtau42 -> Draw("same hist");
    //bkgtttau42 -> Draw("same hist");
    legendtau42 -> Draw("same hist");

    tau42canvas -> SaveAs("tau42_stack.png");

// thrust
    bkgzhccthrust -> SetLineColor(kMagenta + 2);
    bkgzhccthrust -> SetLineWidth(4);
    bkgzhccthrust -> Scale(1.0/bkgzhccthrust->Integral());
    bkgttthrust -> SetLineColor(kGray + 3);
    bkgttthrust -> SetLineWidth(4);
    bkgttthrust -> Scale(1.0/bkgttthrust->Integral());
    bkgzzqqthrust -> SetLineColor(kViolet - 6);
    bkgzzqqthrust -> SetLineWidth(4);
    bkgzzqqthrust -> Scale(1.0/bkgzzqqthrust->Integral());
    bkgzwqqthrust -> SetLineColor(kOrange + 7);
    bkgzwqqthrust -> SetLineWidth(4);
    bkgzwqqthrust -> Scale(1.0/bkgzwqqthrust->Integral());
    bkgzh4qthrust -> SetLineColor(kTeal + 3);
    bkgzh4qthrust -> SetLineWidth(4);
    bkgzh4qthrust -> Scale(1.0/bkgzh4qthrust->Integral());
    bkgzhbbthrust -> SetLineColor(kAzure + 2);
    bkgzhbbthrust -> SetLineWidth(4);
    bkgzhbbthrust -> Scale(1.0/bkgzhbbthrust->Integral());
    bkgzjetsthrust -> SetLineColor(kOrange - 2);
    bkgzjetsthrust -> SetLineWidth(4);
    bkgzjetsthrust -> Scale(1.0/bkgzjetsthrust->Integral());
    sigthrust -> SetLineColor(kPink - 1);
    sigthrust -> SetLineWidth(4);
    sigthrust -> Scale(1.0/sigthrust->Integral());

    TCanvas * thrustcanvas = new TCanvas("thrustcanvas", "Canvas", 1400, 1400, 1400, 1400);
    thrustcanvas -> SetWindowSize(1248, 1228);
    thrustcanvas -> SetCanvasSize(1200, 1200);
    thrustcanvas -> SetLogy(0);
    thrustcanvas -> cd();
    TLegend * legendthrust = new TLegend(0.15, 0.65, 0.35, 0.85);
    legendthrust -> AddEntry(sigthrust, "signal", "l");
    //legendthrust -> AddEntry(bkgzjetsthrust, "bkg Z+jets", "l");
    //legendthrust -> AddEntry(bkgzhbbthrust, "bkg Z+bb", "l");
    legendthrust -> AddEntry(bkgzh4qthrust, "bkg Z+qqqq", "l");
    //legendthrust -> AddEntry(bkgzhccthrust, "bkg Z+cc", "l");
    //legendthrust -> AddEntry(bkgzzqqthrust, "bkg ZZ", "l");
    //legendthrust -> AddEntry(bkgzwqqthrust, "bkg ZW", "l");
    //legendthrust -> AddEntry(bkgttthrust, "bkg tt", "l");
    sigthrust -> Draw("hist");
    sigthrust -> SetStats(0);
    sigthrust ->GetXaxis()->SetTitle("Thrust");
    // bkgzhccthrust -> Draw("hist");
    // bkgzhccthrust -> SetStats(0);
    // bkgzhccthrust ->GetXaxis()->SetTitle("Thrust");
    bkgzh4qthrust-> Draw("same hist");
    // bkgzjetsthrust -> Draw("same hist");
    // bkgzhbbthrust -> Draw("same hist");
    // bkgzh4qthrust -> Draw("same hist");
    // bkgzzqqthrust -> Draw("same hist");
    // bkgzwqqthrust -> Draw("same hist");
    // bkgttthrust -> Draw("same hist");
    legendthrust -> Draw("same hist");

    thrustcanvas -> SaveAs("thrust_stack.png");

// thrust_minor
    bkgzhccthrust_minor -> SetLineColor(kMagenta + 2);
    bkgzhccthrust_minor -> SetLineWidth(4);
    bkgzhccthrust_minor -> Scale(1.0/bkgzhccthrust_minor->Integral());
    bkgttthrust_minor -> SetLineColor(kGray + 3);
    bkgttthrust_minor -> SetLineWidth(4);
    bkgttthrust_minor -> Scale(1.0/bkgttthrust_minor->Integral());
    bkgzzqqthrust_minor -> SetLineColor(kViolet - 6);
    bkgzzqqthrust_minor -> SetLineWidth(4);
    bkgzzqqthrust_minor -> Scale(1.0/bkgzzqqthrust_minor->Integral());
    bkgzwqqthrust_minor -> SetLineColor(kOrange + 7);
    bkgzwqqthrust_minor -> SetLineWidth(4);
    bkgzwqqthrust_minor -> Scale(1.0/bkgzwqqthrust_minor->Integral());
    bkgzh4qthrust_minor -> SetLineColor(kTeal + 3);
    bkgzh4qthrust_minor -> SetLineWidth(4);
    bkgzh4qthrust_minor -> Scale(1.0/bkgzh4qthrust_minor->Integral());
    bkgzhbbthrust_minor -> SetLineColor(kAzure + 2);
    bkgzhbbthrust_minor -> SetLineWidth(4);
    bkgzhbbthrust_minor -> Scale(1.0/bkgzhbbthrust_minor->Integral());
    bkgzjetsthrust_minor -> SetLineColor(kOrange - 2);
    bkgzjetsthrust_minor -> SetLineWidth(4);
    bkgzjetsthrust_minor -> Scale(1.0/bkgzjetsthrust_minor->Integral());
    sigthrust_minor -> SetLineColor(kPink - 1);
    sigthrust_minor -> SetLineWidth(4);
    sigthrust_minor -> Scale(1.0/sigthrust_minor->Integral());

    TCanvas * thrust_minorcanvas = new TCanvas("thrust_minorcanvas", "Canvas", 1400, 1400, 1400, 1400);
    thrust_minorcanvas -> SetWindowSize(1248, 1228);
    thrust_minorcanvas -> SetCanvasSize(1200, 1200);
    thrust_minorcanvas -> SetLogy(0);
    thrust_minorcanvas -> cd();
    TLegend * legendthrust_minor = new TLegend(0.15, 0.65, 0.35, 0.85);
    legendthrust_minor -> AddEntry(sigthrust_minor, "signal", "l");
    //legendthrust_minor -> AddEntry(bkgzjetsthrust_minor, "bkg Z+jets", "l");
    //legendthrust_minor -> AddEntry(bkgzhbbthrust_minor, "bkg Z+bb", "l");
    legendthrust_minor -> AddEntry(bkgzh4qthrust_minor, "bkg Z+qqqq", "l");
    //legendthrust_minor -> AddEntry(bkgzhccthrust_minor, "bkg Z+cc", "l");
    //legendthrust_minor -> AddEntry(bkgzzqqthrust_minor, "bkg ZZ", "l");
    //legendthrust_minor -> AddEntry(bkgzwqqthrust_minor, "bkg ZW", "l");
    //legendthrust_minor -> AddEntry(bkgttthrust_minor, "bkg tt", "l");
    bkgzh4qthrust_minor -> Draw("hist");
    bkgzh4qthrust_minor -> SetStats(0);
    bkgzh4qthrust_minor ->GetXaxis()->SetTitle("thrust_minor");
    //bkgzhccthrust_minor -> Draw("same hist");
    sigthrust_minor -> Draw("same hist");
    //bkgzjetsthrust_minor -> Draw("same hist");
    //bkgzhbbthrust_minor -> Draw("same hist");
    //bkgzzqqthrust_minor -> Draw("same hist");
    //bkgzwqqthrust_minor -> Draw("same hist");
    //bkgttthrust_minor -> Draw("same hist");
    legendthrust_minor -> Draw("same hist");

    thrust_minorcanvas -> SaveAs("thrust_minor_stack.png");

// truth_total
    bkgzhcctruth_total -> SetLineColor(kMagenta + 2);
    bkgzhcctruth_total -> SetLineWidth(4);
    bkgzhcctruth_total -> Scale(1.0/bkgzhcctruth_total->Integral());
    bkgzh4qtruth_total -> SetLineColor(kTeal + 3);
    bkgzh4qtruth_total -> SetLineWidth(4);
    bkgzh4qtruth_total -> Scale(1.0/bkgzh4qtruth_total->Integral());
    bkgtttruth_total -> SetLineColor(kGray + 3);
    bkgtttruth_total -> SetLineWidth(4);
    bkgtttruth_total -> Scale(1.0/bkgtttruth_total->Integral());
    bkgzzqqtruth_total -> SetLineColor(kViolet - 6);
    bkgzzqqtruth_total -> SetLineWidth(4);
    bkgzzqqtruth_total -> Scale(1.0/bkgzzqqtruth_total->Integral());
    bkgzwqqtruth_total -> SetLineColor(kOrange + 7);
    bkgzwqqtruth_total -> SetLineWidth(4);
    bkgzwqqtruth_total -> Scale(1.0/bkgzwqqtruth_total->Integral());
    bkgzhbbtruth_total -> SetLineColor(kAzure + 2);
    bkgzhbbtruth_total -> SetLineWidth(4);
    bkgzhbbtruth_total -> Scale(1.0/bkgzhbbtruth_total->Integral());
    bkgzjetstruth_total -> SetLineColor(kOrange - 2);
    bkgzjetstruth_total -> SetLineWidth(4);
    bkgzjetstruth_total -> Scale(1.0/bkgzjetstruth_total->Integral());
    sigtruth_total -> SetLineColor(kPink - 1);
    sigtruth_total -> SetLineWidth(4);
    sigtruth_total -> Scale(1.0/sigtruth_total->Integral());

    TCanvas * truth_totalcanvas = new TCanvas("truth_totalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    truth_totalcanvas -> SetWindowSize(1248, 1228);
    truth_totalcanvas -> SetCanvasSize(1200, 1200);
    truth_totalcanvas -> SetLogy(0);
    truth_totalcanvas -> cd();
    TLegend * legendtruth_total = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendtruth_total -> AddEntry(sigtruth_total, "signal", "l");
    legendtruth_total -> AddEntry(bkgzjetstruth_total, "bkg Z+jets", "l");
    //legendtruth_total -> AddEntry(bkgzhbbtruth_total, "bkg Z+bb", "l");
    legendtruth_total -> AddEntry(bkgzh4qtruth_total, "bkg Z+qqqq", "l");
    //legendtruth_total -> AddEntry(bkgzhcctruth_total, "bkg Z+cc", "l");
    legendtruth_total -> AddEntry(bkgzzqqtruth_total, "bkg ZZ", "l");
    legendtruth_total -> AddEntry(bkgzwqqtruth_total, "bkg ZW", "l");
    //legendtruth_total -> AddEntry(bkgtttruth_total, "bkg tt", "l");
    bkgzh4qtruth_total -> Draw("hist");
    bkgzh4qtruth_total -> SetStats(0);
    bkgzh4qtruth_total ->GetXaxis()->SetTitle("truth_total");
    //bkgzhcctruth_total -> Draw("same hist");
    sigtruth_total -> Draw("same hist");
    bkgzjetstruth_total -> Draw("same hist");
    //bkgzhbbtruth_total -> Draw("same hist");
    bkgzzqqtruth_total -> Draw("same hist");
    bkgzwqqtruth_total -> Draw("same hist");
    //bkgtttruth_total -> Draw("same hist");
    legendtruth_total -> Draw("same hist");

    truth_totalcanvas -> SaveAs("truth_total_stack.png");

// ntotal
    bkgzhccntotal -> SetLineColor(kMagenta + 2);
    bkgzhccntotal -> SetLineWidth(4);
    bkgzhccntotal -> Scale(1.0/bkgzhccntotal->Integral());
    bkgzh4qntotal -> SetLineColor(kTeal + 3);
    bkgzh4qntotal -> SetLineWidth(4);
    bkgzh4qntotal -> Scale(1.0/bkgzh4qntotal->Integral());
    bkgttntotal -> SetLineColor(kGray + 3);
    bkgttntotal -> SetLineWidth(4);
    bkgttntotal -> Scale(1.0/bkgttntotal->Integral());
    // bkgzzqqntotal -> SetLineColor(kViolet - 6);
    bkgzzqqntotal -> SetLineColor(kGray + 3);
    bkgzzqqntotal -> SetLineWidth(4);
    //bkgzzqqntotal -> SetLineStyle(kDotted);
    bkgzzqqntotal -> Scale(1.0/bkgzzqqntotal->Integral());
    // bkgzwqqntotal -> SetLineColor(kOrange + 7);
    bkgzwqqntotal -> SetLineColor(kGray + 3);
    bkgzwqqntotal -> SetLineWidth(4);
    //bkgzwqqntotal -> SetLineStyle(kDotted);
    bkgzwqqntotal -> Scale(1.0/bkgzwqqntotal->Integral());
    bkgzhbbntotal -> SetLineColor(kAzure + 2);
    bkgzhbbntotal -> SetLineWidth(4);
    bkgzhbbntotal -> Scale(1.0/bkgzhbbntotal->Integral());
    // bkgzjetsntotal -> SetLineColor(kOrange - 2);
    bkgzjetsntotal -> SetLineColor(kBlack);
    bkgzjetsntotal -> SetLineWidth(4);
    bkgzjetsntotal -> Scale(1.0/bkgzjetsntotal->Integral());
    signtotal -> SetLineColor(kPink - 1);
    signtotal -> SetLineWidth(4);
    signtotal -> Scale(1.0/signtotal->Integral());

    TCanvas * ntotalcanvas = new TCanvas("ntotalcanvas", "Canvas", 1400, 1400, 1400, 1400);
    ntotalcanvas -> SetWindowSize(1248, 1228);
    ntotalcanvas -> SetCanvasSize(1200, 1200);
    ntotalcanvas -> SetLogy(0);
    ntotalcanvas -> cd();

    
    TH1D *bkg2bosonntotal = (TH1D*)bkgzzqqntotal->Clone("bkg2bosonntotal");
    bkg2bosonntotal -> Add(bkgzwqqntotal);
    bkg2bosonntotal -> SetLineColor(kGray + 3);
    bkg2bosonntotal -> SetLineWidth(4);
    bkg2bosonntotal -> SetLineStyle(kDotted);
    bkg2bosonntotal -> Scale(1.0/bkg2bosonntotal->Integral());

    TLegend * legendntotal = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendntotal -> AddEntry(signtotal, "signal", "l");
    legendntotal -> AddEntry(bkgzjetsntotal, "bkg Z+jets", "l");
    legendntotal -> AddEntry(bkg2bosonntotal, "bkg ZZ/ZW", "l");
    // legendntotal -> AddEntry(bkgzzqqntotal, "bkg ZZ", "l");
    // legendntotal -> AddEntry(bkgzwqqntotal, "bkg ZW", "l");
    // legendntotal -> AddEntry(bkgzhbbntotal, "bkg Z+bb", "l");
    legendntotal -> AddEntry(bkgzh4qntotal, "bkg Z+H->4q", "l");
    // legendntotal -> AddEntry(bkgzhccntotal, "bkg Z+cc", "l");
    // bkgzh4qntotal -> Draw("hist");
    // bkgzh4qntotal -> SetStats(0);
    // bkgzh4qntotal ->GetXaxis()->SetTitle("ntotal");

    bkg2bosonntotal -> Draw("hist");
    bkg2bosonntotal -> SetStats(0);
    bkg2bosonntotal ->GetXaxis()->SetTitle("ntotal");
    // bkgzhccntotal -> Draw("same hist");
    signtotal -> Draw("same hist");
    bkgzjetsntotal -> Draw("same hist");
    bkgzh4qntotal -> Draw("same hist");
    // bkgzwqqntotal -> Draw("same hist");
    legendntotal -> Draw("same hist");

    ntotalcanvas -> SaveAs("ntotal_stack.png");    

// pt
    bkgzhccpt -> SetLineColor(kMagenta + 2);
    bkgzhccpt -> SetLineWidth(4);
    bkgzhccpt -> Scale(1.0/bkgzhccpt->Integral());
    bkgzh4qpt -> SetLineColor(kTeal + 3);
    bkgzh4qpt -> SetLineWidth(4);
    bkgzh4qpt -> Scale(1.0/bkgzh4qpt->Integral());

    bkgttpt -> SetLineColor(kGray + 3);
    bkgttpt -> SetLineWidth(4);
    bkgttpt -> Scale(1.0/bkgttpt->Integral());
    bkgzzqqpt -> SetLineColor(kViolet - 6);
    bkgzzqqpt -> SetLineWidth(4);
    bkgzzqqpt -> Scale(1.0/bkgzzqqpt->Integral());
    bkgzwqqpt -> SetLineColor(kOrange + 7);
    bkgzwqqpt -> SetLineWidth(4);
    bkgzwqqpt -> Scale(1.0/bkgzwqqpt->Integral());

    bkgzhbbpt -> SetLineColor(kAzure + 2);
    bkgzhbbpt -> SetLineWidth(4);
    bkgzhbbpt -> Scale(1.0/bkgzhbbpt->Integral());
    bkgzjetspt -> SetLineColor(kOrange + 7);
    bkgzjetspt -> SetLineWidth(4);
    bkgzjetspt -> Scale(1.0/bkgzjetspt->Integral());
    sigpt -> SetLineColor(kPink - 1);
    sigpt -> SetLineWidth(4);
    sigpt -> Scale(1.0/sigpt->Integral());

    TCanvas * ptcanvas = new TCanvas("ptcanvas", "Canvas", 1400, 1400, 1400, 1400);
    ptcanvas -> SetWindowSize(1248, 1228);
    ptcanvas -> SetCanvasSize(1200, 1200);
    ptcanvas -> SetLogy(0);
    ptcanvas -> cd();
    TLegend * legendpt = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendpt -> AddEntry(sigpt, "signal", "l");
    legendpt -> AddEntry(bkgzjetspt, "bkg Z+jets", "l");
    legendpt -> AddEntry(bkgzhbbpt, "bkg Z+bb", "l");
    legendpt -> AddEntry(bkgzh4qpt, "bkg Z+qqqq", "l");
    legendpt -> AddEntry(bkgzhccpt, "bkg Z+cc", "l");
    bkgzh4qpt -> Draw("hist");
    bkgzh4qpt -> SetStats(0);
    bkgzh4qpt ->GetXaxis()->SetTitle("pt");
    bkgzhccpt -> Draw("same hist");
    sigpt -> Draw("same hist");
    bkgzjetspt -> Draw("same hist");
    bkgzhbbpt -> Draw("same hist");
    bkgzzqqpt -> Draw("same hist");
    bkgzwqqpt -> Draw("same hist");
    bkgttpt -> Draw("same hist");
    legendpt -> Draw("same hist");

    ptcanvas -> SaveAs("pt_stack.png");

// ratio
    bkgzhccratio -> SetLineColor(kMagenta + 2);
    bkgzhccratio -> SetLineWidth(4);
    bkgzhccratio -> Scale(1.0/bkgzhccratio->Integral());
    bkgzh4qratio -> SetLineColor(kTeal + 3);
    bkgzh4qratio -> SetLineWidth(4);
    bkgzh4qratio -> Scale(1.0/bkgzh4qratio->Integral());
    bkgttratio -> SetLineColor(kGray + 3);
    bkgttratio -> SetLineWidth(4);
    bkgttratio -> Scale(1.0/bkgttratio->Integral());
    bkgzzqqratio -> SetLineColor(kViolet - 6);
    bkgzzqqratio -> SetLineWidth(4);
    bkgzzqqratio -> Scale(1.0/bkgzzqqratio->Integral());
    bkgzwqqratio -> SetLineColor(kOrange + 7);
    bkgzwqqratio -> SetLineWidth(4);
    bkgzwqqratio -> Scale(1.0/bkgzwqqratio->Integral());
    bkgzhbbratio -> SetLineColor(kAzure + 2);
    bkgzhbbratio -> SetLineWidth(4);
    bkgzhbbratio -> Scale(1.0/bkgzhbbratio->Integral());
    bkgzjetsratio -> SetLineColor(kOrange - 2);
    bkgzjetsratio -> SetLineWidth(4);
    bkgzjetsratio -> Scale(1.0/bkgzjetsratio->Integral());
    sigratio -> SetLineColor(kPink - 1);
    sigratio -> SetLineWidth(4);
    sigratio -> Scale(1.0/sigratio->Integral());

    TCanvas * ratiocanvas = new TCanvas("ratiocanvas", "Canvas", 1400, 1400, 1400, 1400);
    ratiocanvas -> SetWindowSize(1248, 1228);
    ratiocanvas -> SetCanvasSize(1200, 1200);
    ratiocanvas -> SetLogy(0);
    ratiocanvas -> SetLogx(0);
    ratiocanvas -> cd();
    TLegend * legendratio = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendratio -> AddEntry(sigratio, "signal", "l");
    legendratio -> AddEntry(bkgzjetsratio, "bkg Z+jets", "l");
    //legendratio -> AddEntry(bkgzhbbratio, "bkg Z+bb", "l");
    //legendratio -> AddEntry(bkgzh4qratio, "bkg Z+qqqq", "l");
    legendratio -> AddEntry(bkgzhccratio, "bkg Z+cc", "l");
    bkgzjetsratio -> Draw("hist");
    bkgzjetsratio -> SetStats(0);
    bkgzjetsratio ->GetXaxis()->SetTitle("ratio");
    sigratio -> Draw("same hist");
    //bkgzhbbratio -> Draw("same hist");
    //bkgzhccratio -> Draw("same hist");
    bkgzwqqratio -> Draw("same hist");
    bkgzzqqratio -> Draw("same hist");
    legendratio -> Draw("same hist");

    ratiocanvas -> SaveAs("ratio_stack.png");

// asym
    bkgzhccasym -> SetLineColor(kMagenta + 2);
    bkgzhccasym -> SetLineWidth(4);
    bkgzhccasym -> Scale(1.0/bkgzhccasym->Integral());
    bkgzh4qasym -> SetLineColor(kTeal + 3);
    bkgzh4qasym -> SetLineWidth(4);
    bkgzh4qasym -> Scale(1.0/bkgzh4qasym->Integral());
    bkgttasym -> SetLineColor(kGray + 3);
    bkgttasym -> SetLineWidth(4);
    bkgttasym -> Scale(1.0/bkgttasym->Integral());
    bkgzzqqasym -> SetLineColor(kViolet - 6);
    bkgzzqqasym -> SetLineWidth(4);
    bkgzzqqasym -> Scale(1.0/bkgzzqqasym->Integral());
    bkgzwqqasym -> SetLineColor(kOrange + 7);
    bkgzwqqasym -> SetLineWidth(4);
    bkgzwqqasym -> Scale(1.0/bkgzwqqasym->Integral());
    bkgzhbbasym -> SetLineColor(kAzure + 2);
    bkgzhbbasym -> SetLineWidth(4);
    bkgzhbbasym -> Scale(1.0/bkgzhbbasym->Integral());
    bkgzjetsasym -> SetLineColor(kOrange - 2);
    bkgzjetsasym -> SetLineWidth(4);
    bkgzjetsasym -> Scale(1.0/bkgzjetsasym->Integral());
    sigasym -> SetLineColor(kPink - 1);
    sigasym -> SetLineWidth(4);
    sigasym -> Scale(1.0/sigasym->Integral());

    TCanvas * asymcanvas = new TCanvas("asymcanvas", "Canvas", 1400, 1400, 1400, 1400);
    asymcanvas -> SetWindowSize(1248, 1228);
    asymcanvas -> SetCanvasSize(1200, 1200);
    asymcanvas -> SetLogy(0);
    asymcanvas -> SetLogx(0);
    asymcanvas -> cd();
    TLegend * legendasym = new TLegend(0.15, 0.65, 0.35, 0.85);
    legendasym -> AddEntry(sigasym, "signal", "l");
    legendasym -> AddEntry(bkgzjetsasym, "bkg Z+jets", "l");
    //legendasym -> AddEntry(bkgzhbbasym, "bkg Z+bb", "l");
    // -> AddEntry(bkgzh4qasym, "bkg Z+qqqq", "l");
    //legendasym -> AddEntry(bkgzhccasym, "bkg Z+cc", "l");
    bkgzjetsasym -> Draw("hist");
    bkgzjetsasym -> SetStats(0);
    bkgzjetsasym -> GetXaxis()->SetTitle("asym");
    sigasym -> Draw("same hist");
    //bkgzhbbasym -> Draw("same hist");
    //bkgzhccasym -> Draw("same hist");
    //bkgzzqqasym -> Draw("same hist");
    //bkgzwqqasym -> Draw("same hist");
    legendasym -> Draw("same hist");

    asymcanvas -> SaveAs("asym_stack.png");

    // Rjj
    bkgzhccRjj -> SetLineColor(kMagenta + 2);
    bkgzhccRjj -> SetLineWidth(4);
    bkgzhccRjj -> Scale(1.0/bkgzhccRjj->Integral());
    bkgzh4qRjj -> SetLineColor(kTeal + 3);
    bkgzh4qRjj -> SetLineWidth(4);
    bkgzh4qRjj -> Scale(1.0/bkgzh4qRjj->Integral());

    bkgttRjj -> SetLineColor(kGray + 3);
    bkgttRjj -> SetLineWidth(4);
    bkgttRjj -> Scale(1.0/bkgttRjj->Integral());
    bkgzzqqRjj -> SetLineColor(kViolet - 6);
    bkgzzqqRjj -> SetLineWidth(4);
    bkgzzqqRjj -> Scale(1.0/bkgzzqqRjj->Integral());
    bkgzwqqRjj -> SetLineColor(kOrange + 7);
    bkgzwqqRjj -> SetLineWidth(4);
    bkgzwqqRjj -> Scale(1.0/bkgzwqqRjj->Integral());

    bkgzhbbRjj -> SetLineColor(kAzure + 2);
    bkgzhbbRjj -> SetLineWidth(4);
    bkgzhbbRjj -> Scale(1.0/bkgzhbbRjj->Integral());
    bkgzjetsRjj -> SetLineColor(kOrange - 2);
    bkgzjetsRjj -> SetLineWidth(4);
    bkgzjetsRjj -> Scale(1.0/bkgzjetsRjj->Integral());
    sigRjj -> SetLineColor(kPink - 1);
    sigRjj -> SetLineWidth(4);
    sigRjj -> Scale(1.0/sigRjj->Integral());

    TCanvas * Rjjcanvas = new TCanvas("Rjjcanvas", "Canvas", 1400, 1400, 1400, 1400);
    Rjjcanvas -> SetWindowSize(1248, 1228);
    Rjjcanvas -> SetCanvasSize(1200, 1200);
    Rjjcanvas -> SetLogy(0);
    Rjjcanvas -> SetLogx(0);
    Rjjcanvas -> cd();
    
    TH1D *bkg2bosonRjj = (TH1D*)bkgzzqqRjj->Clone("bkg2bosonRjj");
    bkg2bosonRjj -> Add(bkgzwqqRjj);
    bkg2bosonRjj -> SetLineColor(kGray + 3);
    //bkg2bosonRjj -> SetLineWidth(4);
    //bkg2bosonRjj -> SetLineStyle(kDotted);
    bkg2bosonRjj -> Scale(1.0/bkg2bosonRjj->Integral());

    TLegend * legendRjj = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendRjj -> AddEntry(sigRjj, "signal", "l");
    legendRjj -> AddEntry(bkgzjetsRjj, "bkg Z+jets", "l");
    legendRjj -> AddEntry(bkg2bosonRjj, "bkg ZZ/ZW", "l");
    // legendRjj -> AddEntry(bkgzzqqRjj, "bkg ZZ", "l");
    // legendRjj -> AddEntry(bkgzwqqRjj, "bkg ZW", "l");
    // legendRjj -> AddEntry(bkgzhbbRjj, "bkg Z+bb", "l");
    // legendRjj -> AddEntry(bkgzh4qRjj, "bkg Z+qqqq", "l");

    bkgzjetsRjj -> Draw("hist");
    bkgzjetsRjj -> SetStats(0);
    bkgzjetsRjj ->GetXaxis()->SetTitle("Rjj");
    sigRjj -> Draw("same hist");
    //bkgzhbbRjj -> Draw("same hist");
    //bkgzhccRjj -> Draw("same hist");
    //bkgzzqqRjj -> Draw("same hist");
    bkg2bosonRjj -> Draw("same hist");
    legendRjj -> Draw("same hist");

    Rjjcanvas -> SaveAs("Rjj_stack.png");

// subjetwidth
    bkgzhccsubjetwidth -> SetLineColor(kMagenta + 2);
    bkgzhccsubjetwidth -> SetLineWidth(4);
    bkgzhccsubjetwidth -> Scale(1.0/bkgzhccsubjetwidth->Integral());
    bkgzh4qsubjetwidth -> SetLineColor(kTeal + 3);
    bkgzh4qsubjetwidth -> SetLineWidth(4);
    bkgzh4qsubjetwidth -> Scale(1.0/bkgzh4qsubjetwidth->Integral());

    bkgttsubjetwidth -> SetLineColor(kGray + 3);
    bkgttsubjetwidth -> SetLineWidth(4);
    bkgttsubjetwidth -> Scale(1.0/bkgttsubjetwidth->Integral());
    bkgzzqqsubjetwidth -> SetLineColor(kViolet - 6);
    bkgzzqqsubjetwidth -> SetLineWidth(4);
    bkgzzqqsubjetwidth -> Scale(1.0/bkgzzqqsubjetwidth->Integral());
    bkgzwqqsubjetwidth -> SetLineColor(kOrange + 7);
    bkgzwqqsubjetwidth -> SetLineWidth(4);
    bkgzwqqsubjetwidth -> Scale(1.0/bkgzwqqsubjetwidth->Integral());

    bkgzhbbsubjetwidth -> SetLineColor(kAzure + 2);
    bkgzhbbsubjetwidth -> SetLineWidth(4);
    bkgzhbbsubjetwidth -> Scale(1.0/bkgzhbbsubjetwidth->Integral());
    bkgzjetssubjetwidth -> SetLineColor(kOrange - 2);
    bkgzjetssubjetwidth -> SetLineWidth(4);
    bkgzjetssubjetwidth -> Scale(1.0/bkgzjetssubjetwidth->Integral());
    sigsubjetwidth -> SetLineColor(kPink - 1);
    sigsubjetwidth -> SetLineWidth(4);
    sigsubjetwidth -> Scale(1.0/sigsubjetwidth->Integral());

    TCanvas * subjetwidthcanvas = new TCanvas("subjetwidthcanvas", "Canvas", 1400, 1400, 1400, 1400);
    subjetwidthcanvas -> SetWindowSize(1248, 1228);
    subjetwidthcanvas -> SetCanvasSize(1200, 1200);
    subjetwidthcanvas -> SetLogy(0);
    subjetwidthcanvas -> SetLogx(0);
    subjetwidthcanvas -> cd();
    TLegend * legendsubjetwidth = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendsubjetwidth -> AddEntry(sigsubjetwidth, "signal", "l");
    legendsubjetwidth -> AddEntry(bkgzjetssubjetwidth, "bkg Z+jets", "l");
    //legendsubjetwidth -> AddEntry(bkgzhbbsubjetwidth, "bkg Z+bb", "l");
    //legendsubjetwidth -> AddEntry(bkgzh4qsubjetwidth, "bkg Z+qqqq", "l");
    //legendsubjetwidth -> AddEntry(bkgzhccsubjetwidth, "bkg Z+cc", "l");
    bkgzjetssubjetwidth -> Draw("hist");
    bkgzjetssubjetwidth -> SetStats(0);
    bkgzjetssubjetwidth ->GetXaxis()->SetTitle("subjetwidth");
    sigsubjetwidth -> Draw("same hist");
    //bkgzhbbsubjetwidth -> Draw("same hist");
    //bkgzhccsubjetwidth -> Draw("same hist");
    //bkgzzqqsubjetwidth -> Draw("same hist");
    //bkgzwqqsubjetwidth -> Draw("same hist");
    legendsubjetwidth -> Draw("same hist");

    subjetwidthcanvas -> SaveAs("subjetwidth_stack.png");

//end    
    output -> cd();
    sigHmass -> Write();
    bkgzjetsHmass -> Write();
    bkgzhbbHmass -> Write();
    bkgzhccHmass -> Write();
    bkgzh4qHmass -> Write();
    bkgzzqqHmass -> Write();
    bkgzwqqHmass -> Write();
    bkgttHmass -> Write();
    Hmass -> Write();
    tree_output -> Write();
    output -> Close();
}
