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
                            TH1D * pt) {

    Double_t ptcut = 300;
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
        if (h_gg.M() < 100 or h_gg.M() > 160) {
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
        for (int j = 0; j < nsmallJet; ++j) {
            smalljet = (Jet * ) branchsmallJet -> At(j);
            if (h_gg.DeltaR(smalljet -> P4()) < 1.0) {
                nsmallB++;
                if (smalljet -> BTag >= 1) {
                    smallBtag = true;
                }
            }
        }

        if (smallBtag == true) {
            continue;
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
                const char * bkgzh4qFile) {
    TFile * output = new TFile("ZH_analysis.root", "RECREATE");
    TTree * tree_output = new TTree("tree_output", "Delphes");

    int histmin = 0;
    int histmax = 250;
    int bins = 50;
    int taubins = 50;
    int thrustbins = 50;
    int ntotalbins = 50;
    int ptbins = 50;

    TH1D * sigHmass = new TH1D("sigHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzjetsHmass = new TH1D("bkgzjetsHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzhbbHmass = new TH1D("bkgzhbbHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzhccHmass = new TH1D("bkgzhccHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);
    TH1D * bkgzh4qHmass = new TH1D("bkgzh4qHmass", "Reconstucted Higgs invariant mass", bins, histmin, histmax);

    TH1D * sigtau21 = new TH1D("sigtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzjetstau21 = new TH1D("bkgzjetstau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzhbbtau21 = new TH1D("bkgzhbbtau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzhcctau21 = new TH1D("bkgzhcctau21", "Reconstucted Higgs tau21", taubins, 0, 1);
    TH1D * bkgzh4qtau21 = new TH1D("bkgzh4qtau21", "Reconstucted Higgs tau21", taubins, 0, 1);

    TH1D * sigtau42 = new TH1D("sigtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzjetstau42 = new TH1D("bkgzjetstau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzhbbtau42 = new TH1D("bkgzhbbtau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzhcctau42 = new TH1D("bkgzhcctau42", "Reconstucted Higgs tau42", taubins, 0, 1);
    TH1D * bkgzh4qtau42 = new TH1D("bkgzh4qtau42", "Reconstucted Higgs tau42", taubins, 0, 1);

    TH1D * sigthrust = new TH1D("sigthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzjetsthrust = new TH1D("bkgzjetsthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzhbbthrust = new TH1D("bkgzhbbthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzhccthrust = new TH1D("bkgzhccthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);
    TH1D * bkgzh4qthrust = new TH1D("bkgzh4qthrust", "Reconstucted Higgs thrust", thrustbins, 0.5, 1);

    TH1D * sigthrust_minor = new TH1D("sigthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzjetsthrust_minor = new TH1D("bkgzjetsthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzhbbthrust_minor = new TH1D("bkgzhbbthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzhccthrust_minor = new TH1D("bkgzhccthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);
    TH1D * bkgzh4qthrust_minor = new TH1D("bkgzh4qthrust_minor", "Reconstucted Higgs thrust_minor", thrustbins, 0, 0.8);

    TH1D * signtotal = new TH1D("signtotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzjetsntotal = new TH1D("bkgzjetsntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzhbbntotal = new TH1D("bkgzhbbntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzhccntotal = new TH1D("bkgzhccntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);
    TH1D * bkgzh4qntotal = new TH1D("bkgzh4qntotal", "Reconstucted Higgs number of constituents", ntotalbins, 0, 140);

    TH1D * sigpt = new TH1D("signpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzjetspt = new TH1D("bkgzjetsnpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzhbbpt = new TH1D("bkgzhbbpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzhccpt = new TH1D("bkgzhccpt", "Reconstucted Higgs pt", ptbins, 200, 1000);
    TH1D * bkgzh4qpt = new TH1D("bkgzh4qpt", "Reconstucted Higgs pt", ptbins, 200, 1000);

    llgg(sigFile, sigHmass, sigtau21, sigtau42, sigthrust, sigthrust_minor, signtotal, sigpt);
    llgg(bkgzjetsFile, bkgzjetsHmass, bkgzjetstau21, bkgzjetstau42, bkgzjetsthrust, bkgzjetsthrust_minor, bkgzjetsntotal, bkgzjetspt);
    llgg(bkgzhbbFile, bkgzhbbHmass, bkgzhbbtau21, bkgzhbbtau42, bkgzhbbthrust, bkgzhbbthrust_minor, bkgzhbbntotal, bkgzhbbpt);
    llgg(bkgzhccFile, bkgzhccHmass, bkgzhcctau21, bkgzhcctau42, bkgzhccthrust, bkgzhccthrust_minor, bkgzhccntotal, bkgzhccpt);
    llgg(bkgzh4qFile, bkgzh4qHmass, bkgzh4qtau21, bkgzh4qtau42, bkgzh4qthrust, bkgzh4qthrust_minor, bkgzh4qntotal, bkgzh4qpt);

    THStack * Hmass = new THStack("Hmass", ";Reconstructed Higgs Invariant Mass [GeV]; Events/5 GeV");

    TH1 * sigHmass_clone = new TH1D( * sigHmass);

    Double_t sig_w = 0.0031132;
    Double_t sig_w_enh = 1000 * sig_w;
    Double_t bkgzjets_w = 92.79567345; //556.6953864;
    Double_t bkgzhbb_w = 0.06324864;
    Double_t bkgzhcc_w = 0.003139626;
    Double_t bkgzh4q_w = 0.01182654;

    Int_t sig = sigHmass -> GetEntries();
    Int_t bkgzjets = bkgzjetsHmass -> GetEntries();
    Int_t bkgzhbb = bkgzhbbHmass -> GetEntries();
    Int_t bkgzhcc = bkgzhccHmass -> GetEntries();
    Int_t bkgzh4q = bkgzh4qHmass -> GetEntries();
    Double_t significance = sig * sig_w / TMath::Sqrt(bkgzjets * bkgzjets_w + bkgzhbb * bkgzhbb_w + bkgzhcc * bkgzhcc_w + bkgzh4q * bkgzh4q_w);
    std::cout << "signal expected" << sig * sig_w << std::endl;
    std::cout << "zjets expected" << bkgzjets * bkgzjets_w << std::endl;
    std::cout << "zhbb expected" << bkgzhbb * bkgzhbb_w << std::endl;
    std::cout << "zhcc expected" << bkgzhcc * bkgzhcc_w << std::endl;
    std::cout << "zh4q expected" << bkgzh4q * bkgzh4q_w << std::endl;
    std::cout << "Expected significance is" << significance << std::endl;
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
    legend -> AddEntry(bkgzh4qHmass, "bkg Z+qqqq", "f");
    legend -> AddEntry(bkgzhccHmass, "bkg Z+cc", "f");
    Hmass -> Draw("HIST");
    Hmass -> GetHistogram() -> GetYaxis() -> SetTitleOffset(1.5);
    gPad -> Update();
    sigHmass_clone -> Draw("same");
    legend -> Draw("same");

    significance = 0.12;
    
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
    bkgzh4qtau21 -> SetLineColor(kTeal + 3);
    bkgzh4qtau21 -> SetLineWidth(4);
    bkgzh4qtau21 -> Scale(1.0/bkgzh4qtau21->Integral());
    bkgzhbbtau21 -> SetLineColor(kAzure + 2);
    bkgzhbbtau21 -> SetLineWidth(4);
    bkgzhbbtau21 -> Scale(1.0/bkgzhbbtau21->Integral());
    bkgzjetstau21 -> SetLineColor(kOrange + 7);
    bkgzjetstau21 -> SetLineWidth(4);
    bkgzjetstau21 -> Scale(1.0/bkgzjetstau21->Integral());
    sigtau21 -> SetLineColor(kPink - 1);
    sigtau21 -> SetLineWidth(4);
    sigtau21 -> Scale(1.0/sigtau21->Integral());

    TCanvas * tau21canvas = new TCanvas("tau21canvas", "Canvas", 1400, 1400, 1400, 1400);
    tau21canvas -> cd();
    tau21canvas -> SetWindowSize(1248, 1228);
    tau21canvas -> SetCanvasSize(1200, 1200);
    tau21canvas -> SetLogy(0);
    TLegend * legendtau21 = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendtau21 -> AddEntry(sigtau21, "signal", "l");
    legendtau21 -> AddEntry(bkgzjetstau21, "bkg Z+jets", "l");
    legendtau21 -> AddEntry(bkgzh4qtau21, "bkg Z+qqqq", "l");
    legendtau21 -> AddEntry(bkgzhbbtau21, "bkg Z+bb", "l");
    legendtau21 -> AddEntry(bkgzhcctau21, "bkg Z+cc", "l");
    bkgzh4qtau21 -> Draw("hist");
    bkgzh4qtau21 -> SetStats(0);
    bkgzh4qtau21->GetXaxis()->SetTitle("#tau_{21}");
    sigtau21 -> Draw("same hist");
    bkgzjetstau21 -> Draw("same hist");
    bkgzhbbtau21 -> Draw("same hist");
    bkgzhcctau21 -> Draw("same hist");
    bkgzh4qtau21 -> Draw("same hist");
    legendtau21 -> Draw("same hist");

    tau21canvas -> SaveAs("tau21_stack.png");

// tau42
    bkgzhcctau42 -> SetLineColor(kMagenta + 2);
    bkgzhcctau42 -> SetLineWidth(4);
    bkgzhcctau42 -> Scale(1.0/bkgzhcctau42->Integral());
    bkgzh4qtau42 -> SetLineColor(kTeal + 3);
    bkgzh4qtau42 -> SetLineWidth(4);
    bkgzh4qtau42 -> Scale(1.0/bkgzh4qtau42->Integral());
    bkgzhbbtau42 -> SetLineColor(kAzure + 2);
    bkgzhbbtau42 -> SetLineWidth(4);
    bkgzhbbtau42 -> Scale(1.0/bkgzhbbtau42->Integral());
    bkgzjetstau42 -> SetLineColor(kOrange + 7);
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
    legendtau42 -> AddEntry(bkgzjetstau42, "bkg Z+jets", "l");
    legendtau42 -> AddEntry(bkgzhbbtau42, "bkg Z+bb", "l");
    legendtau42 -> AddEntry(bkgzh4qtau42, "bkg Z+qqqq", "l");
    legendtau42 -> AddEntry(bkgzhcctau42, "bkg Z+cc", "l");
    sigtau42 -> Draw("hist");
    sigtau42 -> SetStats(0);
    sigtau42 ->GetXaxis()->SetTitle("#tau_{42}");
    bkgzjetstau42 -> Draw("same hist");
    bkgzhbbtau42 -> Draw("same hist");
    bkgzhcctau42 -> Draw("same hist");
    bkgzh4qtau42 -> Draw("same hist");
    legendtau42 -> Draw("same hist");

    tau42canvas -> SaveAs("tau42_stack.png");

// thrust
    bkgzhccthrust -> SetLineColor(kMagenta + 2);
    bkgzhccthrust -> SetLineWidth(4);
    bkgzhccthrust -> Scale(1.0/bkgzhccthrust->Integral());
    bkgzh4qthrust -> SetLineColor(kTeal + 3);
    bkgzh4qthrust -> SetLineWidth(4);
    bkgzh4qthrust -> Scale(1.0/bkgzh4qthrust->Integral());
    bkgzhbbthrust -> SetLineColor(kAzure + 2);
    bkgzhbbthrust -> SetLineWidth(4);
    bkgzhbbthrust -> Scale(1.0/bkgzhbbthrust->Integral());
    bkgzjetsthrust -> SetLineColor(kOrange + 7);
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
    legendthrust -> AddEntry(bkgzjetsthrust, "bkg Z+jets", "l");
    legendthrust -> AddEntry(bkgzhbbthrust, "bkg Z+bb", "l");
    legendthrust -> AddEntry(bkgzh4qthrust, "bkg Z+qqqq", "l");
    legendthrust -> AddEntry(bkgzhccthrust, "bkg Z+cc", "l");
    bkgzhccthrust -> Draw("hist");
    bkgzhccthrust -> SetStats(0);
    bkgzhccthrust ->GetXaxis()->SetTitle("Thrust");
    sigthrust -> Draw("same hist");
    bkgzjetsthrust -> Draw("same hist");
    bkgzhbbthrust -> Draw("same hist");
    bkgzh4qthrust -> Draw("same hist");
    legendthrust -> Draw("same hist");

    thrustcanvas -> SaveAs("thrust_stack.png");
// thrust_minor
    bkgzhccthrust_minor -> SetLineColor(kMagenta + 2);
    bkgzhccthrust_minor -> SetLineWidth(4);
    bkgzhccthrust_minor -> Scale(1.0/bkgzhccthrust_minor->Integral());
    bkgzh4qthrust_minor -> SetLineColor(kTeal + 3);
    bkgzh4qthrust_minor -> SetLineWidth(4);
    bkgzh4qthrust_minor -> Scale(1.0/bkgzh4qthrust_minor->Integral());
    bkgzhbbthrust_minor -> SetLineColor(kAzure + 2);
    bkgzhbbthrust_minor -> SetLineWidth(4);
    bkgzhbbthrust_minor -> Scale(1.0/bkgzhbbthrust_minor->Integral());
    bkgzjetsthrust_minor -> SetLineColor(kOrange + 7);
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
    legendthrust_minor -> AddEntry(bkgzjetsthrust_minor, "bkg Z+jets", "l");
    legendthrust_minor -> AddEntry(bkgzhbbthrust_minor, "bkg Z+bb", "l");
    legendthrust_minor -> AddEntry(bkgzh4qthrust_minor, "bkg Z+qqqq", "l");
    legendthrust_minor -> AddEntry(bkgzhccthrust_minor, "bkg Z+cc", "l");
    bkgzh4qthrust_minor -> Draw("hist");
    bkgzh4qthrust_minor -> SetStats(0);
    bkgzh4qthrust_minor ->GetXaxis()->SetTitle("thrust_minor");
    bkgzhccthrust_minor -> Draw("same hist");
    sigthrust_minor -> Draw("same hist");
    bkgzjetsthrust_minor -> Draw("same hist");
    bkgzhbbthrust_minor -> Draw("same hist");
    legendthrust_minor -> Draw("same hist");

    thrust_minorcanvas -> SaveAs("thrust_minor_stack.png");
// ntotal
    bkgzhccntotal -> SetLineColor(kMagenta + 2);
    bkgzhccntotal -> SetLineWidth(4);
    bkgzhccntotal -> Scale(1.0/bkgzhccntotal->Integral());
    bkgzh4qntotal -> SetLineColor(kTeal + 3);
    bkgzh4qntotal -> SetLineWidth(4);
    bkgzh4qntotal -> Scale(1.0/bkgzh4qntotal->Integral());
    bkgzhbbntotal -> SetLineColor(kAzure + 2);
    bkgzhbbntotal -> SetLineWidth(4);
    bkgzhbbntotal -> Scale(1.0/bkgzhbbntotal->Integral());
    bkgzjetsntotal -> SetLineColor(kOrange + 7);
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
    TLegend * legendntotal = new TLegend(0.65, 0.65, 0.85, 0.85);
    legendntotal -> AddEntry(signtotal, "signal", "l");
    legendntotal -> AddEntry(bkgzjetsntotal, "bkg Z+jets", "l");
    legendntotal -> AddEntry(bkgzhbbntotal, "bkg Z+bb", "l");
    legendntotal -> AddEntry(bkgzh4qntotal, "bkg Z+qqqq", "l");
    legendntotal -> AddEntry(bkgzhccntotal, "bkg Z+cc", "l");
    bkgzh4qntotal -> Draw("hist");
    bkgzh4qntotal -> SetStats(0);
    bkgzh4qntotal ->GetXaxis()->SetTitle("ntotal");
    bkgzhccntotal -> Draw("same hist");
    signtotal -> Draw("same hist");
    bkgzjetsntotal -> Draw("same hist");
    bkgzhbbntotal -> Draw("same hist");
    legendntotal -> Draw("same hist");

    ntotalcanvas -> SaveAs("ntotal_stack.png");

// pt
    bkgzhccpt -> SetLineColor(kMagenta + 2);
    bkgzhccpt -> SetLineWidth(4);
    bkgzhccpt -> Scale(1.0/bkgzhccpt->Integral());
    bkgzh4qpt -> SetLineColor(kTeal + 3);
    bkgzh4qpt -> SetLineWidth(4);
    bkgzh4qpt -> Scale(1.0/bkgzh4qpt->Integral());
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
    legendpt -> Draw("same hist");

    ptcanvas -> SaveAs("pt_stack.png");

//end    
    output -> cd();
    sigHmass -> Write();
    bkgzjetsHmass -> Write();
    bkgzhbbHmass -> Write();
    bkgzhccHmass -> Write();
    bkgzh4qHmass -> Write();
    Hmass -> Write();
    tree_output -> Write();
    output -> Close();
}
