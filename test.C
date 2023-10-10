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
	
void test(const char *inputFile) {

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
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchEvent = treeReader->UseBranch("Event");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    std::cout << "Start loop ever tree" << std::endl;
    for (Long64_t entry = 0; entry < nEntries; entry++) {
	treeReader->ReadEntry(entry);
	Int_t nJet = branchJet->GetEntries();
	TObject *object;
	Jet *jet;
	for (Int_t jentry = 0; jentry < nJet; jentry++) {
	    jet = (Jet *) branchJet -> At(jentry);
	    std::cout << jet->Constituents.GetEntriesFast() << "constituents in jet" << std::endl;
	    for(int i = 0; i < jet->Constituents.GetEntriesFast(); ++i) {
	    	object = jet->Constituents.At(i);
	    	if(object == 0) {
		    std::cout << "Null" << std::endl;
		    continue;
	    	}
	    	if(object->IsA() == Tower::Class()) {
		    Tower *tower;
		    tower = (Tower*) object;
		    std::cout << "is tower with Eem = " << tower->Eem << "and Ehad = " <<tower->Ehad << std::endl;
		}
	    	if(object->IsA() == GenParticle::Class()) {
		    std::cout << "is genpar" << std::endl;
	    	} 
	    	if(object->IsA() == Track::Class()) {
		    Track *track;
		    track = (Track*) object;
		    std::cout << "is track" << std::endl;
	    	}	
	    	if(object->IsA() == ParticleFlowCandidate::Class()) {
		    std::cout << "is candidate" << std::endl;
	    	}
	    }    
	}
    }
}
