

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom.h"
//#include "vector"
#include <vector>
#endif

//------------------------------------------------------------------------------

void gammaGammaHHggbbEventDisplay(const char *inputFile, int NthEvent, bool Back)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

int njets=0;
int ntracks=0;
int nentries=0;
int nparticles=0, nElectrons=0, nMuons=0;
int contEventsHH=0, contB=0, contG=0;
int dispar[10];
int guia=0;


  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis  (DSiDi)
  TClonesArray *branchJet = treeReader->UseBranch("PFJet10");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet10");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowPhoton");
  
  /*//Pointers to branches (example3)
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");*/
  
  
  TObject *object;
  GenParticle *particleConst;
  Track *trackConst;
  Tower *towerConst;
  
  
  // Book histograms
  TH1 *histJetPhi = new TH1F("jet_phi", "jet Phi", 100, -4.0, 4.0);
  TH1 *histJetEta = new TH1F("jet_eta", "jet Eta", 100, -4.0, 4.0);
  TH2F *histJetPhiVsEta = new TH2F("jet_phi_vs_eta", "jet phi vs. eta", 100, -4.0, 4.0, 100, -4.0, 4.0);
  TH2F *histJetPhiVsEtaTracks = new TH2F("jet_phi_vs_eta_tracks", "jet phi vs. eta", 100, -4.0, 4.0, 100, -4.0, 4.0);
  
  TCanvas *c1 = new TCanvas();
  TCanvas *c2 = new TCanvas();
  //TCanvas *c3 = new TCanvas();
  
  histJetPhiVsEta->SetStats(0);
  histJetPhiVsEta->GetXaxis()->SetTitle("eta");
  histJetPhiVsEta->GetYaxis()->SetTitle("phi");
  //histJetPhiVsEtaTracks->SetMarkerColor(kPink);
  gStyle->SetPalette(1);
  histJetPhiVsEta->SetFillColor(0);
  gStyle->SetOptStat(0);
  
  //Filling 2d histogram with jets
  c1->cd();
  histJetPhiVsEta->Draw("col");
  //histJetPhiVsEtaTracks->Draw("same");
  //int NthEvent=14069;  //set which event you want to plot
  treeReader->ReadEntry(NthEvent);
  
  
  //////Drawing ellipse`s jets
  if(branchJet->GetEntries() > 0) njets =  branchJet->GetEntries();
  else njets=0;
  
  int indexb1, indexb2, indexnb1, indexnb2, contBJets=0, contNBJets=0;
  
  cout<<"For event "<<NthEvent<<", number of jets: "<<njets<<endl<<endl;
  for(int i=0;i<njets;i++)
  {
	      Jet *jet = (Jet*) branchJet->At(i); // Take ith jet
	 
	      double valPhi = jet->Phi;
	      double valEta = jet->Eta;
	      //histJetPhiVsEta->Fill(valEta,valPhi);
	      cout<<"Loop for jet: "<<i+1<<endl;
	      cout<<"Phi value: "<<valPhi<<endl;
	      cout<<"Eta value: "<<valEta<<endl;
	      cout<<"Jet flavor: "<<jet->Flavor<<endl;
	      cout<<"b Tag: "<<jet->BTag<<endl<<endl;
	      TEllipse* r = new TEllipse(valEta,valPhi,0.5,0.5); //.5 because parameter R =.5 for the antikt used in SiDi
              if(jet->BTag==1) r->SetLineStyle(2);
              else if(jet->BTag==0) r->SetLineStyle(1);
              r->SetFillStyle(4000);  //makes it transparent
              r->SetLineColor(2);
   	      r->SetLineWidth(1);
	      r->Draw();
	      
	      
	      if(jet->BTag==1)
		{
			contBJets++;
			if(contBJets==1)
			{
				indexb1=i;
			}	
			else if(contBJets==2)
			{
				indexb2=i;
			}
		}
		else
		{
			contNBJets++;
			if(contNBJets==1)
			{
				indexnb1=i;
			}
			else if(contNBJets==2)
			{
				indexnb2=i;
			}
		}
  }
  /////
  
  //////Filling with Higgs Bosons
  
  if(branchParticle->GetEntries() > 0) nparticles =  branchParticle->GetEntries();
  else nparticles=0;
  
  cout<<"For event "<<NthEvent<<", number of particles: "<<nparticles<<endl<<endl;
  
  for(int i=0;i<nparticles;i++)
  {
	GenParticle *particle = (GenParticle*) branchParticle->At(i); 
  	if(particle->PID==25 || particle->PID==36 || Back==true) ////36 is for A0 (BSM Higgs)
  	{
	      double valPhi = particle->Phi;
	      double valEta = particle->Eta;
	      //histJetPhiVsEta->Fill(valEta,valPhi);
	      if(Back==false)
	      {
		      cout<<endl<<"Loop for Higgs: "<<i+1<<endl;
		      cout<<"Phi value: "<<valPhi<<endl;
		      cout<<"Eta value: "<<valEta<<endl;
	      }
	      
	      ///////Ellipse
	      /*
	      TEllipse* r = new TEllipse(valEta,valPhi,0.1,0.1); 
              r->SetLineStyle(1);
              r->SetFillStyle(1001); 
              r->SetLineColor(4);
   	      r->SetLineWidth(1);
   	      r->SetFillColor(3);
	      r->Draw();
	      */
	      ////////
	      
	      //////Line
	      /*
	      TLine* lineX = new TLine(valEta-.2,valPhi,valEta+.2,valPhi);
	      lineX->SetLineColor(3);  //verde
	      lineX->SetLineWidth(2);
	      lineX->Draw("same");
	      TLine* lineY = new TLine(valEta,valPhi-.2,valEta,valPhi+.2);
	      lineY->SetLineColor(3);  //verde
	      lineY->SetLineWidth(2);
	      lineY->Draw("same");
	      */
	      //////
	      
	      int d1Index=particle->D1;
  	      int d2Index=particle->D2;
  	      int d1=0;
  	      int d2=0;
  	      if (d1Index >= 0 && d2Index >= 0)
  	      {
	  	d1=abs(static_cast<GenParticle*>(branchParticle->At(d1Index))->PID);
	  	if(d1==5) contB++;
	  	else if(d1==21) contG++;
	  	d2=abs(static_cast<GenParticle*>(branchParticle->At(d2Index))->PID);
	  	if(d2==5) contB++;
	  	else if(d2==21) contG++;
  	      }
  	      cout<<"D1 PID: "<<d1<<endl;
	      cout<<"D2 PID: "<<d2<<endl;
  	      if(d1==5 || d2==5 || d1==21 || d2==21)
  	      {
		      
		      
		      ///////Ellipse
		      /*
		      TEllipse* rd1 = new TEllipse(valEta,valPhi,0.1,0.1); 
		      rd1->SetLineStyle(1);
		      rd1->SetFillStyle(1001);  
		      rd1->SetLineColor(4);
	   	      rd1->SetLineWidth(1);
	   	      if(d1==21) rd1->SetFillColor(4);
	   	      else if(d1==5) rd1->SetFillColor(2);
		      rd1->Draw();
		      */
		      /////////
		      
		      //////Line
		      if(d1==5)
		      {
			      valPhi=static_cast<GenParticle*>(branchParticle->At(d1Index))->Phi;
			      valEta=static_cast<GenParticle*>(branchParticle->At(d1Index))->Eta;
			      //histJetPhiVsEta->Fill(valEta,valPhi);
			      //cout<<"Loop for D1: "<<i+1<<endl;
			      //cout<<"D1 PID: "<<d1<<endl;
			      cout<<"Entraaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"<<endl;
			      cout<<"Phi value: "<<valPhi<<endl;
			      cout<<"Eta value: "<<valEta<<endl;
			      TLine* lineX1 = new TLine(valEta-.2,valPhi,valEta+.2,valPhi);
			      if(d1==21) lineX1->SetLineColor(kGreen);  //g: green
		   	      else if(d1==5) lineX1->SetLineColor(kBlue);  //b: blue
			      lineX1->SetLineWidth(2);
			      lineX1->Draw("same");
			      TLine* lineY1 = new TLine(valEta,valPhi-.2,valEta,valPhi+.2);
			      if(d1==21) lineY1->SetLineColor(kGreen); //g: green
		   	      else if(d1==5) lineY1->SetLineColor(kBlue);  //b: blue
			      lineY1->SetLineWidth(2);
			      lineY1->Draw("same");
		      
		      ////////
		      
			     
		      //////////
		      }
		      
		      //////Line
		      if(d2==5)
		      {
			      valPhi=static_cast<GenParticle*>(branchParticle->At(d2Index))->Phi;
			      valEta=static_cast<GenParticle*>(branchParticle->At(d2Index))->Eta;
			      //histJetPhiVsEta->Fill(valEta,valPhi);
			      //cout<<"Loop for D2: "<<i+1<<endl;
			      //cout<<"D2 PID: "<<d2<<endl;
			      cout<<"Entraaaaaaa 22222222222222222"<<endl;
			      cout<<"Phi value: "<<valPhi<<endl;
			      cout<<"Eta value: "<<valEta<<endl;
			      cout<<"B counter: "<<contB<<endl<<"G counnter: "<<contG<<endl<<endl;		      
			      ////////Ellipse
			      /*
			      TEllipse* rd2 = new TEllipse(valEta,valPhi,0.1,0.1); 
			      rd2->SetLineStyle(1);
			      rd2->SetFillStyle(1001);  
			      rd2->SetLineColor(4);
		   	      rd2->SetLineWidth(1);
		   	      if(d2==21) rd2->SetFillColor(4);
		   	      else if(d2==5) rd2->SetFillColor(2);
			      rd2->Draw();
		      */
			      TLine* lineX2 = new TLine(valEta-.2,valPhi,valEta+.2,valPhi);
			      if(d2==21) lineX2->SetLineColor(kGreen);  //g: green
		   	      else if(d2==5) lineX2->SetLineColor(kBlue);  //b: blue
			      lineX2->SetLineWidth(2);
			      lineX2->Draw("same");
			      TLine* lineY2 = new TLine(valEta,valPhi-.2,valEta,valPhi+.2);
			      if(d2==21) lineY2->SetLineColor(kGreen);  //g: green
		   	      else if(d1==2) lineY2->SetLineColor(kBlue);  //b: blue
			      lineY2->SetLineWidth(2);
			      lineY2->Draw("same");
		      }
		      ////////
		      
		      if(d1!=d2)
		      {
		      	dispar[guia]=i+1;
		      	guia++;
		      }
	      }
	 }
	 else cout<<i+1<<", not Higgs"<<endl;
  }
  cout<<endl;
  /////
  
  
  ///adding electron
  if(branchElectron->GetEntries() > 0) nElectrons =  branchElectron->GetEntries();
  else nElectrons=0;
  Electron *electron = (Electron*) branchElectron->At(0);
  if(nElectrons>0)
  {
	  double elecEta=electron->Eta;
	  double elecPhi=electron->Phi;
	  TLine* lineElecX = new TLine(elecEta-.2,elecPhi,elecEta+.2,elecPhi);
				      lineElecX->SetLineColor(kViolet);  
				      lineElecX->SetLineWidth(2);
				      lineElecX->Draw("same");
				      TLine* lineElecY = new TLine(elecEta,elecPhi-.2,elecEta,elecPhi+.2);
				      lineElecY->SetLineColor(kViolet); 
				      lineElecY->SetLineWidth(2);
				      lineElecY->Draw("same");
				      
	  cout<<endl<<"Electron eta: "<<elecEta<<endl<<"Electron phi: "<<elecPhi<<endl<<endl;
  }				      
  /////
  
  ///////Filling with tracks
  if(branchEFlowTrack->GetEntries() > 0) ntracks =  branchEFlowTrack->GetEntries();
  else ntracks=0;
  
  cout<<"For event "<<NthEvent<<", number of tracks: "<<ntracks<<endl<<endl;
  
  for(int i=0;i<ntracks;i++)
  {
	      Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
	 
	      double valPhi = track->Phi;
	      double valEta = track->Eta;
	      TEllipse* eTrack = new TEllipse(valEta,valPhi,0.01,0.01); 
              eTrack->SetLineStyle(1);
              eTrack->SetFillColor(kPink);
              eTrack->SetLineColor(kPink);
   	      eTrack->SetLineWidth(1);
	      eTrack->Draw();
	      //cout<<"Loop for track: "<<i+1<<endl;
	      //cout<<"Phi value: "<<valPhi<<endl;
	      //cout<<"Eta value: "<<valEta<<endl<<endl;
  }
  //////
  
  /*
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    
    //HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    //LHEFEvent *event = (LHEFEvent*) branchEvent -> At(0);
    //Float_t weight = event->Weight;

    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
    
      njets =  branchJet->GetEntries();
      }
     else
     {
      njets=0;
      }
      
      //cout<<"Jets per event: "<<njets<<endl;
      for(int i=0;i<njets;i++)
      {
	      // Take ith jet
	      Jet *jet = (Jet*) branchJet->At(i);
	      // Plot jet transverse momentum
	      histJetPhi->Fill(jet->Phi);
	      histJetEta->Fill(jet->Eta);
	      // Print jet transverse momentum
	      //cout << "Jet pt: "<<jet->PT << endl;
    }
    nentries+=njets;
  }
  
  */
  
  
  
  
  
  
  
  cout<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<"Starts Px and Py"<<endl<<endl;
  
  
   treeReader->ReadEntry(NthEvent);
  
  
  
  double ptMax=0;
  /*if(branchEFlowTrack->GetEntries() > 0)
  {
    
      ntracks =  branchEFlowTrack->GetEntries();
  }
  else
  {
      ntracks=0;
  }
  //cout<<"Number of tracks: "<<ntracks<<endl;
  for(int i=0;i<ntracks;i++)
  {
	      Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
	      if(ptMax<(track->PT)) ptMax=track->PT;
  }
  //cout<<"For tracks, max pt: "<<ptMax<<endl<<endl;
  */
  
  
  TLorentzVector jetMomentum;
  int njetsg;
  
 
  if(branchJet->GetEntries() > 0)
  {
    
      njets =  branchJet->GetEntries();
  }
  else
  {
      njets=0;
  }
  if(branchGenJet->GetEntries() > 0)
  {
    
      njetsg =  branchGenJet->GetEntries();
  }
  else
  {
      njetsg=0;
  }
  
  
  for(int i=0;i<njetsg;i++)
  {
	      Jet *genjet = (Jet*) branchGenJet->At(i); // Take ith jet
	
	      
	      // Loop over all jet's constituents
	      for(int j = 0; j < genjet->Constituents.GetEntriesFast(); ++j)
	      {
		object = genjet->Constituents.At(j);
		
		
		if(object == 0) continue;

		if(object->IsA() == GenParticle::Class())
		{
		  particleConst = (GenParticle*) object;
		  double ptSen=particleConst->PT;
		  if(ptSen>ptMax) ptMax=ptSen;
		}
		
	
	      }
  }
  
  
  
  if(ntracks>0)
  {
  // Book histograms
  TH1 *histJetPx = new TH1F("jet_px", "jet px", 100, -100.0, 100.0);
  TH1 *histJetPy = new TH1F("jet_py", "jet py", 100, 0.0, 300.0);
  TH2F *histJetPxAndPy = new TH2F("jet_px_vs_py", "jet py vs. px", 100, -ptMax*1.1, ptMax*1.1, 100, -ptMax*1.1, ptMax*1.1);
  
  c2->cd();
  histJetPxAndPy->Draw("col");
  
  histJetPxAndPy->SetStats(0);
  histJetPxAndPy->GetXaxis()->SetTitle("y");
  histJetPxAndPy->GetYaxis()->SetTitle("x");
  gStyle->SetPalette(1);
  histJetPxAndPy->SetFillColor(0);
  gStyle->SetOptStat(0);
  
  
  TEllipse* detector = new TEllipse(0.0,0.0,ptMax,ptMax);
              detector->SetLineStyle(1);
              detector->SetFillStyle(4000);  //makes it transparent
              detector->SetLineColor(kBlack);
   	      detector->SetLineWidth(1);
	      detector->Draw();
	
  TLine* xAxis = new TLine(-ptMax,0.0,ptMax,0.0);
	      xAxis->SetLineColor(1);
	      xAxis->SetLineWidth(1);
	      xAxis->Draw("same");
	      
  TLine* yAxis = new TLine(0.0,-ptMax,0.0,ptMax);
	      yAxis->SetLineColor(1);
	      yAxis->SetLineWidth(1);
	      yAxis->Draw("same"); 	
  	       	
  
  //Filling 2d histogram with jets
  //histJetPhiVsEta->Draw("col");
  
  
  
  cout<<"For event "<<NthEvent<<", number of jets: "<<njets<<endl<<endl;
  cout<<"For event "<<NthEvent<<", number of gen jets: "<<njetsg<<endl<<endl;
  
  int indexgenb1=0, indexgenb2=0, indexgennb1=0, indexgennb2=0, contB1DeltaR04=0, contNB1DeltaR04=0, contB2DeltaR04=0, contNB2DeltaR04=0;
  
  for(int i=0;i<njetsg;i++)
  {
	      Jet *genjet = (Jet*) branchGenJet->At(i); // Take ith jet
	      jetMomentum = genjet->P4();
	      
	      Jet *jetb1 = (Jet*) branchJet->At(indexb1);
		Jet *jetb2 = (Jet*) branchJet->At(indexb2);
		Jet *jetnb1 = (Jet*) branchJet->At(indexnb1);
		Jet *jetnb2 = (Jet*) branchJet->At(indexnb2);
	      //double length = 1;
	      //histJetPxAndPy->Fill(sin(jet->Phi)*ptMax,cos(jet->Phi)*ptMax);
	      //TLine* line = new TLine(0.0,0.0,length*jetMomentum.Py(),length*jetMomentum.Px());
	      //if (jetMomentum.Pt()>93) line->SetLineColor(2);
	      //else line->SetLineColor(15);
	      //line->SetLineWidth(2);
	      //line->Draw("same");
	      double x=sin(genjet->Phi)*ptMax*1.05;
	      double y=cos(genjet->Phi)*ptMax*1.05;
	      TEllipse* r = new TEllipse(x,y,ptMax*0.09,ptMax*0.09);
              r->SetLineStyle(2);
              //r->SetFillStyle(0);
              r->SetFillColor(2);
              //r->SetFillStyle(4000);  //makes it transparent
              r->SetLineColor(kBlue);
   	      r->SetLineWidth(1);
	      r->Draw();
	      
	      cout<<"GenJet const amount: "<<genjet->Constituents.GetEntriesFast()<<endl<<endl;
	      double constX, constY;
	      
	      // Loop over all jet's constituents
	      for(int j = 0; j < genjet->Constituents.GetEntriesFast(); ++j)
	      {
		object = genjet->Constituents.At(j);
		
		cout<<"gen constituent: "<<j<<endl;

		// Check if the constituent is accessible
		if(object == 0) cout<<"Not accesible"<<endl;
        	if(object != 0) cout<<"Accesible"<<endl;
		
		if(object == 0) continue;

		if(object->IsA() == GenParticle::Class())
		{
		  particleConst = (GenParticle*) object;
		  cout << "    GenPart pt: " << particleConst->PT << ", eta: " << particleConst->Eta << ", phi: " << particleConst->Phi<<", PID: "<<particleConst->PID<<endl;
		  constX = (particleConst->PT) * sin(particleConst->Phi);
		  constY = (particleConst->PT) * cos(particleConst->Phi);
		  int d1Index=particleConst->D1;
		  int d2Index=particleConst->D2;
		  int m1Index=particleConst->M1;
		  int m2Index=particleConst->M2;
		  if (d1Index>=0)
		  {
			int d1=abs(static_cast<GenParticle*>(branchParticle->At(d1Index))->PID);
			cout<<"d1 PID: "<<d1;
  		  }
  		  if (d2Index>=0)
		  {
			int d2=abs(static_cast<GenParticle*>(branchParticle->At(d2Index))->PID);
			cout<<"d2 PID: "<<d2;
  		  }
  		  if (m1Index>=0)
		  {
			int m1=abs(static_cast<GenParticle*>(branchParticle->At(m1Index))->PID);
			cout<<"m1 PID: "<<m1;
  		  }
  		  if (m2Index>=0)
		  {
			int m2=abs(static_cast<GenParticle*>(branchParticle->At(m2Index))->PID);
			cout<<"m2 PID: "<<m2;
  		  }
  		  cout<<endl<<endl<<endl;
		}
		else cout<<"Not GenParticle"<<endl;
		
		if(object->IsA() == Track::Class())
		{
		  trackConst = (Track*) object;
		  cout << "    Track pt: " << trackConst->PT << ", eta: " << trackConst->Eta << ", phi: " << trackConst->Phi << endl<<endl<<endl<<endl;
		  constX = (trackConst->PT) * sin(trackConst->Phi);
		  constY = (trackConst->PT) * cos(trackConst->Phi);
		  //momentum += track->P4();
		}
		else cout<<"Not Track"<<endl;
		
		if(object->IsA() == Tower::Class())
		{
		  towerConst = (Tower*) object;
		  cout << "    Tower pt: " << towerConst->ET << ", eta: " << towerConst->Eta << ", phi: " << towerConst->Phi << endl<<endl<<endl<<endl;
		  constX = (towerConst->ET) * sin(towerConst->Phi);
		  constY = (towerConst->ET) * cos(towerConst->Phi);
		  //momentum += tower->P4();
		}
		else cout<<"Not Tower"<<endl;
		
						double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2, deltaRMin;
						
						double genjetEta = genjet->Eta;
				     		double genjetPhi = genjet->Phi;
				     		
				     		double jetB1Eta = jetb1->Eta;
			     			double jetB1Phi = jetb1->Phi;
			     			deltaRB1 = sqrt(pow(genjetEta-jetB1Eta, 2) + pow(genjetPhi-jetB1Phi, 2));
			     			
			     			double jetB2Eta = jetb2->Eta;
			     			double jetB2Phi = jetb2->Phi;
			     			deltaRB2 = sqrt(pow(genjetEta-jetB2Eta, 2) + pow(genjetPhi-jetB2Phi, 2));
			     			
			     			double jetNB1Eta = jetnb1->Eta;
			     			double jetNB1Phi = jetnb1->Phi;
			     			deltaRNB1 = sqrt(pow(genjetEta-jetNB1Eta, 2) + pow(genjetPhi-jetNB1Phi, 2));
			     			
			     			double jetNB2Eta = jetnb2->Eta;
			     			double jetNB2Phi = jetnb2->Phi;
			     			deltaRNB2 = sqrt(pow(genjetEta-jetNB2Eta, 2) + pow(genjetPhi-jetNB2Phi, 2));
			     			
			     			if(deltaRB1<deltaRB2 && deltaRB1<deltaRNB1 && deltaRB1<deltaRNB2)
			     			{
			     				TLine* line = new TLine(0.0,0.0,constX,constY);
							line->SetLineColor(kBlue);
							line->SetLineWidth(2);
					      	        line->Draw("same");
					      	        //if(indexgenb1!=i && indexgenb1!=0)  contindb1++;
					      	        indexgenb1=i;
			     			}
			     			else if(deltaRB2<deltaRB1 && deltaRB2<deltaRNB1 && deltaRB2<deltaRNB2)
			     			{
			     				TLine* line = new TLine(0.0,0.0,constX,constY);
							line->SetLineColor(kCyan);
							line->SetLineWidth(2);
					      	        line->Draw("same");
					      	        //if(indexgenb2!=i && indexgenb2!=0)  contindb2++; 
					      	        indexgenb2=i;
			     			}
			     			else if(deltaRNB1<deltaRB2 && deltaRNB1<deltaRB1 && deltaRNB1<deltaRNB2)
			     			{
							TLine* line = new TLine(0.0,0.0,constX,constY);
							line->SetLineColor(kGreen+2);
							line->SetLineWidth(2);
					      	        line->Draw("same");
					      	        //if(indexgennb1!=i && indexgennb1!=0)  contindnb1++;
					      	        indexgennb1=i; 
			     			}
			     			else if(deltaRNB2<deltaRB2 && deltaRNB2<deltaRB1 && deltaRNB2<deltaRNB1)
			     			{
			     				TLine* line = new TLine(0.0,0.0,constX,constY);
							line->SetLineColor(kGreen);
							line->SetLineWidth(2);
					      	        line->Draw("same");
					      	        //if(indexgennb2!=i && indexgennb2!=0)  contindnb2++; 
					      	        indexgennb2=i;
			     			}
		
	
	      }
  }
  /////
  
  int contB=0, contNB=0;
  
  for(int i=0;i<njets;i++)
  {
	      Jet *jet = (Jet*) branchJet->At(i); // Take ith jet
	      jetMomentum = jet->P4();
	      histJetPx->Fill(jetMomentum.Px());
	      histJetPy->Fill(jetMomentum.Py()); 
	      //double length = 1;
	      //histJetPxAndPy->Fill(sin(jet->Phi)*ptMax,cos(jet->Phi)*ptMax);
	      //TLine* line = new TLine(0.0,0.0,length*jetMomentum.Py(),length*jetMomentum.Px());
	      //if (jetMomentum.Pt()>93) line->SetLineColor(2);
	      //else line->SetLineColor(15);
	      //line->SetLineWidth(2);
	      //line->Draw("same");
	      double x=sin(jet->Phi)*ptMax*1.05;
	      double y=cos(jet->Phi)*ptMax*1.05;
	      /*TEllipse* r = new TEllipse(x,y,ptMax*0.09,ptMax*0.09);
              r->SetLineStyle(1);
              r->SetFillColor(2);
              //r->SetFillStyle(4000);  //makes it transparent
              r->SetLineColor(2);
   	      r->SetLineWidth(1);
	      r->Draw();*/
	      int jetFlavor=jet->Flavor;
	      int jetBTag=jet->BTag;
	      TString jetName = Form("");
	      TString jetB1Name = Form("jetB1");
	      TString jetB2Name = Form("jetB2");
	      TString jetNB1Name = Form("jetNB1");
	      TString jetNB2Name = Form("jetNB2");
	      if(jetBTag==1)
	      {
	      	contB++;
	      	if(contB==1) jetName=jetB1Name;
	      	else if(contB==2) jetName=jetB2Name;
	      }
	      else
	      {
	      	contNB++;
	      	if(contNB==1) jetName=jetNB1Name;
	      	else if(contNB==2) jetName=jetNB2Name;
	      }
	      TString jetFlavorText = Form("Flavor: %d", jetFlavor);
	      TString jetBTagText = Form(";    b-tag: %d", jetBTag);
	      TString jetPTText = Form(";   PT: %.3f", jet->PT);
	      TString jetMassText = Form(";    Mass: %.3f", jet->Mass);
	      jetFlavorText=jetFlavorText+jetBTagText;
	      jetPTText+=jetMassText;
	      TPaveText *pave = new TPaveText(x-1, y, x+1, y+2, "");	
              pave->AddText(jetFlavorText);
              pave->SetBorderSize(0);  // optionally set border size
              pave->SetFillColor(0);   // set fill color (0 is transparent)
              pave->SetTextFont(1);   // optional: set the font
              pave->SetTextSize(0.03); // optional: set the text size
              pave->Draw();
              TPaveText *pave2 = new TPaveText(x-1, y-2, x+1, y, "");	
              pave2->AddText(jetPTText);
              pave2->SetBorderSize(0);  // optionally set border size
              pave2->SetFillColor(0);   // set fill color (0 is transparent)
              pave2->SetTextFont(1);   // optional: set the font
              pave2->SetTextSize(0.03); // optional: set the text size
              pave2->Draw();
              TPaveText *pave3 = new TPaveText(x-1, y+2, x+1, y+4, "");	
              pave3->AddText(jetName);
              pave3->SetBorderSize(0);  // optionally set border size
              pave3->SetFillColor(0);   // set fill color (0 is transparent)
              pave3->SetTextFont(1);   // optional: set the font
              pave3->SetTextSize(0.03); // optional: set the text size
              pave3->Draw();
	      //cout<<"For jet "<<i+1<<", length: "<<length<<endl;
	      cout<<"For jet "<<i+1<<", mass: "<<jet->Mass<<endl;
	      cout<<"For jet "<<i+1<<", pt: "<<jetMomentum.Pt()<<endl;
	      cout<<"For jet "<<i+1<<", px: "<<jetMomentum.Px()<<endl;
	      cout<<"For jet "<<i+1<<", py: "<<jetMomentum.Py()<<endl;
	      cout<<"for jet "<<i+1<<", phi: "<<jet->Phi<<endl<<endl;
	      cout<<"for jet "<<i+1<<", flavor: "<<jet->Flavor<<endl<<endl;
	      cout<<"for jet "<<i+1<<", b-tag: "<<jet->BTag<<endl<<endl<<endl<<endl;
	      
  }
  
  
  cout<<endl<<endl<<endl<<"Now particles: "<<endl<<endl;
  
  int nParticles=0;
  if(branchParticle->GetEntries() > 0) nParticles =  branchParticle->GetEntries();
    		else nParticles=0;
    		for(int j=0;j<nParticles;j++)
    		{
  			Jet *genjetb1 = (Jet*) branchGenJet->At(indexgenb1);
			Jet *genjetb2 = (Jet*) branchGenJet->At(indexgenb2);
			Jet *genjetnb1 = (Jet*) branchGenJet->At(indexgennb1);
			Jet *genjetnb2 = (Jet*) branchGenJet->At(indexgennb2);
  			GenParticle *particle = (GenParticle*) branchParticle->At(j); 
  			int pid=abs(particle->PID);
  			if(pid==4 || pid==5) cout<<j<<" particle is b"<<endl;
		  	int d1Index=particle->D1;
		  int d2Index=particle->D2;
		  int m1Index=particle->M1;
		  int m2Index=particle->M2;
		  if (d1Index>=0)
		  {
			int d1=abs(static_cast<GenParticle*>(branchParticle->At(d1Index))->PID);
			if((d1==5 || d1==4) && pid==21)
			{ 
				cout<<j<<"DAAUG particle`s d1 PID: "<<d1<<endl;
				double constX, constY;
				constX = (particle->PT) * sin(particle->Phi);
		  		constY = (particle->PT) * cos(particle->Phi);
				TLine* line = new TLine(0.0,0.0,constX,constY);
				line->SetLineColor(kOrange+10);
	        		line->SetLineWidth(2);
      	        		line->Draw("same"); 
			}
  		  }
  		  if (d2Index>=0)
		  {
			int d2=abs(static_cast<GenParticle*>(branchParticle->At(d2Index))->PID);
			if((d2==5 || d2==4) && pid==21)
			{ 
				cout<<j<<"DAAUG particle`s d2 PID: "<<d2<<endl;
				double constX, constY;
				constX = (particle->PT) * sin(particle->Phi);
		  		constY = (particle->PT) * cos(particle->Phi);
				TLine* line = new TLine(0.0,0.0,constX,constY);
				line->SetLineColor(kOrange+10);
	        		line->SetLineWidth(2);
      	        		line->Draw("same"); 
			}
  		  }
  		if (m1Index>=0)
		  {
			int m1=abs(static_cast<GenParticle*>(branchParticle->At(m1Index))->PID);
			if(m1==21 && (pid==5 || pid==4))
			{ 
				cout<<j<<": particle with PID "<<pid<<" and m1 PID: "<<m1<<endl;
				double constX, constY;
				constX = (particle->PT) * sin(particle->Phi);
		  		constY = (particle->PT) * cos(particle->Phi);
				TLine* line = new TLine(0.0,0.0,constX,constY);
				line->SetLineColor(kOrange+10);
	        		line->SetLineWidth(2);
      	        		line->Draw("same"); 
			}
  		  }
  		  if (m2Index>=0)
		  {
			int m2=abs(static_cast<GenParticle*>(branchParticle->At(m2Index))->PID);
			if(m2==21 && (pid==5 || pid==4))
			{ 
				cout<<j<<": particle with PID "<<pid<<" and m2 PID: "<<m2<<endl;
				double constX, constY;
				constX = (particle->PT) * sin(particle->Phi);
		  		constY = (particle->PT) * cos(particle->Phi);
				TLine* line = new TLine(0.0,0.0,constX,constY);
				line->SetLineColor(kOrange+10);
	        		line->SetLineWidth(2);
      	        		line->Draw("same"); 
			}
  		  }
  		  if (pid==5 || pid==4)
		  {
						double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2;
						
						double partEta = particle->Eta;
				     		double partPhi = particle->Phi;
				     		
				     		double genjetB1Eta = genjetb1->Eta;
			     			double genjetB1Phi = genjetb1->Phi;
			     			deltaRB1 = sqrt(pow(partEta-genjetB1Eta, 2) + pow(partPhi-genjetB1Phi, 2));
			     			
			     			double genjetB2Eta = genjetb2->Eta;
			     			double genjetB2Phi = genjetb2->Phi;
			     			deltaRB2 = sqrt(pow(partEta-genjetB2Eta, 2) + pow(partPhi-genjetB2Phi, 2));
			     			
			     			double genjetNB1Eta = genjetnb1->Eta;
			     			double genjetNB1Phi = genjetnb1->Phi;
			     			deltaRNB1 = sqrt(pow(partEta-genjetNB1Eta, 2) + pow(partPhi-genjetNB1Phi, 2));
			     			
			     			double genjetNB2Eta = genjetnb2->Eta;
			     			double genjetNB2Phi = genjetnb2->Phi;
			     			deltaRNB2 = sqrt(pow(partEta-genjetNB2Eta, 2) + pow(partPhi-genjetNB2Phi, 2));
			     			
			     			if (deltaRB1<0.4) contB1DeltaR04++;
			     			if (deltaRB2<0.4) contB2DeltaR04++;
			     			if (deltaRNB1<0.4) contNB1DeltaR04++;
			     			if (deltaRNB2<0.4) contNB2DeltaR04++;
  		  }
  		  
  		  
     		}
  
  cout<<endl<<endl;
  cout<<"contB1DeltaR04: "<<contB1DeltaR04<<endl;
  cout<<"contB2DeltaR04: "<<contB2DeltaR04<<endl;
  cout<<"contNB1DeltaR04: "<<contNB1DeltaR04<<endl;
  cout<<"contNB2DeltaR04: "<<contNB2DeltaR04<<endl;
  
  //Filling with tracks
  for(int i=0;i<ntracks;i++)
  {
	      Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
	 
	      double trackPhi = track->Phi;
	      double trackPt = track->PT;
	      double trackPx = trackPt*cos(trackPhi);
	      double trackPy = trackPt*sin(trackPhi);
	      //histJetPxAndPy->Fill(trackPy,trackPx);
	      
	      
	      /*
	      TLine* line = new TLine(0.0,0.0,trackPy,trackPx);
	      line->SetLineColor(kMagenta);
	      line->SetLineWidth(2);
	      line->Draw("same"); 
	      */
	      
	      
	      //cout<<"For track "<<i+1<<", pt: "<<trackPt<<endl;
	      //cout<<"For track "<<i+1<<", px: "<<trackPx<<endl;
	      //cout<<"For track "<<i+1<<", py: "<<trackPy<<endl;
	        // cout<<"for track "<<i+1<<", phi: "<<track->Phi<<endl<<endl;
  }
  //////
  
  // Show resulting histograms
  //histJetPx->Draw();
  //histJetPy->Draw();
  
  }
  
   
  
  
  /*
  
  
  
  cout<<endl<<endl<<endl<<endl<<"Starts Pz and Pt"<<endl<<endl;
  
  
  
  
  
  
  double pyMax=0, pzMax=0;
  if(branchEFlowTrack->GetEntries() > 0)
  {
    
      ntracks =  branchEFlowTrack->GetEntries();
  }
  else
  {
      ntracks=0;
  }
  cout<<"Number of tracks: "<<ntracks<<endl;
  for(int i=0;i<ntracks;i++)
  {
	      Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
	 
	      double eta = track->Eta;
	      double eEta = pow(TMath::E(), -eta);
	      double theta = 2*TMath::ATan(eEta);
	      double cosTheta = TMath::Cos(theta);
	      double trackPhi = track->Phi;
	      
	      double trackPt = track->PT;
    	      double trackP = track->P;
    	      
    	      double trackPtS=trackPt*trackPt;
    	      double trackPS=trackP*trackP;
    	      
	      double trackPz = trackP*cosTheta;
	      double trackPy = trackPt*sin(trackPhi);
	      
	      if(pyMax<(abs(trackPy))) pyMax=abs(trackPy);
	      if(pzMax<(abs(trackPz))) pzMax=abs(trackPz);
  }
  cout<<"For tracks, max py: "<<pyMax<<endl;
  cout<<"For tracks, max pz: "<<pzMax<<endl<<endl;
  
  if(ntracks>0)
  {
  // Book histograms
  TH1 *histJetPt = new TH1F("jet_pt", "jet pt", 100, -100.0, 100.0);
  TH1 *histJetPz = new TH1F("jet_pz", "jet pz", 100, 0.0, 300.0);
  TH2F *histJetPtAndPz = new TH2F("jet_pt_vs_pz", "jet pz vs. pt", 100, -pzMax*1.1, pzMax*1.1, 100, -pyMax*1.1, pyMax*1.1);
  
  c3->cd();
  histJetPtAndPz->Draw("col");
  
  histJetPtAndPz->SetStats(0);
  histJetPtAndPz->GetXaxis()->SetTitle("z");
  histJetPtAndPz->GetYaxis()->SetTitle("y");
  gStyle->SetPalette(1);
  histJetPtAndPz->SetFillColor(0);
  gStyle->SetOptStat(0);
  
  
  //TEllipse* detectorPtPz = new TEllipse(0.0,0.0,pMax,pMax);
    //          detectorPtPz->SetLineStyle(1);
      //        detectorPtPz->SetFillStyle(4000);  //makes it transparent
        //      detectorPtPz->SetLineColor(4);
   	  //    detectorPtPz->SetLineWidth(1);
	    //  detectorPtPz->Draw();
	      
  TLine* detectorTop = new TLine(-pzMax*1.1,pyMax,pzMax*1.1,pyMax);
	      detectorTop->SetLineColor(4);
	      detectorTop->SetLineWidth(2);
	      detectorTop->Draw("same");
	      
  TLine* detectorBottom = new TLine(-pzMax*1.1,-pyMax,pzMax*1.1,-pyMax);
	      detectorBottom->SetLineColor(4);
	      detectorBottom->SetLineWidth(2);
	      detectorBottom->Draw("same");	
	      
  TLine* detectorLeft = new TLine(-pzMax,-pyMax,-pzMax,pyMax);
	      detectorLeft->SetLineColor(4);
	      detectorLeft->SetLineWidth(2);
	      detectorLeft->Draw("same");    
	      
  TLine* detectorRight = new TLine(pzMax,-pyMax,pzMax,pyMax);
	      detectorRight->SetLineColor(4);
	      detectorRight->SetLineWidth(2);
	      detectorRight->Draw("same"); 	  	      
	
  TLine* zAxisPtPz = new TLine(-pzMax*1.1,0.0,pzMax,0.0);
	      zAxisPtPz->SetLineColor(1);
	      zAxisPtPz->SetLineWidth(1);
	      zAxisPtPz->Draw("same");
	      
  TLine* yAxisPtPz = new TLine(0.0,-pyMax*1.1,0.0,pyMax*1.1);
	      yAxisPtPz->SetLineColor(1);
	      yAxisPtPz->SetLineWidth(1);
	      yAxisPtPz->Draw("same"); 	
  	       	
  
  
  treeReader->ReadEntry(NthEvent);
  
  //Filling with tracks
  for(int i=0;i<ntracks;i++)
  {
	      Track *track = (Track*) branchEFlowTrack->At(i); // Take ith track
	 
	      double eta = track->Eta;
	      double eEta = pow(TMath::E(), -eta);
	      double theta = 2*TMath::ATan(eEta);
	      double cosTheta = TMath::Cos(theta);
	      double trackPhi = track->Phi;
	      
	      double trackPt = track->PT;
    	      double trackP = track->P;
    	      
    	      double trackPtS=trackPt*trackPt;
    	      double trackPS=trackP*trackP;
    	      
	      double trackPz = trackP*cosTheta;
	      double trackPy = trackPt*sin(trackPhi);
	      
	      TLine* line = new TLine(0.0,0.0,trackPz,trackPy);
	      line->SetLineColor(7);
	      line->SetLineWidth(2);
	      line->Draw("same");
	      cout<<"For track "<<i+1<<", p: "<<trackP<<endl; 
	      cout<<"For track "<<i+1<<", pt: "<<trackPt<<endl;
	      cout<<"For track "<<i+1<<", pz (w/theta): "<<trackPz<<endl<<endl;
  }

}*/
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
//cout<<"events: "<<numberOfEntries<<endl;
//cout<<"entries: "<<nentries<<endl;
  // Show resulting histograms
  //histJetPhi->Draw();
  //histJetEta->Draw();
  //histJetPhiVsEta->Draw("colz");
  
  cout<<"B: "<<contB<<endl<<"G: "<<contG<<endl<<endl;
  for(int i=0;i<guia;i++) cout<<dispar[i]<<",";
  cout<<endl;
}
