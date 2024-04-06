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

void trainDNN(const char* sigFilename, const char* bkgFilename, const char* outputFilename) {

    gRandom -> SetSeed(0);
    // Load the signal and background ROOT files
    TFile* sigFile = TFile::Open(sigFilename);
    TFile* bkgFile = TFile::Open(bkgFilename);

    // Retrieve the "HistData" TTrees from the files
    TTree* sigTree = dynamic_cast<TTree*>(sigFile->Get("HistData"));
    TTree* bkgTree = dynamic_cast<TTree*>(bkgFile->Get("HistData"));

    // Create a ROOT output file for TMVA
    TFile* outputFile = TFile::Open(outputFilename, "RECREATE");

    // Create the TMVA factory
    TMVA::Factory factory("TMVAClassification", outputFile,
                          "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

    // Create the TMVA dataloader
    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

    // Add the input variables to the dataloader
    dataloader->AddVariable("hpt", 'F');
    dataloader->AddVariable("htau21", 'F');
    dataloader->AddVariable("htau42", 'F');
    dataloader->AddVariable("hthrust", 'F');
    dataloader->AddVariable("hthrust_minor", 'F');
    dataloader->AddVariable("hsubratio", 'F');
    dataloader->AddVariable("hasym", 'F');
    dataloader->AddVariable("hRjj", 'F');
    dataloader->AddVariable("hsubjetwidth", 'F');
    dataloader->AddVariable("hD2", 'F');
    dataloader->AddVariable("hEntropy", 'F');
    dataloader->AddVariable("hPE", 'F');
    dataloader->AddVariable("hmass", 'F');
    dataloader->AddVariable("heta", 'F');
    dataloader->AddVariable("hphi", 'F');
    dataloader->AddVariable("hsub1eta", 'F');
    dataloader->AddVariable("hsub1phi", 'F');
    dataloader->AddVariable("hsub1pt", 'F');
    dataloader->AddVariable("hsub1mass", 'F');
    dataloader->AddVariable("hsub2eta", 'F');
    dataloader->AddVariable("hsub2phi", 'F');
    dataloader->AddVariable("hsub2pt", 'F');
    dataloader->AddVariable("hsub2mass", 'F');
    dataloader->AddVariable("lep1eta", 'F');
    dataloader->AddVariable("lep1phi", 'F');
    dataloader->AddVariable("lep1pt", 'F');
    dataloader->AddVariable("lep1mass", 'F');
    dataloader->AddVariable("lep2eta", 'F');
    dataloader->AddVariable("lep2phi", 'F');
    dataloader->AddVariable("lep2pt", 'F');
    dataloader->AddVariable("lep2mass", 'F');
    dataloader->AddVariable("hsub1charge", 'F');
    dataloader->AddVariable("hsub2charge", 'F');
    //dataloader->AddVariable("qg_tag", 'I');
    //dataloader->AddVariable("htruth_total", 'I');
    dataloader->AddVariable("hntotal", 'I');

    // Apply the mass cut as an event selection
    TCut massCut = "hmass > 100 && hmass < 150";
    TCut nanCut = "hsubjetwidth == hsubjetwidth && hsubjetwidth != TMath::Infinity() && hsubjetwidth != -TMath::Infinity()";
    TCut combinedCut = massCut && nanCut;

    // Add the signal and background trees to the dataloader
    dataloader->AddTree(sigTree, "Signal", 1.0, combinedCut);
    dataloader->AddTree(bkgTree, "Background", 1.0, combinedCut);

    // Prepare the training and testing data
    dataloader->PrepareTrainingAndTestTree(combinedCut, combinedCut,
                                            "SplitMode=Random:NormMode=NumEvents:V:nTrain_Signal=6180:nTest_Signal=2060:nTrain_Background=959:nTest_Background=320");

    // Book the BDT method
    TString methodName = "DNN";
    TString methodOptions = "!H:!V:VarTransform=N:Layout=RELU|64,RELU|128,RELU|16,SIGMOID:ErrorStrategy=CROSSENTROPY:TrainingStrategy=LearningRate=1e-3,Momentum=0.9,Repetitions=1,ConvergenceSteps=20,BatchSize=64,TestRepetitions=1,WeightDecay=1e-4,Regularization=None,Optimizer=Adam,DropConfig=0.0,DropRepetitions=0,MinimumLoss=0.0001,MaxFallbacksAfterConvergence=5:Architecture=CPU";
    factory.BookMethod(dataloader, TMVA::Types::kDL, methodName, methodOptions);

    // Train, test, and evaluate the DNN
    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    // Close the output file
    outputFile->Close();

    // Clean up
    sigFile->Close();
    bkgFile->Close();

    // Retrieve the trained DNN
    TMVA::MethodDNN* dnn = dynamic_cast<TMVA::MethodDNN*>(factory.GetMethod(dataloader->GetName(), methodName));

    // Create a canvas for the ROC curve
    TCanvas* rocCanvas = new TCanvas("rocCanvas", "ROC Curve", 800, 600);
    rocCanvas->SetGrid();

    // Get the ROC curve
    TGraph* rocGraph = (TGraph*)factory.GetROCCurve(dataloader->GetName(), methodName);
    rocGraph->SetTitle("ROC Curve");
    rocGraph->GetXaxis()->SetTitle("False Positive Rate");
    rocGraph->GetYaxis()->SetTitle("True Positive Rate");
    rocGraph->Draw("AL");

    // Save the ROC curve figure
    rocCanvas->SaveAs("DNN_roc_curve.png");
}