#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"

#include "TruthJetsAnalysis.h"
#include "TruthJetsTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// --------------------------------------------------------------
// Fastjet Selectors ----------------------------------------------
fastjet::Selector SelectorHSTracks() {
    return new SelectorWorkerHSTracks();
}

fastjet::Selector SelectorPUTracks() {
    return new SelectorWorkerPUTracks();
}

fastjet::Selector SelectorTracks() {
    return new SelectorWorkerTracks();
}

fastjet::Selector SelectorHS() {
    return new SelectorWorkerHS();
}

// -------------------------------------------------------------------
// Constructor 
TruthJetsAnalysis::TruthJetsAnalysis(){
    if(fDebug) cout << "TruthJetsAnalysis::TruthJetsAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new TruthJetsTools();

    if(fDebug) cout << "TruthJetsAnalysis::TruthJetsAnalysis End " << endl;
}

// Destructor 
TruthJetsAnalysis::~TruthJetsAnalysis(){
    delete tool;
}

// Begin method
void TruthJetsAnalysis::Begin(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for TruthJets");
    
   DeclareBranches();
   ResetBranches();
   

   return;
}

// End
void TruthJetsAnalysis::End(){
    
    tT->Write();
    tF->Close();
    return;
}

// GetParticles from pythia
bool TruthJetsAnalysis::GetPythiaParticles(Pythia8::Pythia* pythia8, event_type etype, int nPU){

    int n_runs = 1;
    if(etype == pileup){ n_runs = nPU;}

    for(int irun=0; irun<n_runs; ++irun){
        // next event 
        if (!pythia8->next()) return false;

        for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){

            int origin    = -1;
            if      (etype == hardscatter) origin = 0;
            else if (etype == pileup     ) origin = irun +1;

            fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() );
            p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.); // massless particles 

            // tracks: charged particles with pt>0.5 GeV, |eta|<2.4 
            if(! (pythia8->event[ip].isFinal()    && 
                  fabs(pythia8->event[ip].id())  !=12  && 
                  fabs(pythia8->event[ip].id())  !=13  && 
                  fabs(pythia8->event[ip].id())  !=14  && 
                  fabs(pythia8->event[ip].id())  !=16  && 
                  pythia8->event[ip].pT() > 0.5         ) ) continue; 

            // fill containers ------
            bool is_track   = (fabs(p.eta())<2.4  && pythia8->event[ip].isCharged());
            bool is_pu      = (origin >0?true:false);
            bool is_charged = pythia8->event[ip].isCharged();

            p.set_user_info(new MyUserInfo(is_pu, is_track, is_charged  ));
            particles.push_back(p); 



        }
    } // end irun loop

    return true;
}

// Analyze
void TruthJetsAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia8PU, int nPU){
    if(fDebug) cout << "TruthJetsAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if(fDebug) cout << "TruthJetsAnalysis::AnalyzeEvent Event Number " << ievt << endl;
    
    // reset branches 
    ResetBranches();
    

    // new event-----------------------
    fTEventNumber = ievt;
    fTNPV         = nPU;
    particles.clear();
    
    // get new event and return if failed
    bool ok(true);
    ok = GetPythiaParticles(pythia8,   hardscatter);
    ok = GetPythiaParticles(pythia8PU, pileup,     nPU);
    if (!ok) return;
    


    // get pu track and hs tracks 
    fastjet::Selector selector_hs       =  SelectorHS();
    fastjet::Selector selector_tracks   =  SelectorTracks();
    fastjet::Selector selector_hstracks =  selector_hs && selector_tracks;
    fastjet::Selector selector_putracks = !selector_hs && selector_tracks;
    vector<fastjet::PseudoJet> putracks, hstracks, tracks, dummy;
    selector_hstracks.sift(particles, hstracks, dummy);
    selector_putracks.sift(particles, putracks, dummy);
    selector_tracks  .sift(particles, tracks, dummy);
    
    /// Calculate moments 
    for(int ip=0; ip<particles.size(); ++ip){
        // pT moments 
        float pt_moment_hs_0p1_0p3  = tool->PtMoment(particles[ip], hstracks,  0.1, 0.3); 
        float pt_moment_pu_0p1_0p3  = tool->PtMoment(particles[ip], putracks,  0.1, 0.3); 
        float pt_moment_All_0p1_0p3 = tool->PtMoment(particles[ip], particles, 0.1, 0.3); 
        float pt_moment_hs_0p2_0p3  = tool->PtMoment(particles[ip], hstracks,  0.2, 0.3); 
        float pt_moment_pu_0p2_0p3  = tool->PtMoment(particles[ip], putracks,  0.2, 0.3); 
        float pt_moment_All_0p2_0p3 = tool->PtMoment(particles[ip], particles, 0.2, 0.3); 
        float pt_moment_hs_0p3_0p3  = tool->PtMoment(particles[ip], hstracks,  0.3, 0.3); 
        float pt_moment_pu_0p3_0p3  = tool->PtMoment(particles[ip], putracks,  0.3, 0.3); 
        float pt_moment_All_0p3_0p3 = tool->PtMoment(particles[ip], particles, 0.3, 0.3); 
        float pt_moment_hs_0p3_0p5  = tool->PtMoment(particles[ip], hstracks,  0.3, 0.5); 
        float pt_moment_pu_0p3_0p5  = tool->PtMoment(particles[ip], putracks,  0.3, 0.5); 
        float pt_moment_All_0p3_0p5 = tool->PtMoment(particles[ip], particles, 0.3, 0.5); 

        // JVF moments 
        float JVF_moment_0p1_0p3       = tool->JVFMomentCalculator(particles[ip], tracks,     false ,0.1, 0.3); 
        float JVF_moment_0p2_0p3       = tool->JVFMomentCalculator(particles[ip], tracks,     false ,0.2, 0.3); 
        float JVF_moment_0p3_0p3       = tool->JVFMomentCalculator(particles[ip], tracks,     false ,0.3, 0.3); 
        float JVF_moment_0p3_0p5       = tool->JVFMomentCalculator(particles[ip], tracks,     false ,0.3, 0.5); 
        float JVF_moment_0p1_0p3_gaus  = tool->JVFMomentCalculator(particles[ip], tracks,     true  ,0.1, 0.3); 
        float JVF_moment_0p2_0p3_gaus  = tool->JVFMomentCalculator(particles[ip], tracks,     true  ,0.2, 0.3); 
        float JVF_moment_0p3_0p3_gaus  = tool->JVFMomentCalculator(particles[ip], tracks,     true  ,0.3, 0.3); 
        float JVF_moment_0p3_0p5_gaus  = tool->JVFMomentCalculator(particles[ip], tracks,     true  ,0.3, 0.5); 

        // CMS alpha_F 
        float alpha_F = tool->AlphaMoment(particles[ip], particles, 0.3);
        float alpha_C = tool->AlphaMoment(particles[ip], hstracks,  0.3);

        // filling  ------------
        if(fTNParticlesFilled == MaxNParticles) continue;
        fTParticlePt                [fTNParticlesFilled] = particles[ip].pt();
        fTParticleEta               [fTNParticlesFilled] = particles[ip].eta();
        fTParticlePhi               [fTNParticlesFilled] = particles[ip].phi();
        fTParticleIsHS              [fTNParticlesFilled] = (particles[ip].user_info<MyUserInfo>().is_pileup()?0:1);

        fTParticlePtMomentPU01_03   [fTNParticlesFilled] = pt_moment_pu_0p1_0p3;
        fTParticlePtMomentHS01_03   [fTNParticlesFilled] = pt_moment_hs_0p1_0p3;
        fTParticlePtMomentAll01_03  [fTNParticlesFilled] = pt_moment_All_0p1_0p3;
        fTParticlePtMomentPU02_03   [fTNParticlesFilled] = pt_moment_pu_0p2_0p3;
        fTParticlePtMomentHS02_03   [fTNParticlesFilled] = pt_moment_hs_0p2_0p3;
        fTParticlePtMomentAll02_03  [fTNParticlesFilled] = pt_moment_All_0p2_0p3;
        fTParticlePtMomentPU03_03   [fTNParticlesFilled] = pt_moment_pu_0p3_0p3;
        fTParticlePtMomentHS03_03   [fTNParticlesFilled] = pt_moment_hs_0p3_0p3;
        fTParticlePtMomentAll03_03  [fTNParticlesFilled] = pt_moment_All_0p3_0p3;
        fTParticlePtMomentPU03_05   [fTNParticlesFilled] = pt_moment_pu_0p3_0p5;
        fTParticlePtMomentHS03_05   [fTNParticlesFilled] = pt_moment_hs_0p3_0p5;
        fTParticlePtMomentAll03_05  [fTNParticlesFilled] = pt_moment_All_0p3_0p5;

        fTParticleJVFMoment01_03        [fTNParticlesFilled] = JVF_moment_0p1_0p3;
        fTParticleJVFMoment02_03        [fTNParticlesFilled] = JVF_moment_0p2_0p3;
        fTParticleJVFMoment03_03        [fTNParticlesFilled] = JVF_moment_0p3_0p3;
        fTParticleJVFMoment03_05        [fTNParticlesFilled] = JVF_moment_0p3_0p5;
        fTParticleJVFMoment01_03_gaus   [fTNParticlesFilled] = JVF_moment_0p1_0p3_gaus;
        fTParticleJVFMoment02_03_gaus   [fTNParticlesFilled] = JVF_moment_0p2_0p3_gaus;
        fTParticleJVFMoment03_03_gaus   [fTNParticlesFilled] = JVF_moment_0p3_0p3_gaus;
        fTParticleJVFMoment03_05_gaus   [fTNParticlesFilled] = JVF_moment_0p3_0p5_gaus;

        fTParticleAlphaF                [fTNParticlesFilled] = alpha_F;
        fTParticleAlphaC                [fTNParticlesFilled] = alpha_C;

        fTNParticlesFilled++;

    }



    // Fill
    tT->Fill();

    if(fDebug) cout << "TruthJetsAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void TruthJetsAnalysis::DeclareBranches(){
   
   // Event Properties 
   tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
   tT->Branch("NPV",                       &fTNPV,                    "NPV/I");

   // smAllR jets
   tT->Branch("NParticlesFilled",          &fTNParticlesFilled,        "NParticlesFilled/I");
   tT->Branch("Pt",                &fTParticlePt,              "Pt           [NParticlesFilled]/F");
   tT->Branch("Eta",               &fTParticleEta,             "Eta          [NParticlesFilled]/F");
   tT->Branch("Phi",               &fTParticlePhi,             "Phi          [NParticlesFilled]/F");
   tT->Branch("IsHS",              &fTParticleIsHS,            "IsHS         [NParticlesFilled]/F");
   tT->Branch("PtMomentPU01_03",   &fTParticlePtMomentPU01_03, "PtMomentPU01_03[NParticlesFilled]/F");
   tT->Branch("PtMomentHS01_03",   &fTParticlePtMomentHS01_03, "PtMomentHS01_03[NParticlesFilled]/F");
   tT->Branch("PtMomentAll01_03",  &fTParticlePtMomentAll01_03,"PtMomentAll01_03[NParticlesFilled]/F");
   tT->Branch("PtMomentPU02_03",   &fTParticlePtMomentPU02_03, "PtMomentPU02_03[NParticlesFilled]/F");
   tT->Branch("PtMomentHS02_03",   &fTParticlePtMomentHS02_03, "PtMomentHS02_03[NParticlesFilled]/F");
   tT->Branch("PtMomentAll02_03",  &fTParticlePtMomentAll02_03,"PtMomentAll02_03[NParticlesFilled]/F");
   tT->Branch("PtMomentPU03_03",   &fTParticlePtMomentPU03_03, "PtMomentPU03_03[NParticlesFilled]/F");
   tT->Branch("PtMomentHS03_03",   &fTParticlePtMomentHS03_03, "PtMomentHS03_03[NParticlesFilled]/F");
   tT->Branch("PtMomentAll03_03",  &fTParticlePtMomentAll03_03,"PtMomentAll03_03[NParticlesFilled]/F");
   tT->Branch("PtMomentPU03_05",   &fTParticlePtMomentPU03_05, "PtMomentPU03_05[NParticlesFilled]/F");
   tT->Branch("PtMomentHS03_05",   &fTParticlePtMomentHS03_05, "PtMomentHS03_05[NParticlesFilled]/F");
   tT->Branch("PtMomentAll03_05",  &fTParticlePtMomentAll03_05,"PtMomentAll03_05[NParticlesFilled]/F");

   tT->Branch("JVFMoment01_03",         &fTParticleJVFMoment01_03,       "JVFMoment01_03[NParticlesFilled]/F");
   tT->Branch("JVFMoment02_03",         &fTParticleJVFMoment02_03,       "JVFMoment02_03[NParticlesFilled]/F");
   tT->Branch("JVFMoment03_03",         &fTParticleJVFMoment03_03,       "JVFMoment03_03[NParticlesFilled]/F");
   tT->Branch("JVFMoment03_05",         &fTParticleJVFMoment03_05,       "JVFMoment03_05[NParticlesFilled]/F");
   tT->Branch("JVFMoment01_03_gaus",    &fTParticleJVFMoment01_03_gaus,  "JVFMoment01_03_gaus[NParticlesFilled]/F");
   tT->Branch("JVFMoment02_03_gaus",    &fTParticleJVFMoment02_03_gaus,  "JVFMoment02_03_gaus[NParticlesFilled]/F");
   tT->Branch("JVFMoment03_03_gaus",    &fTParticleJVFMoment03_03_gaus,  "JVFMoment03_03_gaus[NParticlesFilled]/F");
   tT->Branch("JVFMoment03_05_gaus",    &fTParticleJVFMoment03_05_gaus,  "JVFMoment03_05_gaus[NParticlesFilled]/F");
   
   tT->Branch("ParticleAlphaC",    &fTParticleAlphaC,  "ParticleAlphaC[NParticlesFilled]/F");
   tT->Branch("ParticleAlphaF",    &fTParticleAlphaF,  "ParticleAlphaF[NParticlesFilled]/F");
   tT->GetListOfBranches()->ls();
    
   return;
}


// resets vars
void TruthJetsAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNPV                         = -999;
      fTNParticlesFilled            = 0;

      for (int iP=0; iP < MaxNParticles; ++iP){
          fTParticlePt               [iP]= -999;
          fTParticlePhi              [iP]= -999;
          fTParticleEta              [iP]= -999;
          fTParticleIsHS             [iP]= -999;

          fTParticlePtMomentPU01_03  [iP]= -999;
          fTParticlePtMomentHS01_03  [iP]= -999;
          fTParticlePtMomentAll01_03 [iP]= -999;
          fTParticlePtMomentPU02_03  [iP]= -999;
          fTParticlePtMomentHS02_03  [iP]= -999;
          fTParticlePtMomentAll02_03 [iP]= -999;
          fTParticlePtMomentPU03_03  [iP]= -999;
          fTParticlePtMomentHS03_03  [iP]= -999;
          fTParticlePtMomentAll03_03 [iP]= -999;
          fTParticlePtMomentPU03_05  [iP]= -999;
          fTParticlePtMomentHS03_05  [iP]= -999;
          fTParticlePtMomentAll03_05 [iP]= -999;
          
          fTParticleJVFMoment01_03      [iP]= -999;
          fTParticleJVFMoment02_03      [iP]= -999;
          fTParticleJVFMoment03_03      [iP]= -999;
          fTParticleJVFMoment03_05      [iP]= -999;
          fTParticleJVFMoment01_03_gaus [iP]= -999;
          fTParticleJVFMoment02_03_gaus [iP]= -999;
          fTParticleJVFMoment03_03_gaus [iP]= -999;
          fTParticleJVFMoment03_05_gaus [iP]= -999;

          fTParticleAlphaC              [iP]=-999;
          fTParticleAlphaF              [iP]=-999;

      }
}



