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
                  fabs(pythia8->event[ip].id())  !=11  && 
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
    vector<fastjet::PseudoJet> putracks, hstracks, dummy;
    selector_hstracks.sift(particles, hstracks, dummy);
    selector_putracks.sift(particles, putracks, dummy);
    
    /// Calculate moments 
    for(int ip=0; ip<particles.size(); ++ip){
        float pt_moment_hs_0p1 = tool->PtMoment(particles[ip], hstracks, 0.1); 
        float pt_moment_pu_0p1 = tool->PtMoment(particles[ip], putracks, 0.1); 
//        cout << particles[ip].pt() << " pu " << particles[ip].user_info<MyUserInfo>().is_pileup() << " " << pt_moment_hs_0p1 << " "  << pt_moment_pu_0p1  << endl;

        // filling  ------------
        if(fTNParticlesFilled == MaxNParticles) continue;
        cout << fTNParticlesFilled << " " << MaxNParticles << endl;
        fTParticlePt             [fTNParticlesFilled] = particles[ip].pt();
        fTParticleEta            [fTNParticlesFilled] = particles[ip].eta();
        fTParticlePhi            [fTNParticlesFilled] = particles[ip].phi();
        fTParticlePtMomentPU01   [fTNParticlesFilled] = pt_moment_pu_0p1;
        fTParticlePtMomentHS01   [fTNParticlesFilled] = pt_moment_hs_0p1;
        fTParticleIsHS           [fTNParticlesFilled] = (particles[ip].user_info<MyUserInfo>().is_pileup()?0:1);
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

   // smallR jets
   tT->Branch("NParticlesFilled",          &fTNParticlesFilled,        "NParticlesFilled/I");
   tT->Branch("ParticlePt",                &fTParticlePt,              "ParticlePt           [NParticlesFilled]/F");
   tT->Branch("ParticleEta",               &fTParticleEta,             "ParticleEta          [NParticlesFilled]/F");
   tT->Branch("ParticlePhi",               &fTParticlePhi,             "ParticlePhi          [NParticlesFilled]/F");
   tT->Branch("ParticlePtMomentPU01",      &fTParticlePtMomentPU01,    "ParticlePtMomentPU01 [NParticlesFilled]/F");
   tT->Branch("ParticlePtMomentHS01",      &fTParticlePtMomentHS01,    "ParticlePtMomentHS01 [NParticlesFilled]/F");
   tT->Branch("ParticleIsHS",              &fTParticleIsHS,            "ParticleIsHS         [NParticlesFilled]/F");
   
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
          fTParticlePt           [iP]= -999;
          fTParticlePhi          [iP]= -999;
          fTParticleEta          [iP]= -999;
          fTParticleIsHS         [iP]= -999;
          fTParticlePtMomentPU01 [iP]= -999;
          fTParticlePtMomentHS01 [iP]= -999;
      }

}



