#ifndef  TruthJetsAnalysis_H
#define  TruthJetsAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "TruthJetsTools.h"
#include "myFastJetBase.h"
#include "Pythia8/Pythia.h"

using namespace std;

class TruthJetsAnalysis{
    public:
        TruthJetsAnalysis ();
        ~TruthJetsAnalysis ();
        
        void Begin();
        void AnalyzeEvent(int iEvt, Pythia8::Pythia *pythia8, Pythia8::Pythia* pythia8PU, int nPU);
        void End();
        void Debug(int debug){
            fDebug = debug;
        }
        void SetOutName(string outname){
            fOutName = outname;
        }

        enum event_type {
            hardscatter,
            pileup
        } event_type_;

    private:
        void            DeclareBranches();
        void            ResetBranches();
        bool            GetPythiaParticles(Pythia8::Pythia* pythia8, event_type etype, int nPU=0);

        int  ftest;
        int  fDebug;
        int  nPU_;
        string fOutName;
        
        TruthJetsTools *tool;

        std::vector <fastjet::PseudoJet>           particles;

        TFile *tF;
        TTree *tT;

        // Tree Vars ---------------------------------------
        int              fTEventNumber;
        int              fTNPV;

        static const int MaxNParticles = 5000;
    
        int              fTNParticlesFilled;
        float            fTParticlePt              [MaxNParticles];
        float            fTParticleEta             [MaxNParticles];
        float            fTParticlePhi             [MaxNParticles];
        float            fTParticleIsHS            [MaxNParticles];

        float            fTParticlePtMomentPU01_03  [MaxNParticles];  
        float            fTParticlePtMomentHS01_03  [MaxNParticles]; 
        float            fTParticlePtMomentAll01_03 [MaxNParticles]; 
        float            fTParticlePtMomentPU02_03  [MaxNParticles]; 
        float            fTParticlePtMomentHS02_03  [MaxNParticles]; 
        float            fTParticlePtMomentAll02_03 [MaxNParticles]; 
        float            fTParticlePtMomentPU03_03  [MaxNParticles]; 
        float            fTParticlePtMomentHS03_03  [MaxNParticles]; 
        float            fTParticlePtMomentAll03_03 [MaxNParticles]; 
        float            fTParticlePtMomentPU03_05  [MaxNParticles]; 
        float            fTParticlePtMomentHS03_05  [MaxNParticles]; 
        float            fTParticlePtMomentAll03_05 [MaxNParticles]; 
        float            fTParticleJVFMoment01_03      [MaxNParticles];      
        float            fTParticleJVFMoment02_03      [MaxNParticles];  
        float            fTParticleJVFMoment03_03      [MaxNParticles];  
        float            fTParticleJVFMoment03_05      [MaxNParticles];          
        float            fTParticleJVFMoment01_03_gaus [MaxNParticles];   
        float            fTParticleJVFMoment02_03_gaus [MaxNParticles];  
        float            fTParticleJVFMoment03_03_gaus [MaxNParticles];  
        float            fTParticleJVFMoment03_05_gaus [MaxNParticles];  

        float            fTParticleAlphaC [MaxNParticles];
        float            fTParticleAlphaF [MaxNParticles];

};
#endif

