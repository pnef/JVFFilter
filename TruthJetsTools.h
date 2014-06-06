#ifndef TruthJetsTOOLS_H
#define TruthJetsTOOLS_H 

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  

#include "Pythia8/Pythia.h"

#include "myFastJetBase.h"

using namespace std;

class TruthJetsTools {
    private:

    public:
      TruthJetsTools();
        
      // methods
      float  JetCharge(fastjet::PseudoJet jet,float kappa);
      bool   BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, float jetrad, int BosonID );
      bool   TruthMatchDR(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float dR, float &truthmatchpt);
      float  dRTruth(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float ptThreshold);

      // moment calculator 
      float PtMoment(fastjet::PseudoJet particle, vector<fastjet::PseudoJet> &otherParticles, float sigma);
};

#endif

