#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/SubjetJVF.hh"

#include "TruthJetsTools.h"
#include "myFastJetBase.h"

#include "TRandom3.h"
#include "TError.h"
#include "TVector3.h"
#include "TMath.h"

using namespace std;

// Constructor 
TruthJetsTools::TruthJetsTools(){
}

// JVF Moment calculator 
float TruthJetsTools::JVFMomentCalculator(fastjet::PseudoJet particle, vector<fastjet::PseudoJet> &otherParticles,  bool gaussianWeight, float sigma, float maxdR){
    // JVF = HS / (PU + HS), HS = sum_{hs tracks} ( p_T * Gaus(dR) * H(dR - maxdR)), PU = sum_{pu tracks} ( p_T * Gaus(dR) * H(dR - maxdR)) 
    // Gausian weight is ignored if bool gaussianWeight=false.
    float HS=0;
    float PU=0;
    for(int ip=0; ip<otherParticles.size(); ++ip){
        if(particle == otherParticles[ip])               continue;
        if(particle.delta_R(otherParticles[ip]) > maxdR) continue;
        float weight = (gaussianWeight) ? TMath::Gaus(otherParticles[ip].delta_R(particle), 0, sigma, false):1; 
        if   (otherParticles[ip].user_info<MyUserInfo>().is_pileup()) PU += otherParticles[ip].pt() * weight; 
        else                                                          HS += otherParticles[ip].pt() * weight;
    }

    float result = -1;
    if( HS + PU >0) result = HS / (HS + PU);
    return result;
}

// Pt Moment calculator
float TruthJetsTools::PtMoment(fastjet::PseudoJet particle, vector<fastjet::PseudoJet> &otherParticles, float sigma, float maxdR){
    // sum p_T * Gaus(delta_R, sigma)*H(delta_R - maxdR), H=heavyside function
    float sum = 0;
    for(int ip=0; ip<otherParticles.size(); ++ip){
        if(particle == otherParticles[ip])               continue; 
        if(particle.delta_R(otherParticles[ip]) > maxdR) continue;
        sum  += otherParticles[ip].pt() * TMath::Gaus(otherParticles[ip].delta_R(particle), 0, sigma, false);
    }

    return sum;
}

// matching to truth jets: return true if dR to closest truth jet < dR
bool TruthJetsTools::TruthMatchDR(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float dR, float& truthmatchpt){
    float highestmatchpt = -999;
    for(unsigned int i=0; i<TruthJets.size(); ++i){
        if(jet.delta_R(TruthJets[i]) < dR) {
            if (TruthJets[i].pt() > highestmatchpt) highestmatchpt = TruthJets[i].pt();
        }
    }
    if(highestmatchpt  > 0) {truthmatchpt = highestmatchpt; return true;}
    else                    {truthmatchpt = -1;             return false;} 
}

// dR to closest truth jet above pT threshold
float TruthJetsTools::dRTruth(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float ptThreshold){
    float dR = 999.;
    for(unsigned int i=0; i<TruthJets.size(); ++i){
        if (TruthJets[i].pt() < ptThreshold ) continue;
        if (jet.delta_R(TruthJets[i]) < dR  ) dR = jet.delta_R(TruthJets[i]);
    }
    return dR;
}

