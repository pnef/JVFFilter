#ifndef MYFASTJETBASE_H
#define MYFASTJETBASE_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"



class MyUserInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
    MyUserInfo(const bool &is_pileup, const bool &is_track, const bool &is_charged): _is_pileup(is_pileup), _is_track(is_track), _is_charged(is_charged){}
    bool is_pileup()  const { return _is_pileup;}
    bool is_track()   const { return _is_track;}
    bool is_charged() const { return _is_charged;}
 protected:
    bool _is_pileup;  // true if ghost particle
    bool _is_track;   // true if track 
    bool _is_charged; // true if charged  
};


// Selector
class SelectorWorkerPUTracks : public fastjet::SelectorWorker {
    public:
    virtual bool pass(const fastjet::PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<MyUserInfo>()  
            &&  particle.user_info<MyUserInfo>().is_pileup() == true
            &&  particle.user_info<MyUserInfo>().is_track()  == true );
    }   

    virtual string description() const {return "track from pileup interaction";}
};

class SelectorWorkerHSTracks : public fastjet::SelectorWorker {
    public:
    virtual bool pass(const fastjet::PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<MyUserInfo>()  
            &&  particle.user_info<MyUserInfo>().is_pileup() == false
            &&  particle.user_info<MyUserInfo>().is_track()  == true );
    }   

    virtual string description() const {return "track from hs interaction";}
};

class SelectorWorkerTracks : public fastjet::SelectorWorker {
    public:
    virtual bool pass(const fastjet::PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<MyUserInfo>()  
            &&  particle.user_info<MyUserInfo>().is_track()  == true );
    }   

    virtual string description() const {return "tracks";}
};

class SelectorWorkerHS : public fastjet::SelectorWorker {
    public:
    virtual bool pass(const fastjet::PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<MyUserInfo>()  
            &&  particle.user_info<MyUserInfo>().is_pileup() == false);
    }   

    virtual string description() const {return "hard-scatter";}
};

#endif

