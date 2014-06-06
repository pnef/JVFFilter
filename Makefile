# --------------------------------------------- #
# Makefile for TruthJets code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall 
LINKLIBS := $(FASTJETLOCATION)/lib/libScJet.a

.PHONY: clean debug all

all: TruthJets

TruthJets:  TruthJets.so TruthJetsTools.so TruthJetsAnalysis.so $(LINKLIBS)
	$(CXX) TruthJets.so TruthJetsTools.so TruthJetsAnalysis.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lSubjetJVF -lVertexJets -lJetCleanser \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

TruthJets.so: TruthJets.C    TruthJetsTools.so TruthJetsAnalysis.so $(LINKLIBS)
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets -lJetCleanser \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

TruthJetsTools.so : TruthJetsTools.cc TruthJetsTools.h $(LINKLIBS)
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets -lJetCleanser \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs`

TruthJetsAnalysis.so : TruthJetsAnalysis.cc TruthJetsAnalysis.h $(LINKLIBS)
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets  -lJetCleanser \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs`   

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~
