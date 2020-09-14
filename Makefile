OPT=-O3
ROOT=`root-config --cflags --glibs`  -lRooFit -lRooFitCore
CFLAGS=-I${CMSSW_BASE}/src
OBIN=./bin

.PHONY: all test

all: mt-2016 mt-2017 mt-2018 et-2016 et-2017 et-2018 ac-reweight create-fakes

mt-2016: plugins/mt_analyzer2016.cc
	g++ $(OPT) plugins/mt_analyzer2016.cc $(ROOT) $(CFLAGS) -o $(OBIN)/analyze2016_mt

mt-2017: plugins/mt_analyzer2016.cc
	g++ $(OPT) plugins/mt_analyzer2017.cc $(ROOT) $(CFLAGS) -o $(OBIN)/analyze2017_mt

mt-2018: plugins/mt_analyzer2016.cc
	g++ $(OPT) plugins/mt_analyzer2018.cc $(ROOT) $(CFLAGS) -o $(OBIN)/analyze2018_mt

et-2016: plugins/et_analyzer2016.cc
	g++ $(OPT) plugins/et_analyzer2016.cc $(ROOT) $(CFLAGS) -o $(OBIN)/analyze2016_et

et-2017: plugins/et_analyzer2016.cc
	g++ $(OPT) plugins/et_analyzer2017.cc $(ROOT) $(CFLAGS) -o $(OBIN)/analyze2017_et

et-2018: plugins/et_analyzer2016.cc
	g++ $(OPT) plugins/et_analyzer2018.cc $(ROOT) $(CFLAGS) -o $(OBIN)/analyze2018_et

ac-reweight:  plugins/ac_reweighter.cc
	g++ $(OPT) plugins/ac_reweighter.cc $(ROOT) $(CFLAGS) -o $(OBIN)/ac-reweight

create-fakes: plugins/fake_creater.cc
	g++ $(OPT) plugins/fake_creater.cc $(ROOT) $(CFLAGS) -o $(OBIN)/create-fakes

# Clean binaries
clean:
	rm $(OBIN)/*