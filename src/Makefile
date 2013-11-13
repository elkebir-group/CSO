SYSTEM        = x86-64_sles10_4.1
#SYSTEM        = x86-64_darwin
LIBFORMAT     = static_pic
ILOGDIR       = $(HOME)/ILOG/CPLEX
CPLEXDIR      = $(ILOGDIR)/cplex
CONCERTDIR    = $(ILOGDIR)/concert
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
BINDIR        = .

ifeq ($(DEBUG),yes)
	CXXFLAGS = -Wall -g3 $(SVNDEV) -DDEBUG -DVERBOSE -Ilemon/include/ -I./ -Icommon/ -Ianalysis/ -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include
	LDFLAGS := -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -lilocplex -lcplex -lconcert -lm -pthread #-framework CoreFoundation -framework IOKit
else
	CXXFLAGS = -Wall -O3 $(SVNDEV) -ffast-math -fcaller-saves -finline-functions -DNDEBUG -I./ -Ilemon/include/ -Icommon/ -Ianalysis/ -DIL_STD -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include
	LDFLAGS := -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -lilocplex -lcplex -lconcert -lm -pthread #-framework CoreFoundation -framework IOKit
endif

all: csoilp2 processDAG computeLB

#csoilp2: ilpsolverheuristic.o crossingschedule.o ilpsolver.o ilpsolverext.o main.o lemon/arg_parser.o analysis/linkageanalysis.o analysis/lowerbounds.o dp/dpitem.o dp/dptable.o dp/dptablehashmap.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o
csoilp2: ilpsolverheuristic.o crossingschedule.o ilpsolver.o ilpsolverext.o main.o lemon/arg_parser.o analysis/linkageanalysis.o analysis/lowerbounds.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS)

csoilp: schedule3.o dp/dpitem.o dp/dptable.o dp/dptablehashmap.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS) -lboost_program_options

#processDAG: crossingschedule.o processDAG.o lemon/arg_parser.o dp/dpitem.o dp/dptable.o dp/dptablehashmap.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o
processDAG: crossingschedule.o processDAG.o lemon/arg_parser.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS)

#computeLB: computeLB.o lemon/arg_parser.o analysis/linkageanalysis.o analysis/lowerbounds.o dp/dpitem.o dp/dptable.o dp/dptablehashmap.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o 
computeLB: computeLB.o lemon/arg_parser.o analysis/linkageanalysis.o analysis/lowerbounds.o common/stopwatch.o common/data.o common/cso.o tinyxml/tinyxml.o tinyxml/tinyxmlparser.o tinyxml/tinyxmlerror.o tinyxml/tinystr.o common/genotype.o 
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS)

csogenerator: common/cso.o common/data.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o analysis/linkageanalysis.o generator/inputgen.o analysis/linkageanalysis.o
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS) -lboost_program_options

clean:
	rm -f *.o tinyxml/*.o common/*.o analysis/*.o generator/*.o dp/*.o csoilp2 csoilp processDAG computeLB csogenerator
