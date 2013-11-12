
ifeq ($(CXX),)
	CXX = g++
endif

ifeq ($(DEBUG),yes)
	LDFLAGS = -lm -lboost_program_options
	CXXFLAGS = -Wall -g3 -DDEBUG -DVERBOSE -Icommon -Itinyxml -Iga -Idp
else
	LDFLAGS = -lm -Wl,-Bstatic -lboost_program_options-mt -Wl,-Bdynamic
	CXXFLAGS = -Wall -DNDEBUG -ffast-math -fcaller-saves -finline-functions -O3 -Icommon -Itinyxml -Iga -Idp -Ianalysis
	#CXXFLAGS = -Wall -fcaller-saves -finline-functions -O3 -Icommon -Itinyxml -Iga
endif

all: csodp csoga csora csogenerator recipe2dag csoanalyse

csodp: common/cso.o common/data.o common/stopwatch.o common/genotype.o common/genotypegamete.o dp/dpitem.o dp/dptablearray.o dp/dptable.o dp/dptablehashmap.o dp/dpmain.o dp/solver.o dp/heuristicsolver.o dp/heuristicsolverfast.o dp/heuristicpurelinkage.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o analysis/linkageanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

csoga: common/cso.o common/data.o common/stopwatch.o common/genotype.o common/genotypegamete.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o analysis/linkageanalysis.o ga/genotypetable.o ga/genotypetablearray.o ga/genotypetablehashmap.o ga/node.o ga/leafnode.o ga/innernode.o ga/gaindividual.o ga/gasolver.o ga/gamain.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

csora: common/cso.o common/data.o common/stopwatch.o common/genotype.o common/genotypegamete.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o ga/genotypetable.o ga/genotypetablearray.o ga/genotypetablehashmap.o ga/node.o ga/leafnode.o ga/innernode.o ga/gaindividual.o ra/rasolver.o ra/ramain.o analysis/coveranalysis.o analysis/linkageanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

csogenerator: common/cso.o common/data.o common/genotype.o common/genotypegamete.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o analysis/linkageanalysis.o generator/inputgen.o analysis/linkageanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

recipe2dag: common/cso.o common/data.o common/genotype.o common/genotypegamete.o dp/dpitem.o dp/dptablearray.o dp/dptable.o dp/dptablehashmap.o dp/solver.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o recipe/recipemain.o analysis/linkageanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

csoanalyse: common/cso.o common/data.o common/genotype.o common/genotypegamete.o tinyxml/tinystr.o tinyxml/tinyxml.o tinyxml/tinyxmlerror.o tinyxml/tinyxmlparser.o analysis/linkageanalysis.o analysis/coveranalysis.o analysis/analysemain.o 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

Makefile.depend: *.h *.cpp Makefile
	$(CC) -M *.cpp > Makefile.depend

clean:
	rm -f *.o csoga csodp csora csogenerator recipe2dag csoanalyse common/*.o ga/*.o dp/*.o tinyxml/*.o ra/*.o recipe/*.o inputgen/*.o analysis/*.o generator/*.o

-include Makefile.depend
