DEBUG ?= 0
PARALLEL ?= 0

# External Libraries
BOOST_ROOT ?= /g/solexa/bin/software/boost_1_53_0/
BAMTOOLS_ROOT ?= /g/solexa/bin/software/bamtools-2.3.0/
SEQTK_ROOT ?= /g/solexa/bin/software/kseq/

# Flags
CXX=g++
CXXFLAGS += -isystem ${BOOST_ROOT}/include -isystem ${BAMTOOLS_ROOT}/include -isystem ${SEQTK_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -L${BOOST_ROOT}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -L${BAMTOOLS_ROOT}/lib -lbamtools -lz -Wl,-rpath,${BAMTOOLS_ROOT}/lib,-rpath,${BOOST_ROOT}/lib

# Additional flags for release/debug
ifeq (${PARALLEL}, 1)
	CXXFLAGS += -fopenmp -DOPENMP
else
	CXXFLAGS += -DNOPENMP
endif

# Additional flags for release/debug
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O9 -DNDEBUG
	#LDFLAGS += --static
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BAMTOOLSSOURCES = $(wildcard src/bamtools/src/api/*.h) $(wildcard src/bamtools/src/api/*.cpp)
DELLYSOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
TARGETS = src/delly src/extract src/cov src/iover src/stats

all: $(TARGETS)

src/delly: $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/extract: $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/cov: $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/iover: $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/stats: $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGETS) $(TARGETS:=.o) 
