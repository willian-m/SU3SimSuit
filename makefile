PWD=$(shell pwd)
BIN=$(PWD)/bin
SRC=$(PWD)/src
MODULES=$(SRC)/modules
FC=gfortran

MODE=D
#MODE can be R (Release, with optimazation flags), D (Debug) or P (Profiling)

ifeq ($(FC),ifort)

ifeq ($(MODE),D)
FFLAGS=-heap-arrays 4096 -O0 -g -traceback -check all -ftrapuv -gen-interfaces -warn interfaces `$(SRC)/lib/FoX/FoX-config` 
endif

ifeq ($(MODE),R)
FFLAGS=-heap-arrays 4096 -O3 -xHost -fp-model precise `$(SRC)/lib/FoX/FoX-config` #fp-model precise may slow down the code, but insures precision.
endif

ifeq ($(MODE),P)
FFLAGS=-p -heap-arrays 4096 `$(SRC)/lib/FoX/FoX-config`
endif

MKL_LINK=-mkl
endif


ifeq ($(FC),gfortran)
MKL_LINK=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl
FFLAGS=-ffree-line-length-none `$(SRC)/lib/FoX/FoX-config`
endif

all: gen_lat_conf.run avrg_plaquette.run tmunu.run tensor_correlator.run green_function_tmunu.run


OBJ_LAT_CONF=$(BIN)/ziggurat.o $(BIN)/types_params.o $(BIN)/math.o $(BIN)/IO.o $(BIN)/lattice.o $(BIN)/objects.o $(BIN)/heat_bath.o $(BIN)/xml_parser.o
gen_lat_conf.run: dir $(OBJ_LAT_CONF) $(SRC)/gen_lat_conf.f90
	$(FC) $(FFLAGS) -I$(BIN) $(OBJ_LAT_CONF) `$(SRC)/lib/FoX/FoX-config` $(SRC)/gen_lat_conf.f90 -o $(BIN)/$@

OBJ_AVRG_PLAQUETTE=$(BIN)/types_params.o $(BIN)/ziggurat.o $(BIN)/math.o $(BIN)/IO.o $(BIN)/lattice.o $(BIN)/objects.o
avrg_plaquette.run: dir $(OBJ_AVRG_PLAQUETTE) $(SRC)/avrg_plaquette.f90
	$(FC) $(FFLAGS) -I$(BIN) $(OBJ_AVRG_PLAQUETTE) $(SRC)/avrg_plaquette.f90 -o $(BIN)/$@

OBJ_TMUNU=$(OBJ_AVRG_PLAQUETTE)
tmunu.run: dir $(OBJ_TMUNU) $(SRC)/tmunu.f90
	$(FC) $(FFLAGS) -I$(BIN) $(OBJ_TMUNU) $(SRC)/tmunu.f90 -o $(BIN)/$@

OBJ_TENSOR_CORRELATOR=$(BIN)/ziggurat.o $(BIN)/types_params.o $(BIN)/math.o $(BIN)/IO.o $(BIN)/lattice.o $(BIN)/xml_parser.o $(BIN)/statistic.o
tensor_correlator.run: dir $(OBJ_TENSOR_CORRELATOR) $(SRC)/tmunu.f90
	$(FC) $(FFLAGS) -I$(BIN) -I$(MKLROOT)/include  $(OBJ_TENSOR_CORRELATOR) `$(SRC)/lib/FoX/FoX-config` $(SRC)/tensor_correlator.f90 $(MKL_LINK) -o $(BIN)/$@

OBJ_GREEN_FUNCTION=$(BIN)/ziggurat.o $(BIN)/types_params.o $(BIN)/math.o $(BIN)/IO.o $(BIN)/lattice.o $(BIN)/xml_parser.o $(BIN)/statistic.o
green_function_tmunu.run: dir $(OBJ_TENSOR_CORRELATOR) $(SRC)/tmunu.f90
	$(FC) $(FFLAGS) -I$(BIN) -I$(MKLROOT)/include  $(OBJ_TENSOR_CORRELATOR) `$(SRC)/lib/FoX/FoX-config` $(SRC)/green_function_tmunu.f90 $(MKL_LINK) -o $(BIN)/$@

dir: 
	mkdir -p $(BIN)

$(BIN)/%.o: $(MODULES)/%.f90
	ln -sf $(MKLROOT)/include/mkl_dfti.f90 $(MODULES)/mkl_dfti.f90 
	if [ $(FC) = ifort ]; then $(FC) -I=$(MODULES) $(MKL_LINK) $(FFLAGS) -c -module $(BIN) -o $@ $<; elif [ $(FC) = gfortran ]; then $(FC) $(MKL_LINK) $(FFLAGS) -c -J$(BIN) -o $@ $<; fi

clean:
	rm -f $(BIN)/*.o $(BIN)/*.mod $(BIN)/*.run
	rmdir $(BIN)
#	rm -rf output

sandwich:
	@USER="$(id -u)"
	@if [ $(USER) != "root" ]; then echo "Nice try, but one should not search for make commands in an internet comic strip."; fi;
	@if [ $(USER) = "root" ]; then echo "Ok. Here is your sandwich"; echo "                                                           _"; echo "                                                          //"; echo "                                                         //"; echo "                                         _______________//__"; echo "                                       .(______________//___)."; echo "                                       |              /      |"; echo "                                       |. . . . . . . / . . .|"; echo "                                       \ . . . . . ./. . . . /"; echo "                                        |           / ___   |"; echo "                   _.---._              |::......./../...\.:|"; echo "                _.-~       ~-._         |::::/::\::/:\::::::|"; echo "            _.-~               ~-._     |::::\::/::::::X:/::|"; echo "        _.-~                       ~---.;:::::::/::\::/:::::|"; echo "    _.-~                                 ~\::::::n::::::::::|"; echo " .-~                                    _.;::/::::a::::::::/"; echo " :-._                               _.-~ ./::::::::d:::::::|"; echo " \`-._~-._                   _..__.-~ _.-~|::/::::::::::::::|"; echo "  /  ~-._~-._              / .__..--~----.YWWWWWWWWWWWWWWWP'"; echo " \_____(_;-._\.        _.-~_/       ~).. . \ "; echo "    /(_____  \`--...--~_.-~______..-+_______)"; echo "  .(_________/\`--...--~/    _/           /\ "; echo " /-._     \_     (___./_..-~__.....__..-~./"; echo " \`-._~-._   ~\--------~  .-~_..__.-~ _.-~"; echo "     ~-._~-._ ~---------'  / .__..--~"; echo "         ~-._\.        _.-~_/"; echo "             \`--...--~_.-~"; echo "              \`--...--~"; fi;
