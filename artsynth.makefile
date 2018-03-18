# Makefile of the library "artsynth"
# Paul Boersma, 8 August 2017

# include ../makefile.defs

OPTIONS = -I ./artsynth -I ./kar -I ./sys -I ./fon -I ./stat -I ./dwsys -I ./dwtools -I ./external/gsl -I ./LPC 
OPTIONS += -D NOMP3=1 -D NOFLAC=1 -D NOGRAPHICS=1 -D FLOATTYPE=double
OPTIONS += -Wlogical-op-parentheses -Wlogical-op-parentheses -Wdeprecated-register -Wshift-op-parentheses
OPTIONS += -O3
CFLAGS = $(OPTIONS)
CXXFLAGS = -std=c++11 $(OPTIONS)
LIBS=-lgsl -lsndfile

OBJECTS = main_artsynth.o artsynth/Artword_Speaker_to_Sound.o artsynth/Speaker.o artsynth/Articulation.o artsynth/Artword.o \
    artsynth/Delta.o artsynth/Speaker_to_Delta.o artsynth/Art_Speaker_Delta.o artsynth/Art_Speaker.o \
	sys/melder.o sys/melder_readtext.o sys/melder_error.o sys/melder_strings.o sys/melder_writetext.o sys/melder_files.o \
	sys/melder_debug.o sys/melder_textencoding.o sys/melder_alloc.o sys/melder_ftoa.o sys/melder_info.o \
	sys/melder_audiofiles.o kar/longchar.o sys/melder_atof.o sys/melder_time.o \
	sys/Thing.o sys/Data.o sys/abcio.o sys/tensor.o \
	dwsys/SVD.o dwsys/NUMclapack.o dwsys/NUMcblas.o dwsys/Eigen.o  \
	num/NUM.o num/NUMarrays.o num/NUMrandom.o num/NUMsort.o dwsys/NUM2.o dwsys/NUMf2c.o dwsys/NUMmachar.o dwsys/NUMsort2.o \
	dwsys/regularExp.o sys/Preferences.o sys/Collection.o sys/Simple.o \
	fon/Sound.o fon/Matrix.o fon/Vector.o fon/Sampled.o fon/SampledXY.o fon/Function.o fon/Sound_files.o

#	../external/gsl/gsl_cdf__fdist.o ../external/gsl/gsl_cdf__fdistinv.o ../external/gsl/gsl_cdf__beta.o ../external/gsl/gsl_cdf__betainv.o ../external/gsl/gsl_cdf__lognormal.o ../external/gsl/gsl_cdf__gauss.o \
#	../external/gsl/gsl_specfunc__gamma_inc.o ../external/gsl/gsl_specfunc__beta.o ../external/gsl/gsl_specfunc__gamma.o ../external/gsl/gsl_specfunc__exp.o ../external/gsl/gsl_specfunc__bessel_In.o \
#	../external/gsl/gsl_specfunc__bessel_I0.o

#  ../fon/praat_Matrix.o ../dwtools/Sound_extensions.o

#     Art_Speaker.o Art_Speaker_to_VocalTract.o Artword_Speaker.o Artword_Speaker_Sound.o \
#     Artword_Speaker_to_Sound.o Artword_to_Art.o \
#     ArtwordEditor.o praat_Artsynth.o manual_Artsynth.o

.PHONY: all clean

#all: libartsynth.a
all: run_artsynth

clean:
	$(RM) $(OBJECTS)
#	$(RM) libartsynth.a
	$(RM) run_artsynth

libartsynth.a: $(OBJECTS)
	touch libartsynth.a
	rm libartsynth.a
	$(AR) cq libartsynth.a $(OBJECTS)
	$(RANLIB) libartsynth.a

run_artsynth: $(OBJECTS)
	$(CXX) -o run_artsynth $(CPPFLAGS) $(LIBS) $(OBJECTS)

$(OBJECTS): ./artsynth/*.h ./kar/*.h ./sys/*.h ./fon/*.h ./stat/*.h ./dwsys/*.h ./dwtools/*.h
