#include "Artword_Speaker_to_Sound.h"
#include "Speaker.h"

#include <vector>
#include <iostream>

#include <sndfile.h>

#define SAMPLERATE 22050
#define BLOCKLEN 1
#define OVERSAMPLING 25
#define DURATION 1.
#define MASSES 2  // 1, 2, 10


inline void sf_write(SNDFILE *sndfile, const float *ptr, sf_count_t frames)
{
    sf_writef_float(sndfile, ptr, frames);
}

inline void sf_write(SNDFILE *sndfile, const double *ptr, sf_count_t frames)
{
    sf_writef_double(sndfile, ptr, frames);
}

static const char *getCmdOption(int argc, const char *const *argv, const char *option)
{
    const size_t optlen = strlen(option);
    for(int i = 1; i < argc; ++i)
        if(!strncmp(argv[i], option, optlen)) {
            const char *res = argv[i]+optlen;
            while(res && isspace(*res)) ++res; // ignore whitespace
            if(*res == '=') {
                ++res; // ignore '='
                while(res && isspace(*res)) ++res;  // ignore whitespace after =
            }
            return res;
        }
    return NULL;
}

template <typename FLOAT>
struct Pair:
std::pair<double, FLOAT>
{
    Pair(double t=-1, FLOAT v=0): std::pair<double, FLOAT>(t, v) {}
};

template <typename FLOAT>
static const char *getPair(const char *option, Pair<FLOAT> &time_value)
{
    char *endptr;
    while(option && isspace(*option)) ++option;
    if(!option || *option != '(') throw std::runtime_error("At parsing (time,value) tuple: opening '(' missing");
    ++option; // skip (
    double time = strtod(option, &endptr);
    option = endptr;
    while(option && isspace(*option)) ++option;
    if(!option || *option != ',') throw std::runtime_error("At parsing (time,value) tuple: comma ',' missing");
    ++option; // skip ,
    while(option && isspace(*option)) ++option;
    double value = strtod(option, &endptr);
    option = endptr;
    if(!option || *option != ')') throw std::runtime_error("At parsing (time,value) tuple: closing ')' missing ");
    ++option; // skip )
    while(option && isspace(*option)) ++option;
    
    time_value.first = time;
    time_value.second = value;
    return option;
}

template <typename FLOAT>
static std::vector<Pair<FLOAT> > *getPairList(const char *option)
{
    std::vector<Pair<FLOAT> > *pairs = new std::vector<Pair<FLOAT> >;
    try {
        while(option && *option) {
            Pair<FLOAT> pair;
            option = getPair(option, pair);
            if(pair.first >= 0)
                pairs->push_back(pair);
            while(option && isspace(*option)) ++option;
            if(!option || !*option)
                break;
            if(*option != ',')
                throw std::runtime_error("At parsing (time,value) tuples: use ',' as separator ");
            option++; // skip ,
            while(option && isspace(*option)) ++option;
        }
    }
    catch(std::exception &x) {
        std::cerr << "Error: " << x.what() << std::endl;
    }
    return pairs;
}

const char *muscle_names[] = {
    "_",
    "Lungs",
    "Interarytenoid",
    "Cricothyroid",
    "Vocalis",
    "Thyroarytenoid",
    "PosteriorCricoarytenoid",
    "LateralCricoarytenoid",
    "Stylohyoid",
    "Sternohyoid",
    "Thyropharyngeus",
    "LowerConstrictor",
    "MiddleConstrictor",
    "UpperConstrictor",
    "Sphincter",
    "Hyoglossus",
    "Styloglossus",
    "Genioglossus",
    "UpperTongue",
    "LowerTongue",
    "TransverseTongue",
    "VerticalTongue",
    "Risorius",
    "OrbicularisOris",
    "LevatorPalatini",
    "TensorPalatini",
    "Masseter",
    "Mylohyoid",
    "LateralPterygoid",
    "Buccinator"
};

int main(int argc, const char *const *argv)
{
    NUMrandom_init();
    
    typedef float sample_t;
    
    const char *help_option = getCmdOption(argc, argv, "--help");
    if(help_option) {
        std::cerr << argv[0] << " does articulatory synthesis derived from PRAAT." << std::endl;
        std::cerr << "Available options:" << std::endl;
        std::cerr << "--help: Get this help" << std::endl;
        std::cerr << "--out: Set output file name" << std::endl;
        std::cerr << "--sr: Set audio sample rate (default=" << SAMPLERATE << ")" << std::endl;
        std::cerr << "--blk: Set block length (default=" << BLOCKLEN << ")" << std::endl;
        std::cerr << "--ovs: Set oversampling factor (default=" << OVERSAMPLING << ")" << std::endl;
        std::cerr << "--dur: Set articulation duration in seconds (default=" << DURATION << ")" << std::endl;
        std::cerr << "--spk: Set kind of speaker: Female (default), Male, Child" << std::endl;
        return -1;
    }
    
    const char *out_file = getCmdOption(argc, argv, "--out");
    
    int samplerate = SAMPLERATE;
    const char *sr_option = getCmdOption(argc, argv, "--sr");
    if(sr_option)
        samplerate = atoi(sr_option);
    
    int oversampling = OVERSAMPLING;
    const char *ovs_option = getCmdOption(argc, argv, "--ovs");
    if(ovs_option)
        oversampling = atoi(ovs_option);
    
    int blocklen = BLOCKLEN;
    const char *blk_option = getCmdOption(argc, argv, "--blk");
    if(blk_option)
        blocklen = atoi(blk_option);
    
    float duration = DURATION;
    const char *dur_option = getCmdOption(argc, argv, "--dur");
    if(dur_option)
        duration = atof(dur_option);
    
    const char32 *spk_choice;
    const char *spk_option = getCmdOption(argc, argv, "--spk");
    if(!spk_option || !strcmp(spk_option, "Female"))
        spk_choice = U"Female";
    else if(!strcmp(spk_option, "Male"))
        spk_choice = U"Male";
    else if(!strcmp(spk_option, "Child"))
        spk_choice = U"Child";
    else
        throw std::runtime_error("--spk option unknown, must be Male, Female, or Child.");
    
    autoSpeaker speaker = Speaker_create(spk_choice, MASSES);

    autoArtword artword = Artword_create(duration);
    
    //  muscle movement options
    for(int mi = (int)kArt_muscle::MIN; mi <= (int)kArt_muscle::MAX; ++mi) {
        std::string muscle("--");
        muscle += muscle_names[mi];
        const char *muscle_option = getCmdOption(argc, argv, muscle.c_str());
        if(muscle_option) {
            typedef std::vector<Pair<sample_t> > PairVec;
            PairVec *muscle = getPairList<sample_t>(muscle_option);
            // set Artword muscle with parameters
            for(PairVec::const_iterator it = muscle->begin(); it != muscle->end(); ++it)
                Artword_setTarget(artword.get(), kArt_muscle(mi), it->first, it->second);
            delete muscle;
        }
    }

    autoSound res = Artword_Speaker_to_Sound(artword.get(), speaker.get(), samplerate, oversampling,
                                             NULL, 0, NULL, 0, NULL, 0,
                                             NULL, 0, NULL, 0, NULL, 0,
                                             NULL, 0, NULL, 0, NULL, 0);
    
    if(out_file) {
        SF_INFO sfi;
        sfi.samplerate = samplerate;
        sfi.channels = 1;
        sfi.format = SF_FORMAT_WAV|SF_FORMAT_PCM_24;
        SNDFILE *sf = sf_open(out_file, SFM_WRITE, &sfi);
        sf_write(sf, res->z[1], res->nx);
    }
    
    return 0;
}
