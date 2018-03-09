#include "Artword_Speaker_to_Sound.h"
#include "Speaker.h"

#define DURATION 1.0
#define SPEAKER U"Male"
#define MASSES 2  // 1, 2, 10
#define FSAMP 22050
#define OVERSAMPLING 25

int main() 
{
    autoArtword artword = Artword_create(DURATION);
    autoSpeaker speaker = Speaker_create(SPEAKER, MASSES);
    autoSound res = Artword_Speaker_to_Sound(artword.get(), speaker.get(), FSAMP, OVERSAMPLING,
	                    NULL, 0, NULL, 0, NULL, 0,
	                    NULL, 0, NULL, 0, NULL, 0,
                        NULL, 0, NULL, 0, NULL, 0);
    
    MelderFile fp;
    str32cpy(fp->path, U"/tmp/out.aiff");
    MelderFile file = MelderFile_create(fp); 
    Sound_saveAsAudioFile(res.get(), file, Melder_AIFF, 24);
    return 0;
}
