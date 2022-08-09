#define framesize 320
#define MaxFrameCount 400
#define Ampscale 5000
#define SpeechLength 600000
#define IgnoreSamples 80
#define samplingrate 16000
#define InitFrames 15
#define testw "test.wav" //store the signal in wave form in the file mentioned here
#define testt "test.txt" //store the signal in text form in the file mentioned here
#define duration 3 //time given for recording in seconds
#define recmod "C:\\RecordingModule.exe" //path of recording module provided
#define path "C:\\universe\\" //path for universe of vowel sound collection
#define p 12
#define precision 10 //precision upto these many number of decimal places
#define delim ' ' //character which you used to sepstral coefficients while storing in file
#define mapdelim '\t' //delimiter used to store the mapping from word(vowel) to vector end indices
#define universe "cep.txt" //universe name. oldname--->filename
#define codebook "codebook.txt"
#define cbsize 32 // code book size
#define cepsize 12 //number of cepstral coefficients
#define wordcount 5 //number of words which are used for estimation
#define thresholdist 0.1 //threshold distortion in percentage w.r.t previous distortion
#define N 5 //number of states in system
#define M 32 //number of dustinct observation symbols
#define ModelIterate 25 //number of re-estimations preferred