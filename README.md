# HMMforSpeech
HMM based Speech Recognition


All the following files and graphs are in "HMM/HMM/" directory

config.h---> configuration file to set various things such as delimiter, filename, thresold distortion allowed, precision of the code vectors, etc.

hmm---> implements hmm on given codebook ("codebook.txt") and desired test speech ("test.txt").

set the values in config.h and run the project.

"codebook.txt" was built for vowels, so this hmm model is used for vowel detection

Then immediately it will ask "Do you want to record a speech?(y/n): ", 

case 1: if you already have a pre-recorded(offline) test signal then type 'n'(without quotes),

	then it will try to apply the hmm on the speech in file "test.txt" file after displaying the message 

	"we will try to recognize the speech in test.txt file".

case 2: if you want to record a speech, then press 'y'(without quotes) and record your speech when given time is over,

	then it displays the message "your speech has been succesfully recorded"

Then it will ask "Did you remove the non speech part of the signal?(y/n): "

case 1: if you had removed non speech part by yourself then press 'y'(without quotes)

case 2: if you had not removed non speech part by yourself then press 'n'(without quotes)

	then the program automatically tries to remove non-speech part of the signal, but this is not as accurate as
	
	you can remove manually

Now program executes, and re-estimates the model lambda for the "ModelIterate"(as mentioned in "config.h" file) number of times.


Output of final model lambda(pi, A, B), probability of an utterance being from the model, and The estimated state sequence

will be stored in "output.txt"

Intermediate processing results are stored in "log.txt" file

Note:
Extra method uniform_model() is provided in the program which you may use based on your interest


P.S:

For sample log files which are generated at the time of development check "log.txt", "log1.txt", and "log2.txt"

"cep.txt"---> universe used for generating codebook

"test.txt"---> test speech in .txt format

"test.wav"---> test speech in .wav format

we also included various other files which helped in generating the codebook.
