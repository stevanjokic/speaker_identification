# Speaker Identification
### MFCC based identification with additional features

After building Prepoznavanje.cpp into Prepoznavanje.exe, Prepoznavanje.exe is run from the command prompt.
Prepoznavanje.cpp is used in experimental purposes, therefore command line have a lot of parameters (arguments), this version have 9 parameters. Command line for Prepoznavanje.exe:
Prepoznavanje -working_mode- -file_for_training- -name_of_model- -file_for_testing- -govornici.txt- -threshold- -first_additional_feature- -second_additional_feature- -third_additional_feature-

Files for training and testing should be in wav format, 16 bit resolution, frequency sampling of fs=44100 Hz (also in code can be defined fs=22050 Hz).

Feature vector contains 21 MFCCs and additional features.
If additional feature is active then this field in command line should be filled with 1, otherwise it should be filled with 0.

Example of command line for making model named frf03 for training file frf03_s08_solo.wav, feature vector of 21 MFCCs and first second and third additional feature:
Prepoznavanje 3 frf03_s08_solo.wav frf03 xy.wav govornici.txt 1.245 1 1 1
(in this example it is not of interest the name of testing file, also the value of threshold is only of interest if recognizer is tested for speaker verification). Created model is written in file modeli.txt and name of model is written in file govornici.txt.

If we have created models of speakers by using command line for each training wave file as in previous example (for feature vector of 21 MFCCs and all three additional features for example), models are written in modeli.txt, then if we want to test speaker recognizer in speaker identification for test wav file, for example frf01_s07_solo.wav, command line should be:
Prepoznavanje xy.wav xy frf01_s07_solo.wav govornici.txt 1.245 1 1 1

Training on speech database can be done for working_mode=11, if we choose feature vector of 21 MFCCs and two additional features then command line should be: Prepoznavanje 11 xy.wav xy xy.wav govornici.txt 1.245 1 1 0

Testing of speaker identification on speech database can be done for working_mode=6, if we choose feature vector of 21 MFCCs and two additional features then command line should be: Prepoznavanje 6 xy.wav xy xy.wav govornici.txt 1.245 1 1 0

Training for working_mode=11 and testing for working_mode=6 are set in Prepoznavanje.cpp for speech database recorded in the project Speaker/style adaptation for digital voice assistants based on image processing methods (acronym: S-ADAPT), with little modifications in code these working modes can be used for training and testing on speech database CHAracterizing INdividual Speakers (CHAINS).
