# Speech-Recognition
Developing and evaluating a keyword speech recognizer.
- Here, we will developed and evaluated a keyword speech recognizer. 
- we Recorded a multiple instances of spoken words, to create a datbase of multiple wordd. 
- Implemented the MFCC algorithm for feature extraction and extract the following features:
1. 12 MFCC (mel frequency cepstral coefficients)
2. 1 energy feature
3. 12 delta MFCC features 
4.  12 double-delta MFCC features
5.  1 delta energy feature
6. 1 double-delta energy feature
- Process: Useed a frame size of 25 ms with a 10 ms overlap. This will result in a different number of frames for each word (since each will likely be of slightly different length). Uniformly sample the frames in the speech so that the same number of frames exist for each of the five words. The result will be a data set of N words. For each word you will have a total of 39*M features, where M is the number of frames in each word. 
