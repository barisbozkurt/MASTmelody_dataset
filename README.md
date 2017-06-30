[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-ff69b4.svg)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

# MASTmelody_dataset

This repository contains a dataset of f0-series data computed from signing and piano samples and <a href="https://github.com/barisbozkurt/MASTmelody_dataset/tree/master/baseline">a baseline system for assessment</a> (comparing a singing sample to a reference piano recording), both of which are described in the scientific paper:
> Bozkurt, B., Baysal, O., Yuret, D. A Dataset and Baseline System for Singing Voice Assessment, 13th Int. Symposium on Computer Music Multidisciplinary Research, Porto, Sept. 25-28, 2017.
```
@inproceedings{inproceedings,
  author={Bozkurt, B., Baysal, O., Yuret, D.},
  title={A Dataset and Baseline System for Singing Voice Assessment},
  year={2017},
  booktitle={13th Int. Symposium on Computer Music Multidisciplinary Research, CMMR 2017}
}
```

Please cite the publication if you use this dataset and/or the baseline system in your work.

<a name="Introduction"></a>Introduction
--------------------
The MASTmelody dataset is designed and shared to facilitate comparison of algorithms in the field of automatic music performance assessment.

The dataset includes pitch(f0) data extracted from audio data recorded during conservatory entrance examinations. Audio data could not be included due to the difficulties involved in completely anonymizing audio files (recognition of the singer is possible via listening). Only a few samples are provided in the 'wavSamples' directory.   

There are two broad categories for the data: f0 data of the reference recording files (melodies played on the piano as reference) and f0 data of the performance recording files (recordings of the candidate in singing the melodies)  

Each recording was subject to f0-detection using <a href="https://github.com/sertansenturk/predominantmelodymakam"> a variant </a> of <a href="http://mtg.upf.edu/technologies/melodia"> Melodia Melody Extraction tool</a>. The results were saved as text files containing two columns: time-stamp and estimated f0 information in Hz.  

The performances have been graded by three jury members, who are teaching staff members of the conservatory. Grades are binary: pass, fail. This dataset includes only the samples for which all jury members agreed in grading with the same score. Hence, there are basically two categories for the performance files: i)performances which was graded as 'fail' by all the jury members , i)performances which was graded as 'pass' by all the jury members.

#### Naming convention:
The dataset is just composed of a list of text files. All other information is coded in the file names:
'ref': reference recording on the piano
'per': performance recording
'fail': performance graded as 'fail'
'pass': performance graded as 'pass'

There are basically 40 distinct melodies performed. The id for the melody makes up the first part of the file name. Examples:

'51_mel1_per101559_fail.f0s.txt': Melody with ID: '51_mel1' and this is a performance file graded as fail

'55_mel2_ref280758.f0s.txt': Melody with ID: '55_mel2' and this is a reference file

<a name="License"></a>License
--------------------
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

<a name="Acknowledgement"></a>Acknowledgement
--------------------
This dataset has been curated within the TUBITAK(The Scientific and Technological Research Council of Turkey) funded research project 1001-215K017 targeting development of automatic assessment tools for music performances.
