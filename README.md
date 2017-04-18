[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-ff69b4.svg)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

# MASTmelody_dataset
A dataset of pitch curves for singing assessment

[IMPORTANT NOTICE: a publication has been submitted announcing this dataset to The 13th International Symposium on Computer Music Multidisciplinary Research (CMMR). The dataset will be partially available(10%) until the publication is accepted]

<a name="Introduction"></a>Introduction
--------------------
The MASTmelody dataset is designed and shared to facilitate comparison of algorithms in the field of automatic music performance assessment.

The dataset includes pitch(f0) data extracted from audio data recorded during conservatory entrance examinations. Audio data could not be included due to the difficulties involved in completely anonymizing audio files (recognition of the singer is possible via listening). Only a few samples are provided in the 'wavSamples' directory.   

There are two broad categories for the data: f0 data of the reference recording files (melodies played on the piano as reference) and f0 data of the performance recording files (recordings of the candidate in singing the melodies)  

Each recording was subject to f0-detection using <a href="https://github.com/sertansenturk/predominantmelodymakam"> a variant </a> of <a href="http://mtg.upf.edu/technologies/melodia"> Melodia Melody Extraction tool</a>. The results were saved as text files containing two columns: time-stamp and estimated f0 information in Hz.  

The performances have been graded by three jury members, who are teaching staff members of the conservatory. Grades are binary: pass, fail. This dataset includes only the samples for which all jury members agreed in grading with the same score. Hence, there are basically two categories for the performance files: i)performances which was graded as 'fail' by all the jury members , i)performances which was graded as 'pass' by all the jury members.

Naming convention:
The dataset is just composed of a list of text files. All other information is coded in the file names:
'ref': reference recording on the piano
'per': performance recording
'fail': performance graded as 'fail'
'pass': performance graded as 'pass'

There are basically 40 distinct melodies performed. The id for the melody makes up the first part of the file name. Examples:
'51_mel1_per101559_fail.f0s.txt': Melody with ID: '51_mel1' and this is a performance file graded as fail
'55_mel2_ref280758.f0s.txt': Melody with ID: '55_mel2' and this is a reference file

Please cite the following publication if you use this dataset in your work:
> Bozkurt, Baysal, Yuret. Reference goes here if paper accepted, 2017.

<a name="License"></a>License
--------------------
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

<a name="Acknowledgement"></a>Acknowledgement
--------------------
This dataset has been curated within the TUBITAK(The Scientific and Technological Research Council of Turkey) funded research project 1001-215K017 targeting development of automatic assessment tools for music performances.
