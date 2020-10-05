###### NZE Dataset

- Background: please refer to the paper (file Gubian_et_al2019_Tracking.pdf)

Gubian, M., Harrington, J., Stevens, M., Schiel, F., & Warren, P. (2019). Tracking the New Zealand English NEAR/SQUARE Merger Using Functional Principal Components Analysis. In INTERSPEECH (pp. 296-300).

- Dataset: NZE.csv

The table contains formant tracks and meta-data from 4222 diphthong tokens.
Each token was sampled at 11 equidistant time points. 
The table columns are:

tokenId		integer from 1 to 4222
Diphthong	I@, e@	
time_ms		time in ms starting from 0
Region		Hamilton or Wellington
Age		3 levels age groups (not used in the paper)
Age2		2 levels age groups (used in the paper)
Phrase		position in phrase, final or not
Sex		M, F
Speaker		speaker Id
Word		word containing the target diphthong
F1, F2		normalised values of the first two formants
 
- Notes:

The F1 and F2 columns were originally in semitones. The mean of each formant of each token
was subtracted, i.e. each F1 or F2 track is centered on zero. 

In the paper, the figures showing formants were obtained by adding/subtracting
0.5 to/from F2/F1. This was done for convenience of representation, in order
to fit F1 and F2 in the same plot while keeping them separated from each
other. Since the tracks were centered, they would in general cross each other,
making the plots harder to interpret. The constant shift was not used in
any of the models (of course). 

- Suggested project ideas:

1) [baseline] Reproduce the main results in the paper using the same
techniques, i.e. FPCA + LMER. 

The results reported in Sec. 2 can be reproduced by following the description
in the paper (please ignore Sec. 3, for which no data are provided) and using the script ex2D.R as template.

Question: Table 1 in the paper reports which random slope terms where included
in the full LMER model. The remaining terms, e.g. (Diphthong|Word), were not
included because they don't make sense. Why? 

2) [advanced] Include total duration in the FPCA + LMER model

Use the time warping representation technique to include duration in the
analysis. Each token becomes a 3D trajectory: F1, F2, log rate. As there are
no landmarks (boundaries) between track start and end, log rate is a constant line (and h(t) is a straight
line), which you can construct without using landmarkreg.nocurve().
You can use the script exLand.R as template (with appropriate modifications).

3) [baseline] Apply GAMMs separately on F1 and F2 and compare the results with the ones in the paper

Apply GAMMs following the steps indicated in the tutorial by M. Wieling
(excluding Sec. 4.9).
Construct two independent GAMMs, one on F1 and one on F2. 
You don't need to use all the fixed and random factors listed in Table 1 of
the NZE paper, rather limit to those that allow you to get to a solution
without convergence problems and/or too long computation (this applies to all
GAMM projects). You may also try and subsample the dataset in order to speed
up the estimation time (e.g. eliminate the MidAge class).

4) [advanced] Experiment applying GAMMs on F1 and F2 together (multi-response
GAMM)

There are two suggested approaches to this task

a) Using a 'trick'. Predict a generic formant track and add a factor to the
RHS of the model controlling for the formant. In other words, just
like you add a factor for Sex (M,F) or AgeGroup (Young,Old), you add
another one for Formant (F1,F2). This approach has been used e.g. in:
https://icphs2019.org/icphs2019-fullpapers/pdf/full-paper_56.pdf

Why is this a sort of trick or shortcut? What information is missing and how
could you recover it? (if possible}

b) The proper way. Use the example provided in the excerpt
notes/MultiResponseGAMM.pdf (Sec. 7.10). 
What are the limitations of the bam/gam command? (if you find any)
(cf. Sec. 4.9 in this paper:
https://www.frontiersin.org/articles/10.3389/fevo.2018.00149/full)
