############# EMA dataset

- Background:  please refer to the paper (file GubianPastaetterPouplier2019.pdf)

Gubian, M., Pastätter, M., & Pouplier, M. (2019). Zooming in on Spatiotemporal V-to-C Coarticulation with Functional PCA. In INTERSPEECH (pp. 889-893).

- Dataset: EMA.csv

The table contains 73 EMA trajectories from one sensor and relevant meta-data. The
sampling rate is 250 Hz.
The table columns are:

ID, idx		unique identifiers for trajectory (either is sufficient)
CONDITION	ka, ki
SPEAKER		speaker id
ACC		accentuation condition (not used in the paper)
REPETITION	repetition from the same speaker (not used in the paper)
TRIAL		(not used in the paper)
INTERVAL	interval i is delimited by t_{i-1} and t_i in Fig. 1 in the
paper
ONSET, OFFSET	sample indices delimiting an interval (not starting from zero)
INTSAMPLE	sample index within an interval
TOTSAMPLE	sample index within a trajectory (utterance)
POSX, POSY	position of the sensor in the sagittal plane (mm)


- Suggested project ideas:

1) [advanced] Reproduce the main results in the paper using the same
techniques, i.e. landmark reg. + FPCA + LMER.

To run the analysis you can use exLand.R as template, with the appropriate
modification and data prep. (e.g. you have to construct the land table
yourselves from the information in EMA.csv).
Please follow the instructions in Sec. 2.2 in the paper to carry out all
necessary normalisations and zero-offsetting (e.g. POSX,Y are neither
speaker-normalised nor zero-offsetted). z-score normalisation means
subtracting the mean and dividing by the standard deviation.

2) [advanced] Apply GAMMs on these data

The suggested approach is the one used in this paper:
https://research-information.bris.ac.uk/en/publications/investigating-dialectal-differences-using-articulography
which corresponds to the one denoted as 'trick' in the NZE project. 
Note that this paper uses linear time normalisation. You may want to apply
landmark registration instead (but without using the log rate curves further
in the analysis).
