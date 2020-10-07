######## FPCA-based resynthesis - Example sound files

- Background: please refer to presentations/fda_resynthesis.pdf

- Filename structure:

S1, S2, S3 = narrow focus on 1st, 2nd, 3rd lexical stress, i.e. on Subject,
Verb or Prop. Phrase (cf. presentations/FPCA-phonetics.pdf)
(Note: capital "S", not to be confused with s1, s2, which are the PC scores)

A1A, A1T = speakers (male, female)

S1_A1A.wav and similar = original recordings

S1_A1A_XXX....wav = PSOLA resynthesised sounds based on S1_A1A recording

FPCA-based manipulation is based on the first 2 PCs (s1, s2) from FPCA applied to pairs of (f0, log rate) curves on 21 speakers.

S1_A1A_null_f0_logvel.wav = resynth based on the original PC scores, i.e. no
transformation, only FPCA-based reconstruction.

S1_A1A_to_midP_f0_logvel.wav = resynth using (s1, s2) coordinates of the
mid point (cf. slides)

S1_A1A_to_S2_A1A_f0_logvel.wav = resynth using (s1, s2) of condition S2, i.e.
transform S1 into S2.

- Notes: The sentence is always "Danilo vola da Roma". For some reason I don't
know, most of the files are cut off at the end.

What was manipulated was f0 and log rate, i.e. segment durations. Probably
manipulating intensity would improve the quality, as the vowel carrying a narrow
focus is also louder. Clearly there's also something missing in the vowel
quality, as vowels not carrying any accent tend to be centralised
(schwa-like). When one such centralised vowel is artificially 'promoted' to carry the accent, its quality does not change, and the result is a bit unnatural. 

