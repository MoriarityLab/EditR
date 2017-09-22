## Test suite of files

These files provide test cases to use to validate that EditR is working correctly


**Testing: how it handles a long tail of noise**

file: PD-1_2.ab1

guide: TGCAGATCCCACAGGCGCCC

guide matched 132--151

Editing is at position 8, with 40.9% percent area for T, false positive at 5 for A with 8.2%

Long tail of noise that doesn't get trimmed using our method of filtering bases that have a 10th of the average peak area

manually edited to 20 -- 433

**Testing: guide is a reverse complement**

file: 4d.ab1

guide: CGTGCATCAGATGCTTCACC

guide matched 319--338

Editing is at position 13 of the guide, 19.9 for A


**Just a regular file**

file: BE MAFB5.ab1

guide: AGCCGGCTGGCTGCAGGCGT


Editing at position 16, 46.2% for A

manually trimmed 20 to 318
