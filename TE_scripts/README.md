### README

#### !!! IMPORTANT !!!
All scripts assume a RepeatModeler/RepeatMasker formatting style. Es : Lapu_rnd1-family4_64#DNA

---

TE scripts usefull in TE annotation and curation

 1. BEE_v2.0.py : Python script to produce a set of alignments with new raw consensus sequences that MUST be manually curated (FOR MANUAL CURATION ONLY)
 2. AutomaticBEE.py : Python script to produce a set of alignments and automatically extended consensus sequences (FOR AUTOMATIC EXTENSION ONLY)
 3. TE_AnnoSum.py : Python script usefull to calculate some summary statistics about a TE library and genome annotation. Usefull to compare the quality of different libraries.
 4. Merge_TELibraries.py : Merge multiple consensus libraries with redundant entries (*i.e* sequences with the same identifier before the classification) in order to keep only the "best one" (i.e. the classified and/or the longest one). It is possibile to provide a manually curated library that will be prioritize (*i.e.* all its sequences will be included in the final library).
 
#### Possible workflow

**a**. Identify raw consensus sequences that must be prioritized in the manual curation process.  
**b.** Run the BEE_v2.0.py script.  
**c.** Run on the same raw consensus the AutomaticBEE.py script.  
**d.** Annotate the genome indipendenlty with all three libraries (Manually curated, aumatically curated, raw).  
**e.** Compare the three libraries and annotations with TE_AnnoSum.py script.  
**f.** If necessary try to change Trimal parameters in AutomaticBEE.py script or to perform additional rounds of automatic curation.  
**g.** When results of automatic curation are satisfactory (*e.g* similar to manual curation results), use the same parameters to automatically improve all consensus sequences not included in manual curation.  
**e.** Classify automatically refinied conensus using *e.g* RepeatClassifier from the RepeatMasker package.
**e.** Combine manually and automatically curated libraries using Merge_TELibraries.py  
**f.** Remove redundancy following 80 - 80 rule.  
**g.** Compare the new annotations with the one obtained with the raw library.  

#### To Do :
 - Automatic removal of low complexity repeats and host genes from raw consensus sequences
 - Script to help identify raw consensus sequences that must be prioritized in the manual curation process (e.g longest one, with more hits, with proteins similarities).
 - Identify autonomous insertions 
