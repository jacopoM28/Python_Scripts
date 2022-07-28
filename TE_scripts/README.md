### README

TE scripts usefull in TE annotation and curation

 1. BEE_v2.0.py : Python script to produce a set of alignments with new raw consensus sequences tha can be manually curated (FOR MANUAL CURATION ONLY)
 2. AutomaticBEE.py : Python script to produce a set of alignments and automatically extended consensus sequences (FOR AUTOMATIC EXTENSION ONLY)
 3. TE_AnnoSum.py : Python script usefull to calculate some summary statistics about a TE library and genome annotation. Usefull to compare annotation coming from different libraries.
 
#### Suggested workflow

**a**. Identify raw consensus sequences that must be prioritized in the manual curation process.  
**b.** Run the BEE_v2.0.py script on that.  
**c.** Run on the same consensus the AutomaticBEE.py script.  
**d.** Annotate the genome indipendenlty with all three libraries (Manually curate, aumatically curated, raw consensus).  
**e.** Compare the three libraries and annotations with TE_AnnoSum.py script.  
**f.** If necessary try to change Trimal parameters in AutomaticBEE.py script or perform additional round of automatic curation.  
**g.** When results of automatic curation are satisfactory (*e.g* similar to manual curation results), use the same parameters to automatically improve all consensus sequences not included in manual curation.  
**e.** Classify automatically refinied conensus using *e.g* RepeatClassifier from the RepeatMasker package.  
**e.** Combine manually and automatically curated libraries, remove redundancy and perform genome annotation of repeats.  
**f.** Compare the new annotations with the one obtained with the raw library.  

#### To Do :
 - Script to help identify raw consensus sequences that must be prioritized in the manual curation process (e.g longest one, with more hits, with proteins similarities)
