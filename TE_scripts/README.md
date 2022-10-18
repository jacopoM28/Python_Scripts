### README

#### !!! IMPORTANT !!!
All scripts assume a RepeatModeler/RepeatMasker formatting style. Es : Lapu_rnd1-family4_64#DNA

---

TE scripts useful in TE annotation and curation

 1. **BEE_v3.0.py** : Python script to produce a set of alignments with new raw consensus sequences that MUST be manually curated (FOR MANUAL CURATION ONLY)
 2. **AutomaticBEE_v2.0.py** : Python script to produce a set of alignments and automatically extended consensus sequences (FOR AUTOMATIC EXTENSION ONLY)
 3. **TE_AnnoSum.py** : Python script useful to calculate some summary statistics about a TE library and genome annotation (*e.g* to compare the quality of different libraries).
 4. **Merge_TELibraries.py** : Merge multiple consensus libraries with redundant entries (*i.e* sequences with the same identifier before the classification) in order to keep only the "best one" (i.e. the classified and/or the longest one). It is possibile to provide a manually curated library that will be prioritize (*i.e.* all its sequences will be included in the final library).
 5. **BlastX_TE_Anno.py** : Annotation of consensus sequences based on Blastx of the 20 longer instertions for each one. If all significant hits are against the same TE (*e.g.* DNA/hAT; LINE/R2; DNA/TcMar-Trigger) the annotation is transferred to the element. If not the script will check the class level.
 
#### Possible workflow

**a**. Identify raw consensus sequences that must be prioritized in the manual curation process.  
**b.** Run the BEE_v3.0.py script.  
**d.** Run on the same raw consensus the AutomaticBEE_v2.0.py script.  
**e.** Annotate the genome indipendenlty with all three libraries (Manually curated, aumatically curated, raw).  
**f.** Compare the three libraries and annotations with TE_AnnoSum.py script. This script could also help to check the manually curated library and eventually add anotations by manually screening the Blastx results.  
**g.** Classify automatically refinied conensus using *e.g* RepeatClassifier from the RepeatMasker package.  
**h.** Combine manually and automatically curated libraries using Merge_TELibraries.py.  
**i.** Classify still unknown sequences using both blastx on N longest insertions (BlastX_TE_Anno.py) and similarity to known consensus (*e.g* previously classified; To Do).  
**l.** Remove redundancy following 80 - 80 rule.  
**m.** Compare the new annotations with the one obtained with the raw library.  

#### To Do :
 - Automatic removal of low complexity repeats and host genes from raw consensus sequences.
 - Script to help identify raw consensus sequences that must be prioritized in the manual curation process (e.g longest one, with more hits, with proteins similarities).
 - Identify autonomous insertions.  
