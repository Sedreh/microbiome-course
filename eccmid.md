# ESGMD Postgraduate Technical Workshop- metagenomics

## 1. Enhancing Pathogen Detection in Mono and Polymicrobial Samples with 16S Amplicon-based Metagenomics


After physicians confirm the presence of an infection, it is typically necessary to determine the specific microorganism responsible for causing the infection. Various microorganisms have the potential to cause a particular infection and accurate identification of pathogens is crucial for effective diagnosis and treatment[^1]. While standard diagnostic methods like quantitative PCR (qPCR) or broad-range eubacterial PCR followed by Sanger (bPCR) sequencing are sensitive, reliable, and fast, they are limited in their ability to determin the specific composition of polymicrobial infections[^2].


## Case study

A hospitalized patient with symptoms such as fever, indicating a potential infection, for example related to periodontal disease. The initial diagnostic approach will be involved performing routine diagnostic tests, including bPCR. However, in the presence of several bacteria harboring a different target gene sequence, overlapping peaks forming mixed sequences and complex chromatograms will complicate the interpretation of results and the distinction of present microorganisms in the sample[^3]

> **Figure 1: Broad-range PCR + Sanger sequencing**
>![assembly graph](ESGMD_tutorial_doc/docs/img/PCR.png){: width="560px" width:"100%" .center}


Recognizing the limitations of conventional tests and the need for a complementary method to identify all pathogens present in the sample, one potential approach is to employ high-throughput sequencing techniques, enabling recovering all bacterial DNA directly from clinical specimens[^4]. 


### 1.1. *In silico* validation of 16S taxonomy assignment

16s amplicon-based metagenomics is a method focused on sequencing specific target regions of *16S rRNA* genes (amplicons). Since *16S rRNA* genes are ubiquitous in prokaryotes, hence in bacteria, they can be targeted to identify the members of a bacterial community. Indeed, each  16S rRNA gene has an approximate lenght of 1600 base pairs and includes nine hypervariable regions of varying conservation (V1-V9) interspersed by conserved sequences. The variable regions are taxa-specific and they provide sufficient sequence diversity to differentiate between microbial *Genera*, and often even at the *Species* level[^5]. Thus, through the construction of primers shaped on the constant sequences, only variable regions of interest can be amplified, in our case V3 and V4 regions. This strategy rather than using the full length of *16S rRNA* gene maximizes cost efficiency and minimizes technical challenges associated with amplification and analysis. 

The taxonomical assignment is a challenging process. Indeed, the close similarity of the 16s gene of two closely related species make  difficult and sometimes impossibile, describe an isolate at the *Species* level. Moreover, from 1 to 7 copies copies of the *16S rRNA* gene are usually present in a single genome (4.2 copies on average) and their classification can be discordant[^6]. Therefore, the selection of a reference database for taxonomical annotation is a key factor that impact the quality of the prediction and the taxonomical resolution[^7].

> **Figure 2: 16S amplicon-based metagenomics workflow for pathogen detection**
>![assembly graph](ESGMD_tutorial_doc/docs/img/16SWorkflow.png){: width="560px" width:"100%" .center}



To evaluate the ability of our approach in correclty classifying the species of interest we can run an *in silico* analysis, that consist of the following steps:

*provide the bacteria genome of known species from which the sequence of V3-V4 region(s) is extracted
*classify the identified region(s) using a reference database
*evaluate if the obtained classification corresponds to what expected
   
??? info "Precomputed results"
Here you can find the results of an  *in silico* taxonomy assignment using two widely used 16S databases, Silva and EzBioCloud:
*Pseudomonas aeruginosa* genomes
*Staphylococcus aureus* genomes
*Salmonella enterica* genomes   
- [**EzBioCloud**](ESGMD_tutorial_doc/docs/img/InSilico_results_EzBioCLoud.html)   
- [**Silva**](ESGMD_tutorial_doc/docs/img/InSilico_results_Silva.html)


Questions: 
!!! question "Question 1.1"
        How many *16S rRNA* gene copies (V3-V4 regions) are identified in each genome?
!!! question "Question 1.2"
        Are all the copies classified to the same taxa in a single genome? If, not which species have discordant predictions?
!!! question "Question 1.3"
        Does the classification reach the *Species* level? Why? (to do for the why questions: allign V3-V4 of the two species wrongly classified)

??? info "Answer"
-The number of  16S rRNA genes detected per genome varies from 4 to 7.
-While all V3-V4 regions of Pseudomonas aeruginosa and Staphylococcus aureus genomes are grouped in the same amplicon variant, the V3-V4 regions of Salmonella enterica genomes, are assigned to 2 and 6 amplicon variants, respectively.
- Pseudomonas aeruginosa and Staphylococcus aureus genomes are correctly classified at species level, and even the classification at the strain level is provided. Instead, Salmonella enterica genomes are wrongly misclassified as Enterobacter cloacae, thus not only the species, but also the Genus is wrongly assigned.
-The reason of such misclassification can be identified in the high similarity of 16S rRNA V3-V4 regions of Salmonella enterica and Enterobacter cloacae [alignment](https://pages.github.com/)(see alignment_](https://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=clustalo-I20230608-134158-0534-15996657-p1m)

!!! question "Question 1.4" 
      Is the classification different? Link to table


??? info "Answer"
Silva database, similarly to EzbioCLoud, correctly classifies Pseudomonas aeruginosa and Staphylococcus aureus genomes. However, a different classification is achive for Salmonella enterica genomes. The two amplicon variants of the first one are both assigned to the correct Genus, Salmonella, but no assignment at Species level is provided. While among the 6 amplicon variants of the second genome, 3 are assigned to Samonella genus and no species classification, but three others are wrongly assigned to  Enterobacter cloacae.


### 1.2. *In vitro* validation of 16S with MOCK community

Microbial mock community is a defined mixture of microbial nucleic acids or cells that are commonly used in micbiota study as a positive control or reference for evaluating the performance and accuracy of experimental and bioinformatic procedures. 

For example, the standard ATCC-2002 mock community consist of 20 bacterial species evenly distributed, mimicking polymicrobial metagenomic samples.


> **Figure 3: The 20 bacterial taxa of ATCC mock community**
>![assembly graph](ESGMD_tutorial_doc/docs/img/atcc.png){: width="560px" width:"100%" .center}

1) The mock community can be used as a positive control for sequencing runs. The composition of one positive control is represented as a KRONA plot, an interactive pie chart.

[Click here](https://rpubs.com/ychoi/mock_krona)

!!! question "Question 1.5"
      Are the expected 20 bacterial members of the MOCK community well identified?
      
??? info "Answer"
No, not every species.

!!! question "Question 1.6"      
Is there any unexpected bacterial species?

??? info "Answer"
Yes, *Cutibacterium acnes*.

!!! question "Question 1.7"       
      Are the relative abundances as expected? Why or why not?
      
??? info "Answer"
No, while the relative abundance is expected to be 5 % for each species, some species are more abundant than others.

2) The results of 16S sequencing can vary according to the targeted variable region.


> **Figure 4: The relative abundance of mock community members based on the variable region (modified from [^8])**
>![assembly graph](ESGMD_tutorial_doc/docs/img/atcc2.png){: width="560px" width:"100%" .center}


!!! question "Question 1.8"       
      How much % of expected species did they detect? Different from our expectation? Why?

??? info "Answer"
75%. We expect the result to be 100% but we are missing some taxa. It might be due to the low copy number or harsh extraction method for some species.

!!! question "Question 1.9"       
      59%. They might have been introduced by contamination. 
      
### 1.3. Diagnostic performance of 16S rRNA amplicon-based metagenomics compared to gold standards (bPCR)

To compare the performance of targeted metagenomics with broad-range PCR followed by Sanger sequencing, several bPCR positive samples were selected for retrospective analysis. These samples originate from usually sterile body sites including peritoneal fluid, bone fragments, cervical abscesses, biopsie. After sequencing V3-V4 regions of 16S rRNA gene, we processed raw reads using [zAMP](https://pages.github.com/) bioinformatics pipeline and Amplicon Sequence Variants (ASVs) were classified by the [RDP classifier](https://sourceforge.net/projects/rdp-classifier/) using the [EzBioCloud reference database (release 2018.05)](https://www.ezbiocloud.net/). To evaluate the event of
contaminations along the workflow, several no-template controls have been included.

#### 1.3.1. Single pathogen PCR positives (mono-microbial samples)

> **Figure 5: Broad range PCR positive samples with single pathogen)**
>![assembly graph](ESGMD_tutorial_doc/docs/img/EPCR_mono.png){: width="560px" width:"100%" .center}


!!! question "Question 1.10"       
      Did targeted metagenomics yield comparable accuracy to the pathogen identified by bPCR?

??? info "Answer"
      Yes, same expected pathogen identified in the PCR positive samples.

!!! question "Question 1.11"       
      How does the abundance of the identified pathogen using targeted metagenomics affect its accuracy? 

??? info "Answer"
      If targeted metagenomics identifies the same pathogen as bPCR but in lower abundance, it can still be considered a valid detection. However, in such cases, it is important to interpret the results with caution and it would be beneficial to validate the findings using additional methods such as expert evaluation. 


#### 1.3.2. Single pathogen PCR positives with additional taxa identified by targeted metagenomics (mono- or polymicrobial samples?)

> **Figure 6: Mono microbial PCR positive samples with multiple pathogen detected by targeted metagenomics)**
>![assembly graph](ESGMD_tutorial_doc/docs/img/ePCR_strange.png){: width="560px" width:"100%" .center}


!!! question "Question 1.12"       
      How do you interpret the results when multiple pathogens are identified in samples that were supposed to be mono-microbial?

??? info "Answer"
When multiple pathogens are identified in samples that were expected to be mono-microbial, there would be several possibilities such as:
- Mixed infections: Mixed infections involve the presence of multiple pathogens in a single sample, with each pathogen contributing to the disease leading to complex clinical symptoms.
- Contamination: The presence of multiple pathogens could be due to contamination during sample collection, handling, or laboratory procedures. For example sample 1015435733 is Mandibular bone biopsyand all of the detected bacteria are oral microbiota microorganisms (contamination from sample collection).


!!! question "Bonus question1"       
      What is cross-contamination and how it can happen in this context?
      
!!! question "Bonus question2"       
      How do you interpret the results in terms of the abundance of the identified bacteria?

#### 1.3.3. Suspected poly-microbial samples

> **Figure 7: poly-microbial samples. Suspected poly-microbial PCR positive samples with single or no pathogen reported**
>![assembly graph](ESGMD_tutorial_doc/docs/img/polymicrobial.png){: width="560px" width:"100%" .center}

!!! question "Question 1.13"       
      Is targeted metagenomics more sensitive than bPCR? why?

??? info "Answer"
- Broad detection capability: Targeted metagenomics can detect a wide range of pathogens simultaneously, as it involves sequencing and analyzing the entire DNA present in the sample. Based on the results, this approach not only recovered expected pathogens (reported by bPCR) but also identified additional species in suspected polymicrobial samples.  
- Identification of novel or emerging pathogens: It allows for the detection of unknown or unexpected pathogens, enhancing sensitivity in situations where the causative agent is not well characterized.
- Quantification of abundance: It provides information about the relative abundance of different pathogens within a sample. This quantitative aspect enables the assessment of the pathogen's presence even at low abundance.

However the sensitivity of targeted metagenomics can still vary depending on factors such as sequencing depth, data analysis methods, and the specific assay design. 




## References
[^1]: 
    Kabiraz, Meera & Majumdar, Priyanka & Mahmud, Chayan & Bhowmik, Shuva & Ali, Mohammad. (2023). Conventional and advanced detection techniques of foodborne pathogens: A comprehensive review. Heliyon. 9. 10.1016/j.heliyon.2023.e15482.

[^2]: 
    O. Kommedal, K. Kvello, R. Skjåstad, N. Langeland, and H. Wiker. Direct 16s rrna gene sequencing from clinical specimens, with special focus on polybacterial samples and interpretation of mixed dna chromatograms. Journal of clinical microbiology, 47:3562–8, 10 2009. doi: 10.1128/JCM.00973-09.
[^3]: 
    K. Mongkolrattanothai and J. Dien Bard. The utility of direct specimen detection by sanger sequencing in hospitalized pediatric patients. Diagnostic Microbiology and Infectious Disease, 87(2):100–102, 2017. ISSN 0732-8893. doi: https://doi.org/10.1016/j.diagmicrobio. 2016.10.024.
[^4]: 
    E. E. Hilt and P. Ferrieri. Next generation and other sequencing technologies in diagnostic microbiology and infectious diseases. Genes, 13(9):1566, Aug 2022. ISSN 2073-4425. doi: 10.3390/genes13091566.
[^5]: 
    Kamble, Asmita & Sawant, Shriya & Singh, Harinder. (2020). 16S ribosomal RNA gene-based metagenomics: A review. Biomedical Research Journal. 7. 5. 10.4103/BMRJ.BMRJ_4_20.
[^6]: 
    Vetrovsky, Tomas & Baldrian, Petr. (2013). The Variability of the 16S rRNA Gene in Bacterial Genomes and Its Consequences for Bacterial Community Analyses. PloS one. 8. e57923. 10.1371/journal.pone.0057923.
[^7]:
    Park, Sang-Cheol & Won, Sungho. (2018). Evaluation of 16S rRNA Databases for Taxonomic Assignments Using Mock Community. Genomics & Informatics. 16. e24. 10.5808/GI.2018.16.4.e24.
[^8]:
    Fouhy, Fiona & Clooney, Adam & Stanton, Catherine & Claesson, Marcus & Cotter, Paul. (2016). 16S rRNA gene sequencing of mock microbial populations- impact of DNA extraction method, primer choice and sequencing platform. BMC Microbiology. 16. 10.1186/s12866-016-0738-z.
    
