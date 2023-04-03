# IgPPI

Cells detect changes of external environments or communicate with each other through proteins on their surfaces. These cell surface proteins form a complicated network of interactions in order to fulfill their functions. The interactions between cell surface proteins are highly dynamic and thus are challenging to be detected by traditional experimental techniques. Here we developed a new computational framework to tackle this challenge. The primary focus of the framework is to identify interactions between domains in immunoglobulin (Ig) fold, which is the most abundant domain family in cell surface proteins. In practice, we collected all structural data of Ig domain interactions and transformed them into an interface fragment pair library. A high dimensional profile can be then constructed from the library for a given pair of query protein sequences. Multiple machine learning models were used to read this profile, so that the probability of interaction between the query proteins can be predicted. We tested our models to an experimentally derived dataset which contains 564 cell surface proteins in human. The cross-validation results show that we can achieve higher than 70% accuracy in identifying the PPIs within this dataset. We then applied this method to a group of 46 cell surface proteins in C elegans. We screened every possible interaction between these proteins. Many interactions recognized by our machine learning classifiers have been experimentally confirmed in the literatures. In conclusion, our computational platform serves a useful tool to help identifying potential new interactions between cell surface proteins in addition to current state-of-the-art experimental techniques. Moreover, the general framework of the machine learning classification can also be extended to study interactions of proteins in other domain superfamilies.


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
The details of the non-redundant structural database for Ig domain interactions are provided in the following list:
IgDDI_NR90_list.dat
The file contains the PDB and chain information of 831 pairs of domain-domain interactions.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
The benchmark datasets of an experimentally-derived human interactome of cell surface proteins are listed as follows:
HmSurfaceInteractome_Info.dat
This file provide the detailed information for all human cell surface proteins in the benchmark
HmSurfaceInteractome_PPI.dat
This file provide the detailed information for experimentally identified interactions among all human cell surface proteins in the benchmark
HmSurfaceInteractome_fasta.dat
This file provide the fasta sequences for all human cell surface proteins in the benchmark
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
The input datasets for our machine-learning-based prediction program are listed as follows:
HSIPPI_Profile_positiveset.dat gives the positive training data derived from the human cell surface proteins interactome. 
HSIPPI_Profile_negativeset.dat gives the negative training data derived from the human cell surface proteins interactome. 
CElegantsIgPPI_Profile_Final.dat contains the profiles of ranked alignment scores (PRAS) for all tested 46 proteins in C. elegans.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
The source codes that uses the multiple machine-learning algorithms to predict the interactions of cell surface proteins in C. elegans are in the Fortran77 format. They are listed as follows:
The source codes in CElegantsIgPPI_MLtest_main.f is the main program for machine-learning-based prediction using the datasets listed above as inputs.
The source codes in sub_BP_learn.f is the subroutine for the back-propagation neural network training process.
The source codes in sub_BP_recall.f is the subroutine for the back-propagation testing after neural network training.
The source codes in sub_SVM.f is the subroutine for the support vector machine classification process.
The source codes in sub_RandForest.f is the subroutine for the random forest classification process.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
The outputs from our machine-learning-based prediction program are listed as follows:
CElegantIgPPI_MLtest_Result.dat give the predicted probability for all 1081 possible interactions among 46 cell surface proteins in C. elegans.
CElegantIgPPI_MLtest_Result_ranked.dat contains all 1081 possible interactions among 46 cell surface proteins in C. elegans ranked by their overall predicted probability.
CElegansIgPPI_MLtestResult_Ranked.xlsx is the Excel format file of all 1081 possible interactions among 46 cell surface proteins in C. elegans ranked by their overall predicted probability.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
