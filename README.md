# IgPPI

Cells detect changes of external environments or communicate with each other through proteins on their surfaces. These cell surface proteins form a complicated network of interactions in order to fulfill their functions. The interactions between cell surface proteins are highly dynamic and thus are challenging to be detected by traditional experimental techniques. Here we developed a new computational framework to tackle this challenge. The primary focus of the framework is to identify interactions between domains in immunoglobulin (Ig) fold, which is the most abundant domain family in cell surface proteins. In practice, we collected all structural data of Ig domain interactions and transformed them into an interface fragment pair library. A high dimensional profile can be then constructed from the library for a given pair of query protein sequences. Multiple machine learning models were used to read this profile, so that the probability of interaction between the query proteins can be predicted. We tested our models to an experimentally derived dataset which contains 564 cell surface proteins in human. The cross-validation results show that we can achieve higher than 70% accuracy in identifying the PPIs within this dataset. We then applied this method to a group of 46 cell surface proteins in C elegans. We screened every possible interaction between these proteins. Many interactions recognized by our machine learning classifiers have been experimentally confirmed in the literatures. In conclusion, our computational platform serves a useful tool to help identifying potential new interactions between cell surface proteins in addition to current state-of-the-art experimental techniques. Moreover, the general framework of the machine learning classification can also be extended to study interactions of proteins in other domain superfamilies.
