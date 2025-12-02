

IMPUTATION documentation

- loading data, initialize module
- log parameters
- replace zeros with NaNs in data, if switched on
- if no values are missing, skip imputation
- if hierarchical_imputation is activated (hierarchical_imputation==1), load (if imputation_table_fp existent and overwrite_imputation_procedures==0) or create and save the imputation procedures table (imputation_table_fp not existent or overwrite_imputation_procedures==1) specifying imputation strategies over hierarchies:
	- one can specify different hierarchies (e.g. superbatch, batch, plate), sorted from highest to lowest hierarchy and/or missingness thresholds
    - hierarchy instances need to be unique,  e.g. there cannot be plate "1" in several batches, or name them as "batch1_plate1" and "batch2_plate1" etc. by concatenating hierarchies
    - an imputation method can be specified for every missingness type (MAR, MNAR, ProtAbsence, ProtAbsence/MAR, ProtAbsence/MNAR), missingness threshold and hierarchy combination
    - e.g. select as hierarchies superbatch;superbatch_batch;superbatch_batch_plate and missingness_thresh 1_0.33_0.25_0.3 to impute for every plate in every batch in every superbatch with missingness up to 25% with one method (e.g. knn with k=15) and with missingness up to 0.3 with a different (or the same) method (e.g knn with k=5), impute over the whole batch with missingness up to 0.33 with a different method and in every superbatch all remaining missing values, again with a different method. Within every imputation step, one can additionally distinguish different imputation methods for all the different missingness types. Also, the parameters for the imputation methods can be adjusted in the imputation procedures table 
    - the table can be created before calling the module_imputation and manually adjusted with different parameters for advanced application, e.g. setting different missingness thresholds within a hierarchy (like batch1 with a differnt threshold than batch2) or e.g. different parameter settings like k-values for knn. The parameters missingness_imputation_thresholds, imputation_sample_hierarchies and imputation_method_* are to generate a first template and in most cases sufficient for the table generation. 
	 
- if eval_misingness_on_protein_presence create a mask specifying whether any precursor of a protein is present in a sample or not
- if eval_missingness_by_detection_proportion identifies hierarchy to evaluate missingness (e.g. over whole dataset or per batch) evaluates missingness (MAR or MNAR) by modelling precursor completeness over abundance and stores missingness type in another mask
- combines mask specifying missingness by protein presence with mask specifying missingness by detection proportion (two masks from previous steps); a precursor in a sample can fall into several categories and thus be assigned with several missingness types (e.g. missing by protein presence and MAR); if save_intermediate_output_for_QA this mask is saved
- if missingness types (by Protein presence, MAR, MNAR) should be saved into filter mask, assigns their values
- if hierarchical imputation, do hierarchical imputation:
	- hierarchies could be something like batch - superbatch or different conditions (healthy/diseased, growing conditions) or covariates (age, sex)
    - loop over all hierarchies (select samples of interest)
        - log number of samples in hierarchy level
        - loop over all missingness thresholds (select precursors of interest)
            - log number of precursors in hierarchy and missingness level
            - extract all missingness types present in precursor and sample subset
            - loop over missingness types
                - extract imputation method and imputation parameters for this hierarchy, missingness threshold and missingness type
                - remember missing values
                - apply imputation on precursor and sample subset, precursors and samples in given subset without any measured values not imputed
                - update data with imputed values for given missingness type, missingness threshold and hierachy samples
                - count number of still missing values
                - update mask, log number of imputed values

- if no hierarchical imputation or imputation over whole dataset (e.g. after hierarchical imputation) is switched on, apply imputation over complete dataset:
    - extract missingness types
    - loop over missingness types
        - define imputation method from parameters dict
        - remember which and how many values were missing before imputation
        - apply imputation
        - remember which and how many values are missing after imputation
        - log them and update mask
- count final number of missing and imputed values and log them
- close module
