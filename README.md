# BONZ_ALICE_BDT
All scripts and files used for my bachelor thesis project at University Utrecht, done in 2018/2019  
  
Workflow for these files is:  
First create the training files with create_training_file.C  
After this train.C can be run to train the BDT  
With apply_BDT.C, you can apply the trained BDT  
Optionally, regular_select.C can be run to create selections from the regular ALICE method  
Next use fit.C to fit the invariant mass spectra  
Finally, with eff_correction.C any obtained fits can be corrected for efficiency  
  
There are some redundant (and legacy) files and methods here and there, so read and understand all code before executing!  
