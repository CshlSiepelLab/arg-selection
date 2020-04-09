### WORKFLOW TO RUN RELATE AND EXTRACT FEATURES ###

(To be completed)

1) Partition pickle files to enable parallel processing

IN:  <pkl_path>/<pkl_pref>_*.pkl
OUT: <pkl_pref>/<pkl_pref>_pgv_*.pkl

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2) Run relate and extract features

IN:  <handle>/<handle>_pgv_*.pkl
OUT: <handle>/<handle>_<TAG>_inf_fea_*.pickle

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### DEPRECATED ###
### 3) IF any thread fails due to memory issues, you will need to re-run the specific thread ###

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3) Combine the extracted features

IN:  <handle>/<handle>_<TAG>_inf_fea_*.pickle
OUT: <handle>_inf_fea_<oTAG>.pkl