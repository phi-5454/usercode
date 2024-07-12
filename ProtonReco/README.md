# Use

e.g. to analyze the file `root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM21/270000/F8A4502E-0D50-E911-A832-0025904B302E.root`, and output it to the scripts directory, run (from this directory):

```
# Change to correct architecture
cmssw-el6

# Set up CMS environment
cmsenv

# Build the code
scram b -j 8

# Permission for files
voms-proxy-init

cd scripts
./run_reco.sh TOTEM21/270000/F8A4502E-0D50-E911-A832-0025904B302E.root ./output.root
```
