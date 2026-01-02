

## Database Setup



## Defense-finder Setup
Due to a bug in the current macsyfinder version (see https://github.com/mdmparis/defense-finder/issues/91), defense-finder models must be set up manually before running. Here's how:
```bash
# 1. Create a new conda environment
conda create -n defensefinder -c bioconda -c conda-forge python=3.10 hmmer
conda activate defensefinder

# 2. install defense-finder 2.0.1
pip install mdmparis-defense-finder==2.0.1

# 3. update defense-finder models
defense-finder update

# 4. Manually install CasFinder 3.1.0
cd ~/.macsyfinder/models && rm -rf CasFinder 
git clone https://github.com/macsy-models/CasFinder.git && cd CasFinder && git checkout 3.1.0

# DO NOT RUN defense-finder update again after this step
```
