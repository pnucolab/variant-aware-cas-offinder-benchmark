
- First create a conda env and the install requirements.txt after activating the environment
 ```
conda create -n crispr
```
- Make sure all files have all permissions
- Change this based on your requirements from benchmark.py.

    - num_guides,
    - guide_length,
    - pam, and
    - mismatches

- To execute
  ```
  conda activate crispr
  ```
  - And then
  
  ```
 ./vcf-cas-offinder.py -i chr1.vcf.gz -r /home/abyot/training/crispritz/genome-file/genome-file-test/gzchr1.fa -t generated_sequences.txt -d G1
```
