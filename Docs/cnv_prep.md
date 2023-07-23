# cmd prep

### Prep gatk env:

You can update choose one of those options:

1. create gatk env and use it while runnig the CNVworkflow”:
    
    ```bash
    cd workflow/env
    mamba env create -n gatk -f gatkcondaenv.yml
    ```
    
2. update current env with required packages (not recommended):
    
    ```bash
    cd workflow/env
    conda env update --file gatkcondaenv.yml --name wes_gatk --prune
    ```
    
3. to update Python env only (not recommended):
    
    ```bash
    pip install workflow/env/gatkPythonPackageArchive.zip
    ```
    

### To check install is successful

run:

```bash
conda activate gatk # or wes_gatk if you updated this env 
conda list 
```

You should find `gatkpythonpackages` in the packages list

you can also test the python required package as follows:

```bash
conda activate gatk # or wes_gatk if you updated this env 
python -c "import gcnvkernel"
```

it should run without any errors

### Some troubleshooting:

if you face this error `ModuleNotFoundError: No module named 'numpy.testing.decorators'`  you can try the solution posted in this post: [https://gatk.broadinstitute.org/hc/en-us/community/posts/360056743551-Problem-Installing-GATK-python-environment-SOLUTION-POSTED-](https://gatk.broadinstitute.org/hc/en-us/community/posts/360056743551-Problem-Installing-GATK-python-environment-SOLUTION-POSTED-) 

or this video: [https://www.youtube.com/watch?v=je5iskPlNkg&ab_channel=Sky'sBioinformatics](https://www.youtube.com/watch?v=je5iskPlNkg&ab_channel=Sky%27sBioinformatics)
