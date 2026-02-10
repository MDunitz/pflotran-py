## Installation 
### linux
See https://documentation.pflotran.org/user_guide/how_to/installation/linux.html#linux-install for more details
```
git clone https://gitlab.com/petsc/petsc petsc
cd petsc
git checkout v3.21.5
```
If running on python=3.13….
```
conda create -n petscpy311 python=3.11 -y
conda activate petscpy311
cd ~/petsc
```
Then configure petsc
```
./configure \
  --with-python-exec=$(which python) \
  --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
  --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3' \
  --with-debugging=0 \
  --download-mpich \
  --download-hdf5 --download-hdf5-fortran-bindings \
  --download-fblaslapack \
  --download-metis --download-parmetis \
  --download-hdf5-configure-arguments="--with-zlib=yes"
```
and follow along with the prompts
(yours will be slightly different based on where petsc is set)
```
   make PETSC_DIR=/home/mdunitz/2025/pflotran/petsc PETSC_ARCH=arch-linux-c-opt all
```
check if the installation is successful (use cmd given by program, this one is just for reference)
```
make PETSC_DIR=/home/mdunitz/2025/pflotran/petsc PETSC_ARCH=arch-linux-c-opt check
```

Set vars
```
export PETSC_DIR=$PWD
export PETSC_ARCH=$(ls -d arch-* | head -n1)
```

Exit petsc dir
```
cd ..
```

Clone pflotran
```
git clone https://bitbucket.org/pflotran/pflotran
```
Compile pflotran
```
cd pflotran/src/pflotran
make pflotran
```
## Running your first test file <3
First, make sure you're in the pflotran directory. Then...
```
module load mpi/openmpi-x86_64
```

Also, you must do your equivalent of this:
```
export PFLOTRAN_DIR=/home/mdunitz/2025//pflotran
```
You're ready to run your first file! Make sure that it's somewhere within the pflotran directory (not necessary, but for organization purposes). Replace filenamehere with your PFLOTRAN input file name (don't include the .in)

```
mpirun -n 1 $PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix filenamehere
```
I suggest that you use one of the test files in our directory, e.g. 9_add_nitrogen.in. In order to run it, you will need to also download/upload the hanford.dat file that we have (it includes a few species that are not in the default file that comes with PFLOTRAN.) Make sure to change the relevant line in the .in file.

Our hanford.dat is located in pflotran_testing.


## Regenerating our 'main' input file, with different constants
Feel free to ignore this section if you don't want to change any of the reaction rates/etc. etc. (this will apply for most things)

However, if you do want to, navigate to saltyBiomass/software_module/master_input_file_generator.py

Change the constants you want, using the dictionaries at the top.

Also, it is default going to output .tec files. Follow the instructions if you want hdf5 output (recommended if you want to test multiple salinity conditions.)


## Generating files with different salt concentrations
Navigate to software_module/reproducible_varconc

Copy in the 9_addnitrogen.in generated earlier, when you changed any variables. If you did not change anything, what it is already in the folder should work.

Run create_modified_files.py (in the same directory). By default, it will use 9_addnitrogen.in. You can drag and drop another .in file into the directory, and change the variable name here: 
```
def create_modified_pflotran_files(source_file="9_addnitrogen.in"):
```

In addition, the default is to create 20 files, from 1 mol to 20 mol of sea water. You may change the dictionary input as necessary.

First four major constituents of sea water were calculated from here, in order to generate the different concentrations: https://docs.google.com/spreadsheets/d/1iNVlg_OOcvQkkKXAuV_2iWS9l-c619CPPVjG7pkcaQE/edit?pli=1&gid=0#gid=0

## Running multiple files at the same time
1) If you're running multiple files, first copy in saltyBiomass/software_module/reproducible_varconc/run_pflotran_batch.sh
2) Make it executable
```
chmod +x run_pflotran_batch.sh
```
3) Run it, from the directory with your pflotran .in files:
```
./run_pflotran_batch.sh
```

* It will run all .in files in the specified directory
* However, note that it runs in sequence, not parallel, so the only speedup you will get is not having to manually run all the input files in the directory
* The script will skip any input files that fail for whatever reason (e.g. concentrations were too high, couldn't simulate), so make sure you keep an eye out for those. Error messages will be printed if so.

## Visualizing PFLOTRAN (single files)

### For plotting CO2/CH4 flu
1. Navigate to saltyBiomass/software_module/pflotran_visualization/step_orchestra
2. Follow the instructions at the top of the file in order to customize outputs/allow it to run.
* E.g. point to the right directory, hdf or tec input, etc. etc.
* You can also choose to visualize something other than CH4 or CO2.
* Note that if you want to do debugging, you can go through the different steps, and set verbose=True in the functions. 
3. Run the orchestra file.
```
python step_orchestra.py
```
4. Open the three different generated html files.

### To generate a reaction diagram (for presentations)
1) Navigate here: saltyBiomass/software_module/pflotran_testing/graphing/graph_diagram.py
2) Copy any reactions from the PFLOTRAN input file into a .txt file into the folder, and name it reactions.txt
3) It will generate an image :) You'll probably want to clean it up on LucidChart.

## Visualizing PFLOTRAN (multiple conditions - default varying water activity)
This will plot CH4/CO2 flux/output as a function of water activity, and also over time.
1) Move all your .tec files into a known location.
2) Run this notebook: saltyBiomass/software_module/pflotran_compare/comparing_aw.ipynb
* Remember to change this as necessary:
working_directory = "saltyBiomass/software_module/pflotran_vis15/modified_pflotran_files_test1"

## Running the reaction sandbox
1) Copy a WORKING file -e.g. an example already in the bitbucket, or another that you already know works. Make sure it is in the same location as the other reaction sandboxes (this matters).
2) Make changes to the file as you wish (look at the other examples). Note that the following steps are for reaction sandbox file with the unique name 'awinhibit'.
3) You must add these lines to the reaction_sandbox.F90, where you can see other instances of them.
```
use Reaction_Sandbox_AWInhibit_class
```

```
 case('AWINHIBIT')
    new_sandbox => AWInhibitCreate()
```

4) Edit pflotran_object_files.txt (put this in the chem_obj section):
```
${common_src}reaction_sandbox_awinhibit.o \
```

5) Edit pflotran_dependencies.txt; put these edits where you see the rest.
```
reaction_sandbox_awinhibit.o : \
  reaction_sandbox_base.o \
  reactive_transport_aux.o \
  global_aux.o \
  reaction_aux.o
```
NOTE you also have to add it in dependencies in the reaction sandbox, e.g.
```
reaction_sandbox.o : \
  global_aux.o \
  input_aux.o \
  material_aux.o \
  option.o \
  output_aux.o \
  pflotran_constants.o \
  reaction_aux.o \
  reaction_sandbox_base.o \
  reaction_sandbox_bioTH.o \
  reaction_sandbox_biohill.o \
  reaction_sandbox_calcite.o \
  reaction_sandbox_chromium.o \
  reaction_sandbox_clm_cn.o \
  reaction_sandbox_equilibrate.o \
  reaction_sandbox_example.o \
  reaction_sandbox_flexbiohill.o \
  reaction_sandbox_gas.o \
  reaction_sandbox_pnnl_cyber.o \
  reaction_sandbox_pnnl_lambda.o \
  reaction_sandbox_radon.o \
  reaction_sandbox_simple.o \
  reaction_sandbox_ufd_wp.o \
  reaction_sandbox_awinhibit.o\
  reaction_sandbox_awinhibitacetate.o\
  reaction_sandbox_awinhibitmethyl.o\
  reactive_transport_aux.o \
  string.o \
  utility.o

```

6) Special note: to copy our results, repeat the steps above with our three unique water activity files, located in software_module/pflotran_sandbox.

7) Also presumably, you can change constants in the reaction sandboxes without needing to remake PFLOTRAN.

## General PFLOTRAN tips:
* If you're adding any new reactions, make sure they're balanced, and check out the doccumentation before doing anything
* Also recall that any new species you add MUST be in the hanford.dat file. If it is not, you will need to edit it in.
