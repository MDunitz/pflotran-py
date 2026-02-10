# Step by step guide for making a new reaction sandbox:
1) **copy a WORKING file** -e.g. an example already in the bitbucket, or another that you already know works.
2) change it as you desire. make sure the classnames and everything is consistent. when in doubt, do command-f for the old file name. e.g. if you were copying reaction_sandbox_awinhibit.F90 and wanted to change it to reaction_sandbox_specialcase, you would change EVERY MENTION of awinhibit. All upper/lowercase mentions. Also note that the file name should match. e.g. reaction_sandbox_specialcase.F90
3) **edit reaction_sandbox.F90.** you have two changes to make. First, near the top:
use Reaction_Sandbox_AWInhibit_class
(add another line that is: Reaction_Sandbox_Specialcase_class)

Then, under where you see this:
 case('AWINHIBIT')
    new_sandbox => AWInhibitCreate()
You need to add this:
 case('SPECIALCASE')
    new_sandbox => SpecialcaseCreate()

4) **Edit pflotran_object_files.txt**
        ${common_src}reaction_sandbox_awinhibit.o \
      (put this in the chem_obj section)
5) **edit pflotran_dependencies.txt**

<pre>
reaction_sandbox_awinhibit.o : \
  reaction_sandbox_base.o \
  reactive_transport_aux.o \
  global_aux.o \
  reaction_aux.o
reaction_sandbox_awinhibitacetate.o : \
  reaction_sandbox_base.o \
  reactive_transport_aux.o \
  global_aux.o \
  reaction_aux.o
</pre>

NOTE you also have to add it in dependencies in the reaction sandbox, e.g.
<pre>
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
</pre>

6) **After you’re done all of these steps, you MUST rebuild. You MUST make clean first. **
Note that rebuilding will delete any/all files you have, so if you make any changes to hanford.dat or any other pflotran files, make sure you have them saved elsewhere
7) **Things to keep in mind:**
* apparently it crashed because my error warning was 38 characters. Don’t do that, keep it under 32
