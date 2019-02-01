#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Bernadou Basis p3
#----------------------------
cd Validation

echo "Running Bernadou Basis p3 validation "
mkdir RESLT
../multiple_element_tester > OUTPUT_multiple_element_tester
echo "done"
echo " " >> validation.log
echo "Bernadou Basis p3 validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > multiple_element_tester_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/multiple_element_tester_result.dat.gz   \
  multiple_element_tester_results.dat  >> validation.log
fi

# Validation for Bernadou Basis p5
#----------------------------
echo "Running Bernadou Basis p5 validation "
../multiple_element_tester --higher_order > OUTPUT_multiple_element_tester
echo "done"
echo " " >> validation.log
echo "Bernadou Basis p5 validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > multiple_element_tester_p5_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/multiple_element_tester_p5_result.dat.gz   \
  multiple_element_tester_p5_results.dat  >> validation.log
fi

# Validation for Bell Basis
#----------------------------
echo "Running Bell Basis validation "
../bell_element_tester > OUTPUT_bell_element_tester
echo "done"
echo " " >> validation.log
echo "Bell basis validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace_bell.dat > bell_tester_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/bell_tester_result.dat.gz  \
  bell_tester_results.dat  >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
 . $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
 exit 10
