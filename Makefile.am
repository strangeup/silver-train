# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS =  single_element_tester multiple_element_tester bell_element_tester integral_tester
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Local sources that Jack's Poisson code depends on:
single_element_tester_SOURCES = single_element_tester.cc 
multiple_element_tester_SOURCES = multiple_element_tester.cc 
bell_element_tester_SOURCES = bell_element_checker.cc 
integral_tester_SOURCES = integral_test.cc 
# C1_curved_elements.cc C1_curved_association_matrix.cc C1_curved_basic_basis.cc
#---------------------------------------------------------------------------

single_element_tester_LDADD = -L@libdir@ -lgeneric -lC1_basis $(EXTERNAL_LIBS) $(FLIBS)
multiple_element_tester_LDADD = -L@libdir@ -lgeneric -lC1_basis $(EXTERNAL_LIBS) $(FLIBS)
bell_element_tester_LDADD = -L@libdir@ -lgeneric -lC1_basis $(EXTERNAL_LIBS) $(FLIBS)
integral_tester_LDADD = -L@libdir@ -lgeneric -lC1_basis $(EXTERNAL_LIBS) $(FLIBS)

# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
#jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
# jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.

AM_CPPFLAGS += -I@includedir@ -Wall -g -DRANGE_CHECKING -DPARANOID #-fpermissive#-static# -DPURIFY -Wextra -Wuninitialized -Wpedantic #-DRANGE_CHECKING -DPARANOID
