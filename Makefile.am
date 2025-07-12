# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS = ks_fd_testing \
                  triangle_element \
				  basis_function_check_fvk \
				  basis_function_check_ks \
				  internal_kink_fvk \
				  internal_kink_ks
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Local sources that each code depends on:
ks_fd_testing_SOURCES = \
 ks_fd_testing.cc
triangle_element_SOURCES = \
 triangle_element.cc
basis_function_check_fvk_SOURCES = \
 basis_function_check_fvk.cc
basis_function_check_ks_SOURCES = \
 basis_function_check_ks.cc
internal_kink_fvk_SOURCES = \
 internal_kink_fvk.cc
internal_kink_ks_SOURCES = \
 internal_kink_ks.cc
#---------------------------------------------------------------------------

ks_fd_testing_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
triangle_element_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
basis_function_check_fvk_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
basis_function_check_ks_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
internal_kink_fvk_LDADD = -L@libdir@ -lc1_foeppl_von_karman \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
internal_kink_ks_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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

AM_CPPFLAGS += -I@includedir@ -Wall -ggdb3 -O0
#EXTRA_DIST = Circle1.1.ele Circle1.1.node Circle1.1.poly
