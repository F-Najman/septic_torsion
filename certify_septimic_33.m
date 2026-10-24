/* Standalone complete exact septimic X1(33) class certificate. */
load "probe_33_class2.m";
load "verify_counts_33.m";
load "class_33_degree7_continuation.m";
assert #unknown eq 0;
load "verify_known_33.m";
print "C33_SEPTIMIC_COMPLETENESS_DONE classes=10 unknown=0 full_torsion=33 complete=true",Cputime(t0);
quit;
