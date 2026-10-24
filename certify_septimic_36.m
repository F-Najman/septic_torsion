/* Pending certificate until its final completion marker is observed. */
load "probe_36_class5_v2.m";
load "prepare_36_degree7_orbits.m";
load "quotient_36_to_2_18.m";
load "class_36_quotient_filter.m";
load "class_36_known_lifts_orbits.m";
assert #unknown eq 0 and #knownclasses eq 12;
load "verify_known_36.m";
load "export_36_quotient_certificate.m";
print "C36_SEPTIMIC_COMPLETENESS_DONE known_orbits=2 known_classes=12 unknown=0 full_torsion=36 complete=true",Cputime(t0);
quit;
