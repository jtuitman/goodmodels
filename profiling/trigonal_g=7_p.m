SetOutputFile("./profiling/data/trigonal_g=7_p.out");

load "goodmodels_p.m";
load "trigonal.m";

g:=7;

print "trigonal, p, g=", g;

N:=5;

for p in [11,67,521,4099,32771] do
  T:=Cputime();
  failures:=0;
  for i:=1 to N do
    f:=generic_trigonal(p,g);
    Q,W0,Winf:=lift_trigonal(f);
    if Q ne 0 then
      chi:=num_zeta(Q,p:W0:=W0,Winf:=Winf);
    else
      failures:=failures+1;
    end if;
  end for;
  print p,Cputime(T)/(N-failures),failures;
  GetMaximumMemoryUsage() div 1024^2;
end for;


