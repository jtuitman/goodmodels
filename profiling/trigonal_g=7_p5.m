SetOutputFile("./profiling/data/trigonal_g=7_p5.out");

load "goodmodels_q.m";
load "trigonal.m";

g:=7;

print "trigonal, p^5, g=", g;

N:=5;

for q in [3^5,7^5,17^5,37^5,79^5] do
  T:=Cputime();
  failures:=0;
  for i:=1 to N do
    f:=generic_trigonal(q,g);
    Q,W0,Winf:=lift_trigonal(f);
    if Q ne 0 then
      chi:=num_zeta(Q,q:W0:=W0,Winf:=Winf);
    else
      failures:=failures+1;
    end if;
  end for;
  print q,Cputime(T)/(N-failures),failures;
  GetMaximumMemoryUsage() div 1024^2;
end for;


