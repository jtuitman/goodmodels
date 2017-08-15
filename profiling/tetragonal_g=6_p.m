SetOutputFile("./profiling/data/tetragonal_g=6_p.out");

load "goodmodels_p.m";
load "tetragonal.m";

g:=6;

print "tetragonal, p, g=", g;

N:=5;

for p in [11,67,521,4099] do
  T:=Cputime();
  failures:=0;
  for i:=1 to N do
    f1,f2:=generic_tetragonal(p,g);
    Q,W0,Winf:=lift_tetragonal(f1,f2);
    if Q ne 0 then
      chi:=num_zeta(Q,p:W0:=W0,Winf:=Winf);
    else
      failures:=failures+1;
    end if;
  end for;
  print p,Cputime(T)/(N-failures),failures;
  GetMaximumMemoryUsage() div 1024^2;
end for;


