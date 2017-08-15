SetOutputFile("./profiling/data/tetragonal_g=7_q.out");

load "goodmodels_q.m";
load "tetragonal.m";

g:=7;

print "trigonal, q, g=", g;

N:=5;

for q in [3^5,7^5,3^10,7^10] do
  T:=Cputime();
  failures:=0;
  for i:=1 to N do
    f1,f2:=generic_tetragonal(q,g);
    Q,W0,Winf:=lift_tetragonal(f1,f2);
    if Q ne 0 then
      chi:=num_zeta(Q,q:W0:=W0,Winf:=Winf);
    else
      failures:=failures+1;
    end if;
  end for;
  print q,Cputime(T)/(N-failures),failures;
  GetMaximumMemoryUsage() div 1024^2;
end for;


