SetOutputFile("./profiling/data/g=4_chi=-1_p.out");

load "goodmodels_p.m";

print "g=4, chi=-1, p";

for p in [11,67,521,4099,32771] do

l:=[**];

T:=Cputime();

failures:=0;

for i:=1 to 1000 do

  S2,S3:=random_genus4(p,-1);
  Q,W0,Winf:=optimal_model_genus4(S2,S3);

  if Q ne 0 then
    l:=Append(l,[*Q,W0,Winf*]);
  else
    failures:=failures+1;
  end if;

end for;

time_lift:=Cputime(T)/1000;
T:=Cputime();

for i:=1 to 10 do
  chi:=num_zeta(l[i][1],p:W0:=l[i][2],Winf:=l[i][3]);
end for;

time_pcc:=Cputime(T)/10;

total_mem:=GetMaximumMemoryUsage() div 1024^2;
ResetMaximumMemoryUsage();

print "p =", p, "failures =", failures, "time lift =", time_lift, "time pcc =", time_pcc, "memory =", total_mem;

end for;

