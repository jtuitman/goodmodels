SetOutputFile("./profiling/data/g=5_trig_p.out");

load "goodmodels_p.m";

print "g=5, trigonal, p";

for p in [11,67,521,4099,32771] do

L:=[];

for i:=1 to 1000 do
  S21,S22,S23,S31,S32:=random_genus5_trigonal(p);
  L:=Append(L,[S21,S22,S23,S31,S32]);
end for;

T:=Cputime();

l:=[**];
failures:=0;

for i:=1 to 1000 do

  Q,W0,Winf:=optimal_model_genus5_trigonal(L[i][1],L[i][2],L[i][3],L[i][4],L[i][5]);

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

