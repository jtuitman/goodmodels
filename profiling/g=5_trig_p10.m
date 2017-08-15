SetOutputFile("./profiling/data/g=5_trig_p10.out");

load "goodmodels_q.m";

print "g=5, trigonal, p^10";

for q in [3^10,7^10,17^10,37^10,79^10] do

L:=[];

for i:=1 to 1000 do
  S21,S22,S23,S31,S32:=random_genus5_trigonal(q);
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
  chi:=num_zeta(l[i][1],q:W0:=l[i][2],Winf:=l[i][3]);
end for;

time_pcc:=Cputime(T)/10;

total_mem:=GetMaximumMemoryUsage() div 1024^2;
ResetMaximumMemoryUsage();

print "q =", q, "failures =", failures, "time lift =", time_lift, "time pcc =", time_pcc, "memory =", total_mem;

end for;

