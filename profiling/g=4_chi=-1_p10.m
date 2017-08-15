SetOutputFile("./profiling/data/g=4_chi=-1_p10.out");

load "goodmodels_q.m";

print "g=4, chi=-1, p^10";

for q in [3^10,7^10,17^10,37^10,79^10] do

l:=[**];

T:=Cputime();

failures:=0;

for i:=1 to 1000 do

  S2,S3:=random_genus4(q,-1);
  Q,W0,Winf:=optimal_model_genus4(S2,S3:alternative_Winf:=true);

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

