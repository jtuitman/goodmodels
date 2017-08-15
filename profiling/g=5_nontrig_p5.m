SetOutputFile("./profiling/data/g=5_nontrig_p5.out");

load "goodmodels_q.m";

print "g=5, non trigonal, p^5";

for q in [3^5,7^5,17^5,37^5,79^5] do

l:=[**];

T:=Cputime();

failures:=0;

for i:=1 to 1000 do

  Q1,Q2,Q3:=random_genus5_nontrigonal(q);
  Q,W0,Winf:=optimal_model_genus5_nontrigonal(Q1,Q2,Q3);

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

