SetOutputFile("./profiling/data/g=3_quartic_p10.out");

load "goodmodels_q.m";

print "g=3, p^10, quartic model";

for q in [3^10,7^10,17^10,37^10,79^10] do

l:=[];

T:=Cputime();

failures:=0;

for i:=1 to 1000 do

  f:=random_genus3(q);

  if Evaluate(f,[0,1,0]) ne 0 then
    Q:=lift_genus3(f);
  else
    Q:=0;
  end if;

  if Q ne 0 then
    l:=Append(l,Q);
  else
    failures:=failures+1;  
  end if;

end for;

T:=Cputime();

for i:=1 to 10 do
  chi:=num_zeta(l[i],q);
end for;

time_pcc:=Cputime(T)/10;

total_mem:=GetMaximumMemoryUsage() div 1024^2;
ResetMaximumMemoryUsage();

print "q =", q, "time pcc =", time_pcc, "memory =", total_mem;

end for;

