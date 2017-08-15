naive_count_genus3:=function(f,k)

  // Takes as input f in Fq[x,y,z] smooth of degree 4 and positive 
  // integer k. Returns the number of points of f=0 in P^2 over F_{q^k}.

  Fq:=BaseRing(Parent(f));
  Fqk:=ext<Fq|k>;
  P2:=ProjectiveSpace(Fq,2);
  X:=Scheme(P2,f);

  return #Points(X,Fqk);

end function;


naive_num_zeta_genus3:=function(f)

  // Takes as input f in Fq[x,y,z] smooth of degree 4. Returns the 
  // numerator of the zeta function of f, computed by naive counting.

  Fq:=BaseRing(Parent(f));
  q:=#Fq;

  cnts:=[];

  for k:=1 to 3 do
    cnts:=Append(cnts,naive_count_genus3(f,k));
  end for;

  chi:=[RationalField()|];

  chi[1]:=1;
  chi[2]:=cnts[1]-(1+q);
  chi[3]:=cnts[2]/2+(chi[2]^2)/2-(1+q^2)/2;
  chi[4]:=cnts[3]/3+chi[3]*chi[2]-(chi[2]^3)/3-(1+q^3)/3;

  for i:=4 to 6 do
    chi[i+1]:=chi[7-i]*q^(i-3);
  end for;

  Z:=IntegerRing();
  ZT<T>:=PolynomialRing(Z);

  return ZT!chi;

end function;


naive_count_genus4:=function(S2,S3,k)

  // Takes as input S2,S3 in Fq[x,y,z,w] homogeneous of degree 2 and 3 defining 
  // a smooth curve of genus 4 and positive integer k. Returns the number of points 
  // of this curve over F_{q^k}.

  Fq:=BaseRing(Parent(S2));
  Fqk:=ext<Fq|k>;
  P3:=ProjectiveSpace(Fq,3);
  X:=Scheme(P3,[S2,Parent(S2)!S3]);

  return #Points(X,Fqk);

end function;


naive_num_zeta_genus4:=function(S2,S3)

  // Takes as input S2,S3 in Fq[x,y,z,w] homogeneous of degree 2 and 3 defining 
  // a smooth curve of genus 4. Returns the numerator of the zeta function of this
  // curve, computed by naive counting.

  Fq:=BaseRing(Parent(S2));
  q:=#Fq;

  cnts:=[];

  for k:=1 to 4 do
    cnts:=Append(cnts,naive_count_genus4(S2,S3,k));
  end for;

  chi:=[RationalField()|];

  chi[1]:=1;
  chi[2]:=cnts[1]-(1+q);
  chi[3]:=cnts[2]/2+chi[2]^2/2-(1+q^2)/2;
  chi[4]:=cnts[3]/3+chi[2]*chi[3]-chi[2]^3/3-(1+q^3)/3;
  chi[5]:=cnts[4]/4+chi[2]*chi[4]+chi[3]^2/2-chi[3]*chi[2]^2+chi[2]^4/4-(1+q^4)/4;

  for i:=5 to 8 do
    chi[i+1]:=chi[9-i]*q^(i-4);
  end for;  

  Z:=IntegerRing();
  ZT<T>:=PolynomialRing(Z);

  return ZT!chi;

end function;


naive_count_genus5_trigonal:=function(S21,S22,S23,S31,S32,k)

  // Takes as input quadrics S21,S22,S23 and cubics S31,S3 in 5 variables
  // over a finite field Fq defining a trigonal curve of genus 5 and a 
  // positive integer k. Returns the number of points of this curve over 
  // F_{q^k}, computed by naive counting.

  Fq:=BaseRing(Parent(S21));
  Fqk:=ext<Fq|k>;
  P4:=ProjectiveSpace(Fq,4);
  X:=Scheme(P4,[S21,S22,S23,S31,S32]);  

  return #Points(X,Fqk);

end function;


naive_num_zeta_genus5_trigonal:=function(S21,S22,S23,S31,S32);

  // Takes as input quadrics S21,S22,S23 and cubics S31,S3 in 5 variables
  // over a finite field Fq defining a trigonal curve of genus 5. Returns 
  // the zeta function of this curve, computed by naive counting.

  Fq:=BaseRing(Parent(S21));
  q:=#Fq;

  cnts:=[];

  for k:=1 to 5 do
    cnts:=Append(cnts,naive_count_genus5_trigonal(S21,S22,S23,S31,S32,k));
  end for;

  chi:=[RationalField()|];

  chi[1]:=1;
  chi[2]:=cnts[1]-(1+q);
  chi[3]:=cnts[2]/2+chi[2]^2/2-(1+q^2)/2;
  chi[4]:=cnts[3]/3+chi[2]*chi[3]-chi[2]^3/3-(1+q^3)/3;
  chi[5]:=cnts[4]/4+chi[2]*chi[4]+chi[3]^2/2-chi[3]*chi[2]^2+chi[2]^4/4-(1+q^4)/4;
  chi[6]:=cnts[5]/5+chi[2]*chi[5]+chi[3]*chi[4]-chi[2]*chi[3]^2-chi[2]^2*chi[4]+chi[2]^3*chi[3]-chi[2]^5/5-(1+q^5)/5;

  for i:=6 to 10 do
    chi[i+1]:=chi[11-i]*q^(i-5);
  end for;  

  Z:=IntegerRing();
  ZT<T>:=PolynomialRing(Z);

  return ZT!chi;

end function;


naive_count_genus5_nontrigonal:=function(S1,S2,S3,k)

  // Takes as input S1,S2,S3 in Fq[x,y,z,w,v] homogeneous of degree 2 defining
  // a smooth curve of genus 5 and positive integer k. Returns the number of 
  // points of this curve over F_{q^k}, computed by naive counting.

  Fq:=BaseRing(Parent(S1));
  Fqk:=ext<Fq|k>;
  P3:=ProjectiveSpace(Fq,4);
  X:=Scheme(P3,[S1,Parent(S1)!S2,Parent(S1)!S3]);

  return #Points(X,Fqk);

end function;


naive_num_zeta_genus5_nontrigonal:=function(S1,S2,S3);

  // Takes as input S1,S2,S3 in Fq[x,y,z,w,v] homogeneous of degree 2 defining 
  // a smooth curve of genus 5. Returns the zeta function of this curve, computed 
  // by naive counting.

  Fq:=BaseRing(Parent(S1));
  q:=#Fq;

  cnts:=[];

  for k:=1 to 5 do
    cnts:=Append(cnts,naive_count_genus5_nontrigonal(S1,S2,S3,k));
  end for;

  chi:=[RationalField()|];

  chi[1]:=1;
  chi[2]:=cnts[1]-(1+q);
  chi[3]:=cnts[2]/2+chi[2]^2/2-(1+q^2)/2;
  chi[4]:=cnts[3]/3+chi[2]*chi[3]-chi[2]^3/3-(1+q^3)/3;
  chi[5]:=cnts[4]/4+chi[2]*chi[4]+chi[3]^2/2-chi[3]*chi[2]^2+chi[2]^4/4-(1+q^4)/4;
  chi[6]:=cnts[5]/5+chi[2]*chi[5]+chi[3]*chi[4]-chi[2]*chi[3]^2-chi[2]^2*chi[4]+chi[2]^3*chi[3]-chi[2]^5/5-(1+q^5)/5;

  for i:=6 to 10 do
    chi[i+1]:=chi[11-i]*q^(i-5);
  end for;  

  Z:=IntegerRing();
  ZT<T>:=PolynomialRing(Z);

  return ZT!chi;

end function;


