print "";
print "************************************";
print "*  goodmodels Magma library v1.17  *";
print "*                                  *";
print "* by Wouter Castryck & Jan Tuitman *";
print "************************************";
print "";

load "auxpolys.m";
load "quadraticforms.m";
load "convertscroll.m";
load "genus3model.m";
load "genus4model.m";
load "genus5model.m";
SetPath("./pcc_p");
load "pcc_p.m";

zeta_genus3:=function(f)

  // Starting from a plane quartic in 3 variables over a finite field Fp with p prime, 
  // this function tries to compute the zeta function of the corresponding curve.

  Fp:=BaseRing(Parent(f));
  p:=#Fp;

  Q,W0,Winf:=optimal_model_genus3(f);

  if Q ne 0 then
    chi:=num_zeta(Q,p:W0:=W0,Winf:=Winf);
  else
    return 0;
  end if;  

  ZT<T>:=PolynomialRing(IntegerRing());
  zeta:=(ZT!chi)/((1-T)*(1-p*T));
  return zeta;

end function;


zeta_genus4:=function(S2,S3)

  // S2 is a quadric and S3 a cubic, both in 4 variables over a finite field Fp 
  // with p prime, that define a smooth curve of genus 4 in P^3. This function 
  // tries to compute the zeta function of this curve.

  Fp:=BaseRing(Parent(S2));
  p:=#Fp;

  Q,W0,Winf:=optimal_model_genus4(S2,S3);

  if Q ne 0 then
    chi:=num_zeta(Q,p:W0:=W0,Winf:=Winf);
  else
    return 0;
  end if;  

  ZT<T>:=PolynomialRing(IntegerRing());
  zeta:=(ZT!chi)/((1-T)*(1-p*T));
  return zeta;

end function;


zeta_genus5_trigonal:=function(S21,S22,S23,S31,S32)

  // S21,S22,S23 are quadrics and S31,S32 are cubics in 5 variables
  // over a finite field Fp with p prime defining a trigonal curve of genus 5.  
  // This function tries to compute the zeta function of this curve.

  Fp:=BaseRing(Parent(S21));
  p:=#Fp;

  Q,W0,Winf:=optimal_model_genus5_trigonal(S21,S22,S23,S31,S32);

  if Q ne 0 then
    chi:=num_zeta(Q,p:W0:=W0,Winf:=Winf);
  else
    return 0;
  end if;  

  ZT<T>:=PolynomialRing(IntegerRing());
  zeta:=(ZT!chi)/((1-T)*(1-p*T));
  return zeta;

end function;


zeta_genus5_nontrigonal:=function(Q1,Q2,Q3)

  // Q1,Q2,Q3 are quadrics in 5 variables over a finite field Fp with p prime,
  // that define a smooth curve of genus 5 in P^4. This function tries to compute
  // the zeta function of this curve.

  Fp:=BaseRing(Parent(Q1));
  p:=#Fp;

  Q,W0,Winf:=optimal_model_genus5_nontrigonal(Q1,Q2,Q3);

  if Q ne 0 then
    chi:=num_zeta(Q,p:W0:=W0,Winf:=Winf);
  else
    return 0;
  end if;  

  ZT<T>:=PolynomialRing(IntegerRing());
  zeta:=(ZT!chi)/((1-T)*(1-p*T));
  return zeta;

end function;
