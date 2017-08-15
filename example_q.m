//////////////////////////////////////////////////////// 
// This code is part of the goodmodels MAGMA library  //
//                                                    //
// copyright (c) 2016 Jan Tuitman and Wouter Castryck //
////////////////////////////////////////////////////////

load "goodmodels_q.m";

q:=3^8;

print "//////////////";
print "// genus 3: //";
print "//////////////";
print "";

S:=random_genus3(q);
Q,W0,Winf:=optimal_model_genus3(S);
if Q ne 0 then
  chi:=num_zeta(Q,q:verbose:=true,W0:=W0,Winf:=Winf);
end if;

print "//////////////////////////////";
print "// genus 4, discriminant 0: //";
print "//////////////////////////////";
print "";

S2,S3:=random_genus4(q,0);
Q,W0,Winf:=optimal_model_genus4(S2,S3);
if Q ne 0 then
  chi:=num_zeta(Q,q:verbose:=true,W0:=W0,Winf:=Winf);
end if; 

print "///////////////////////////////////////////";
print "// genus 4, discriminant nonzero square: //";
print "///////////////////////////////////////////";
print "";

S2,S3:=random_genus4(q,1);
Q,W0,Winf:=optimal_model_genus4(S2,S3);
if Q ne 0 then
  chi:=num_zeta(Q,q:verbose:=true,W0:=W0,Winf:=Winf);
end if; 

print "///////////////////////////////////////";
print "// genus 4, discriminant non-square: //";
print "///////////////////////////////////////";
print "";

S2,S3:=random_genus4(q,-1);
n:=Degree(FiniteField(q));
if n lt 8 then
  Q,W0,Winf:=optimal_model_genus4(S2,S3);
else
  Q,W0,Winf:=optimal_model_genus4(S2,S3:alternative_Winf:=true);
end if;
if Q ne 0 then
  chi:=num_zeta(Q,q:verbose:=true,W0:=W0,Winf:=Winf);
end if; 

print "////////////////////////";
print "// genus 5, trigonal: //";
print "////////////////////////";
print "";

S21,S22,S23,S31,S32:=random_genus5_trigonal(q);
Q,W0,Winf:=optimal_model_genus5_trigonal(S21,S22,S23,S31,S32);
if Q ne 0 then
  chi:=num_zeta(Q,q:verbose:=true,W0:=W0,Winf:=Winf);
end if;

print "////////////////////////////";
print "// genus 5, non-trigonal: //";
print "////////////////////////////";
print "";

S21,S22,S23:=random_genus5_nontrigonal(q);
n:=Degree(FiniteField(q));
if n lt 6 then
  Q,W0,Winf:=optimal_model_genus5_nontrigonal(S21,S22,S23);
else
  Q,W0,Winf:=optimal_model_genus5_nontrigonal(S21,S22,S23:alternative_Winf:=true);
end if;
if Q ne 0 then
  chi:=num_zeta(Q,q:verbose:=true,W0:=W0,Winf:=Winf);
end if;



