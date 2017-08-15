//////////////////////////////////////////////////////// 
// This code is part of the goodmodels MAGMA library  //
//                                                    //
// copyright (c) 2016 Jan Tuitman and Wouter Castryck //
////////////////////////////////////////////////////////

load "goodmodels_p.m";

p:=1009;

print "//////////////";
print "// genus 3: //";
print "//////////////";
print "";

S:=random_genus3(p);
Q,W0,Winf:=optimal_model_genus3(S);
if Q ne 0 then
  chi:=num_zeta(Q,p:verbose:=true,W0:=W0,Winf:=Winf);
end if;

print "//////////////////////////////";
print "// genus 4, discriminant 0: //";
print "//////////////////////////////";
print "";

S2,S3:=random_genus4(p,0);
Q,W0,Winf:=optimal_model_genus4(S2,S3);
if Q ne 0 then
  chi:=num_zeta(Q,p:verbose:=true,W0:=W0,Winf:=Winf);
end if; 

print "///////////////////////////////////////////";
print "// genus 4, discriminant nonzero square: //";
print "///////////////////////////////////////////";
print "";

S2,S3:=random_genus4(p,1);
Q,W0,Winf:=optimal_model_genus4(S2,S3);
if Q ne 0 then
  chi:=num_zeta(Q,p:verbose:=true,W0:=W0,Winf:=Winf);
end if; 

print "///////////////////////////////////////";
print "// genus 4, discriminant non-square: //";
print "///////////////////////////////////////";
print "";

S2,S3:=random_genus4(p,-1);
Q,W0,Winf:=optimal_model_genus4(S2,S3);
if Q ne 0 then
  chi:=num_zeta(Q,p:verbose:=true,W0:=W0,Winf:=Winf);
end if; 

print "////////////////////////";
print "// genus 5, trigonal: //";
print "////////////////////////";
print "";

S21,S22,S23,S31,S32:=random_genus5_trigonal(p);
Q,W0,Winf:=optimal_model_genus5_trigonal(S21,S22,S23,S31,S32);
if Q ne 0 then
  chi:=num_zeta(Q,p:verbose:=true,W0:=W0,Winf:=Winf);
end if;

print "////////////////////////////";
print "// genus 5, non-trigonal: //";
print "////////////////////////////";
print "";

S21,S22,S23:=random_genus5_nontrigonal(p);
Q,W0,Winf:=optimal_model_genus5_nontrigonal(S21,S22,S23);
if Q ne 0 then
  chi:=num_zeta(Q,p:verbose:=true,W0:=W0,Winf:=Winf);
end if;

