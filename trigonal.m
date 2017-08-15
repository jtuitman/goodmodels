function generic_trigonal(q,genus)
  
  Fq:=FiniteField(q);
  R<x,y> := PolynomialRing(Fq,2);
  b := (genus - 1) div 2; // GENERIC MARONI INVARIANTS
  a := genus - 2 - b;
  AA := AffinePlane(Fq);
  repeat
    repeat
      repeat
      f := R ! 0;
      for i in [0..(2*b + 2 - a)] do
        f +:= Random(Fq)*x^i;
      end for;
      for i in [0..(2*b + 2 - a - (b - a))] do
        f +:= Random(Fq)*x^i*y;
      end for;
      for i in [0..(2*b + 2 - a - 2*(b-a))] do
        f +:= Random(Fq)*x^i*y^2;
      end for;
      for i in [0..(2*b + 2 - a - 3*(b-a))] do
        f +:= Random(Fq)*x^i*y^3;
      end for;
      until IsIrreducible(f);
      C := Curve(AA,f);
    until IsAbsolutelyIrreducible(C);
    gen := Genus(C);
  until genus eq gen;
  return f;

end function;


lift_trigonal:=function(f)
  
  Z:=IntegerRing();
  Fq:=BaseRing(Parent(f));
  q:=#Fq;

  Fqx:=PolynomialRing(Fq);
  Fqxy:=PolynomialRing(Fqx);
  f:=Evaluate(f,[Fqx.1,Fqxy.1]);

  if IsPrime(q) then
    S<x>:=PolynomialRing(Z);
    R<y>:=PolynomialRing(S);
    g:=R!0;
    for i:=0 to Degree(f) do
      poly:=Coefficient(f,i);
      for j:=0 to Degree(poly) do
        g:=g+(Z!Coefficient(poly,j))*y^i*x^j;
      end for;
    end for;
    lc:=LeadingCoefficient(g);
    g:=Numerator(lc^(Degree(g)-1)*Evaluate(g,y/lc));
    succes,W0,Winf:=check_model_p(g,q);
    if succes then
      return g,W0,Winf;
    end if;
  else // q is not prime
    n:=Degree(Fq);
    p:=Characteristic(Fq);
    conway:=PolynomialRing(IntegerRing())!DefiningPolynomial(Fq);
    K<a>:=NumberField(conway);
    S<x>:=PolynomialRing(K);
    R<y>:=PolynomialRing(S);
    h:=R!0;
    for i:=0 to Degree(f) do
      poly:=Coefficient(f,i);
      for j:=0 to Degree(poly) do
        coefs:=Eltseq(Coefficient(poly,j));
        for k:=1 to n do
          h:=h+(Z!coefs[k])*a^(k-1)*x^j*y^i;
        end for;
      end for;
    end for;
    lc:=LeadingCoefficient(h);
    h:=Numerator(lc^(Degree(h)-1)*Evaluate(h,y/lc));
    O<a>:=PolynomialRing(Z);
    S<x>:=PolynomialRing(O);
    R<y>:=PolynomialRing(S); // required by pcc_q
    g:=R!0;
    for i:=0 to Degree(h) do
      poly:=Coefficient(h,i);
      for j:=0 to Degree(poly) do
        coefs:=Eltseq(Coefficient(poly,j));
        for k:=1 to n do
          g:=g+(Z!coefs[k])*a^(k-1)*x^j*y^i;
        end for;
      end for;
    end for;
    succes,W0,Winf:=check_model_q(g,p,n);
    if succes then
      return g,W0,Winf;
    end if;
  end if;

  return 0,0,0;
end function;
