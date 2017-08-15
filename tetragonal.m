function generic_tetragonal(q,genus)
  
  Fq:=FiniteField(q);
  R<x,y,z> := PolynomialRing(Fq,3);
  b1 := (genus - 4) div 2; // GENERIC SCHREYER INVARIANTS
  b2 := genus - 5 - b1;
  c := (genus - 1) div 3; // GENERIC MARONI INVARIANTS
  b := (genus - 2) div 3;
  a := genus - 3 - b - c;
  AA := AffineSpace(Fq,3);
  repeat
    repeat
      repeat
        f := R ! 0;
        g := R ! 0;
        for i in [0..(2*b - b1)] do
          f +:= Random(Fq)*x^i*y^2;
        end for;
        for i in [0..(2*b - b2)] do
          g +:= Random(Fq)*x^i*y^2;
        end for;
        for i in [0..(2*b - b1 + (c-b))] do
          f +:= Random(Fq)*x^i*y;
        end for;
        for i in [0..(2*b - b2 + (c-b))] do
          g +:= Random(Fq)*x^i*y;
        end for;
        for i in [0..(2*b - b1 + 2*(c-b))] do
          f +:= Random(Fq)*x^i;
        end for;
        for i in [0..(2*b - b2 + 2*(c-b))] do
          g +:= Random(Fq)*x^i;
        end for;
        for i in [0..(2*c - b1 - (c-a))] do
          f +:= Random(Fq)*x^i*z;
        end for;
        for i in [0..(2*c - b2 - (c-a))] do
          g +:= Random(Fq)*x^i*z;
        end for;
        for i in [0..(2*c - b1 - 2*(c-a))] do
          f +:= Random(Fq)*x^i*z^2;
        end for;
        for i in [0..(2*c - b2 - 2*(c-a))] do
          g +:= Random(Fq)*x^i*z^2;
        end for;
        for i in [0..(2*b - b1 - (b-a))] do
          f +:= Random(Fq)*x^i*y*z;
        end for;
        for i in [0..(2*b - b2 - (b-a))] do
          g +:= Random(Fq)*x^i*y*z;
        end for;
        S := Scheme(AA,[f,g]);
      until IsIrreducible(S) and Dimension(S) eq 1;
      C := Curve(S);
    until IsAbsolutelyIrreducible(C);
    gen := Genus(C);
  until genus eq gen;
  
  return f,g;
end function;


lift_tetragonal:=function(f,g);

  Z:=IntegerRing();
  Fq:=BaseRing(Parent(f));
  q:=#Fq;

  if IsPrime(q) then
    R:=PolynomialRing(IntegerRing(),3);
    S<x,y,z>:=PolynomialRing(RationalField(),3);
    
    fQ:=S!(R!f);
    gQ:=S!(R!g);
    
    F:=Resultant(fQ,gQ,3);

    C:=Coefficients(F); 
    denom:=1;
    for i:=1 to #C do
      denom:=LCM(denom,Denominator(C[i]));
    end for;
    F:=F*denom;

    S<x>:=PolynomialRing(RationalField());
    R<y>:=PolynomialRing(S);
    h:=Evaluate(F,[x,y,0]);
    
    lc:=LeadingCoefficient(h);
    h:=Numerator(lc^(Degree(h)-1)*Evaluate(h,y/LeadingCoefficient(h)));

    S<x>:=PolynomialRing(Z);
    R<y>:=PolynomialRing(S);
    g:=R!0;
    for i:=0 to Degree(h) do
      poly:=Coefficient(h,i);
      for j:=0 to Degree(poly) do
        g:=g+(Z!Coefficient(poly,j))*y^i*x^j;
      end for;
    end for;
    succes,W0,Winf:=check_model_p(g,q);
    if succes then
      return g,W0,Winf;
    end if;
  else // q is not prime 
    n:=Degree(Fq);
    p:=Characteristic(Fq);
    conway:=PolynomialRing(IntegerRing())!DefiningPolynomial(Fq);
    K<a>:=NumberField(conway);
    R<x,y,z>:=PolynomialRing(K,3);    

    C:=Coefficients(f);
    terms:=Terms(f);
    fK:=R!0;
    for i:=1 to #C do
      E:=Exponents(terms[i]);
      for j:=1 to n do
        fK:=fK+(Z!Eltseq(C[i])[j])*a^(j-1)*R.1^E[1]*R.2^E[2]*R.3^E[3];
      end for;
    end for;

    C:=Coefficients(g);
    terms:=Terms(g);
    gK:=R!0;
    for i:=1 to #C do
      E:=Exponents(terms[i]);
      for j:=1 to n do
        gK:=gK+(Z!Eltseq(C[i])[j])*a^(j-1)*R.1^E[1]*R.2^E[2]*R.3^E[3];
      end for;
    end for;

    F:=Resultant(fK,gK,3);

    C:=Coefficients(F); 
    denom:=1;
    for i:=1 to #C do
      for j:=1 to n do
        denom:=LCM(denom,Denominator(C[i][j]));
      end for;
    end for;
    F:=F*denom;

    S<x>:=PolynomialRing(K);
    R<y>:=PolynomialRing(S);
    h:=Evaluate(F,[x,y,0]);

    lc:=LeadingCoefficient(h);
    h:=Numerator(lc^(Degree(h)-1)*Evaluate(h,y/LeadingCoefficient(h)));

    O<a>:=PolynomialRing(Z);
    S<x>:=PolynomialRing(O);
    R<y>:=PolynomialRing(S); // required by pcc_q
    g:=R!0;
    for i:=0 to Degree(h) do
      poly:=Coefficient(h,i);
      for j:=0 to Degree(poly) do
        for k:=1 to n do
          g:=g+(Z!Coefficient(poly,j)[k])*a^(k-1)*y^i*x^j;
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
