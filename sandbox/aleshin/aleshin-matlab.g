a := function(n,i)
  if n=1 then return 1;
  elif i > n/2 then return b(n/2,i-n/2);
  else return n/2+c(n/2,i);
  fi;
end;
b := function(n,i)
  if n=1 then return 1;
  elif i > n/2 then return c(n/2,i-n/2);
  else return n/2+b(n/2,i);
  fi;
end;
c := function(n,i)
  if n=1 then return 1;
  elif i > n/2 then return n/2+a(n/2,i-n/2);
  else return a(n/2,i);
  fi;
end;
A := function(n,i)
  if n=1 then return 1;
  elif i > n/2 then return C(n/2,i-n/2);
  else return n/2+B(n/2,i);
  fi;
end;
B := function(n,i)
  if n=1 then return 1;
  elif i > n/2 then return B(n/2,i-n/2);
  else return n/2+C(n/2,i);
  fi;
end;
C := function(n,i)
  if n=1 then return 1;
  elif i > n/2 then return n/2+A(n/2,i-n/2);
  else return A(n/2,i);
  fi;
end;

n := 4096;
PrintTo("matlab.mat","% aleshin group sparse matrix for n=",n,"\n");
PI := function(x)AppendTo("matlab.mat",x," ");end;
Pln := function()AppendTo("matlab.mat","\n");end;

for j in [1..6] do for i in [1..n] do PI(i); od; od; Pln();
for j in [a,b,c,A,B,C] do for i in [1..n] do PI(j(n,i)); od; od; Pln();
for i in [1..6*n] do PI(1); od; Pln();
