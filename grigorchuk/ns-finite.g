G := MakeGGroup(6);
a := G.1; b := G.2; c := G.3; d := G.4; x := Comm(a,b);
NC := function(x)
  return NormalClosure(G,Subgroup(G,x));
end;

G01 := G; SetName(G01,"G");
G0 := [G01];

G11 := NC([a,b]); SetName(G11,"a,b");
G12 := NC([a,c]); SetName(G12,"a,c");
G13 := NC([a,d]); SetName(G13,"a,d");
G14 := NC([a*b,c]); SetName(G14,"ab,c");
G15 := NC([a*c,d]); SetName(G15,"ac,d");
G16 := NC([a*d,b]); SetName(G16,"ad,b");
G17 := NC([b,c]); SetName(G17,"H");
G1 := [G11,G12,G13,G14,G15,G16,G17];

G21 := NC([Comm(a,c),b]); SetName(G21,"G',b");
G22 := NC([Comm(a,b),c]); SetName(G22,"G',c");
G23 := NC([Comm(a,b),d]); SetName(G23,"G',d");
G24 := NC([Comm(a,c),a*b]); SetName(G24,"G',ab");
G25 := NC([Comm(a,b),a*c]); SetName(G25,"G',ac");
G26 := NC([Comm(a,b),a*d]); SetName(G26,"G',ad");
G27 := NC([Comm(a,b),a]); SetName(G27,"G',a");
G2 := [G21,G22,G23,G24,G25,G26,G27];

G31 := NC([Comm(a,c),d^a*b]); SetName(G31,"[a,c],Db");
G32 := NC([c]); SetName(G32,"c");
G33 := NC([x,c^a*d]); SetName(G33,"x,Cd");
G34 := NC([Comm(a,c),x]); SetName(G34,"G'");
G35 := NC([b]); SetName(G35,"b");
G36 := NC([Comm(a,d),b^a*c]); SetName(G36,"[a,d],Bc");
G37 := NC([x^2,d]); SetName(G37,"x^2,d");
G3 := [G31,G32,G33,G34,G35,G36,G37];

G41 := NC([Comm(a,c)]); SetName(G41,"[a,c]");
G42 := NC([x]); SetName(G42,"K");
G43 := NC([Comm(a,d),x^2]); SetName(G43,"[a,d],x^2");
G44 := NC([d]); SetName(G44,"d");
G45 := NC([Comm(a,d),x^2*d]); SetName(G45,"[a,d],x^2d");
G4 := [G41,G42,G43,G44,G45];    
