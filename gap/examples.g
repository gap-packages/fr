#Testelements

#Grig
b := FRElement(["b","c","d","a"],[[[4],[2]],[[4],[3]],[[],[1]],[[],[]]],[(),(),(),(1,2)],[1]);
c := State(b,2);
d := State(c,2);
a := State(b,1);


adding_machine := FRElement([[[],[1]]],[(1,2)],[1]);
finitary := FRElement(["fin"],[[[a],[]]],[(1,2)],[1]);
nonfin := FRElement(["nin"],[[[1],[a]]],[()],[1]); 

aut_inf_order := FRElement([[[1,1],[1]]],[(1,2)],[1]);#Be carefull this are not bounded...
afin := FRElement(["c","d","e","f","id",],[[[2],[4]],[[],[3]],[[],[3]],[[],[]],[[],[]]],[(1,2),(),(1,2),(1,2),()],[1]);
z := FRElement(["z"],[[[1],[1]]],[(1,2)],[1]);
degtwo := FRElement(["d2","d1"],[[[1],[2]],[[2],[]]],[(),(1,2)],[1]);

add_10 := FRElement([[[],[],[],[],[],[],[],[],[],[1]]],[(1,2,3,4,5,6,7,8,9,10)],[1]);
add_4 := FRElement(["a4"],[[[],[],[],[1]]],[(1,2,3,4)],[1]);
bounded_4 := FRElement(["b4","c4"],[[[],[],[],[2]],[[2],[],[],[]]],[(),(1,2,3,4)],[1]);

Elms := [[b,"b"],[c,"c"],[d,"d"],[a,"a"],[adding_machine,"adding_machine"],[finitary,"finitary"],[nonfin,"nonfin"],[add_10,"add_10"],[add_4,"add_4"],[bounded_4,"bounded_4"]];
#Groups
Grig := Group(b,c,d,a);
SetName(Grig,"Grig");

RAut_bin := FullSCGroup([1,2],IsFRElement);
FAut_bin := FullSCGroup([1,2],IsFiniteStateFRElement);
Poly0_bin := FullSCGroup([1,2],IsBoundedFRElement);
FINAut_bin := FullSCGroup([1,2],IsFinitaryFRElement);
RAut_4 := FullSCGroup([1,2,3,4],IsFRElement);

Groups := [[RAut_bin,"RAut_bin"],[FAut_bin,"FAut_bin"],[Poly0_bin,"Poly0_bin"],[FINAut_bin,"FINAut_bin"],[RAut_4,"RAut_4"]];

TEST_Grig:= function()
	local testelements, f, g;
	testelements:=[a,b,c,d,a^(b*a*d*a*c),a^(d*a*d*a*c*a*c),b^(d*a*c*a),b^((d*a*c*a)^13*(b*a)^5),c^(a*d)];
	for f in testelements do
		for g in testelements do
			Print("Conjugacy of pair(",f![2],",",g![2],"): ",IsConjugate(Grig,f,g),"\n");
			Print("                  by: ",RepresentativeActionOp(Grig,f,g),"\n");
		od;
	od;
end;
		
