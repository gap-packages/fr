if START_CP_BRANCH@ <> [] then
	InstallValue(START_CP_BRANCH@,[]);
fi;
InstallMethod(InitConjugateForBranchGroups,"Branch Groups",[IsFRGroup,IsList], function(G,L)
		if not HasName(G) then 
			SetName(G,Concatenation("BranchGroup_",String(Size(START_CP_BRANCH@)+1)));
		fi;		
		Add(START_CP_BRANCH@,[Name(G),L]);
	end);
MAKE_READ_WRITE_GLOBAL("CONJUGATORS_BRANCH@");
UNBIND_GLOBAL("CONJUGATORS_BRANCH@");
BindGlobal("CONJUGATORS_BRANCH@",function(G,g,h)
	local Start, Gen, K, R,RR, Conjugators_branch_rek,l,k,InStart,GetStart;
	InStart := function(N)
		local H;
		for H in START_CP_BRANCH@ do
			if H[1] = N then
				return true;
			fi;
		od;
		return false;
	end;
	GetStart := function(G)
		local H;
		for H in START_CP_BRANCH@ do
			if H[1] = G then
				return H[2];
			fi;
		od;
		return fail;
	end;
	if Alphabet(G) = [1,2] and HasName(G) and InStart(Name(G)) then
		Start := GetStart(Name(G));
	else
		return fail;
	fi;
	if g = h then
		return One(g);
	fi;
	Gen := GeneratorsOfGroup(G);
	#Branching Subgroup
	K := BranchingSubgroup(G);
	#Representative System
	RR := RightCosets(G,K);
	R := List(RR,x->Representative(x));
	
	#First for binary alphabet.
	Conjugators_branch_rek := function(g,h)
		local L1,L2,c,CT,i;
		Info(InfoFRCP,3,"Computing g,h=",g,",",h,"\n");
		if IsOne(g) or IsOne(h) then
			if g = h then
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     g=h=1 So return R\n");
				return R;
			else
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     g,h is One but the other not. So return []\n");
				return [];
			fi;
		fi;
		if Activity(g) <> Activity(h) then
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Act(g)<>Act(h). So return []\n");
			return [];
		fi;
		if g in Gen and h in Gen then
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     g,h are Generators. So return Start[g][h]\n");
			return Start[Position(Gen,g)][Position(Gen,h)];
		fi;

		
		if IsOne(Activity(g)) then
		#----------Conjugator with trivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Conjugator for elm with trivial activity...\n");
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Try a Conjugator with trivial activity...\n");
			L1 := Conjugators_branch_rek(State(g,1),State(h,1));
			CT := [];
			if Size(L1) > 0 then
				L2 := Conjugators_branch_rek(State(g,2),State(h,2));
				for l in L1 do
					for k in L2 do
						c := FRElement([[[l],[k]]],[()],[1]);
						for i in [1..Size(RR)] do
							if c in RR[i] then
								CT[i] := c;
								break;
							fi;
						od;
					od;
				od;
				if Size(CT)>0 then
					Print("Something found..\n");
					Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Conjugator with trivial activity found. Return it\n");
					return CT;
				fi;
			fi;
		#----------Conjugator with nontrivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Try a Conjugator with non-trivial activity...\n");
			L1 := Conjugators_branch_rek(State(g,1),State(h,2));
			if Size(L1) > 0 then
				L2 := Conjugators_branch_rek(State(g,2),State(h,1));
				for l in L1 do
					for k in L2 do
						c := FRElement([[[l],[k]]],[(1,2)],[1]);
						for i in [1..Size(RR)] do
							if c in RR[i] then
								CT[i] := c;
								break;
							fi;
						od;
					od;
				od;
			fi;
			if Size(CT) > 0 then
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Conjugator with non-trivial activity found. Return it\n");
			else
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     No Conjugator found. Return []\n");
			fi;
			return CT;	
		else
		#----------Conjugator with trivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Conjugator for elm with non-trivial activity...\n");
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Try a Conjugator with trivial activity...\n");
			L1 := Conjugators_branch_rek(State(g,1)*State(g,2),State(h,1)*State(h,2));
			CT := [];
			for l in L1 do
				c := FRElement([[[l],[State(g,1)^-1*l*State(h,1)]]],[()],[1]);
				for i in [1..Size(RR)] do
					if c in RR[i] then
						CT[i] := c;
						break;
					fi;
				od;
			od;
			if Size(CT)>0 then
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Conjugator with trivial activity found. Return it\n");
				return CT;
			fi;		
		#----------Conjugator with nontrivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Try a Conjugator with non-trivial activity...\n");
			L1 := Conjugators_branch_rek(State(g,1)*State(g,2),State(h,2)*State(h,1));
			CT := [];
			for l in L1 do
				c := FRElement([[[l],[State(g,1)^-1*l*State(h,1)]]],[(1,2)],[1]);
				for i in [1..Size(RR)] do
					if c in RR[i] then
						CT[i] := c;
						break;
					fi;
				od;
			od;
			if Size(CT) > 0 then
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     Conjugator with non-trivial activity found. Return it\n");
			else
				Info(InfoFRCP,3,"Computing g,h=",g,",",h,"     No Conjugator found. Return []\n");
			fi;
			return CT;	
		fi;
	end;
	return Conjugators_branch_rek(g,h);
end);

###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%      IsConjugate        %%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%	 RepresentativeActionOp %%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################
InstallOtherMethod(RepresentativeActionOp,
	"Computes a conjugator in the given Branch group ",
	[ IsFRGroup,IsFRElement,IsFRElement], 
	function(G,g,h)
		local con;
		con := CONJUGATORS_BRANCH@(G,g,h);
		if con <> fail then
			if Size(con)>0 then
				return Representative(con);
			fi;
			return fail;
		fi;
		TryNextMethod();
		end);
InstallMethod(IsConjugate,
	"For Branch Groups",
	[ IsFRGroup,IsFRElement,IsFRElement], 
  function(G,g,h)
  	local con;
  	Print("For Branch Method...\n");
  	con := CONJUGATORS_BRANCH@(G,g,h);
		if con <> fail then
			if Size(con)>0 then
				return true;
			fi;
			return false;
		fi;
		Print("Next Method...\n");
		TryNextMethod();
		end);	
	
######################################Example###########################################
AddGrig_toStart := function()
	local G,a,b,c,d,x,S,SL,y,C,RR,K,find_elm;
	G:= GrigorchukGroup;
	a:= G.1; b:=G.2; c:=G.3; d:=G.4;
	K := BranchingSubgroup(G);
	RR := RightCosets(G,K);
	find_elm := function(elm)
		local i;
		for i in [1..Size(RR)] do
			if elm in RR[i] then
				return i;
			fi;
		od;
	end;
	S := [];
	for x in GeneratorsOfGroup(G) do
		SL := [];
		for y in GeneratorsOfGroup(G) do
			if x = y then
				if x=a then
					C:= [];
					C[find_elm(One(a))] := One(a);
					C[find_elm(a)] := a;
					C[find_elm(a*d*a*d)] := a*d*a*d;
					C[find_elm(d*a*d)] := d*a*d;
					Add(SL,C);
				fi;
				if x = d then
					C:= [];
					C[find_elm(One(a))] := One(a);
					C[find_elm(a*d*a)] := a*d*a;
					C[find_elm(a*d*a*d)] := a*d*a*d;
					C[find_elm(d)] := d;
					C[find_elm(b)] := b;
					C[find_elm(b*a*d*a)] := b*a*d*a;
					C[find_elm(b*a*d*a*d)] := b*a*d*a*d;
					C[find_elm(c)] := c;
					Add(SL,C);
				fi;
				if x in [b,c] then
					C:=[];
					C[find_elm(One(a))] := One(a);
					C[find_elm(d)] := d;
					C[find_elm(b)] := b;
					C[find_elm(c)] := c;					
				Add(SL,C);
				fi;
			else
				Add(SL,[]);
			fi;
		od;
		Add(S,SL);
	od;
	InitConjugateForBranchGroups(G,S);
end;


	
	
	
	
	

