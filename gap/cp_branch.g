if not IsList(START_CP_BRANCH@) then
	InstallValue(START_CP_BRANCH@,[]);
fi;

InstallMethod(InitConjugateForBranchGroups,"Branch Groups",[IsFRGroup,IsList], function(G,L)
		if not HasName(G) then 
			SetName(G,Concatenation("BranchGroup_",String(Size(START_CP_BRANCH@)+1)));
		fi;
				
		Add(START_CP_BRANCH@,rec(	name:=Name(G),
															inital_con:=L,
															Branchstructure:=BranchStructure(G),
															RepSystem:=List(~.Branchstructure.group,x->PreImagesRepresentative(~.Branchstructure.quo,x))));
	end);
InstallOtherMethod(InitConjugateForBranchGroups,"Branch Groups",[IsFRGroup], function(G)
	local N,n,m;
	N := Nucleus(G);
	for n in N do
		for m in N do
			#TODO 
		od;
	od;
end);
BindGlobal("CONJUGATORS_BRANCH@",function(G,g,h)
	local CP_init, Start, Gen, B, BS, Conjugators_branch_rek,l,k,InStart,GetStart,rek_count;
	InStart := function(name)
		local H;
		for H in START_CP_BRANCH@ do
			if H.name = name then
				return true;
			fi;
		od;
		return false;
	end;
	GetStart := function(name)
		local H;
		for H in START_CP_BRANCH@ do
			if H.name = name then
				return H;
			fi;
		od;
		return fail;
	end;
	if Alphabet(G) = [1,2] and HasName(G) and InStart(Name(G)) then
		CP_init := GetStart(Name(G));
		BS := CP_init.Branchstructure;
		B := List(BS.group);
		Start := CP_init.inital_con;
	else
		return fail;
	fi;
	#if g = h then
	#	return [One(g)];
	#fi;
	Gen := GeneratorsOfGroup(G);
	rek_count := 1;
	#First for binary alphabet.
	Conjugators_branch_rek := function(g,h)
		local L1,L2,c,CT,i,g1K,h1K;
		if not HasName(g) then
			SetName(g,Concatenation("g_",String(rek_count)));
		fi;
		if not HasName(h) then
			SetName(h,Concatenation("h_",String(rek_count)));
		fi; 
		rek_count := rek_count +1;
		Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"");
		if IsOne(g) or IsOne(h) then
			if g = h then
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     g=h=1 So return B");
				return CP_init.RepSystem;
			else
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     g,h is One but the other not. So return []");
				return [];
			fi;
		fi;
		if Activity(g) <> Activity(h) then
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Act(g)<>Act(h). So return []");
			return [];
		fi;
		if g in Gen and h in Gen then
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     g,h are Generators. So return Start[g][h]");
			return Start[Position(Gen,g)][Position(Gen,h)];
		fi;

		
		if IsOne(Activity(g)) then
		#----------Conjugator with trivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator for elm with trivial activity...");
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Try a Conjugator with trivial activity...");
			L1 := Conjugators_branch_rek(State(g,1),State(h,1));
			CT := [];
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     L1 calculated...");
			if Size(L1) > 0 then
				L2 := Conjugators_branch_rek(State(g,2),State(h,2));
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     L2 calculated...");
				for l in [1..Size(B)] do
					for k in [1..Size(B)] do
						if IsBound(L1[l]) and IsBound(L2[k]) then
							c := (B[l]^Embedding(BS.wreath,1) * B[k]^Embedding(BS.wreath,2))^BS.epi;
							if c <> fail then
								CT[Position(B,c)] := FRElement([[[L1[l]],[L2[k]]]],[()],[1]);
							fi;
						fi;
					od;
				od;
				if Size(CT)>0 then
					Print("Something found..\n");
					Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator with trivial activity found. Return it");
					return CT;
				fi;
			fi;
		#----------Conjugator with nontrivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Try a Conjugator with non-trivial activity...");
			L1 := Conjugators_branch_rek(State(g,1),State(h,2));
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     L1 calculated...");
			if Size(L1) > 0 then
				L2 := Conjugators_branch_rek(State(g,2),State(h,1));
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     L2 calculated...");
				for l in [1..Size(B)] do
					for k in [1..Size(B)] do
						if IsBound(L1[l]) and IsBound(L2[k]) then
							c := (B[l]^Embedding(BS.wreath,1) * B[k]^Embedding(BS.wreath,2))^BS.epi;
							if c <> fail then
								CT[Position(B,c)] := FRElement([[[L1[l]],[L2[k]]]],[(1,2)],[1]);
							fi;
						fi;
					od;
				od;
			fi;
			if Size(CT) > 0 then
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator with non-trivial activity found. Return it");
			else
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     No Conjugator found. Return []");
			fi;
			return CT;	
		else
		#----------Conjugator with trivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator for elm with non-trivial activity...");
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Try a Conjugator with trivial activity...");
			L1 := Conjugators_branch_rek(State(g,1)*State(g,2),State(h,1)*State(h,2));
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     L calculated...");		
			CT := [];
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Computing g_1K and h_1K");
			g1K := (State(g,1)^-1)^BS.quo;
			h1K := (State(h,1))^BS.quo;
			for l in [1..Size(B)] do
				if IsBound(L1[l]) then
					c := (B[l]^Embedding(BS.wreath,1) * (g1K * B[l] * h1K)^Embedding(BS.wreath,2))^BS.epi;
					if c <> fail then
						CT[Position(B,c)] := FRElement([[[L1[l]],[State(g,1)^-1*L1[l]*State(h,1)]]],[()],[1]);
					fi;
				fi;
			od;
			if Size(CT)>0 then
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator with trivial activity found. Return it");
				return CT;
			fi;		
		#----------Conjugator with nontrivial activity-------------------------
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Try a Conjugator with non-trivial activity...");
			L1 := Conjugators_branch_rek(State(g,1)*State(g,2),State(h,2)*State(h,1));
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     L calculated...");
			CT := [];
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Computing h_2K");
			h1K := (State(h,2))^BS.quo; #Change to h2
			for l in [1..Size(B)] do
				if IsBound(L1[l]) then
					c := (B[l]^Embedding(BS.wreath,1) * (g1K * B[l] * h1K)^Embedding(BS.wreath,2))^BS.epi;
					if c <> fail then
						CT[Position(B,c)] := FRElement([[[L1[l]],[State(g,1)^-1*L1[l]*State(h,2)]]],[(1,2)],[1]);
					fi;
				fi;
			od;
			if Size(CT) > 0 then
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator with non-trivial activity found. Return it");
			else
				Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     No Conjugator found. Return []");
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
		Print("For Branch Method...\n");
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
	local G,a,b,c,d,x,S,SL,y,C,B,BS;
	G:= GrigorchukGroup;
	a:= G.1; b:=G.2; c:=G.3; d:=G.4;
	BS := BranchStructure(GrigorchukGroup);
	B := List(BS.group);
	S := [];
	for x in GeneratorsOfGroup(G) do
		SL := [];
		for y in GeneratorsOfGroup(G) do
			if x = y then
				if x=a then
					C:= [];
					C[Position(B,One(BS.group))] := One(a);
					C[Position(B,a^BS.quo)] := a;
					C[Position(B,(a*d*a*d)^BS.quo)] := a*d*a*d;
					C[Position(B,(d*a*d)^BS.quo)] := d*a*d;
					Add(SL,C);
				fi;
				if x = d then
					C:= [];
					C[Position(B,One(BS.group))] := One(a);
					C[Position(B,(a*d*a)^BS.quo)] := a*d*a;
					C[Position(B,(a*d*a*d)^BS.quo)] := a*d*a*d;
					C[Position(B,d^BS.quo)] := d;
					C[Position(B,b^BS.quo)] := b;
					C[Position(B,(b*a*d*a)^BS.quo)] := b*a*d*a;
					C[Position(B,(b*a*d*a*d)^BS.quo)] := b*a*d*a*d;
					C[Position(B,c^BS.quo)] := c;
					Add(SL,C);
				fi;
				if x in [b,c] then
					C:=[];
					C[Position(B,One(BS.group))] := One(a);
					C[Position(B,d^BS.quo)] := d;
					C[Position(B,b^BS.quo)] := b;
					C[Position(B,c^BS.quo)] := c;					
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


	
	
	
	
	

