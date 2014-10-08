CONJUGATORS_BRANCH := function(G,g,h)
	local CP_init, Start, B, BS, saved_quo, quo, Alph, Conjugators_branch_rek,l,k,InStart,GetStart,rek_count;
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
	if HasName(G) and InStart(Name(G)) then
		CP_init := GetStart(Name(G));
		BS := CP_init.Branchstructure;
		B := List(BS.group);
		Start := CP_init.inital_con;
		saved_quo := [];
		quo := function(elm) #Calculate only if asked for.
			local s,q;
			for s in saved_quo do
				if s.selm = elm then
					Info(InfoFRCP,4,"Saved work.");
					return s.squo;
				fi;
			od;
			Info(InfoFRCP,4,"Computing elm^BS.quo. May take some time...");
			q := elm^BS.quo;
			Info(InfoFRCP,4,"Finished");
			Add(saved_quo,rec(selm := elm, squo:=q));
			return q;
		end;
	else
		return fail;
	fi;
	if g = h then
		return [One(g)];
	fi;
	Alph := AlphabetOfFRSemigroup(G);
	rek_count := 1;
	Conjugators_branch_rek := function(g,h)
		local L,LC,C,orbits,orb_repr,p,L_Pos,dep,L_PosC,Pos_Con,c,CT,Con,i,j;
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
		if g in CP_init.inital_set and h in CP_init.inital_set then
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     g,h are inital. So return Start[g][h]");
			return Start[Position(CP_init.inital_set,g)][Position(CP_init.inital_set,h)];
		fi;
		orbits := List(Orbits(Group(g),Alph),SortedList);
		orb_repr := List(orbits,Minimum);
		CT := []; # Resulting Conjugator Tuple
		Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Orbit: ",orbits);
		for p in LEVEL_PERM_CONJ@(g,h,BS.top) do
			Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Try a conjugator with activity ",p);
			L := [];
			L_Pos := []; #Stores the position at which the conjugator tuples are defined.
			dep := []; #Stores the dependencies
			for i in [1..Length(orb_repr)] do
				C := Conjugators_branch_rek(State(g^Length(orbits[i]),orb_repr[i]),State(h^Length(orbits[i]),orb_repr[i]^p));
				if Length(C)=0 then #not a valid conjugator
					L:=[];
					break;
				fi;
				for j in [0..Length(orbits[i])-1] do
					LC := [];
					L_PosC := [];
					for k in [1..Length(C)] do
						if IsBound(C[k]) then
							LC[k] := [State(g^j,orb_repr[i])^-1,C[k],State(h^j,orb_repr[i]^p)];
							L_PosC[k] := k;
						fi;
					od;
					L[orb_repr[i]^(g^j)]:=LC ; 
					L_Pos[orb_repr[i]^(g^j)]:=L_PosC;	
				od;
				Add(dep,orbits[i]);
			od; 
			if Size(L)>0 then
				Con := DEP_CARTESIAN@(L,dep);
				Pos_Con := DEP_CARTESIAN@(L_Pos,dep);
				for i in [1..Size(Con)] do #Now possable Conjugators.
					c:= Product([1..Size(Pos_Con[i])],x->(quo(Con[i][x][1])*B[Pos_Con[i][x]]*quo(Con[i][x][3]))^Embedding(BS.wreath,x));
					c:= (c*p^Embedding(BS.wreath,Size(Alph)+1))^BS.epi;;	
					if c <> fail then #Con is a valid element with representative c;
						Info(InfoFRCP,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator found. Add to conjugator tuple ");
						CT[Position(B,c)] := FRElement([Con[i]],[p],[1]); #Will just be multiplied so no reason for mealy.
					fi;
				od;	
			fi;
		od;
		return CT;				
	end;
	return Conjugators_branch_rek(g,h);
end;
#############################Example##############################
GuptaSidkiConjugateBranchInit := function()
	local G,a,t,x,S,SL,y,C,B,BS;
	G:= GuptaSidkiGroup;
	a:= G.1; t:=G.2;
	BS := BranchStructure(G);
	B := List(BS.group);
	S := [];
	for x in [a,a^2,t,t^2] do
		SL := [];
		for y in [a,a^2,t,t^2] do
			if x = y then
				if x in [a,a^2] then
					C:= [];
					C[Position(B,One(BS.group))] := One(a);
					C[Position(B,a^BS.quo)] := a;
					C[Position(B,(a^2)^BS.quo)] := a^2;
					Add(SL,C);
				fi;
				if x in [t,t^2] then
					C:= [];
					C[Position(B,One(BS.group))] := One(a);
					C[Position(B,t^BS.quo)] := t;
					C[Position(B,(t^2)^BS.quo)] := t^2;
					Add(SL,C);
				fi;
			else
				Add(SL,[]);
			fi;
		od;
		Add(S,SL);
	od;
	return S;
end;
	
	
	
	
	

