#############################################################################
##
#W cp.gi                                                       Thorsten Groth
##
#Y Copyright (C) 2014, Thorsten Groth
##
#############################################################################
##
##  This file implements the conjugacy problem for branch groups
##
#############################################################################

#---------------------------------------------------------------
#------              Dep-Cartesian          --------------------
#--  Calculates a cartesian product of ordered lists with    ---
#-- respect to dependencies which entries belong together    ---
#-- Example: L=[[A,B,X],[C,D],[c,d]], dep=[[1],[2,3]] results --
#--   in [[A,C,c],[A,D,d],[B,C,c],[B,D,d],[X,C,c],[X,D,d]	   ---
#--     The Lists, which are joined by the dependencies      ---
#--            have to be of the same length	               ---	
#---------------------------------------------------------------
BindGlobal("DEP_CARTESIAN@", function(L,dep)
	local res_list, temp_cart, container, al, d, i ,j,a;
	res_list := [];
	temp_cart := [];
	for d in dep do
		container := [];
		for j in [1..Size(L[d[1]])] do
			al := [];
			for i in [1..Size(d)] do
				if IsBound(L[d[i]][j]) then
					al[i]:=L[d[i]][j];
				fi;
			od;
			if al <> [] then container[j]:=al; fi;
		od; 
		Add(temp_cart,container);
	od; 
	temp_cart := Cartesian(temp_cart);
	for i in temp_cart do
		container := [];
		for d in i do
			Append(container,d);
		od;
		Add(res_list,container);
	od;
	return res_list;
end);

#--------------------------------------------------------------
#------             LEVEL_PERM_CONJ                     -------
#------  Takes two FRElements and computes a list of    -------
#------  conjugators of the action on the first level.  -------
#--------------------------------------------------------------
BindGlobal("LEVEL_PERM_CONJ@", function(arg)
	local G, pi_x, pi_y, c;
	if Length(arg) < 3 then
		G:=SymmetricGroup(AlphabetOfFRObject(arg[1]));
	else 
		if not IsFRObject(arg[1]) or not IsFRObject(arg[2]) or not IsPermGroup(arg[3]) then
			Error("Usage: FRelm, FRelm, [PermGroup]");
		fi;
		G:=arg[3];
	fi;
	pi_x := PermList(DecompositionOfFRElement(arg[1])[2]);
 	pi_y := PermList(DecompositionOfFRElement(arg[2])[2]);
 	c := RepresentativeAction(G,pi_x,pi_y);
 	if c= fail then
 		return [];
 	fi;
 	return c*List(Centralizer(G,pi_y));
end);
##################################################################
#````````````````````````````````````````````````````````````````#
#```````````````````    OrbitSignalizer   ```````````````````````#
#````````````````                            ````````````````````#
#````````````````    Guaranteed to stop on    ```````````````````#
#```````````````` BoundedFRElements as input ````````````````````#
#``````  Computes {a^m@v|v∊X*} with m = |Orb_a(v)|   ````````````#
#````````````````````````````````````````````````````````````````#
##################################################################
InstallMethod(OrbitSignalizer,
	"Returns the finite Orbit Signalizer",
	[IsFRElement],
function(a)
	local OS_list,i,OS_unvisited,OS_new,elm,x,new,suc;
	suc := function(state,x)
		return  State(state^Size(ForwardOrbit(state,[x])),[x]);
	end;
	OS_list := [];
	OS_unvisited := [a];
	while Length(OS_unvisited) > 0 do
		OS_new := [];
		for elm in OS_unvisited do
			for x in AlphabetOfFRObject(a) do
				new := suc(elm,x);
				if (not new in OS_list) and (not new in OS_unvisited) and (not new in OS_new) then
					Add(OS_new,new);
				fi;
			od;
			
		od; 
		Append(OS_list,OS_unvisited);
		OS_unvisited := OS_new;
	od;
	return OS_list;
end
);
##################################################################
#````````````````````````````````````````````````````````````````#
#`````````````````````                  `````````````````````````#
#`````````````````````  ConjugatorGraph `````````````````````````#
#`````````````````````     DrawGraph    `````````````````````````#
#`````````````````````                  `````````````````````````#
#````````````````````````````````````````````````````````````````#
##################################################################
BindGlobal("CONJUGATOR_GRAPH@", function(a,b)
	local Alph, Vertices, Edges, c, d, p, v_id, e_id, v, orbits, orb_repr, i, new_con_pair, new_v, w, change, found, e, all_found;
	
	Alph := AlphabetOfFRObject(a);
	Vertices := [];
	Edges := [];
	#--------------------- Save some work, in easy cases--------
	if Size(LEVEL_PERM_CONJ@(a,b))=0 then
		return [[],[]];
	fi;
	#--------------------- Generate the Vertex list ------------
	v_id := 1;
	for c in OrbitSignalizer(a) do
		for d in OrbitSignalizer(b) do
			for p in LEVEL_PERM_CONJ@(c,d) do
				Add(Vertices,rec(	id:= v_id,
													conj_pair := [c,d],
													action := p));
				v_id := v_id+1;
			od;
		od;
	od;
	#Print("Vertexlist generated\n");
	
	#--------------------- Find the Edges  -------------------
	e_id := 1;
	for v in Vertices do
		c := v.conj_pair[1];
		d := v.conj_pair[2];
		orbits := Orbits(Group(c),Alph);
		orb_repr := List(orbits,Minimum);
		for i in [1..Length(orbits)] do
			new_con_pair := [State(c^Length(orbits[i]),orb_repr[i]),State(d^Length(orbits[i]),orb_repr[i]^v.action)];
			for p in LEVEL_PERM_CONJ@(new_con_pair[1],new_con_pair[2]) do
				new_v := 0;
				for w in Vertices do
					if w.conj_pair = new_con_pair and w.action = p then;
						new_v := w;
						break;
					fi;
				od;
				if new_v <> 0 then
					#Print("Add Edge from ",v.id," to ",new_v.id," along ",orb_repr[i],"\n");
					Add(Edges,rec(	from:=v.id,
													to := new_v.id,
													read := orb_repr[i],
													write := orb_repr[i]^v.action,
													id := e_id));
					e_id := e_id +1;	
				else #This case should never happen...
					Error("Error the element is not in the vertex set!\n");
				fi;
			od;
		od;
	od;
	#Print(Size(Edges)," Edges generated\n");
	
	#--------------------- Delete dead Vertices  -------------------
	change:=true;
	while change do
		change := false;
		for v in Vertices do
			orbits := Orbits(Group(v.conj_pair[1]),Alph);
			orb_repr := List(orbits,Minimum);
			found := [];
			for e in Edges do
				if e.from = v.id then
					Add(found,e.read);
				fi;
			od;
			#Are all outgoing edges there?
			all_found := true;
			for i in orb_repr do
				if not i in found then
					all_found := false;
					#Print("At Vertex ",v.id," no Edge ",i," was found. So remove it.\n");
					break;
				fi;
			od;
			if not all_found then
				Unbind(Vertices[v.id]);
				#Delete all Edges from or to the removed vertex
				for e in Edges do
					if e.from = v.id or e.to = v.id then
						Unbind(Edges[e.id]);
					fi;
				od;
				change:=true;
			fi;
		od;
	od;
	#Print("Dead Vertices removed\n");
	return [Vertices,Edges];
end);

BindGlobal("DRAW_GRAPH@",function(Vertices,Edges)
	local v,S,e;
	for v in Vertices do
		Print("ID: ",v.id," Name: (",v.conj_pair[1]![2],",",v.conj_pair[2]![2],")\n");
	od;
	#Draw Conjugacy Graph
	S := "digraph finite_state_machine {\n";
	for e in Edges do
		#Print(Vertices[e.from].conj_pair[1]![2]);
		Append(S,"\"");
		Append(S,String(Vertices[e.from].id));
		Append(S,": (");
		Append(S,String(Vertices[e.from].conj_pair[1]![2]));
		Append(S,",");
		Append(S,String(Vertices[e.from].conj_pair[2]![2]));
		Append(S,")");
		Append(S,String(Vertices[e.from].action));
		Append(S,"\"");
		Append(S," -> ");
		Append(S,"\"");
		Append(S,String(Vertices[e.to].id));
		Append(S,": (");
		Append(S,String(Vertices[e.to].conj_pair[1]![2]));
		Append(S,",");
		Append(S,String(Vertices[e.to].conj_pair[2]![2]));
		Append(S,")");
		Append(S,String(Vertices[e.to].action));
		Append(S,"\"");
		Append(S," [label=\"");
		Append(S,String(e.read));
		Append(S,"|");
		Append(S,String(e.write));
		Append(S,"\"];\n");
	od;
	Append(S,"}");
	Print(S);
	DOT2DISPLAY@(S,"dot");
end);
##################################################################
#````````````````````````````````````````````````````````````````#
#```````````````````  MEALY_FROM_STATES@  ```````````````````````#
#````````````````                            ````````````````````#
#`````````````````  Computes a mealy machine ````````````````````#
#````````````````     with Statesset L and   ````````````````````#
#````````````````        activity act        ````````````````````#
#````````````````````````````````````````````````````````````````#
##################################################################
#Computes a mealy machine with stateset L and given activity act.
BindGlobal("MEALY_FROM_STATES@", function(L,act)
		local m,tran,out,i,j;
		#Force L to contain elements and not lists of elements.
		for i in [1..Size(L)] do
			if IsList(L[i]) then
				L[i] := Product(L[i]);
			fi;
		od;
		#Force act to be a list of output symbols.
		if IsPerm(act) then
			act := List(AlphabetOfFRObject(L[1]),x->x^act);
		fi;
		
		tran := [[]];
		out := [act];
		i := 1;
		j := 2;
		for m in L do
		  Add(tran[1],i+1);
		  Append(tran,List(m!.transitions,x->List(x,y->y+i)));
		  Append(out,m!.output);
		  i := i + Length(m!.transitions);
		  j := j+1;
		od;
		return MealyElement(tran,out,1);
	end);
##################################################################
#````````````````````````````````````````````````````````````````#
#`````````````````````                  `````````````````````````#
#`````````````````````    F.S. Worker   `````````````````````````#
#`````````````````````                  `````````````````````````#
#````````````````````````````````````````````````````````````````#
##################################################################
BindGlobal("CONJUGATORS_FINITE_STATE_WRAPPER@",function(start,CG)
	local v,AS,to_visit, Alph, new_v, i, found, e, Tran, Act, c,d, orbit;
			#--------- Choose one subgraph, as automaton  ---------
			AS := [start.id]; #Contains IDs of vertices, which build the subgraph
			to_visit := [start.id]; 
			Alph := AlphabetOfFRObject(start.conj_pair[1]);
			while Length(to_visit) > 0 do 
				new_v := [];
				for i in to_visit do
					v := CG[1][i];
					found := [];
					for e in CG[2] do
						if e.from = v.id then
							if e.read in found then
								Unbind(CG[2][e.id]);
							else 
								Add(found,e.read);
								if not e.to in AS then
									Add(new_v,CG[1][e.to].id);
									Add(AS,e.to);
								fi;
							fi;
						fi;
					od;
				od;
				to_visit := new_v;
			od;
			#Form an automaton out of the subgraph
			Tran := [];
			Act := [];
			for i in AS do
				Add(Tran,[]);
				Add(Act,CG[1][i].action);
			od;
			for e in CG[2] do
				if e.from in AS and e.to in AS then
					Tran[Position(AS,e.from)][Position(Alph,e.read)] := [Position(AS,e.to)];
					c := CG[1][e.from].conj_pair[1];
					d := CG[1][e.from].conj_pair[2];
					orbit := ForwardOrbit(c,e.read);
					for i in [2..Length(orbit)] do
					#The missing edges...
						Tran[Position(AS,e.from)][Position(Alph,orbit[i])] := [State(c^(i-1),e.read)^(-1),Position(AS,e.to),State(d^(i-1),e.read^(CG[1][e.from].action))];
					od;
				fi;
			od;
			return FRElement(Tran,Act,[1]);
end);
##################################################################
#````````````````````````````````````````````````````````````````#
#`````````````````````                  `````````````````````````#
#`````````````````````  Finitary Worker `````````````````````````#
#`````````````````````                  `````````````````````````#
#````````````````````````````````````````````````````````````````#
##################################################################
BindGlobal("CONJUGATORS_FINITARY_WRAPPER@",function(v,Graph,Seen,Known_vertex_conjugator)
	local CONJUGATORS_FINITARY_REK;
	CONJUGATORS_FINITARY_REK := function(v,Graph,Seen,Known_vertex_conjugator)
		local Vertices,Edges,sons,starts,conj_cand,conjugators_found,son,x,NewSeen,a,b,e,son_conj,son_conjs,tempo_conj, htemp,tempoconj,err,Indices,diction,pos,real_conj,real_conjugators,i,j,ip,son_orbit_size,Alph,action,
		w,Circle,Conjs,con;
	
		#Print("@Vertex ",v,"\n");
		if IsBound(Known_vertex_conjugator[v]) then #Don't do the same work twice.
			return Known_vertex_conjugator[v];
		fi;
		Vertices := Graph[1];
		Edges := Graph[2];
		a := Vertices[v].conj_pair[1];
		b := Vertices[v].conj_pair[2];
		Alph := AlphabetOfFRObject(a);
		action := Vertices[v].action;
	
		if v in Seen then
			Circle:=Seen{[Position(Seen,v)..Size(Seen)]}; #are all one then.
			for w in Circle do
				if not Vertices[w].action = () then
					return [];
				fi;
				if not Vertices[w].conj_pair[1] = Vertices[w].conj_pair[2] then
					return [];
				fi;
			od;
			return [One(a)];
		fi;
	
		sons := [];
		for e in Edges do
			if e.from = v then
				Add(sons,[e.to,e.read]);
			fi;
		od;
		real_conjugators := [];
		conj_cand := EmptyPlist(Size(Alph));
		for x in Alph do
			conj_cand[x] := [];
		od;
		conjugators_found := [];
	
		for son in sons do
			NewSeen := ShallowCopy(Seen);
			Add(NewSeen,v);
			son_conjs := CONJUGATORS_FINITARY_REK(son[1],Graph,NewSeen,Known_vertex_conjugator);
			son_orbit_size := Size(Orbit(Group(a),son[2]));
			for son_conj in son_conjs do
				tempo_conj := [];
				err:=0;
				for j in [1..son_orbit_size-1] do
					htemp := (State(a^j,son[2]))^(-1) * son_conj * State(b^j,son[2]^action);
					if IsFinitaryFRElement(htemp) then
						tempo_conj[son[2]^(a^j)] := htemp;
					else
						err := 1;
						break;
					fi;
				od;
				tempo_conj[son[2]] := son_conj;
				if err = 0 then #tempo_conj is indeed a valid partial_conjugator for (a,b)
					for i in Alph do
						if IsBound(tempo_conj[i]) then
							Add(conj_cand[i],[tempo_conj[i]]);
							conjugators_found[i] := 1;
						fi;
					od;
				fi;
			od;
		od;
		#Test if we have enough partial conjugators
		if IsDenseList(conjugators_found) and Size(conjugators_found) = Size(Alph) then
			#puzzle them together!	
			Conjs := DEP_CARTESIAN@(conj_cand,Orbits(Group(a),Alph));
			for con in Conjs do
				Add(real_conjugators,FRElement([con],[action],[1]));
			od;
		fi;
		Known_vertex_conjugator[v] := real_conjugators;
		return real_conjugators;
	end;
	return CONJUGATORS_FINITARY_REK(v,Graph,Seen,Known_vertex_conjugator);
end);
##################################################################
#````````````````````````````````````````````````````````````````#
#`````````````````````                  `````````````````````````#
#`````````````````````  BoundedWorker   `````````````````````````#
#`````````````````````                  `````````````````````````#
#````````````````````````````````````````````````````````````````#
##################################################################
BindGlobal("CONJUGATORS_BOUNDED_WRAPPER@",function(v,Graph,Seen,readwrite_path,Known_vertex_conjugator)
	local CONJUGATORS_BOUNDED_REK;
	CONJUGATORS_BOUNDED_REK := function(v,Graph,Seen,readwrite_path,Known_vertex_conjugator)
		local Vertices,Edges,sons,starts,conj_cand,conjugators_found,son,x,NewSeen,a,b,e,son_conj,son_conjs,tempo_conj, htemp,tempoconj,err,Indices,diction,pos,real_conj,real_conjugators,conj_cand_aux,i,j,ip,son_orbit_size,Alph,Alph_num,action, alph,beta,Conj_elm,read,write,orb_size,Conj_Tran,Conj_act,m,read_path,write_path,action_path,New_read_path, New_write_path,New_action_path,Conj_Tran_el,Conjs,circle_length,New_Seen,X,con, check_need ;
	
		#Print("@Vertex ",v,"\n");
		if IsBound(Known_vertex_conjugator[v]) then #Don't do the same work twice.
			return Known_vertex_conjugator[v];
		fi;
		#TODO REKURSIONSABBRUCH, wenn start schon gefunden wurde!!!!
		#if IsBound(Known_vertex_conjugator[start]) then
		#	if Size(Known_vertex_conjugator[start]>0) then
		#		Print("Warum denn noch weitermachen... ist doch schon alles klar...\n\n");
		#	fi;
		#fi;
		Vertices := Graph[1];
		Edges := Graph[2];
		a := Vertices[v].conj_pair[1];
		b := Vertices[v].conj_pair[2];
		Alph := AlphabetOfFRObject(a);
		Alph_num := [1..Size(Alph)];
		action := Vertices[v].action;
		read_path := readwrite_path[1];
		write_path := readwrite_path[2];
		action_path := readwrite_path[3];
		if v in Seen then
			m := Size(Orbit(Group(a),read_path));

			Conj_Tran:= [];
			Conj_act := [];
			circle_length:=Size(Seen)-Position(Seen,v)+1;
			check_need := false;
			for i in [Position(Seen,v)..Size(Seen)] do
				alph := Vertices[Seen[i]].conj_pair[1];
				beta := Vertices[Seen[i]].conj_pair[2];
				read:=read_path[i];
				write:=write_path[i];
				orb_size:=	Size(Orbit(Group(alph),read));
				Conj_elm :=[];
				for x in Alph_num do
					Conj_elm[x]:=[];
				od;
				if not orb_size = Size(Alph) then		#Here some extra work, search for a conjugator for the non determined states, there may be more than one.
					x:= read;
					#Get a representative for the orbits
					X:=Difference(Alph,Orbit(Group(alph),x));
					while Size(X)>0 do 
						x := Minimum(X);				
						sons := [];
						for e in Edges do
							if e.from = v and e.read = x then
								Add(sons,e.to);
							fi;
						od;
						New_Seen := Seen{[1..i]};
						son_conjs:= [];
						#Here there is already one circle, so the other states have to be finitary
						for son in sons do
								Conjs :=CONJUGATORS_FINITARY_WRAPPER@(son,Graph,New_Seen,Known_vertex_conjugator);
								Append(son_conjs,Conjs);
						od;		
						for son_conj in son_conjs do
							for j in [0..Size(Orbit(Group(alph),x))-1] do 
								Add(Conj_elm[Position(Alph,x^(alph^j))],[(State(alph^j,x))^-1*son_conj*State(beta^j,x^action_path[i])]); 
							od;
						od;
						X:=Difference(X,Orbit(Group(alph),x)); #next orbit
					od;
				fi;
				Conj_elm[Position(Alph,read)]:= [[(i mod circle_length) +1]];
				for j in [1..orb_size-1] do
					Conj_elm[Position(Alph,read^(alph^j))] := [[(State(alph^j,read))^-1,Conj_elm[Position(Alph,read)][1][1],State(beta^j,write)]];
					#This may be not bounded, so check it later.
					check_need := true;
				od;
				#Remove duplicates in
				for x in Alph_num do
					Conj_elm[x] := Set(Conj_elm[x]);
				od;
				Conj_elm := DEP_CARTESIAN@(Conj_elm,Orbits(Group(alph),Alph)); #puzzle the conjugators together
				Add(Conj_act,action_path[i]);
				Add(Conj_Tran,Conj_elm);
			od;
			Conjs := [];
			for Conj_Tran_el in Cartesian(Conj_Tran) do
				conj_cand := FRElement(Conj_Tran_el,Conj_act,[1]);
				if check_need then
					if IsBoundedFRElement(conj_cand) then
						check_need := false;
					fi;
				fi;
				if not check_need then
					if IsBound(Known_vertex_conjugator[v]) then
						Add(Known_vertex_conjugator[v],conj_cand);
					else
						Known_vertex_conjugator[v] := [conj_cand];
					fi;
					Add(Conjs,conj_cand);
				fi;
			od;
			return Conjs;

		fi;
		sons := [];
		for e in Edges do
			if e.from = v then
				Add(sons,[e.to,e.read]);
			fi;
		od;
		real_conjugators := [];
		conj_cand := EmptyPlist(Size(Alph));
		for x in Alph do
			conj_cand[x] := [];
		od;
		conjugators_found := [];
	
		for son in sons do
			NewSeen := ShallowCopy(Seen);
			Add(NewSeen,v);
			New_read_path := ShallowCopy(read_path);
			Add(New_read_path,son[2]);
			New_write_path := ShallowCopy(write_path);
			Add(New_write_path,son[2]^action);
			New_action_path := ShallowCopy(action_path);
			Add(New_action_path,action);
				
			son_conjs := CONJUGATORS_BOUNDED_REK(son[1],Graph,NewSeen,[New_read_path,New_write_path,New_action_path],Known_vertex_conjugator);
			for son_conj in son_conjs do
			son_orbit_size := Size(Orbit(Group(a),son[2]));
				for j in [0..son_orbit_size-1] do
					Add(conj_cand[son[2]^(a^j)],[(State(a^j,son[2]))^(-1) * son_conj * State(b^j,son[2]^action)]);
					conjugators_found[son[2]^(a^j)] := 1;
				od;
			od;
		od;
		#Test if we have enough partial conjugators
		if IsDenseList(conjugators_found) and Size(conjugators_found) = Size(Alph) then
			#puzzle them together!	
			Conjs := DEP_CARTESIAN@(conj_cand,Orbits(Group(a),Alph));
			for con in Conjs do
				Add(real_conjugators,FRElement([con],[action],[1]));
			od;
		fi;
		Known_vertex_conjugator[v] := real_conjugators;
		return real_conjugators;
	end;
	return CONJUGATORS_BOUNDED_REK(v,Graph,Seen,readwrite_path,Known_vertex_conjugator);
end);
##################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%      IsConjugate        %%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%	 RepresentativeActionOp %%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##################################################################
InstallMethod(IsConjugate,
	"For Aut, RAut, FAut, Poly-1, Poly0",
	#The attribute FullSCVertex charakterizes all FullSCGroups
	[ IsFRGroup and HasFullSCVertex,IsFRElement,IsFRElement], 
  function(G,a,b)
  	local v, Graph, sons, starts;
  	if AlphabetOfFRSemigroup(G) <> AlphabetOfFRObject(a) or AlphabetOfFRSemigroup(G) <> AlphabetOfFRObject(b) then
  		return false;
  	fi;
  	if a = b then #Spare Computing Time in trivial case.
  	 return true; 
  	fi; 
  	Graph := CONJUGATOR_GRAPH@(a,b);
  	if FullSCFilter(G) = IsFRElement or FullSCFilter(G) = IsFiniteStateFRElement then
		#In this cases the conjugacy problems are equivalent,
		#----------------------FiniteState-----------------;
		#------------------FunctionalRecursive-------------;
			for v in Graph[1] do
				if v.conj_pair = [a,b] then;
					return true;
				fi;
			od;
			return false;
		elif FullSCFilter(G) = IsFinitaryFRElement then
		#----------------------Finitary--------------------;
			starts := [];
			for v in Graph[1] do
				if v.conj_pair = [a,b] then
					Add(starts,v.id);
				fi;
			od;
			for v in starts do;
				if Size(CONJUGATORS_FINITARY_WRAPPER@(v,Graph,[],[]))>0 then
					return true;
				fi;
			od;
			return false;	
		elif FullSCFilter(G) = IsBoundedFRElement then
		#----------------------Bounded---------------------;		
			starts := [];
			for v in Graph[1] do
				if v.conj_pair = [a,b] then
					Add(starts,v.id);
				fi;
			od;
			for v in starts do
				if Size(CONJUGATORS_BOUNDED_WRAPPER@(v,Graph,[],[[],[],[]],[]))>0 then
					return true;
				fi;
			od;
			return false;
		else
		#----------------------Else------------------------;	
			TryNextMethod();
		fi;
  end
);

InstallOtherMethod(RepresentativeActionOp,
	"Computes a conjugator in the given FullSCGroup ",
	#The attribute FullSCVertex charakterizes all FullSCGroups
	[ IsFRGroup and HasFullSCVertex,IsFRElement,IsFRElement], 
	function(G,a,b)
  	local CG, v, start, Conjugators;
  	if AlphabetOfFRSemigroup(G) <> AlphabetOfFRObject(a) or AlphabetOfFRSemigroup(G) <> AlphabetOfFRObject(b) then
  		return fail;
  	fi;
  	if a=b then
  	 return One(G); 
  	fi;
  	CG := CONJUGATOR_GRAPH@(a,b);
  	if FullSCFilter(G) = IsFRElement or FullSCFilter(G) = IsFiniteStateFRElement then
		#In this cases the conjugacy problems are equivalent,
		#----------------------FiniteState-----------------;
		#------------------FunctionalRecursive-------------;
			for v in CG[1] do
				if v.conj_pair = [a,b] then;
					start := v;
					break;
				fi;
			od;
			if not IsBound(start) then
				return fail;
			fi;
			return CONJUGATORS_FINITE_STATE_WRAPPER@(start,CG);
		elif FullSCFilter(G) = IsFinitaryFRElement then
		#----------------------Finitary--------------------;
			start := [];
			for v in CG[1] do
				if v.conj_pair = [a,b] then
					Add(start,v.id);
				fi;
			od;
			for v in start do
				Conjugators :=CONJUGATORS_FINITARY_WRAPPER@(v,CG,[],[]);
				if Size(Conjugators)>0 then
					return Conjugators[1];
				fi;
			od;
			return fail;		
		elif FullSCFilter(G) = IsBoundedFRElement then
		#----------------------Bounded---------------------;		
			start := [];
			for v in CG[1] do
				if v.conj_pair = [a,b] then
					Add(start,v.id);
				fi;
			od;
			for v in start do
				Conjugators :=CONJUGATORS_BOUNDED_WRAPPER@(v,CG,[],[[],[],[]],[]);
				if Size(Conjugators)>0 then
					return Conjugators[1];
				fi;
			od;
			return fail;
		else
		#----------------------Else------------------------;	
			TryNextMethod();
		fi;
  end);
#****************************************************************
################################################################*
################################################################*
###############                               ##################*
###############  Algorithm for branch groups  ##################*
###############                               ##################*
################################################################*
################################################################* 
#****************************************************************


#---------------------------------------------------------------
#------      InitConjugateForBranchGroups      -----------------
#--   Sets the Precomputed initial data for the branch			 ---
#--  Algorithm. Stores this data for later computations.     ---  	
#---------------------------------------------------------------

InstallMethod(FRBranchGroupConjugacyData,
	[ IsFRGroup ],	
	 function(G)
	 	local init, N, g, h, b, CT, c, i;
	 	Info(InfoFR, 1, "Init FRBranchGroupConjugacyData");
		init := rec(initial_conj_dic:=NewDictionary([One(G),One(G)],true),
								Branchstructure:=BranchStructure(G),
								RepSystem:=List(~.Branchstructure.group,x->PreImagesRepresentativeNC(~.Branchstructure.quo,x)));
		N := TORSIONNUCLEUS@(G);
		if N = fail then return fail;fi;
		SEARCH@.INIT(G);
		for g in N do
			for h in N do
				#Find one conjugator b
				repeat 
				 	b := SEARCH@.CONJUGATE(G,g,h);
				 	while b=fail and SEARCH@.EXTEND(G)=fail do
            SEARCH@.ERROR(G,"RepresentativeAction");
        	od;
        	Info(InfoFR, 3, "RepresentativeAction: searching at level ",G!.FRData.level," and in sphere of radius ",G!.FRData.radius);
        until b<>fail;
		    CT := []; #The Conjugator tuple
		    if b <> false then
		    	i := 1;
		    	for c in init.Branchstructure.group do
		    		repeat
        			b := SEARCH@.CONJUGATE_COSET(G,c,g,h);
				    	while b=fail and SEARCH@.EXTEND(G)=fail do
				        SEARCH@.ERROR(G,"RepresentativeAction");
				    	od;
				    until b<>fail;
				    if b <> false then
							CT[i] := b;
						fi;
						i := i+1;
					od;
		    fi;
				AddDictionary(init.initial_conj_dic,[g,h],CT);
			od;
		od;
		Info(InfoFR, 1, "Finished Init FRBranchGroupConjugacyData");
		return init;
	 end);
##################################################################
#````````````````````````````````````````````````````````````````#
#`````````````````````                  `````````````````````````#
#`````````````````````   Branch Worker  `````````````````````````#
#`````````````````````                  `````````````````````````#
#````````````````````````````````````````````````````````````````#
##################################################################	
BindGlobal("CONJUGATORS_BRANCH@",function(G,g,h)
	local CP_init, Start, B, BS, Con_dic, saved_quo, quo, Alph, Conjugators_branch_rek,l,k,rek_count;
	CP_init := FRBranchGroupConjugacyData(G);
	if CP_init = fail then
		return fail;
	fi;
	BS := CP_init.Branchstructure;
	B := List(BS.group);
	Con_dic := CP_init.initial_conj_dic;
	saved_quo := NewDictionary(One(G),true);
	quo := function(elm) #Calculate only if asked for.
		local q;
		if not KnowsDictionary(saved_quo,elm) then
			Info(InfoFR,4,"Computing elm^BS.quo. May take some time...");
			q := elm^BS.quo;
			Info(InfoFR,4,"Finished");
			AddDictionary(saved_quo,elm,q);
			return q;
		fi;
		return LookupDictionary(saved_quo,elm);
	end;
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
		Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"");
		if IsOne(g) or IsOne(h) then
			if g = h then
				Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"     g=h=1 So return B");
				return CP_init.RepSystem;
			else
				Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"     g,h is One but the other not. So return []");
				return [];
			fi;
		fi;
		if KnowsDictionary(Con_dic,[g,h]) then
			Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"     g,h are already known. So return them]");
			return LookupDictionary(Con_dic,[g,h]);
		fi;
		orbits := List(Orbits(Group(g),Alph),SortedList);
		orb_repr := List(orbits,Minimum);
		CT := []; # Resulting Conjugator Tuple
		Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"     Orbit: ",orbits);
		for p in LEVEL_PERM_CONJ@(g,h,BS.top) do
			Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"     Try a conjugator with activity ",p);
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
						Info(InfoFR,3,"Computing g,h=",Name(g),",",Name(h),"     Conjugator found. Add to conjugator tuple ");
						CT[Position(B,c)] := MEALY_FROM_STATES@(Con[i],p);
						#CT[Position(B,c)] := FRElement([Con[i]],[p],[1]);
					fi;
				od;	
			fi;
		od;
		AddDictionary(Con_dic,[g,h],CT); #Save work in case is it again asked for a CT for (g,h).
		return CT;				
	end;
	return Conjugators_branch_rek(g,h);
end);

##################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%      IsConjugate         %%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%	RepresentativeActionOp  %%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##################################################################
InstallOtherMethod(RepresentativeActionOp,
	"Computes a conjugator in the given Branch group ",
	[ IsBranched and IsFinitelyGeneratedGroup,IsFRElement,IsFRElement,IsFunction], 
	function(G,g,h,f)
		local con;
		if f <> OnPoints then TryNextMethod(); fi;
		Info(InfoFR,2,"Try method for branch groups.");
		con := CONJUGATORS_BRANCH@(G,g,h);
		if con <> fail then
			if Size(con)>0 then
				return Representative(con);
			fi;
			return fail;
		fi;
		Info(InfoFR,2,"Doesn't work. Try next...");
		TryNextMethod();
		end);
InstallMethod(IsConjugate,
	"For Branch Groups",
	[ IsBranched and IsFinitelyGeneratedGroup,IsFRElement,IsFRElement], 
  function(G,g,h)
  	local con;
		Info(InfoFR,2,"Try method for branch groups.");
  	con := CONJUGATORS_BRANCH@(G,g,h);
		if con <> fail then
			if Size(con)>0 then
				return true;
			fi;
			return false;
		fi;
		Info(InfoFR,2,"Doesn't work. Try next...");
		TryNextMethod();
		end);	


#############################Example##############################
####								Setting the Branch Data										 ###
####			for GrigorchukGroup and GuptaSidkiGroup							 ###
#SetFRBranchGroupConjugacyData(GrigorchukGroup,
#	 rec(	initial_conj_dic:=NewDictionary([One(GrigorchukGroup),One(GrigorchukGroup)],true),
#				Branchstructure:=BranchStructure(GrigorchukGroup),
#				RepSystem:=List(~.Branchstructure.group,x->PreImagesRepresentativeNC(~.Branchstructure.quo,x)))
#	 );
#CallFuncList(function(a,b,c,d) 
#							local G,D,g,h;
#							G:= GrigorchukGroup;
#							D:= FRBranchGroupConjugacyData(G).initial_conj_dic;
#							for g in [a,b,c,d] do
#								for h in [a,b,c,d] do
#									if g<>h then
#										AddDictionary(D,[g,h],[]);
#									fi;
#								od;
#							od;
#							AddDictionary(D,[a,a],[One(G),a,d*a*d,a*d*a*d]);
#							AddDictionary(D,[b,b],[One(G),,,, b,,,,c,,,,d]);
#							AddDictionary(D,[c,c],[One(G),,,, b,,,,c,,,,d]);
#							AddDictionary(D,[d,d],[One(G),,,a*d*a*d,b,,,b*a*d*a*d,c,,,b*a*d*a,d,,,a*d*a]);
#						end,GeneratorsOfGroup(GrigorchukGroup)
#		);
		
#SetFRBranchGroupConjugacyData(GuptaSidkiGroup,
#	 rec(	initial_conj_dic:=NewDictionary([One(GuptaSidkiGroup),One(GuptaSidkiGroup)],true),
#				Branchstructure:=BranchStructure(GuptaSidkiGroup),
#				RepSystem:=List(~.Branchstructure.group,x->PreImagesRepresentativeNC(~.Branchstructure.quo,x)))
#	 );
#CallFuncList(function(a,t) 
#							local G,D,g,h;
#							G:= GuptaSidkiGroup;
#							D:= FRBranchGroupConjugacyData(G).initial_conj_dic;
#							for g in [a,a^2,t,t^2] do
#								for h in [a,a^2,t,t^2] do
#									if g<>h then
#										AddDictionary(D,[g,h],[]);
#									fi;
#								od;
#							od;
#							AddDictionary(D,[a,a],[One(G),a,a^2]);
#							AddDictionary(D,[a^2,a^2],[One(G),a,a^2]);
#							AddDictionary(D,[t,t],[One(G),,,t,,,t^2]);
#							AddDictionary(D,[t^2,t^2],[One(G),,,t,,,t^2]);
#						end,GeneratorsOfGroup(GuptaSidkiGroup)
# 	);
#****************************************************************
################################################################*
################################################################*
###############                               ##################*
###############     Algorithm for the         ##################*
###############       GrigorchukGroup         ##################*
################################################################*
################################################################* 
#****************************************************************
BindGlobal("GRIG_CON@",function(G,g,h)
local f,gw,hw,Gen,a, b, c, d, Fam, aw, dw, ae, be, ce, de, Alph, x_1, x_2, K_repr, K_repr_words, D, ConTup_a, Check, alternating_a_form, shorten_word, compute_conjugates, compute_conjugates_of_word, L_Decomp, Compute_K_rep, L_word_to_Grig, Merge_Ls, conjugators_grig_rek, Res, r, Join_to_first;

############   Spare Computing Time in trivial case.     #########
 	if AlphabetOfFRSemigroup(G) <> AlphabetOfFRObject(g) or AlphabetOfFRSemigroup(G) <> AlphabetOfFRObject(h) then
		return fail;
	fi;
	if g = h then 
	 return One(G); 
	fi;	
############       (Local) GLOBALS           #####################
	f := EpimorphismFromFreeGroup(G);
	gw:=PreImagesRepresentativeNC(f,g);
	hw:=PreImagesRepresentativeNC(f,h);
	
	Gen := GeneratorsOfGroup(G);
	a:= Position(Gen,MealyElement([[4,2],[4,3],[5,1],[5,5],[5,5]],[(),(),(),(1,2),()],4));	
	b:= Position(Gen,MealyElement([[4,2],[4,3],[5,1],[5,5],[5,5]],[(),(),(),(1,2),()],1));	
	c:= Position(Gen,MealyElement([[4,2],[4,3],[5,1],[5,5],[5,5]],[(),(),(),(1,2),()],2));	
	d:= Position(Gen,MealyElement([[4,2],[4,3],[5,1],[5,5],[5,5]],[(),(),(),(1,2),()],3));	
	Fam := FamilyObj(gw);		
##################################################################
	aw :=AssocWordByLetterRep(Fam,[a]);  
	dw :=AssocWordByLetterRep(Fam,[d]);
	ae := ImageElm(f,AssocWordByLetterRep(Fam,[a]));
	be := ImageElm(f,AssocWordByLetterRep(Fam,[b]));
	ce := ImageElm(f,AssocWordByLetterRep(Fam,[c]));
	de := ImageElm(f,AssocWordByLetterRep(Fam,[d]));
	Alph:=AlphabetOfFRObject(g);
	x_1 := Alph[1];
	x_2 := Alph[2];

#Precomputed K-representatives:
#[[],a,ad,ada,adad,adada,adadada,adadada,b,ba,bad,bada,badad,badada,badadad,badadada]
	K_repr := [[],[a],[a,d],[a,d,a],[a,d,a,d],[a,d,a,d,a],[a,d,a,d,a,d],[a,d,a,d,a,d,a],[b],[b,a],[b,a,d],[b,a,d,a],[b,a,d,a,d],[b,a,d,a,d,a],[b,a,d,a,d,a,d],[b,a,d,a,d,a,d,a]];
	K_repr_words := List(K_repr,x->AssocWordByLetterRep(Fam,x));
	
	#Precomputed words, which decompose to the K_repr.: <K_repr[i]·l,f(K_repr[i])·l'> = D[i]·<l,l'>
	D:= List([[],[c],[c,a,c,a],[c,a,c,a,c],[c,a,c,a,c,a,c,a],[c,a,c,a,c,a,c,a,c],[c,a,c,a,c,a,c,a,c,a,c,a],[c,a,c,a,c,a,c,a,c,a,c,a,c],[a,d,a],[a,d,a,c],[a,d,a,c,a,c,a],[a,d,a,c,a,c,a,c],[a,d,a,c,a,c,a,c,a,c,a],[a,d,a,c,a,c,a,c,a,c,a,c],[a,d,a,c,a,c,a,c,a,c,a,c,a,c,a],[a,d,a,c,a,c,a,c,a,c,a,c,a,c,a,c]],x->AssocWordByLetterRep(Fam,x));
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%        Functions       %%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	
	#TeporaryDebug function to locate possable errors.
	Check := function(s,g,h,C)
		local c;
		if InfoLevel(InfoFR)>2 then
			for c in C do
				if g^c <> h then
					Info(InfoFR,2,"Error at ",s);
					Info(InfoFR,3,"Error happened here: g=",g,", and h=",h,", and Conjugator c=",c," number: ",Position(C,c),"in ",C);
					return fail;
				fi;
			od;
		fi;
		return true;
	end;
	
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%     Magic on words     %%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	#Given a word w in Generators of Grig, computes the form w= (a) x1 a x2 a x3… where xi in b,c,d
	alternating_a_form := function(w)
		local i,L,red_L,change,last,last_ind;
		red_L:=List(LetterRepAssocWord(w),AbsInt);
		change := true;
		while (change) do
			change := false;
			last := 5; #not a possable generator
			last_ind := -1;
			for i in [1..Size(red_L)] do
				if IsBound(red_L[i]) then
					if red_L[i]=last then #all generators are of order two
						Unbind(red_L[i]);
						Unbind(red_L[last_ind]);
						last:=5;
						last_ind:=-1;
						change:= true;
					else 
						L := [b,c,d];
						if (not last in [a,5]) and red_L[i] in L then
							Remove(L,Position(L,last));
							Remove(L,Position(L,red_L[i]));
							red_L[last_ind] := L[1]; #bc=cb=d, bd=db=c, cd=dc=b
							last := L[1];
							Unbind(red_L[i]);
							change := true;
						else
							last:=red_L[i];
							last_ind:=i;
						fi;
					fi;
				fi;	
			od;
		od;
		#Fill the gaps
		L:= [];
		for i in red_L do
		 Add(L,i);
		od;
		return AssocWordByLetterRep(FamilyObj(w),L);
	end;
#-----------------------------------------------------------------

	#Shortens a given Letter-word over Generators of L by killing all instances of x,-x and 1,1 and 2,2, -1,-1,-2,-2
	shorten_word := function(w) 
		local change, last_pos, l, new_w;
		change := true;
		last_pos := Size(w)+1;
		w[last_pos] := 0;
		while change do
			last_pos := Size(w);
			change := false;
			for l in [1..Size(w)-1] do 
				if IsBound(w[l]) then
					if w[l] = -1*(w[last_pos]) or (w[l]= w[last_pos] and w[l] in [-2,-1,1,2]) then
						change := true;
						Unbind(w[l]);
						Unbind(w[last_pos]);
						last_pos := Size(w);
					else
						last_pos := l;
					fi;
				fi;
			od;
		od;
		new_w:=[];
		for l in w do		
			Add(new_w,l);
		od;
		return new_w{[1..Size(new_w)-1]};
	end;
#-----------------------------------------------------------------
	#Given a generator gen of L and a Letter-word w in Grig, computes the gen^w in generators of L
	compute_conjugates := function(gen,w) 
		local gen_conjugates, Gen, x, g, L;
		#Precomputed list gen_conjugates[x][y] is x^y as word in L_gen
		#where L_gen = [[b],[a,b,a],[b,a,d,a,b,a,d,a],[a,b,a,d,a,b,a,d]];
		#and y in [b,c,d,a]
		gen_conjugates := [];
		gen_conjugates[1]  := [];
		gen_conjugates[1][a] := [2];
		gen_conjugates[1][b] := [1];
		gen_conjugates[1][c] := [1];
		gen_conjugates[1][d] := [1];
		
		gen_conjugates[2]  := [];
		gen_conjugates[2][a] := [1];
		gen_conjugates[2][b] := [1,2,1];
		gen_conjugates[2][c] := [1,-4,2,1];
		gen_conjugates[2][d] := [-4,2];
		
		gen_conjugates[3]  := [];
		gen_conjugates[3][a] := [4];
		gen_conjugates[3][b] := [1,3,1];
		gen_conjugates[3][c] := [-3];
		gen_conjugates[3][d] := [1,-3,1];
		
		gen_conjugates[4]  := [];
		gen_conjugates[4][a] := [3];
		gen_conjugates[4][b] := [1,4,1];
		gen_conjugates[4][c] := [1,-4,1];
		gen_conjugates[4][d] := [-4];
		
		#gen_conjugates := [[[1],[1],[1],[2]],
		#										[[1,2,1],[1,-4,2,1],[-4,2],[1]],
		#										[[1,3,1],[-3],[1,-3,1],[4]],
		#										[[1,4,1],[1,-4,1],[-4],[3]]];
		Gen := [gen];
		for x in w do
			L:= [];
			for g in Gen do
				if g<0 then
					Append(L,List(Reversed(gen_conjugates[-1*g][x]),y->-1*y));
				else 
					Append(L,gen_conjugates[g][x]);
				fi;
			od;
			Gen := shorten_word(L);
		od;
		return Gen;
	end;
#-----------------------------------------------------------------
	#Given a Letter-word w over G and a Letter-word v over L_Gen returns v^w as word over L_Gen.
	compute_conjugates_of_word := function(v,w)
		local con, x;
		con := [];
		for x in v do
			Append(con,compute_conjugates(x,w));
		od;
		return con;
	end;
#-----------------------------------------------------------------
	#Given a word w in Generators of Grig computes a unique representative of w·L and the corresponding word in Letters
	#of generators of L. 
	#The resulting representative is an element of K_repr_words{[1..8]}
	L_Decomp := function(w)
		local l_elm,l_elm_compl,k,l,i,L,red_L,new_L,change,gen_conjugates;
		w:=alternating_a_form(w);
		l_elm := []; #Will contain tuples [v,w,...] meaning l = ...b^w·b^v
		#L_gen := [[b],[a,b,a],[b,a,d,a,b,a,d,a],[a,b,a,d,a,b,a,d]];
		change := true;
		while change do
			change := false;
			new_L := [];
			red_L:=List(LetterRepAssocWord(w),AbsInt);
			for i in Reversed([1..Size(red_L)]) do
				if red_L[i] = b then
					change := true;
					Add(l_elm,Reversed(new_L));
				elif red_L[i] = c then
					change := true;
					Add(l_elm,Reversed(new_L));
					Add(new_L,d);
				else 
					Add(new_L,red_L[i]);
				fi;
			od;
			new_L := Reversed(new_L);
			w :=	alternating_a_form(AssocWordByLetterRep(Fam,new_L));
		od;
		l_elm_compl := [];
		for l in Reversed(l_elm) do
			Append(l_elm_compl,compute_conjugates(1,l));
		od;
		#Force the form unique beginning with a.
		if Length(w)>7 then 
			w:=Subword(w,1,Length(w) mod 8);
		fi;
		if Length(w)>0 and Subword(w,1,1) = AssocWordByLetterRep(Fam,[d]) then
			w := Subword(AssocWordByLetterRep(Fam,[a,d,a,d,a,d,a,d]),1,8-Length(w));
		fi;
		return [w,shorten_word(l_elm_compl)];
	end;
#-----------------------------------------------------------------
	#Given a word w in Generators of Grig computes a unique represantative of w·K.
	#The result is an element of K_repr_words
	Compute_K_rep := function(w)
		local l,L,red_L,new_L,change,nb,b_exist;
		w:=alternating_a_form(w);
		change := true;
		while change do
			change := false;
			new_L := [];
			#In Grig all generators are selfinverse.
			red_L:=List(LetterRepAssocWord(w),AbsInt);
			nb := 0; #Stores the number of b's occuring
			for l in red_L do
				if l = b then
					nb := nb +1;
				elif l = c then
					nb := nb +1;
					Add(new_L,d);
				else 
					Add(new_L,l);
				fi;
			od;
			w :=	alternating_a_form(AssocWordByLetterRep(Fam,new_L));
			if IsOddInt(nb) then
				w := AssocWordByLetterRep(Fam,[b])*w;
				b_exist := true;
			else 
				b_exist := false;
			fi;
		od;
		#Force the word to begin with a (after the possable b) to gain a unique form.	
		if b_exist then
			w := Subword(w,2,Length(w));
		fi;
		if Length(w)>7 then 
			w:=Subword(w,1,Length(w) mod 8);
		fi;
		if Length(w)>0 and Subword(w,1,1) = AssocWordByLetterRep(Fam,[d]) then
			w := Subword(AssocWordByLetterRep(Fam,[a,d,a,d,a,d,a,d]),1,8-Length(w));
		fi;
		if b_exist then
			w := AssocWordByLetterRep(Fam,[b])*w;
		fi;
		return w;
	end;
#-----------------------------------------------------------------
	#Given a Letter-word w in L_gen computes a Letter-word res in [a,b,c,d] such that <w,1>=res in Grig.
	L_word_to_Grig := function(w)
		local Pre,Res,x;
		#Precomputed set:
		Pre := [[a,d,a],[b,a,d,a,b],[c,b,a,d,a,b,a,c,a,b,a,d,a,b,a,c,a,c],[b,a,d,a,b,a,c,a,b,a,d,a,b,a,c,a]];
		Res := [];
		for x in w do
			if x<0 then
				Append(Res,Reversed(Pre[-1*x]));
			else
				Append(Res,Pre[x]);
			fi;
		od;
		return Res;
	end;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%   Helping Functions    %%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#Computes the conjugator tuple for the pair (g,a): 
	ConTup_a := function (g)
		local g1_modL,l,Allowed_reps,Connected_conjs,con_at_1,con_word,con,Centr_a,Con_tuple;
		if IsOne(Activity(g)) then
			return [];
		fi;
		if not IsOne(State(g,1)*State(g,2)) then
			return [];
		fi;
		#L_gen := [[b],[a,b,a],[b,a,d,a,b,a,d,a],[a,b,a,d,a,b,a,d]];
		g1_modL:=L_Decomp(PreImagesRepresentativeNC(f,State(g,1))); 
		l:=g1_modL[2];
		g1_modL:=LetterRepAssocWord(g1_modL[1]);
		#See Lemma lem:conjugators_of_a for Details
		Allowed_reps:= [[],[a,d],[a,d,a,d,a,d],[a,d,a,d]];
		if not g1_modL in Allowed_reps then
			return [];
		fi;
		#See Lemma lem:conjugators_of_a for the appix conjugator
		Connected_conjs := [[],[c],[a,c],[c,a,c]];
		con := Connected_conjs[Position(Allowed_reps,g1_modL)];
		#resulting conjugator is of the form <l^((g_1modL^-1)),1>·con
		con_at_1 := compute_conjugates_of_word(l,Reversed(g1_modL));
		con_word := L_word_to_Grig(con_at_1);
		Append(con_word,con);
		Info(InfoFR,4,"Conjugator in gen_L: <",con_at_1,",1>",con,"\nConjugator in gen_Grig: ",con_word,"\n");
		#Determine Cosets of K in which the conjugator lies.
		#See Roskov CP Lemma3 for centralizer of a
		Centr_a := List([[],[a],[a,d,a,d],[a,d,a,d,a]],x -> AssocWordByLetterRep(Fam,Concatenation(con_word,x)));
		Con_tuple:= [];
		for con in Centr_a do
			Con_tuple[Position(K_repr,LetterRepAssocWord(Compute_K_rep(con)))] := ImageElm(f,con);
		od;	
		Check("ConTup_a",g,ae,Con_tuple);
		return Con_tuple;
	end;	
	#Finds all Elements <l1,l2> with <l1,l2> in Grig, for l1 in L1, l2 in L2 and return result as Conjugator tuple.
	Merge_Ls := function(L1,L2,with_action)
		local aw_w,aw_t,dw_w,res_Con,i,x;	
		aw_w := One(aw);
		aw_t := ();
		dw_w := One(dw);
		if with_action then
			aw_w := aw;
			aw_t := (x_1,x_2);
		fi;
		Info(InfoFR,4,"Computing ",g,",",h,"  Sub Conjugators: ",L1,"\n");
		Info(InfoFR,4,"Computing ",g,",",h,"  Sub Conjugators: ",L2,"\n");
		#See Lemma 6.16 for <g1,g2<in Grig,  <=> g1=v(a,d)l g2=v(d,a)l
		#So <K_repr[i],K_repr[j]> in Grig  <=> j in [17-x mod 16 +1, 25-x mod 16 +1]
		res_Con := [];
		for i in [1..16] do
			if IsBound(L1[i]) then
				#Find second entry:
				for x in [((17-i) mod 16) +1,((25-i) mod 16) +1] do
					if IsBound(L2[x]) then
						if x>8 then
							dw_w := dw;
						else
							dw_w := One(dw);
						fi;
						Info(InfoFR,4,"Computing ",g,",",h,"  Conjugator found:",i,",",x,"\n");
						if L1[i]=ImageElm(f,K_repr_words[i]) and L2[x]=ImageElm(f,K_repr_words[x]) then
							res_Con[Position(K_repr_words,Compute_K_rep(dw_w*D[i]*aw_w))] := ImageElm(f,dw_w*D[i]*aw_w);
						else #Could always compute the words as generators, but seems uneccassary
							res_Con[Position(K_repr_words,Compute_K_rep(dw_w*D[i]*aw_w))] := MEALY_FROM_STATES@([L1[i],L2[x]],aw_t);
						fi;
					fi;
				od;
			fi;
		od;
		return res_Con;
	end;
	#Joins two Lists in the first with overwriting eventually existing values.
	Join_to_first := function(L,K)
		local i;
		for i in [1..Length(K)] do
			if IsBound(K[i]) then
				L[i] := K[i];
			fi;
		od;
	end;		
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%     Main Computor      %%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	conjugators_grig_rek := function(g,h)
		local Centr_bc, Centr_d, L1, L1_temp, L2, res_Con, g1, h1, x ,y;
		Info(InfoFR,3,"Compute Conjugator for pair: ",g,",",h,".\n");
		if Activity(g) <> Activity(h) then
			return [];
		fi;
		#-#-#-#-#-#-#-#-#-#-#-#-   g = identity   -#-#-#-#-#-#-#-#-#-#
		if IsOne(g) then
			if IsOne(h) then
				return List(K_repr,x -> ImageElm(f,AssocWordByLetterRep(Fam,x)));
			fi;
				return [];
		fi;
		#-#-#-#-#-#-#-#-#-#-#-   g in a,b,c,d    -#-#-#-#-#-#-#-#-#-#
		if g in [be,ce] then
			if h=g then
				Centr_bc := [[],,,,,,,[a,d,a,d,a,d,a],[b],,,,,,,[b,a,d,a,d,a,d,a]];
				return List(Centr_bc,x -> ImageElm(f,AssocWordByLetterRep(Fam,x)));
			fi;
		fi;
		if g = de then
			if h=g then
				Centr_d := [[],,,[a,d,a],[a,d,a,d],,,[a,d,a,d,a,d,a],[b],,,[b,a,d,a],[b,a,d,a,d],,,[b,a,d,a,d,a,d,a]];
				return List(Centr_d,x -> ImageElm(f,AssocWordByLetterRep(Fam,x)));
			fi;
		fi;
		if g=ae then
			return List(ConTup_a(h),x->x^-1);
		fi;
		if g in [be,ce,de] then
			if h in [be,ce,de,One(h),ae] then
				return []; #As g=h already considered in an earlier case
			fi;
			#---------------------    |h|>1     -----------------------
			#Test for Conjugator with trivial Activity
			res_Con := [];
			L1 := conjugators_grig_rek(State(g,x_1),State(h,x_1));
			if Size(L1) > 0 then
				L2 := conjugators_grig_rek(State(g,x_2),State(h,x_2));
				res_Con := Merge_Ls(L1,L2,false);
			fi;		

			#Test for Conjugator with non-trivial Activity
			L1 := conjugators_grig_rek(State(g,x_1),State(h,x_2));
			if Size(L1) = 0 then
				return res_Con;
			fi;
			L2 := conjugators_grig_rek(State(g,x_2),State(h,x_1));
			Join_to_first(res_Con,Merge_Ls(L1,L2,true));
			Check("h>1 nontrivial",g,h,res_Con);
			return res_Con;
		fi;

		#-#-#-#-#-#-#-#-#-   |g| > 1, act(g) = 1    -#-#-#-#-#-#-#
		res_Con := [];
		if IsOne(Activity(g)) then
		#Test for Conjugator with trivial Activity
			L1 := conjugators_grig_rek(State(g,x_1),State(h,x_1));
			if Size(L1) > 0 then
				L2 := conjugators_grig_rek(State(g,x_2),State(h,x_2));
				res_Con := Merge_Ls(L1,L2,false);
			fi;	
			#Test for Conjugator with non-trivial Activity
			L1 := conjugators_grig_rek(State(g,x_1),State(h,x_2));
			L2 := conjugators_grig_rek(State(g,x_2),State(h,x_1));
			Join_to_first(res_Con,Merge_Ls(L1,L2,true));
			Check("g>1, act = 1, non-trivial",g,h,res_Con);
			return res_Con;
		else
		#-#-#-#-#-#-#-#-#-   |g| > 1, act(g) = (1,2)    -#-#-#-#-#
			#Test for Conjugator with trivial Activity
			g1 := Compute_K_rep(PreImagesRepresentativeNC(f,State(g,x_1)^-1));
			h1 := Compute_K_rep(PreImagesRepresentativeNC(f,State(h,x_1)));
			L1 := conjugators_grig_rek(State(g,x_1)*State(g,x_2),State(h,x_1)*State(h,x_2));
			res_Con := [];
			if Size(L1) > 0 then
				for x in L1 do
					#Force that only <x,g1xh1> is checked. #Seems to be a bit too complicated, may be simplified.
					L1_temp := [];
					L1_temp[Position(L1,x)]:=x;
					L2 := [];
					L2[Position(K_repr_words,Compute_K_rep(Compute_K_rep(g1)*K_repr_words[Position(L1,x)]*Compute_K_rep(h1)))] := State(g,x_1)^-1*x*State(h,x_1);
					L2 :=Merge_Ls(L1_temp,L2,false);
					for y in L2 do
						res_Con[Position(L2,y)] := y;
					od;
				od;
			fi;
			#Test for Conjugator with non-trivial Activity
			h1 := Compute_K_rep(PreImagesRepresentativeNC(f,State(h,x_2)));
			L1 := conjugators_grig_rek(State(g,x_1)*State(g,x_2),State(h,x_2)*State(h,x_1));
			if Size(L1) = 0 then
				return res_Con;
			fi;
			for x in L1 do
				#Force that only <x,g1xh1> is checked.
				L1_temp := [];
				L1_temp[Position(L1,x)]:=x;
				L2 := [];
				L2[Position(K_repr_words,Compute_K_rep(Compute_K_rep(g1)*K_repr_words[Position(L1,x)]*Compute_K_rep(h1)))] := State(g,x_1)^-1*x*State(h,x_2);
				L2 :=Merge_Ls(L1_temp,L2,true);
				for y in L2 do
					res_Con[Position(L2,y)] := y;
				od;
			od;
			Check("g>1, act = (1,2), non-trivial",g,h,res_Con);
			return res_Con;
		fi;
	end;
	
	Res:= conjugators_grig_rek(g,h);
	Info(InfoFR,3,"Result of recursive computation: ",Res,"\n");
	if Size(Res) = 0 then
		return fail;
	fi;
	return Representative(Res);
end);
SetFRConjugacyAlgorithm(GrigorchukGroup,GRIG_CON@);
#################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%      IsConjugate        %%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%	 RepresentativeActionOp %%%%%%%%%%%%%%%%%%#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#################################################################	
InstallMethod(IsConjugate,
	" For FR groups with optimized conjugacy algorithm ",
	[ IsFRGroup and HasFRConjugacyAlgorithm,IsFRElement,IsFRElement], 
  function(G,a,b)
  	return FRConjugacyAlgorithm(G)(G,a,b) <> fail;
 	end);
InstallOtherMethod(RepresentativeActionOp,
 " For FR groups with optimized conjugacy algorithm ",
	[ IsFRGroup and HasFRConjugacyAlgorithm,IsFRElement,IsFRElement,IsFunction], 
  function(G,a,b,f)
  	if f <> OnPoints then TryNextMethod(); fi;
  	return FRConjugacyAlgorithm(G)(G,a,b);
 	end);
 	
  















	
	
	
	
	











