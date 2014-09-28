InstallMethod(IteratedOrbit,
	"Computes the Orbit of a word under <first>",
	[ IsFRElement, IsList],
  function(a, L )
		return Orbit(Group(a),L);
  end
);
InstallMethod(IteratedOrbit,
	"Computes the Orbit of the second argument under <first>",
	[ IsFRElement, IsInt],
  function(a, x )
		return Orbit(Group(a),x);
  end
);
#---------------------------------------------------------------
#------              Dep-Cartesian          --------------------
#--  Calculates a cartesian product of ordered lists with    ---
#-- respect to dependencies, which entries belong together   ---
#-- Example: L=[[A,B,X],[C,D],[c,d]], dep=[[1],[2,3]] results --
#--   in [[A,C,c],[A,D,d],[B,C,c],[B,D,d],[X,C,c],[X,D,d]	   ---
#--     The Lists, which are joined by the dependencies      ---
#--            have to be of the same length	               ---	
#---------------------------------------------------------------
MAKE_READ_WRITE_GLOBAL("DEP_CARTESIAN@");
UNBIND_GLOBAL("DEP_CARTESIAN@");
BindGlobal("DEP_CARTESIAN@", function(L,dep)
	local res_list, temp_cart, container, al, d, i ,j,a;
	res_list := [];
	temp_cart := [];
	for d in dep do
		container := [];
		for j in [1..Size(L[d[1]])] do
			al := [];
			for i in [1..Size(d)] do
				Add(al,L[d[i]][j]);
			od;
			Add(container,al);
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
MAKE_READ_WRITE_GLOBAL("LEVEL_PERM_CONJ@");
UNBIND_GLOBAL("LEVEL_PERM_CONJ@");
BindGlobal("LEVEL_PERM_CONJ@", function(x,y)
	local pi_x, pi_y, c;
	if Alphabet(x) <> Alphabet(y) then
		Error("Not a valid argument, they must be defined on the same alphabet.");
		return [];
	fi;
	pi_x := PermList(DecompositionOfFRElement(x)[2]);
 	pi_y := PermList(DecompositionOfFRElement(y)[2]);
 	c := RepresentativeAction(SymmetricGroup(Alphabet(x)),pi_x,pi_y);
 	if c= fail then
 		return [];
 	fi;
 	return c*List(Centralizer(SymmetricGroup(Alphabet(x)),pi_y));
end);
######################################################################
#````````````````````````````````````````````````````````````````````#
#`````````````````````    OrbitSignalizer   `````````````````````````#
#``````````````````                            ``````````````````````#
#``````````````````    Garantied to stop on    ``````````````````````#
#`````````````````` BoundedFRElements as input ``````````````````````#
#````````````````````````````````````````````````````````````````````#
######################################################################
InstallMethod(OrbitSignalizer,
	"Returns the finite Orbit Signalizer",
	[IsFRElement],
function(a)
	local OS_list,i,OS_unvisited,OS_new,elm,x,new,suc;
	suc := function(state,x)
		return  State(state^Size(IteratedOrbit(state,[x])),[x]);
	end;
	OS_list := [];
	OS_unvisited := [a];
	while Length(OS_unvisited) > 0 do
		OS_new := [];
		for elm in OS_unvisited do
			for x in Alphabet(a) do
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
######################################################################
#````````````````````````````````````````````````````````````````````#
#```````````````````````                  ```````````````````````````#
#```````````````````````  ConjugatorGraph ```````````````````````````#
#```````````````````````     DrawGraph    ```````````````````````````#
#```````````````````````                  ```````````````````````````#
#````````````````````````````````````````````````````````````````````#
######################################################################
MAKE_READ_WRITE_GLOBAL("CONJUGATOR_GRAPH@");
UNBIND_GLOBAL("CONJUGATOR_GRAPH@");
BindGlobal("CONJUGATOR_GRAPH@", function(a,b)
	local Alph, Vertices, Edges, c, d, p, v_id, e_id, v, orbits, orb_repr, i, new_con_pair, new_v, w, change, found, e, all_found;
	
	Alph := Alphabet(a);
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
		#Print("Looking at ",v.id,"\n");
		c := v.conj_pair[1];
		d := v.conj_pair[2];
		orbits := Orbits(Group(c),Alph);
		#Print("orbits: ",orbits,"\n");
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

MAKE_READ_WRITE_GLOBAL("DRAW_GRAPH@");
UNBIND_GLOBAL("DRAW_GRAPH@");
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
######################################################################
#````````````````````````````````````````````````````````````````````#
#```````````````````````                  ```````````````````````````#
#```````````````````````    F.S. Worker   ```````````````````````````#
#```````````````````````                  ```````````````````````````#
#````````````````````````````````````````````````````````````````````#
######################################################################
MAKE_READ_WRITE_GLOBAL("CONJUGATORS_FINITE_STATE_WRAPPER@");
UNBIND_GLOBAL("CONJUGATORS_FINITE_STATE_WRAPPER@");
BindGlobal("CONJUGATORS_FINITE_STATE_WRAPPER@",function(start,CG)
	local v,AS,to_visit, Alph, new_v, i, found, e, Tran, Act, c,d, orbit;
			#--------- Choose one subgraph, as automaton  ---------
			AS := [start.id]; #Contains IDs of vertices, which build the subgraph
			to_visit := [start.id]; 
			Alph := Alphabet(start.conj_pair[1]);
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
					orbit := IteratedOrbit(c,e.read);
					for i in [2..Length(orbit)] do
					#The missing edges...
						Tran[Position(AS,e.from)][Position(Alph,orbit[i])] := [State(c^(i-1),e.read)^(-1),Position(AS,e.to),State(d^(i-1),e.read^(CG[1][e.from].action))];
					od;
				fi;
			od;
			return FRElement(Tran,Act,[1]);
end);
######################################################################
#````````````````````````````````````````````````````````````````````#
#```````````````````````                  ```````````````````````````#
#```````````````````````  Finitary Worker ```````````````````````````#
#```````````````````````                  ```````````````````````````#
#````````````````````````````````````````````````````````````````````#
######################################################################
MAKE_READ_WRITE_GLOBAL("CONJUGATORS_FINITARY_WRAPPER@");
UNBIND_GLOBAL("CONJUGATORS_FINITARY_WRAPPER@");
BindGlobal("CONJUGATORS_FINITARY_WRAPPER@",function(v,Graph,Seen,Known_vertex_conjugator)
	local CONJUGATORS_FINITARY_REK;
	CONJUGATORS_FINITARY_REK := function(v,Graph,Seen,Known_vertex_conjugator)
		local Vertices,Edges,sons,starts,conj_cand,conjugators_found,son,x,NewSeen,a,b,e,son_conj,son_conjs,tempo_conj, htemp,tempoconj,err,Indices,diction,pos,real_conj,real_conjugators,i,j,ip,son_orbit_size,Alph,action,
		w,Circle,Conjs,con;
	
		#Print("@Vertex ",v,"\n");
		if IsBound(Known_vertex_conjugator[v]) then #Don't do the same work twice.
			#Print("Save work\n Go up...\n");
			return Known_vertex_conjugator[v];
		fi;
		Vertices := Graph[1];
		Edges := Graph[2];
		a := Vertices[v].conj_pair[1];
		b := Vertices[v].conj_pair[2];
		Alph := Alphabet(a);
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
			#Print("Go Down....\n");
			son_conjs := CONJUGATORS_FINITARY_REK(son[1],Graph,NewSeen,Known_vertex_conjugator);
			#Print("Back....\n");
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
		#Print("conjugators_found: ",conjugators_found,"\n");
		if IsDenseList(conjugators_found) and Size(conjugators_found) = Size(Alph) then
			#puzzle them together!	
			Conjs := DEP_CARTESIAN@(conj_cand,Orbits(Group(a),Alph));
			for con in Conjs do
				Add(real_conjugators,FRElement([con],[action],[1]));
			od;
		fi;
		Known_vertex_conjugator[v] := real_conjugators;
		#Print("Known_vertex_conjugator: ",Known_vertex_conjugator,"\n");
		#Print("Go up...\n");
		return real_conjugators;
	end;
	return CONJUGATORS_FINITARY_REK(v,Graph,Seen,Known_vertex_conjugator);
end);
######################################################################
#````````````````````````````````````````````````````````````````````#
#```````````````````````                  ```````````````````````````#
#```````````````````````  BoundedWorker   ```````````````````````````#
#```````````````````````                  ```````````````````````````#
#````````````````````````````````````````````````````````````````````#
######################################################################
MAKE_READ_WRITE_GLOBAL("CONJUGATORS_BOUNDED_WRAPPER@");
UNBIND_GLOBAL("CONJUGATORS_BOUNDED_WRAPPER@");
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
		Alph := Alphabet(a);
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
				#Print("Running through Circle. Now at vertex ",Seen[i],".\n");
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
			#Print("Finished this vertex (",v,") with the following Conjugators:\n",Conjs,"\n");
			#Print("Aktually known are:",Known_vertex_conjugator,"\n");
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
		#Print("conjugators_found: ",conjugators_found,"\n");
		if IsDenseList(conjugators_found) and Size(conjugators_found) = Size(Alph) then
			#puzzle them together!	
			Conjs := DEP_CARTESIAN@(conj_cand,Orbits(Group(a),Alph));
			for con in Conjs do
				Add(real_conjugators,FRElement([con],[action],[1]));
			od;
		fi;
		Known_vertex_conjugator[v] := real_conjugators;
		#Print("Known_vertex_conjugator: ",Known_vertex_conjugator,"\n");
		#Print("Go up...\n");
		#Print("Finished this vertex (",v,") with the following Conjugators:\n",real_conjugators,"\n");
		return real_conjugators;
	end;
	return CONJUGATORS_BOUNDED_REK(v,Graph,Seen,readwrite_path,Known_vertex_conjugator);
end);

###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%      IsConjugate        %%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%	 RepresentativeActionOp %%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###################################################################
InstallMethod(IsConjugate,
	"For Aut, RAut, FAut, Poly-1, Poly0",
	#The attribute FullSCVertex charakterizes all FullSCGroups
	[ IsFRGroup and HasFullSCVertex,IsFRElement,IsFRElement], 
  function(G,a,b)
  	local v, Graph, sons, starts;
  	if Alphabet(G) <> Alphabet(a) or Alphabet(G) <> Alphabet(b) then
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
  	if Alphabet(G) <> Alphabet(a) or Alphabet(G) <> Alphabet(b) then
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











