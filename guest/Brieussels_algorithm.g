################################################################################################################
#
# This GAP programme implements the algorithm by Brieussel to obtain the pull-back of two given
# words with letters in Grigorchuk's group, see [Bri08].
#
# Brieussel's algorithm is a combination of a modified version of an algorithm by Leonov [Leo68]
# and an algorithm by Grigorchuk [Gri84]. Therefore this programme is split into the functions
# "leomod_algo" and "grigor_algo", implementing the two mentioned algorithms.
# These two functions work with a list [w1,w2] of two reduced words, whereas the reduction of words
# is indirectly implemented with the function "add_reduced". The assumption of a certain word format 
# Brieussel's algorithm works with is implemented in the function "hypothesis".
#
# The complete algorithm was coded in two versions: 
# "bri_algo" rewrites the input words according to the hypothesis, so works for any pair of words. 
# "bri_algo_hypofulfilled" only works for words fulfilling the assumption, but is therefore 
# significantly faster. 
# Both functions work with two words as input in the format (x0*x1*..., y0*y1*...) with xi,yi in {a,b,c,d}.
# In both cases the input words are compared with the output taking the mistake into account with
# the function "check_equality".
#
# "test_with_hypo" tests the algorithm with random words fulfilling the assumption and creates a list
# with the length of the output word. It takes as parameters the path of the data file, two ranges
# of integers constructing a lattice of integer points in IR² and the number of tests for each node.
#
# The exact input format is explained for each function in detail in the comment section before each function.
# This code uses the package "FR" by Laurent Bartholdi [Bar].
#
################################################################################################################
#
# Programmed 2012 by Michael Duppré as part of the bachelor thesis "Lower bounds on the growth of 
# Grigorchuk’s group, according to Leonov and Brieussel", supervised by Laurent Bartholdi at the 
# Georg-August-Universität Göttingen.
#
################################################################################################################
# References:
# [Bar]   http://www.uni-math.gwdg.de/laurent/FR/
# [Bri08] Jérémie Brieussel, Croissance et moyennabilité de certains groupes d’automorphimes d’un arbre enracine,
#         http://www.institut.math.jussieu.fr/theses/2008/brieussel/these-brieussel.pdf, 2008
# [Gri84] R. I. Grigorchuk, Degrees of growth of finitely generated groups and the theory of invariant means, 
#         Izvestiya Akademii Nauk SSSR. Seriya Matematicheskaya. vol.48 (1984), no. 5, pp. 939–985.5, 1984
# [Leo68] Yu. G. Leonov, On a lower bound for the growth of a 3-generator 2-group, Mat.Sb., 192:77-92, 2001
################################################################################################################

G := GrigorchukGroup;

# define epimorphism from free group, call it F and assign generator variables
epi := EpimorphismFromFreeGroup(G : names := ["a","b","c","d"]);
F := Source(epi);
AssignGeneratorVariables(F);

# prompt for user
Print("\nPlease apply the function \"algo()\" on two words in the generators {a,b,c,d}\n(e.g. algo(a*b*a*c*a*d,d*a*c*b*a*d))\n\n");
Print("To perform tests on the length of input and output words use \ntest_with_hypo := function(path,range1,range2,n)\nwhich creates a datafile at path with the length of each input word together with the average, minimal and maximal length for n tries. range1 indicates which length of word1 is given together with the length given in range2, where every combination is computed. The input words are random and calculated according to the hypothesis of having dc-type.");





###################################################################################################
# delete_last
#	deletes last entry of a list x times
###################################################################################################
#	input:		list and integer indicating how many last entries should be deleted
#	input form:	([x1,x2,x3,...,xm],n) with n an integer
#	output:		list with deleted entries
#	output form:	[x1,x2,x3,...,x(m-n)]
####################################################################################################

delete_last := function(l,x)
	local i;
	
	i := x;
	while i > 0 do
		Remove(l);
		i := i-1;
	od;
return l;
end;





####################################################################################################
# add_reduced function
#	adds an entry to a list, with respect to the relations a² = b² = c² = d² = bcd = 1
####################################################################################################
#	input:		list representing a word, letter to be added as integer
#	input form:	([x0, x1, x2, x3, ...],v) with xi,v in {1,2,3,4}
#	output:		list with letter added and possibly reduced again
#	output form:	[1, y1, 1, y2, 1, y3, ...] or [y0, 1, y1, 1, y2, 1, ...] with yi in {2,3,4}
####################################################################################################

add_reduced := function(l,v)
	
	local i;

	for i in [1..Length(v)] do

		# check if l is the empty list
		if Length(l) = 0 then
			Add(l,v[i]);
		else
			# check for double entries, according to a²=b²=c²=d²=1
			if l[Length(l)] = v[i] then
				Remove(l);
			# replace bc=cb=d
			elif ((l[Length(l)] = 2 and v[i] = 3) or (l[Length(l)] = 3 and v[i] = 2)) then
				Remove(l);
				Add(l, 4);
			# replace cd=dc=b
			elif ((l[Length(l)] = 3 and v[i] = 4) or (l[Length(l)] = 4 and v[i] = 3)) then
				Remove(l);
				Add(l, 2);
			# replace db=bd=c
			elif ((l[Length(l)] = 4 and v[i] = 2) or (l[Length(l)] = 2 and v[i] = 4)) then
				Remove(l);
				Add(l, 3);
			else
				Add(l,v[i]);
			fi;
		fi;
	od;

return l;
end;





####################################################################################################
# list_to_double_list function
#	converts a word in list representation to a representation in a double list with the 
#	relations a = (1,1)e, b = (a,c), c = (a,d), d = (1,b) - see psi in thesis
####################################################################################################
#	input:		list representing a word
#	input form:	[x0, x1, x2, ...] with xi in {1,2,3,4}
#	output:		reduced list representing the initial word
#	output form:	[ [y0,y1,...], [z0,z1,...], v ] 
#			with yi, zi in {1,2,3,4}, v boolean, indicating if first level subtrees are permuted or not
####################################################################################################

list_to_double_list := function(l)
	local l1, l2, i, swap;

	# set initial values
	l1 := [];
	l2 := [];
	swap := false;

	# go through every entry of given list
	for i in [1..Length(l)] do
		if l[i] = 1 then
			swap := not swap;
		elif l[i] = 2 then
			if swap = false then
				add_reduced(l1, [1]);
				add_reduced(l2, [3]);
			else
				add_reduced(l1, [3]);
				add_reduced(l2, [1]);
			fi;
		elif l[i] = 3 then
			if swap = false then
				add_reduced(l1, [1]);
				add_reduced(l2, [4]);
			else
				add_reduced(l1, [4]);
				add_reduced(l2, [1]);
			fi;
		elif l[i] = 4 then
			if swap = false then
				add_reduced(l2, [2]);
			else
				add_reduced(l1, [2]);
			fi;
		fi;
	od;

return [l1,l2,swap];
end;





####################################################################################################
# dab function
#	converts dab... or adab... at beginning of word, to fulfill hypothesis
####################################################################################################
#	input:		list representing a reduced word
#	input form:	[x0, x1, x2, x3, ...]
#	output:		list with beginning dab... or adab... replaced
#	output form:	[y0, y1, y2, y3, ...]
####################################################################################################

dab := function(l)
	local start, rev;
	
	# only apply if list is longer than 2
	if Length(l) > 2 then
		start := 0;

		# set start value
		if l[1] = 4 then
			start := 1;
		elif l[2] = 4 then
			start := 2;
		fi;	

		# if word starts with dab... or adab...
		if start <> 0 and Length(l) >= start + 2 and l[start+1] = 1 and l[start+2] = 2 then
			# dab = bacacacacacacad = ba(ca)^6d
			rev := Reversed(l);
			delete_last(rev, 2+start);
			add_reduced(rev, [4,1,3,1,3,1,3,1,3,1,3,1,3,1,2]);
			if start = 2 then
				add_reduced(rev, [1]);
			fi;
			l := Reversed(rev);
		fi;
	fi;
		
return l;
end;





####################################################################################################
# dada function
#	replaces all sequences of the form (da)^n
####################################################################################################
#	input:		list representing a reduced word
#	input form:	[x0, x1, x2, x3, ...]
#	output:		reduced list with all sequences (da)^n replaced using (ad)^4
#	output form:	[y0, y1, y2, y3, ...]
####################################################################################################

dada := function(l)
	local i, d_pos, last_length;

	d_pos := [];
	last_length := -1;

	# concept: only replace dada -> adad
	# enough, when reduce_list after that applied
	# advantage: length of list doesn't change

	# run replacement as long as Length doesn't change anymore OR if length becomes smaller than 4
	while Length(l) <> last_length and Length(l) >= 4 do

		last_length := Length(l);
		for i in [1..Length(l)] do
			if l[i] = 4 then
				Add(d_pos, i);
			fi;
		od;


		# replace d(ada) -> (ada)d
		while Length(d_pos) > 1 do
			# prevent case [...,1,4,1,4] to fail
			if d_pos[1] + 2 = d_pos[2] and d_pos[2] <> Length(l) then
				Add(l, 4, d_pos[2]+2);
				Remove(l, d_pos[1]);
				Remove(d_pos, 1);
				Remove(d_pos, 1);
			# case [...,1,4,1,4]
			elif d_pos[1] + 2 = d_pos[2] and d_pos[2] = Length(l) then
				Add(l, 1);
				Remove(l,d_pos[1]-1);
				Remove(d_pos, 1);
				Remove(d_pos, 1);
			else
				Remove(d_pos, 1);
			fi;
		od;
		l := add_reduced([], l);
	od;

return l;
end;





####################################################################################################
# d_to_c_type function
#	converts d-type to c-type representation
####################################################################################################
#	input:		list representing a reduced word
#	input form:	[x0, x1, x2, x3, ...]
#	output:		list with all d-type subwords replaced by c-type subwords
#	output form:	[y0, y1, y2, y3, ...]
####################################################################################################

d_to_c_type := function(l)
	local i, d_pos, in_subword, start, ending, len, k, lambda, z0, zk, last_pos;

	d_pos := [];
	in_subword := false;
	start := [];
	ending := [];

	# localise subwords of the form
	# s = z0*ada*z1*ada*...*ada*zk

	# store all positions of d's in a list
	for i in [1..Length(l)] do
		if l[i] = 4 then
			Add(d_pos,i);
		fi;
	od;

	# set start and ending of subwords and take care at beginning of word
	while IsEmpty(d_pos) = false do
		len := Length(d_pos);

		# check if more than one entry is in d_pos
		if len > 1 then
			last_pos := d_pos[len];
			# end of subword
			# prevent case e.g. [4,1,3,1,4,1,2,...]
			if in_subword = false and d_pos[len] - 4 = d_pos[len-1] and d_pos[len] - 4 > 2 then
				# set ending of subword 
				Add(ending, d_pos[len]);
				in_subword := true;
				Remove(d_pos);
			# end of subword
			elif in_subword = false and d_pos[len] - 4 = d_pos[len-1] and d_pos[len] - 4 <= 2 then
				# set ending of subword 
				Add(ending, d_pos[len]);
				Add(start, d_pos[len]);
				in_subword := false;
				Remove(d_pos);
			# single d
			elif in_subword = false and d_pos[len] - 4 <> d_pos[len-1] then
				Add(ending, d_pos[len]);
				Add(start, d_pos[len]);
				Remove(d_pos);
			# no start yet of subword
			elif in_subword = true and d_pos[len] - 4 = d_pos[len-1] then
				# still in the subword, don't have to set start
				Remove(d_pos);
			# found start of subword
			elif in_subword = true and d_pos[len] - 4 <> d_pos[len-1] then
				# set start of subword
				Add(start, d_pos[len]);
				Remove(d_pos);
				in_subword := false;
			fi;
		# handle last (first) entry of d_pos
		else
			# prevent case [ 4, 1, 3, 1, 4, 1, 2, 1, 4, 1, 3, 1, 2,...]
			if in_subword = true and d_pos[len] > 2 then
				Add(start, d_pos[len]);
				Remove(d_pos);
			elif in_subword = true and d_pos[len] <= 2 then
				Add(start, last_pos);
				Remove(d_pos);
			elif in_subword = false and d_pos[len] >= 3 then
				Add(start, d_pos[len]);
				Add(ending, d_pos[len]);
				Remove(d_pos);
			else
				Remove(d_pos);
			fi;
		fi;
	od;

	# actually replace the subwords
	while IsEmpty(start) = false do
		# calculate k (see lemma in thesis)
		k := (ending[Length(ending)]-start[Length(start)])/4 + 1;
		# k <= 15, since (adac)^16 = 1
		k := k mod 16;

		# calculate lambda (see lemma in thesis)
		lambda := [];
		# go from z0 to zk in subword z0*ada*z1*ada*...*ada*zk
		z0 := start[Length(start)]-2;
		zk := ending[Length(ending)]+2;
		# only replace, if there is a zk so word doesn't end before
		# only if z0 > 0 to prevent e.g. [1,4,1,...]
		if Length(l) >= zk and z0 > 0 then

			for i in [z0, z0+4..zk] do
				# b = (a,c)
				if l[i] = 2 then
					add_reduced(lambda, [3]);
				# c = (a,d)
				elif l[i] = 3 then
					add_reduced(lambda, [4]);
				# in case there is a d in between (should exclude this before!) do nothing
				else
				fi;
			od;

			# replace the subwords
			# k odd and (a(ba)^k,1) replace by c(d^ac)^k
			if k mod 2 = 1 and lambda = [] then
				l[z0] := 3;
			# k odd and (a(ba)^k,b) replace by b(d^ac)^k
			elif k mod 2 = 1 and lambda = [2] then
				l[z0] := 2;
			# k even and (a(ba)^k,d) replace by c(d^ac)^k
			elif k mod 2 = 0 and lambda = [4] then
				l[z0] := 3;
			# k even and (a(ba)^k,c) replace by b(d^ac)^k
			elif k mod 2 = 0 and lambda = [3] then
				l[z0] := 2;
			fi;
		
			# replace rest
			for i in [z0+4,z0+8..zk] do
				l[i] := 3;
			od;
		fi;

		# remove the used entries from lists
		Remove(ending);
		Remove(start);
	od;

return l;
end;





####################################################################################################
# leomod_algo - modified Leonov algorithm
#	also considering modified sequences of cases
####################################################################################################
# input:	list of two words
# input form: 	[w1, w2] with wi reduced words in list representation
# output:	list of two words with pulled back letters after m steps of algorithm and pull-back W_Leo
# output form: 	[w1_m, w2_m, W_Leo]
####################################################################################################

leomod_algo := function(w)

local i, count_b, last_pos, output, state, input1, input2, len1, len2, ending, modified, first_0d, parity, inbetween_0d, 3_or_4_possible, l1, l2, l_1, reset_modified, swap, modified_output;


reset_modified := function()
	first_0d := true;
	inbetween_0d := false;
	3_or_4_possible := false;
	parity := 0;
end;

output := [];



# only apply algorithm, if l(w1) >= 11 and l(w2) >= 11
if Length(w[1]) < 11 or Length(w[2]) < 11 then
	return [w[1], w[2], output];
else
# for finite state machine first reverse w[1] and w[2]
w := [Reversed(w[1]),Reversed(w[2])];
modified := [];
reset_modified();
state := 1;
swap := false;

# repeat as long as both lists are longer than 23 (case 45d with q = 3, l = 11)
while Length(w[1]) >= 23 and Length(w[2]) >= 23 do


	len1 := Length(w[1]);
	len2 := Length(w[2]);


	# type I
	# w1 = x0*a*x1*a*x2*a*x3*a*x4*...
	# w2 = a*y1*a*y2*a*y3*a*y4*a*...
	#
	# state = 1
	if state = 1 and w[1][len1] <> 1 and w[2][len2] = 1 then
		state := 2;


	# type II
	# w1 = a*x1*a*x2*a*x3*...
	# w2 = y0*a*y1*a*y2*a*...
	#
	# state = 1
	elif state = 1 and w[1][len1] = 1 and w[2][len2] <> 1 then
		# swap words
		w := w{[2,1]};

		# add first a of a*W*a to output
		# (last a will be added before next step, see below)
		swap := true;
		add_reduced(output, [1]);
		state := 2;


	# type III
	# w1 = a*x1*a*x2*a*x3*...
	# w2 = a*y1*a*y2*a*y3*a*...
	#
	# state = 1
	elif state = 1 and  w[1][len1] = 1 and w[2][len2] = 1 then
		# alter (w1,w2) for next step
		Add(w[1], 2);
		# add to output a*d*a (initialising)
		add_reduced(output, [1,4,1]);
		state := 2;


	# type IV
	# w1 = x0*a*x1*a*x2*...
	# w2 = y0*a*y1*a*y2*a*...
	#
	# state = 1
	#
	# x0 = b
	# go to type I (RATHER II, I guess!)
	elif state = 1 and w[1][len1] = 2 and w[2][len2] <> 1 then
		# alter (w1,w2) for next step
		Remove(w[1]);
		# add to output a*d*a (initialising)
		add_reduced(output, [1,4,1]);
		state := 2;
	# x0 = c
	# go to type III
	elif w[1][len1] = 3 and w[2][len2] <> 1 then
		# alter (w1,w2) for next step
		Remove(w[1]);
		Add(w[2], 1);
		# add to output a*b*a (initialising)
		add_reduced(output, [1,2,1]);
		state := 1;
	# x0 = d
	# go to type III
	elif w[1][len1] = 4 and w[2][len2] <> 1 then
		# alter (w1,w2) for next step
		Remove(w[1]);
		Add(w[2], 1);
		# add to output a*c*a (initialising)
		add_reduced(output, [1,3,1]);
		state := 1;




	# case 0d
	# w1 = d*a*x1*a*x2*...
	# w2 = a*y0*a*y1*...
	#
	# state = 2
	elif state = 2 and w[1][len1] = 4 then

		# modified algorithm
		#
		# reset modified algorithm, if after second 0d 
		if 3_or_4_possible = true then
			reset_modified();
		fi;

		# start counting parity
		if first_0d = true then
			# store current words and output (before stepping 0d) to possibly start modified algorithm later
			modified := [ShallowCopy(w[1]),ShallowCopy(w[2])];
			modified_output := ShallowCopy(output);
			inbetween_0d := true;
			first_0d := false;
		# stop counting parity, set first back, set marker to look for case 3 or 4x
		elif first_0d = false then
			inbetween_0d := false;
			first_0d := true;
			3_or_4_possible := true;
		fi;



		# alter (w1,w2) for next step
		Remove(w[1]);
		Remove(w[2]);
		# add to output a*c*a
		add_reduced(output, [1,3,1]);
		state := 1;


	# case 0c
	# w1 = c*a*x1*a*x2*...
	# w2 = a*y0*a*y1*...
	#
	# state = 2
	elif state = 2 and w[1][len1] = 3 then
		# alter (w1,w2) for next step
		Remove(w[1]);
		Remove(w[2]);
		# add to output a*b*a
		add_reduced(output, [1,2,1]);
		state := 1;

		# modified algorithm
		#
		# if this case occurs between two cases of 0d, or after the second 0d, 
		# reset modified algorithm
		if inbetween_0d = true or 3_or_4_possible = true then
			reset_modified();
		fi;


	# case 1
	# w1 = b*a*x1*a*x2*...
	# w2 = a*b*a*y1*...
	#
	# state = 2 	(only consider first letters)
	elif state = 2 and w[1][len1] = 2 and w[2][len2-1] = 2 then
		# alter (w1,w2) for next step
		Remove(w[1]);
		delete_last(w[2], 3);
		# add to output a*b*a*d*a*c*a
		add_reduced(output, [1,2,1,4,1,3,1]);
		state := 1;

		# modified algorithm
		# 
		# if this case occurs between two cases of 0d, or after the second 0d, 
		# reset modified algorithm
		if inbetween_0d = true or 3_or_4_possible = true then
			reset_modified();
		fi;


	# cases 2-41x
	# w1 = b*a*x1*...
	# w2 = a*y1*a*...   y1 <> b
	#
	# set state = 3
	elif state = 2 and w[1][len1] = 2 and w[2][len2-1] <> 2 then
		state := 3;


		# case 2d
		# w1 = b*a*d*a*c*a*...
		# w2 = a*d*a*c*a*y3*...
		#
		# state = 3
		elif state = 3 and w[1][len1-2] = 4 and w[2][len2-1] = 4 then
			# alter (w1,w2) for next step
			delete_last(w[1], 5);
			Add(w[1], 2);
			delete_last(w[2], 4);
			# add to output a*b*a*c*a*c*a*b
			add_reduced(output, [1,2,1,3,1,3,1,2]);
			state := 1;

			# modified algorithm
			#
			# if this case occurs between two cases of 0d, or after the second 0d, 
			# reset modified algorithm
			if inbetween_0d = true or 3_or_4_possible = true then
				reset_modified();
			fi;


		# case 2c
		# w1 = b*a*d*a*c*a*...
		# w2 = a*c*a*y2*...
		#
		# state = 3
		elif state = 3 and w[1][len1-2] = 4 and w[2][len2-1] = 3 then
			# alter (w1,w2) for next step
			delete_last(w[1], 5);
			Add(w[1], 2);
			Add(w[1], 1);
			delete_last(w[2], 3);
			# add to output a*b*a*b*a*c*a
			add_reduced(output, [1,2,1,2,1,3,1]);
			state := 1;

			# modified algorithm
			#
			# if this case occurs between two cases of 0d, or after the second 0d, 
			# reset modified algorithm
			if inbetween_0d = true or 3_or_4_possible = true then
				reset_modified();
			fi;


		# case 3
		# w1 = b*a*c*a*x2*a*...
		# w2 = a*d*a*c*a*y3*...
		#
		# state = 3
		elif state = 3 and w[1][len1-2] = 3 and w[2][len2-1] = 4 then

			# if this case occurs between two cases of 0d, reset modified algorithm
			if inbetween_0d = true then
				reset_modified();
			fi;

			# if occurs after second 0d, check parity and apply modified algorithm
			if 3_or_4_possible = true then
				state := 5;
			# prevent from words being too short for modified algo after applying this case
			else
				# alter (w1,w2) for next step
				delete_last(w[1], 4);
				delete_last(w[2], 4);
				Add(w[2], 2);
				Add(w[2], 1);
				# add to output a*d*a*c*a*b*a*c
				add_reduced(output, [1,4,1,3,1,2,1,3]);
				state := 1;
			fi;


		# case 41d
		# w1 = b*a*c*a*x2...
		# w2 = a*c*a*d*a*c*a*y4...
		#
		# state = 3
		elif state = 3 and w[1][len1-2] = 3 and w[2][len2-1] = 3  and w[2][len2-3] = 4 then
			
			# if this case occurs between two cases of 0d, reset modified algorithm
			if inbetween_0d = true then
				reset_modified();
			fi;

			# if occurs after second 0d, check parity and apply modified algorithm
			if 3_or_4_possible = true then
				state := 5;
			# prevent from words being too short for modified algo after applying this case
			else
				# alter (w1,w2) for next step
				delete_last(w[1], 3);
				Add(w[1], 2);
				delete_last(w[2], 6);
				# add to output a*b*a*b*a*c*a*c*a*c*a*b
				add_reduced(output, [1,2,1,2,1,3,1,3,1,3,1,2]);
				state := 1;
			fi;


		# case 41c
		# w1 = b*a*c*a*x2...
		# w2 = a*c*a*c*a*y3...
		#
		# state = 3
		elif state = 3 and w[1][len1-2] = 3 and w[2][len2-1] = 3  and w[2][len2-3] = 3 then
			
			# if this case occurs between two cases of 0d, reset modified algorithm
			if inbetween_0d = true then
				reset_modified();
			fi;

			# if occurs after second 0d, check parity and apply modified algorithm
			if 3_or_4_possible = true then
				state := 5;
			# prevent from words being too short for modified algo after applying this case
			else
				# alter (w1,w2) for next step
				delete_last(w[1], 3);
				Add(w[1], 2);
				Add(w[1], 1);
				delete_last(w[2], 5);
				# add to output a*b*a*b*a*c*a*b*a*c*a
				add_reduced(output, [1,2,1,2,1,3,1,2,1,3,1]);
				state := 1;
			fi;


		# cases 42-45
		# w1 = b*a*c*a*x2*...
		# w2 = a*c*a*b*y1*...
		#
		# set state = 4
		elif state = 3 and w[1][len1-2] = 3 and w[2][len2-1] = 3  and w[2][len2-3] = 2 then



			# left over cases, discuss on x2*a*x3*a*... which is of form (ba)^k*x*a*...
			# count b's to determine k

			# set count_b according to 5th position
			if w[1][len1-4] = 2 then
				count_b := 1;
			else
				count_b := 0;
			fi;
			# set ending to have the right end for step of two range
			ending := 1;
			if (len1 mod 2) = 0 then
				ending := 2;
			fi;
			last_pos := len1-4;
			# begin at 7th position from back, end at 1 or 2 from front
			for i in [len1-6,len1-8..ending] do
				# check if step-2 entry is b and if there is a b at 5th position
				if w[1][i] = 2 and w[1][len1-4] = 2 then
				count_b := count_b+1;
				# first letter is not b
				else
					break;
				fi;
			od;

			state := 4;


			
			# case 42d
			# w1 = b*a*c*a*d*a*c*a*x4...
			# w2 = a*c*a*b*a*y3...
			elif state = 4 and count_b = 0 and w[1][len1-4] = 4 then
				
				# if this case occurs between two cases of 0d, reset modified algorithm
				if inbetween_0d = true then
					reset_modified();
				fi;

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else
					# alter (w1,w2) for next step
					delete_last(w[1], 7);
					delete_last(w[2], 5);
					# add to output a*d*a*c*a*d*a*c*a*c*a*b*a*c*a*d*a*b*a
					add_reduced(output, [1,4,1,3,1,4,1,3,1,3,1,2,1,3,1,4,1,2,1]);
					state := 1;
				fi;


			# case 42c
			# w1 = b*a*c*a*c*a*x3...
			# w2 = a*c*a*b*a*y3...
			elif state = 4 and count_b = 0 and w[1][len1-4] = 3 then
				
				# count parity if in between two 0d cases
				if inbetween_0d = true then
					parity := parity + 1;
				fi;

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else
					# alter (w1,w2) for next step
					delete_last(w[1], 6);
					delete_last(w[2], 4);
					# add to output .d.c.b.c.b.c.c.d.c.d
					add_reduced(output, [1,4,1,3,1,2,1,3,1,2,1,3,1,3,1,4,1,3,1,4]);
					state := 1;
				fi;


			# case 43d
			# w1 = b*a*c*a*b*a*d*a*c*a*x5...
			# w2 = a*c*a*b*a*y3...
			elif state = 4 and count_b = 1 and w[1][len1-6] = 4 then
				
				# if this case occurs between two cases of 0d, reset modified algorithm
				if inbetween_0d = true then
					reset_modified();
				fi;

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else
					# alter (w1,w2) for next step
					delete_last(w[1], 10);
					delete_last(w[2], 2);
					# add to output .d.c.b.c.b.c.c.c.d.b
					add_reduced(output, [1,4,1,3,1,2,1,3,1,2,1,3,1,3,1,3,1,4,1,2]);
					state := 1;
				fi;



			# case 43c
			# w1 = b*a*c*a*b*a*c*a*x4...
			# w2 = a*c*a*b*a*y3...
			elif state = 4 and count_b = 1 and w[1][len1-6] = 3 then

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else
					# alter (w1,w2) for next step
					delete_last(w[1], 8);
					# add to output .d.c.b.c.d.c.b.c
					add_reduced(output, [1,4,1,3,1,2,1,3,1,4,1,3,1,2,1,3]);
					state := 1;
				fi;


			# case 44d
			# w1 = b*a*c*a*(b*a)^k*d*a*c*a*x...	k = count_b is even >= 2
			# w2 = a*c*a*b*a*y3*a*...
			elif state = 4 and count_b > 1 and (count_b mod 2) = 0 and  w[1][len1-(4 + 2*count_b)] = 4 then
				
				# if this case occurs between two cases of 0d, reset modified algorithm
				if inbetween_0d = true then
					reset_modified();
				fi;

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else
					# alter (w1,w2) for next step
					delete_last(w[1], 7+2*count_b);
					delete_last(w[2], 5);
					# add to output .d.c.b.c(.d.c.d.c)^(count_b/2-1).d.c.b.c.c.b.c.d.b.
					add_reduced(output, [1,4,1,3,1,2,1,3]);
					for i in [1..(count_b/2-1)] do
						add_reduced(output, [1,4,1,3,1,4,1,3]);
					od;
					add_reduced(output, [1,4,1,3,1,2,1,3,1,3,1,2,1,3,1,4,1,2,1]);
					state := 1;
				fi;


			# case 44c
			# w1 = b*a*c*a*(b*a)^k*c*a*x...	k = count_b is even >= 2
			# w2 = a*c*a*b*a*y3*a*...
			elif state = 4 and count_b > 1 and (count_b mod 2) = 0 and  w[1][len1-(4 + 2*count_b)] = 3 then
				
				# count parity if in between two 0d cases
				if inbetween_0d = true then
					parity := parity + 1;
				fi;

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else
					# alter (w1,w2) for next step
					delete_last(w[1], 6+2*count_b);
					delete_last(w[2], 4);
					# add to output .d.c.b.c(.d.c.d.c)^(count_b/2).b.c.c.d.c.d
					add_reduced(output, [1,4,1,3,1,2,1,3]);
					for i in [1..count_b/2] do
						add_reduced(output, [1,4,1,3,1,4,1,3]);
					od;
					add_reduced(output, [1,2,1,3,1,3,1,4,1,3,1,4]);
					state := 1;
				fi;


			# case 45d
			# w1 = b*a*c*a*(b*a)^k*d*a*c*a*x...	k = count_b is odd >= 3
			# w2 = a*c*a*b*a*y3*

			# need length at least 8+2*count_b!?
			elif state = 4 and count_b > 2 and (count_b mod 2) = 1 and Length(w[1]) > 8+2*count_b and  w[1][len1-(4 + 2*count_b)] = 4 then

				# if this case occurs between two cases of 0d, reset modified algorithm
				if inbetween_0d = true then
					reset_modified();
				fi;

				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				# prevent from words being too short for modified algo after applying this case
				else

					# alter (w1,w2) for next step
					delete_last(w[1], 8+2*count_b);
					delete_last(w[2], 2);
					# add to output .d.c.b.c(.d.c.d.c)^(count_b-1)/2.d.c.c.c.b.b
					add_reduced(output, [1,4,1,3,1,2,1,3]);
					for i in [1..(count_b-1)/2] do
						add_reduced(output, [1,4,1,3,1,4,1,3]);
					od;
					add_reduced(output, [1,4,1,3,1,3,1,3,1,2,1,2]);
					state := 1;
				fi;


			# case 45c
			# w1 = b*a*c*a*(b*a)^k*c*a*x...	k = count_b is odd >= 3
			# w2 = a*c*a*b*a*y3*a*...

			# need length at least 6+2*count_b!?
			elif state = 4 and count_b > 2 and (count_b mod 2) = 1 and Length(w[1]) > 6+2*count_b and  w[1][len1-(4 + 2*count_b)] = 3 then
				
				# if occurs after second 0d, check parity and apply modified algorithm
				if 3_or_4_possible = true then
					state := 5;
				else

					# alter (w1,w2) for next step
					delete_last(w[1], 6+2*count_b);
					# add to output .d.c.b.c(.d.c.d.c)^(count_b-1)/2.d.c.b.c
					add_reduced(output, [1,4,1,3,1,2,1,3]);
					for i in [1..(count_b-1)/2] do
						add_reduced(output, [1,4,1,3,1,4,1,3]);
					od;
					add_reduced(output, [1,4,1,3,1,2,1,3]);
					state := 1;
				fi;


		# case 5
		# w1 = b*a*b*a*x2*...
		# w2 = a*y1*a*...
		elif state = 3 and w[1][len1-2] = 2 then
			# alter (w1,w2) for next step
			delete_last(w[1], 4);
			# add to output .d.c.d.c
			add_reduced(output, [1,4,1,3,1,4,1,3]);
			state := 1;


	# modified algorithm
	# sequence of 0d, {42c,43c,44c,45c,5}, 0d, 5,...,5, {3,4}
	elif state = 5 then

		w := [ShallowCopy(modified[1]),ShallowCopy(modified[2])];
		output := ShallowCopy(modified_output);

		# modified algorithm b) - parity is even
		if parity mod 2 = 0 then
			state := 6;
		# modified algorithm c) - parity is odd
		elif parity mod 2 = 1 then
			state := 7;
		fi;



		# modified algorithm b)
		# w1 = d*(a*c*a*b)^l*a*(b*a*b*a)^l_1*b*a*c*a*x*a*...
		# w2 = a*v1*d*a*y*a*...
		# parity even
		elif state = 6 then
			# alter (w1,w2) for next step
			delete_last(w[1], (4 + 4*parity));
			delete_last(w[2], 2);
			# add to output .c.(c.b.c.d.c.b.c.d.)^l1*b.d.c
			add_reduced(output, [1,3,1]);
			for i in [1..parity/2] do
				add_reduced(output, [3,1,2,1,3,1,4,1,3,1,2,1,3,1,4,1]);
			od;
			add_reduced(output, [2,1,4,1,3]);

			parity := 0;
			state := 1;


		# modified algorithm c)
		# w1 = d*(a*c*a*b)^l*a*(b*a*b*a)^l_1*b*a*c*a*x*a*...
		# w2 = a*v1*d*a*y*a*...
		#    = a*(b*a*b*a)^l2*b*a*c*a*y*...
		# parity odd
		elif state = 7 then
			len1 := Length(w[1]);
			len2 := Length(w[2]);

			# compute l_1
			# count b's in w1 from (3 + 4*parity)th position on
			if w[1][len1-(2+4*parity)] = 2 then
				count_b := 1;
			else
				count_b := 0;
			fi;
			# set ending to have the right end for step of range two
			ending := 1;
			if (len1 mod 2) = 0 then
				ending := 2;
			fi;
			# begin at (4 + 4*parity)th position from back, end at 1 or 2 from front
			for i in [len1-(4+4*parity),len1-(6+4*parity)..ending] do
				# check if step-2 entry is b and if there is a b at 5th position
				if w[1][i] = 2 and w[1][len1-(2+4*parity)] = 2 then
				count_b := count_b+1;
				# first letter is not b
				else
					break;
				fi;
			od;

			l_1 := (count_b - 1)/2;


			# compute l2
			# count b's in w2 from 2nd position on
			if w[2][len2-1] = 2 then
				count_b := 1;
			else
				count_b := 0;
			fi;
			# set ending to have the right end for step of range two
			ending := 1;
			if (len2 mod 2) = 1 then
				ending := 2;
			fi;
			# begin at 4th position from back, end at 1 or 2 from front
			for i in [len2-3,len2-5..ending] do
				# check if step-2 entry is b and if there is a b at 5th position
				if w[2][i] = 2 and w[2][len2-1] = 2 then
				count_b := count_b+1;
				# first letter is not b
				else
					break;
				fi;
			od;

			l2 := (count_b - 1)/2;

			# compute l1
			l1 := (parity-1)/2;
	

			# alter (w1,w2) for next step
			delete_last(w[1], (10 + 8*l1 + 4*l_1));
			add_reduced(w[2], [1]);
			for i in [1..l2] do
				add_reduced(w[2], [2,1,2,1]);
			od;
			add_reduced(w[2], [3,1,4,1,4]);

			# delete_last(w[2], (6 + 4*l2));
			# add to output .c.(d.c.d.c.)^l2*(c.b.c.d.)^(2*l1)*b.b.c(.d.c.d.c)^(l_1+1).b.c
			add_reduced(output, [1,3,1]);
			for i in [1..l2] do
				add_reduced(output, [4,1,3,1,4,1,3,1]);
			od;
			for i in [1..2*l1] do
				add_reduced(output, [3,1,2,1,3,1,4,1]);
			od;
			add_reduced(output, [2,1,2,1,3]);
			for i in [1..l_1+1] do
				add_reduced(output, [1,4,1,3,1,4,1,3]);
			od;
			add_reduced(output, [1,2,1,3]);


			parity := 0;
			state := 1;


	# no case applies, terminate algorithm
	# e.g. cases 45, w[1] too short
	else
		break;
	fi;



	# if in step before type II was detected, add last a of a*W*a to output, only if state = 1
	# (first a already added in type II case)
	if state = 1 and swap = true then
		add_reduced(output, [1]);
		# swap words back
		w := w{[2,1]};
		swap := false;
	fi;

od;
fi;

return [Reversed(w[1]), Reversed(w[2]), output];
end;






####################################################################################################
# grigor_algo - Grigorchuk's algorithm
#	should be applied after Leonov's modified algorithm
####################################################################################################
# input:	list of two words, and optional output of modified Leonov algorithm
# input form: 	[w1, w2, {,output}] with wi words in list representation, 
#					optional output to which the output here will be added
# output:	output word and mistake
# output form: 	[W_Gri, z_0]
####################################################################################################

grigor_algo := function(w)

	local i, z0, z1, z1w2, w1z0, output;

	z1 := [];
	z0 := [];

	# if third parameter given, it is output from Leonov's modified algorithm
	if Length(w) = 3 then
		output := w[3];
	else
		output := [];
	fi;


	# check if w[1] is the empty list
	if Length(w[1]) = 0 and Length(w[2]) = 0 then
		return [output,z0];
	else
		# set entries of z1
		#
		# e.g.
		# w0: a * b * a * c * a * d * a
		# W0: b*(ada)*b*(aba)*b*(aca)*b
		# z1: c * 1 * c * a * c * a * c


		for i in [1..Length(w[1])] do
			# current letter is a
			if w[1][i] = 1 then
				# add b to output, c to z1
				add_reduced(output, [2]);
				add_reduced(z1, [3]);
			# current letter is b
			elif w[1][i] = 2 then
				# add a*d*a to output, 1 to z1
				add_reduced(output, [1,4,1]);
				continue;
			# current letter is c
			elif w[1][i] = 3 then
				# add a*b*a to output, a to z1
				add_reduced(output, [1,2,1]);
				add_reduced(z1, [1]);
			# current letter is d
			elif w[1][i] = 4 then
				# add a*c*a to output, a to z1
				add_reduced(output, [1,3,1]);
				add_reduced(z1, [1]);
			fi;
		od;


		# calculate z1^-1*w1
		z1w2 := Reversed(z1);

		for i in [1..Length(w[2])] do
			add_reduced(z1w2, [w[2][i]]);
		od;


		# compute z0


		for i in [1..Length(z1w2)] do
			# current letter is a
			if z1w2[i] = 1 then
				# add a*b*a to output, c to z0
				add_reduced(output, [1,2,1]);
				add_reduced(z0, [3]);
			# current letter is b
			elif z1w2[i] = 2 then
				# add d to output, 1 to z0
				add_reduced(output, [4]);
				continue;
			# current letter is c
			elif z1w2[i] = 3 then
				# add b to output, a to z0
				add_reduced(output, [2]);
				add_reduced(z0, [1]);
			# current letter is d
			elif z1w2[i] = 4 then
				# add c to output, a to z0
				add_reduced(output, [3]);
				add_reduced(z0, [1]);
			fi;
		od;


		# calculate w0*z0
		w1z0 := w[1];

		for i in [1..Length(z0)] do
			add_reduced(w1z0, [z0[i]]);
		od;
	fi;

	# if parameter true is given, add a to output
	if Length(w) >= 3 and w[3] = true then
		add_reduced(output, [1]);
	fi;

return [output, z0];
end;





####################################################################################################
# check_equality
#	checks if two elements are equal in Grigorchuk's group
####################################################################################################
# input:	list with two input words, output word in list representation
# input form: 	([[w0,w1,...],[x0,x1,...]],[[y0,y1,...],[z0,z1,...]]) with wi, xi, yi, zi in {1,2,3,4}
# output:	boolean value
# output form: 	true/false
####################################################################################################

check_equality := function(inword,outword)
local in1, in2,out1,out2, double_out;

	in1 := DecompositionOfFRElement(MappedWord(AssocWordByLetterRep(FamilyObj(a),inword[1]), GeneratorsOfGroup(F), GeneratorsOfGroup(GrigorchukGroup)));
	in2 := DecompositionOfFRElement(MappedWord(AssocWordByLetterRep(FamilyObj(a),inword[2]), GeneratorsOfGroup(F), GeneratorsOfGroup(GrigorchukGroup)));
	out1 := DecompositionOfFRElement(MappedWord(AssocWordByLetterRep(FamilyObj(a), outword[1]), GeneratorsOfGroup(F), GeneratorsOfGroup(GrigorchukGroup)));
	out2 := DecompositionOfFRElement(MappedWord(AssocWordByLetterRep(FamilyObj(a), outword[2]), GeneratorsOfGroup(F), GeneratorsOfGroup(GrigorchukGroup)));
	
return in1 = out1 and in2 = out2;
end;





####################################################################################################
# hypothesis
#	puts together "dab", "dada" and "d_to_c_type"
####################################################################################################
# input:	two words in a,b,c,d
# input form: 	(x0*x1*...,y0*y1*...) with xi,yi in {a,b,c,d}
# output:	word with first Leonov's modified and then Grigorchuk's algorithm applied on and mistake
# output form: 	z0*z1*... with zi in {a,b,c,d}
####################################################################################################

hypothesis := function(l)
	local list, last_length;

	# reduce list
	list := add_reduced([], l);
	last_length := -1;

	while Length(list) <> last_length do
		last_length := Length(list);

		# apply dada function (output is already reduced)
		list := dada(list);

		# apply dab function (don't have to reduce after that, since we replace dab -> ba(ca)^6d
		list := dab(list);

		# apply d_to_c_type (don't have to reduce after that)
		list := d_to_c_type(list);

		# apply dab again, since there is a possible case, without changing length that dab occurs again!
		list := dab(list);
	od;

return list;
end;





####################################################################################################
# Main function bri_algo
#	puts together all the parts above
####################################################################################################
# input:	two words in a,b,c,d
# input form: 	(x0*x1*...,y0*y1*...) with xi,yi in {a,b,c,d}
# output:	word with first Leonov's modified and then Grigorchuk's algorithm applied on and mistake
# output form: 	z0*z1*... with zi in {a,b,c,d}
####################################################################################################

bri_algo := function(w1, w2)
	local l, w, output_word, mistake, before, after, mistake_added, w1_list, w2_list, check_worked;

	check_worked := false;

	# convert the input words to list representation
	w1_list := LetterRepAssocWord(w1);
	w2_list := LetterRepAssocWord(w2);

	# apply hypothesis
	w := [hypothesis(w1_list), hypothesis(w2_list)];
	
	# store w[1] before algorithm applied on it
	mistake_added := ShallowCopy(w[1]);
	before := [mistake_added,ShallowCopy(w[2])];

	# apply Leonov's modified algorithm and then Grigorchuk's algorithm
	l := grigor_algo(leomod_algo(w));
	output_word := list_to_double_list(l[1]);
	mistake := l[2];

	Append(mistake_added,mistake);

	# compare [w[1]+mistake,w[2]] with output in double list

	after := [output_word[1],output_word[2]];
	


	# give promts if check was successful or not
	if check_equality(before,after) = true then
		Print(before,"\n");
		#Print("\nCheck positive - Your input words (w0,w1) give rise to the computed word (w0*z0,w1).\n\noutput below: [(w0*z0,w1),z0] = \n\n");
		check_worked := true;

	else
		#Print("\nCheck negative - Your input words (w0,w1) do not give rise to the computed word (w0*z0,w1).\n\noutput below: [(w0*z0,w1),z0] = \n\n");
	fi;

return [AssocWordByLetterRep(FamilyObj(a),l[1]),AssocWordByLetterRep(FamilyObj(a),mistake),check_worked];
end;





####################################################################################################
# Main function for input words of dc-type
#	puts together all the parts above and does not convert the input form to dc-type
####################################################################################################
# input:	two words in a,b,c,d
# input form: 	(x0*x1*...,y0*y1*...) with xi,yi in {a,b,c,d}
# output:	word with first Leonov's modified and then Grigorchuk's algorithm applied on and mistake
# output form: 	z0*z1*... with zi in {a,b,c,d}
####################################################################################################

bri_algo_hypofulfilled := function(w1, w2)
	local l, w, output_word, mistake, before, after, mistake_added, w1_list, w2_list, check_worked;

	check_worked := false;

	# convert the input words to list representation
	w1_list := LetterRepAssocWord(w1);
	w2_list := LetterRepAssocWord(w2);

	w := [w1_list, w2_list];
	
	# store w[1] before algorithm applied on it
	mistake_added := ShallowCopy(w[1]);
	before := [mistake_added,ShallowCopy(w[2])];

	# apply Leonov's modified algorithm algorithm and then Grigorchuk's algorithm
	l := grigor_algo(leomod_algo(w));
	output_word := list_to_double_list(l[1]);
	mistake := l[2];

	Append(mistake_added,mistake);

	# compare [w[1]+mistake,w[2]] with output in double list

	after := [output_word[1],output_word[2]];
	

	# give promts if check was successful or not
	if check_equality(before,after) = true then
		Print(before,"\n");
		#Print("\nCheck positive - Your input words (w0,w1) give rise to the computed word (w0*z0,w1).\n\noutput below: [(w0*z0,w1),z0] = \n\n");
		check_worked := true;

	else
		#Print("\nCheck negative - Your input words (w0,w1) do not give rise to the computed word (w0*z0,w1).\n\noutput below: [(w0*z0,w1),z0] = \n\n");
	fi;

return [AssocWordByLetterRep(FamilyObj(a),l[1]),AssocWordByLetterRep(FamilyObj(a),mistake),check_worked];
end;





####################################################################################################
# test_with_hypo
#	tests the algorithm with random words of dc-type
####################################################################################################
# input:	path, length of word 1, length of word 2, number of tests
# input form: 	(path,n1,n2,n)
# output:	list of length of output word in the file given at "path"
####################################################################################################

test_with_hypo := function(path,range1,range2,n)
	local i, j, w1, w2, output, a_added1, a_added2, length, k, l, m, sum, average, d_added1, d_added2, ending, max, max_set, min, min_set;

	AppendTo(path, "Length(w1), Length(w2), average Length(output), maximal Length(output), minimal Length(output)\n\n");

for l in range1 do
	for m in range2 do
		k := n;
		AppendTo(path, l," ",m," ");
		sum := 0;
		max_set := false;
		min_set := false;

		while k > 0 do
			w1 := [];
			w2 := [];
			d_added1 := false;
			d_added2 := false;

			# create random word w1
		
			# decide if it starts with a
			if Random([0,1]) = 0 then
				# start with a
				Add(w1, 1);
				a_added1 := true;
				i := 2;
			else 
				a_added1 := false;
				i := 1;
			fi;

			while i <= l do
				# add {b,c,d}
				if i <= l then
					# check if letter added before was d
					if d_added1 = true then
						Add(w1, 3);
						d_added1 := false;
					else 
						Add(w1, Random([2,3,4]));
					fi;
					# note if last one added was d
					if w1[Length(w1)] = 4 then
						d_added1 := true;
					fi;
				fi;
				i := i + 1;

				# add a
				if i <= l then
					Add(w1, 1);
				fi;
				i := i + 1;
			od;


			# create random word w2
		
			# decide if it starts with a
			if Random([0,1]) = 0 then
				# start with a
				Add(w2, 1);
				a_added2 := true;
				i := 2;
			else 
				a_added2 := false;
				i := 1;
			fi;

			while i <= m do
				# add {b,c,d}
				if i <= m then
					# check if letter added before was d
					if d_added2 = true then
						Add(w2, 3);
						d_added2 := false;
					else 
						Add(w2, Random([2,3,4]));
					fi;
					# note if last one added was d
					if w2[Length(w2)] = 4 then
						d_added2 := true;
					fi;
				fi;
				i := i + 1;

				# add a
				if i <= m then
					Add(w2, 1);
				fi;
				i := i + 1;
			od;

			output := algo_hypofulfilled(AssocWordByLetterRep(FamilyObj(a), w1), AssocWordByLetterRep(FamilyObj(a),w2));
			length := Length(LetterRepAssocWord(output[1]));
			sum := sum + length;

			k := k - 1;

			# maximum
			# set first result as maximum
			if max_set = false then
				max := length;
				max_set := true;
			# current result is larger than max
			elif max_set = true and max < length then
				max := length;
			fi;

			# minimum
			# set first result as minimum
			if min_set = false then
				min := length;
				min_set := true;
			# current result is smaller than min
			elif min_set = true and min > length then
				min := length;
			fi;
		od;
		AppendTo(path, Float(sum/n)," ",max," ",min,"\n");
	od;
od;

end;
