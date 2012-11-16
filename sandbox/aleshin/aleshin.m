n = 12
a = sparse([1],[1],[1],1,1);
b = a;
c = a;
for i=1:n,
    olda = a;
    oldb = b;
    oldc = c;
    a = kron([[0,1];[0,0]],oldb) + kron([[0,0];[1,0]],oldc);
    b = kron([[0,1];[0,0]],oldc) + kron([[0,0];[1,0]],oldb);
    c = kron([[1,0];[0,1]],olda);
end
m = a+b+c+transpose(a)+transpose(b)+transpose(c);