
fun2 = @(x) x.^(2-1) .* exp(-x);
fun4 = @(x) x.^(4-1) .* exp(-x);
fun6 = @(x) x.^(6-1) .* exp(-x);
fun8 = @(x) x.^(8-1) .* exp(-x);
fun10 = @(x) x.^(10-1) .* exp(-x);
[gamma2,neval2] = quad(fun2,0,5*2)
[gamma4, neval4] = quad(fun4,0,5*4)
[gamma6, neval6] = quad(fun6,0,5*6)
[gamma8, neval8] = quad(fun8,0,5*8)
[gamma10, neval10] = quad(fun10,0,5*10)

2,57 ; 4,97; 6, 177