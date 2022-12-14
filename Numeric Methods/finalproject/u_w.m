function out = u_w(z, f)

f = f(2:end);
x = floor(z);
p = z - x;
if x==0 % f(0)=0
    out = p*f(x+1) + 0.5*p*(p-1)*(f(x+2)-2*f(x+1));
    return
elseif x==length(f)-3 || x == length(f)-2 || x == length(f)-1
    out = f(x) + p *(f(x)-f(x-1));
    return
end
out = f(x) + p*(f(x+1)-f(x)) + 0.5*p*(p-1)*(f(x+2)-2*f(x+1)+f(x));
end