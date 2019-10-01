%Bessel function is a mathematical function needed to generate the Kaiser window in FIR filter
function [ Ix ] = getBessel(x)% k is considered as 500
sum=1;
for k=1:500
    sum=sum+((1/factorial(k))*(x/2).^k).^2;
end
Ix=sum;
end
