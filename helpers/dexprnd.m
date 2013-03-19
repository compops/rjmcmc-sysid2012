function x = dexprnd(lambda)
% Samples from a discrete exponential distribution. Probably the worst
% possible way to do this, but it works ok...

const = 1/(1-exp(-lambda))-0.5; % Normalising constant
U = const*rand(1);
x = 0;
y = 0.5;
while(U > y)
    x = x + 1;
    y = y + exp(-lambda*x);
end

U = 2*((rand(1)>0.5)-0.5);
x = x*U;
