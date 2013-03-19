function Phi = buildPhi(y,u,na,nb)

Phi = zeros(length(y),na+nb);
for(i = 1:na)
    Phi(i+1:end,i) = -y(1:end-i);
end

for(i = 1:nb)
    Phi(i+1:end, i+na) = u(1:end-i);
end
