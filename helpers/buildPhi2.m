function Phi = buildPhi2(y,u,na,nb,y0,u0)

if(nargin == 4)
    y0 = zeros(na,1);
    u0 = zeros(nb,1);
end
T = length(y);

% Add old values to the data vectors
y = [y0(end-na+1:end) ; y];
u = [u0(end-nb+1:end) ; u];

Phi = zeros(T,na+nb);
for(i = 1:na)
    Phi(:,i) = -y(na-i+1:na-i+T);
end

for(i = 1:nb)
    Phi(:, i+na) = u(nb-i+1:nb-i+T);
end
