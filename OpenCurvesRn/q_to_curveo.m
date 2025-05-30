function p = q_to_curveo(q,transVect,refc)
%
% transVect - is the translation of the expression
%
%

if nargin<2
    transVect = zeros(size(q, 1), 1);
end

[n,T] = size(q);


for i = 1:T
    qnorm(i) = norm(q(:,i),'fro');
end


for i = 1:n
    p(i,:) = [ cumtrapz( q(i,:).*qnorm )/(T) ] ;
end

%size(transVect)
%size(repmat(transVect, 1, T))
% size(p)
% size(transVect)
% T
% pause
p = p + repmat(transVect, 1, T);

error=p-refc;


for i=1:n
    for j=2:T-1
        if(abs(error(i,j))>0.00005)
            p(i,j)=p(i,j)-error(i,j);
        end
    end
end


