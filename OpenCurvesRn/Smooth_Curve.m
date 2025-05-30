function pnew = Smooth_Curve(p)

[n,T] = size(p);
for i = 1:n
    pnew(i,:) = wden(p(i,:),'heursure','s','one',1,'sym8');
end

return;
