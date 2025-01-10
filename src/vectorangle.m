%Compute the angle in degrees between two vectors
function va = vectorangle(vs,ws)
for i=1:size(vs,1)
    v=vs(i,:);
    w=ws(i,:);
    va(i,1)=acosd(dot(v,w)/sqrt(sum(v(:).^2))/sqrt(sum(w(:).^2)));
end
va=real(va);