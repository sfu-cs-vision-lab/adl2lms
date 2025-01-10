% convert an optimal reflectance function given by 
% wavelengths L1 and L2 (columns in WL) to ADL notation
% opt: n*2 matrix
function ADL=OPT2ADL(opt,wavelengths)

start=wavelengths(1);
finish=wavelengths(end);
wlrange=finish-start;
wlmiddle=(finish+start)/2;

ADL=ones(size(opt,1),3);
for i=1:size(opt,1)
    L1=opt(i,1);
    L2=opt(i,2);
    if L1<=L2% type 1
        ADL(i,2)=L2-L1;
        ADL(i,3)=(L1+L2)/2;
    else %type 2
        ADL(i,2)=wlrange-(L1-L2);
        ADL(i,3)=(L1+L2)/2;
        if ADL(i,3)<=wlmiddle
            ADL(i,3)=ADL(i,3)+wlrange/2;
        else
            ADL(i,3)=ADL(i,3)-wlrange/2;
        end
    end
end
        