
%convert optimal reflectance function given by L1 and L2 to lms coors given 
%  illum and cone fundamentals
% Christoph Godau Feb 2010
% Logvinenko's paper Journal of Vision (2009) 9(11):5, 1–23
function LMS = OPT2LMS(WLS, illum, cones, wavelengths)

start = wavelengths(1);
finish = wavelengths(end);
Nwavelengths = size(wavelengths,1);
Nsamples=size(WLS,1);
Ncones=size(cones,2);

% precalculate illumination on cones
cones=cones.*repmat(illum,1,Ncones);

% create interpolation function for fast LMS calculation
% goal: function getsum(wl1,wl2) giving integral from wl1 to wl2
ppcones=spline([1:length(cones)], cones');
xs=linspace(1,size(cones,1),(size(cones,1)-1)*10)'; % .1nm steps
% interpolate values for all these steps
ys=ppval(ppcones,xs)';
% Calculate and interpolate cumulative integrals
cumcones=cumtrapz(xs,ys);
ppcumcones=spline(xs,cumcones');
% now the final function is just a difference
getsum=@(x1, x2) ppval(ppcumcones,x2)'-ppval(ppcumcones,x1)';

LMS =zeros(Nsamples,3);

L1=WLS(:,1)-start+1;
L2=WLS(:,2)-start+1;

T1= find(L1<L2);
T2= find(L1>L2);
if(sum(T1)) LMS(T1,:)= getsum(L1(T1),L2(T1)); end
if(sum(T2)) LMS(T2,:)= getsum(0*L2(T2)+1,L2(T2)) + getsum(L1(T2),0*L2(T2)+Nwavelengths); end