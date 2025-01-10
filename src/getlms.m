function LMS=getlms(ppcumcones,wavelengths,L1,L2)
start = wavelengths(1);
finish = wavelengths(end);
N = numel(wavelengths);


% convert x back to normal "optimal" wavelength range
x=[L1,L2];
x=x-N+1;
x=x+(x<1)*(N-1) - (x>(N))*(N-1)+start-1;
L1=x(:,1);
L2=x(:,2);

% now the final function is just a difference
getsum=@(x1, x2) ppval(ppcumcones,x2)'-ppval(ppcumcones,x1)';

LMS =zeros(size(L1,1),3);

L1=L1-start+1;
L2=L2-start+1;

T1= find(L1<L2);
T2= find(L1>L2);
if(sum(T1)) LMS(T1,:)= getsum(L1(T1),L2(T1)); end
if(sum(T2)) LMS(T2,:)= getsum(0*L2(T2)+1,L2(T2)) + getsum(L1(T2),0*L2(T2)+N); end