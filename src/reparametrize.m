% Calculate reparametrized coordinates from ADL

% Christoph Godau March 2010

% Based on Logvinenko's 2009 paper

function rADL=reparametrize(adl, cmf, wavelengths, ill)

if nargin > 3
    cmf=cmf.*repmat(ill,1,3);
end

% create the reparametrization function
    sigma=cumtrapz(sqrt(sum(cmf.^2,2)));
    sigma=sigma/max(sigma);    
    omegaspline=spline(wavelengths,sigma);
    omega=@(wl) ppval(omegaspline,wl);


start = wavelengths(1);
finish = wavelengths(end);
Nwavelengths = size(wavelengths,1);
Nsamples=size(adl,1);
Ncones=size(cmf,2);

WL=zeros(Nsamples,2);

% determine lambda1 and lambda2 for each sample
for i=1:Nsamples
    rADL(i,1) = adl(i,1);
    
    Delta = adl(i,2);
    Lcentral = adl(i,3);
   
    Lc = Lcentral-start+1;   % lambda central converted to 1 to N range
    if Lc>Nwavelengths || Lc<1
        disp(sprintf('Lcentral out of range: %.4f',Lc));
    end
    if Delta>(finish-start) || Delta < 0
        disp(sprintf('Delta out of range: %.4f',Delta));
    end
    Deltahalf=(Delta/2);
    L1=Lc-Deltahalf;    %Compute begin and end wavelengths of rectangle
    L2=Lc+Deltahalf;
    
    % Do the wraparounds and determine type
    if L1<1
        L1=L1+Nwavelengths-1;
    end
    if L2>Nwavelengths 
        L2=L2-Nwavelengths+1;
    end
    
    om1=omega(L1+start-1);    
    om2=omega(L2+start-1);
    
    if om1<=om2
        delta=abs(om1-om2);
        lambda=(om1+om2)/2;
    else
        delta=1-abs(om1-om2);        
        if om1+om2<1
            lambda=(om1+om2+1)/2;
        else
            lambda=(om1+om2-1)/2;
        end
    end
        
    rADL(i,2:3)=[delta lambda];
end