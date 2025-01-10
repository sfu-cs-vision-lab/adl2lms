%convert adl coors to lms coors given illum and cone fundamentals
%Brian Funt Jan 2010
% modified my Christoph Godau Apr 2010
% Logvinenko's paper Journal of Vision (2009) 9(11):5, 1–23
function LMS = ADL2LMS(adl, illum, cones, wavelengths)

start = wavelengths(1);
finish = wavelengths(end);
Nwavelengths = size(wavelengths,1);
Nsamples=size(adl,1);
Ncones=size(cones,2);

grey=(illum*0.5)'*cones;

WL=zeros(Nsamples,2);
% determine lambda1 and lambda2 for each sample
for i=1:Nsamples
    Alpha = adl(i,1);
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
    
    WL(i,:)=[L1 L2];
end

WL=WL+start-1;

% calculate LMS for optimal reflectance function
oLMS=OPT2LMS(WL, illum, cones, wavelengths); % o for optimal

% calculate final LMS using given alpha

oLMSng=oLMS-repmat(grey,Nsamples,1); % ng = subtract grey
LMS=repmat(grey,Nsamples,1)+repmat(adl(:,1),1,3).*oLMSng;