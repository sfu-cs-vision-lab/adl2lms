% function ADL=LMS2ADL(LMS, illum, cones, wavelengths, maxAngularError)
%
% compute Logvinenko's ADL coordinates for the given LMS, illumination and
% cones by interpolation and optimization
%
% LMS: sensor response values to convert, is size N*3 for N samples
% illum: the spectrum of the illuminant
% cones: the sensor responsitivities, one column per sensor
% wavelengths: the sampling points used for the above
% maxAngularError: the target maximum angular error (origin at grey)
%    optimization is performed for the samples that have an error greater
%    than this after interpolation. This does not guarantee a resulting
%    error smaller than maxAngularError, but usually improves results.
%
%
% Feb. 2010 Christoph Godau & Brian Funt
% Logvinenko's paper Journal of Vision (2009) 9(11):5, 1–23

function ADL=LMS2ADL(LMS, illum, cones, wavelengths, maxAngularError, stepsize)

if nargin < 6
    stepsize = 1;
end

start = wavelengths(1);
finish = wavelengths(end);
N = numel(wavelengths);

% calculate sensor response of the grey center
grey=( (illum*0.5)'*cones );

% substract grey LMS from all input LMS values, moving origin to grey
LMS=LMS-repmat(grey,size(LMS,1),1);

%%%%%%%%%%% Create interpolation functions %%%%%%%%%%%%
    disp('Creating interpolation functions...');
    % Preallocate "enough" space
    numSteps=ceil((finish-start)/stepsize);
    WLS1=zeros(numSteps^2,2);   %% Type 1 Wavelengths
    WLS2=WLS1;  %% Type 2 Wavelengths
    i=0;
    for L1=start:stepsize:finish
        for L2=(L1+stepsize/2):stepsize:finish
            i=i+1;
            WLS1(i,:)=[L1 L2]; % type 1 optimal reflectance function
            WLS2(i,:)=[L2 L1]; % type 2
        end
    end
    % Delete the unused allocated rows
    WLS1(i+1:end,:)=[];
    WLS2(i+1:end,:)=[];
    % Calculate ADL for the optimal reflectance functions of both types
    ADL1=OPT2ADL(WLS1,wavelengths); %% Type 1
    ADL2=OPT2ADL(WLS2,wavelengths); %% Type 2
    
    % calculate XYZ from those results, subtract grey
    LMS1=(OPT2LMS(WLS1,illum, cones, wavelengths))-repmat(grey,size(WLS1,1),1);
    LMS2=(OPT2LMS(WLS2,illum, cones, wavelengths))-repmat(grey,size(WLS2,1),1);
    
    % Convert to spherical coordinates to allow for 2D interpolation
    SPH1=LMS2SPH(LMS1);
    SPH2=LMS2SPH(LMS2);
    
    % Create the interpolation function for lambda and delta
    disp('Generating Interpolation functions...');
    % delete duplicates for interpolation
    [SPH1,rows]=unique(SPH1,'rows');
    ADL1=ADL1(rows,:);
    [SPH2,rows]=unique(SPH2,'rows');
    ADL2=ADL2(rows,:);    
    L1=TriScatteredInterp(SPH1(:,1), SPH1(:,2), ADL1(:,3),'linear');
    D1=TriScatteredInterp(SPH1(:,1), SPH1(:,2), ADL1(:,2),'linear');
    L2=TriScatteredInterp(SPH2(:,1), SPH2(:,2), ADL2(:,3),'linear');
    D2=TriScatteredInterp(SPH2(:,1), SPH2(:,2), ADL2(:,2),'linear');

% try type 1 and 2 interpolation, choose the better one (smaller angular error)...
    display('Interpolating... ...');
    SPH=LMS2SPH(LMS);

    fprintf('\n Calculating Type 1 error...');
        T1ADL=ones(size(SPH,1),3);
        T1ADL(:,3)=L1(SPH); T1ADL(isnan(T1ADL))=500;
        T1ADL(:,2)=D1(SPH); T1ADL(isnan(T1ADL))=0;
        T1LMS=(ADL2LMS(T1ADL,illum, cones, wavelengths))-repmat(grey,size(T1ADL,1),1);
        T1ERROR=vectorangle(T1LMS,LMS);
        
    fprintf('\n Calculating Type 2 error...');
        T2ADL=ones(size(SPH,1),3);
        T2ADL(:,3)=L2(SPH); T2ADL(isnan(T2ADL))=500;
        T2ADL(:,2)=D2(SPH); T2ADL(isnan(T2ADL))=0;
        T2LMS=(ADL2LMS(T2ADL,illum, cones, wavelengths))-repmat(grey,size(T2ADL,1),1);
        T2ERROR=vectorangle(T2LMS,LMS);
    
    % For each sample pick the type with a smaller error
    TYPE1 = repmat((T1ERROR<=T2ERROR),1,3);
    TYPE2 = 1-TYPE1;
    ADL = TYPE1.*T1ADL + TYPE2.*T2ADL;
        
    
        
        
        
%%%%%%%%%%%%%%%%%% now optimization to reduce angular error %%%%%%%%%%%%%%%

% Extend the wavelength range to make everything type 1
illcones=cones.*repmat(illum,1,3);
% econes=[illcones(1:end-1,:);illcones;illcones(2:end,:)]; % e for extended
% 
% % create interpolation function for fast LMS calculation
% % goal: function getelms(wl1,wl2) for extended wl range
%     ppcones=spline(1:length(econes), econes');
%     xs=linspace(1,size(econes,1),(size(econes,1)-1)*10)'; % .1nm steps
%     % interpolate values for all these steps
%     ys=ppval(ppcones,xs)';
%     cumcones=cumtrapz(xs,ys);
%     ppcumcones=spline(xs,cumcones');
%     % and now the final function to calculate LMS
%     getlms=@(x1, x2) (x1<x2)*ppval(ppcumcones,x2)'-ppval(ppcumcones,x1)';
%     
%     
    
    
% create interpolation function for fast LMS calculation
% goal: function getsum(wl1,wl2) giving integral from wl1 to wl2
ppcones=spline([1:length(cones)], illcones');
xs=linspace(1,size(cones,1),(size(cones,1)-1)*10)'; % .1nm steps
% interpolate values for all these steps
ys=ppval(ppcones,xs)';
% Calculate and interpolate cumulative integrals
cumcones=cumtrapz(xs,ys);
ppcumcones=spline(xs,cumcones');
    
    
    


for i=1:size(LMS,1)
    % the objective functions
    iLMS=LMS(i,:);
    
    % get 2 perpendicular vectors
    %p1=cross(rand(1,3),iLMS);
    %p2=cross(p1,iLMS);
    
    
    f=@(WL)vectorangle(iLMS,(getlms(ppcumcones,wavelengths,WL(1),WL(2)))-grey);
    %f=@(WL) dot(p1+p2,getlms(ppcumcones,wavelengths,WL(1),WL(2))-grey).^2;
    
    % get starting values from interpolation results
    x0=[ADL(i,3)-ADL(i,2)/2 ADL(i,3)+ADL(i,2)/2];
    x0=x0- start+1 + (size(cones,1)-1);
    
    % bounds for optimization
    %lb=[1 1];
    %ub=[size(econes,1) size(econes,1)];
    
    if f(x0)>maxAngularError
        fprintf('\n sample %d, error %.4f : not good enough, optimizing...',i,f(x0));
        % set the options for optimization
        options=optimset('MaxFunEvals',50000,'MaxIter',50000);
        x(i,:)=optimize(f, x0, [], [], [], [], [], [], [], [], options);
        nopts=0;
        while f(x(i,:))>maxAngularError && nopts<5
            x0=x0+rand(1,2)*20;
            xn=optimize(f, x0, [], [], [], [], [], [], [], [], options);
            if f(xn)<f(x(i,:))
                x(i,:)=xn;
            end
            nopts=nopts+1;
        end
            
        %x(i,:)=optimize(f,x0);
        if (x(i,2)-x(i,1))>(finish-start) 
            x(i,1)=x(i,1)+(x(i,2)-x(i,1))-(finish-start)+0.0001; 
        end
        if x(i,1)>x(i,2)
            x(i,1)=x(i,2);
        end
        fprintf('\n after optimization the error is: %.4f ',f(x(i,:)));
    else
        fprintf('\n sample %d, error %.4f : good enough, use interpolation results.',i,f(x0));
        x(i,:)=x0;
    end
    % save the angular error
    errors(i,1)=f(x(i,:));
end

% convert x back to normal wavelength range
x=x-N+1;
x=x+(x<1)*(N-1) - (x>(N))*(N-1)+start-1;
ADL=[];
ADL=OPT2ADL(x, wavelengths);

% calculate alphas as ratios between target and optimal lms
normLMS=sqrt(sum(LMS.^2,2));
normOPT=sqrt(sum(((ADL2LMS(ADL,illum,cones,wavelengths))-repmat(grey,size(ADL,1),1)).^2,2));
ADL(:,1)=normLMS./normOPT;



