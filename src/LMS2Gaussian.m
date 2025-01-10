% function ADL=LMS2ADL(LMS, illum, cones, wavelengths, maxAngularError)
%
% compute Logvinenko's Gaussian coordinates for the given LMS, illumination and
% cones by interpolation and optimization
%
% LMS: sensor response values to convert, is size N*3 for N samples
% illum: the spectrum of the illuminant
% cones: the sensor responsitivities, one column per sensor
% wavelengths: the sampling points used for the above
% maxAngularError: the target maximum angular error
%    optimization is performed for the samples that have an error greater
%    than this after interpolation. This does not guarantee a resulting
%    error smaller than maxAngularError, but usually improves results.
%
%
% Hamidreza Mirzaei
% Simon Fraser University
%
%
function GParameters=LMS2Gaussian(LMS, illum, cones, wavelengths, maxAngularError)
%
%
start = wavelengths(1);
finish = wavelengths(end);
N = numel(wavelengths);


%%%%%%%%%%% Create interpolation functions %%%%%%%%%%%%
disp('Creating interpolation functions...');

ThetaStepsize =   10;
% maxAngularError=.01;

MStepsize=1.5;
tic;
i=0;
WLS1=zeros(16800,2);
for theta=.5:1:400%theta=.5:ThetaStepsize:60000
    for M=start:MStepsize:finish
        i=i+1;
        WLS1(i,:)=[1/theta/theta M];%theta=s!!
        %WLS1(i,:)=[1/theta M]; % type 1 optimal reflectance function
    end
end


toc
time1=toc;
G=ones(size(WLS1,1),1);
G(:,2:3)=WLS1;

LMS1=zeros(size(WLS1,1),3);
%
tic
for j=1:size(WLS1,1)
    refl=G2RSamples(1,WLS1(j,1),WLS1(j,2),wavelengths);
    
    
    LMS1(j,:)=(G2LMS(refl,illum, cones, wavelengths));
end
toc
time2=toc;
% % Convert to spherical coordinates to allow for 2D interpolation
SPH1=LMS2SPH(LMS1);


% Create the interpolation function for Theta and M
disp('Generating Interpolation functions...');


%     % delete duplicates for interpolation
[SPH1,rows]=unique(SPH1,'rows');
G=G(rows,:);
L1=TriScatteredInterp(SPH1(:,1), SPH1(:,2),  G(:,2),'linear'); % L1 is a function, output=S
D1=TriScatteredInterp(SPH1(:,1), SPH1(:,2),  G(:,3),'linear');%output=M

display('Interpolating... ...');
SPH=LMS2SPH(LMS);

% fprintf('\n Calculating error...');
G=ones(size(LMS,1),3);
a=L1(SPH(:,1),SPH(:,2));%a= initial point for S
a(isnan(a))=1;

G(:,2)=a;

b=D1(SPH(:,1),SPH(:,2));% initial point for M
b(isnan(b))=500;

G(:,3)=b;

%
%         for j=1:size(G,1)
%     refl=G2RSamples(1,G(j,2),G(j,3),wavelengths);
%
%
%     GLMS(j,:)=(G2LMS(refl,illum, cones, wavelengths));%
% % %     LMS2=(G2LMS(1,WLS2(1),WLS2(2),illum, cones, wavelengths));
% end
%
%         ERROR=vectorangle(GLMS,LMS);
%
%



%%%%%%%%%%%%%%%%%% now optimization to reduce angular error %%%%%%%%%%%%%%%


LMS=double(LMS);
GParameters=zeros(size(LMS,1),4);
tic
for i=1:size(LMS,1)
    % the objective functions
    
    refl=G2RSamples(1,G(i,2),G(i,3),wavelengths);% finds a spectra based on K 1/S/S M
    GLMS=G2LMS(refl,illum, cones, wavelengths);
    
    iLMS=LMS(i,:);
    err=vectorangle(GLMS,iLMS);
    if err>maxAngularError
        fprintf('\n sample %d, error %.4f : not good enough, optimizing...',i,err);
        % set the options for optimization
        % get starting values from interpolation results
        
        init=G(i,2:3);
        
        options=optimset('MaxFunEvals',50000,'MaxIter',50000);
        [out,fval,exitflag] = optimize(@(vars) mycostfunction(vars,illum, cones, wavelengths,iLMS),init,[0 wavelengths(1)],[1000 wavelengths(end)],[],[],[],[],[],options );% nopts=0;
%         [100 wavelengths(end)] --> [1 wavelength(start)] [400 wavelengths(end)]
        nopts=0;
        while fval>maxAngularError && nopts<10
            init(2)=init(2)+70;
            if(init(2)>wavelengths(end))
                init(2)=init(2)-wavelengths(end);
            end
            [out,fval,exitflag] = optimize(@(vars) mycostfunction(vars,illum, cones, wavelengths,iLMS),init,[0 wavelengths(1)],[1000 wavelengths(end)],[],[],[],[],[],options );% nopts=0;
            nopts=nopts+1;
        end
        
        refl=G2RSamples(1,out(1),out(2),wavelengths);
        GLMS=G2LMS(refl,illum, cones, wavelengths);
        % calculate alphas as ratios between target and optimal lms
        normLMS=sqrt(sum(iLMS.^2,2));
        normGLMS=sqrt(sum((GLMS).^2,2));
        GParameters(i,:)=[normLMS./normGLMS,out,fval];
        
        fprintf('\n after optimization the error is: %.4f ',fval);
        
    else
        fprintf('\n sample %d, error %.4f : good enough, use interpolation results.',i,err);
        % calculate alphas as ratios between target and optimal lms
        normLMS=sqrt(sum(iLMS.^2,2));
        normGLMS=sqrt(sum((GLMS).^2,2));
        GParameters(i,:)=[normLMS./normGLMS,G(i,2:3),err];
    end
end


toc
time3=toc;
time1
time2
time3