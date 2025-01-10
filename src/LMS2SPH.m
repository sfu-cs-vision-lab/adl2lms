% function RES=lms2sph(lms)
% Conversion from LMS (3d)  to spherical angles (2d)
% one lms sample per row
% Hamidreza Mirzaei
% Simon Fraser University

function RES=LMS2SPH(lms)

N=sqrt(sum(lms.^2,2));
RES=[acos(lms(:,3)./(N+(N==0)))   atan2(lms(:,2), lms(:,1))];
