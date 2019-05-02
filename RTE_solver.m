%% Main code for the RTE Solver
function y=RTE_solver(c,albedo,start,length_y,length_z,step_y,step_z,aperture,K,g,alph,g1,g2,mu,nn,ord,psf)
b=c*albedo;
I=(length_y/step_y)+1;
M=floor(aperture/step_y)+1;
r(1)=step_y/2;
s=zeros(1,(M-1)/2+1);
s(1,1)=pi*r(1)^2;
for n=2:(M-1)/2+1
r(n)=r(n-1)+step_y;
s(n)=pi*r(n)^2-s(n-1);
end
    J=(length_z/step_z)+1;
    q=zeros(I,J,K);
q(((I-1)/2)+1,1,1)=1; %% src at the middle
%% Non-Uniform discretization
phi=rte_angle_dist_org_pr_biss(g,alph,g1,g2,mu,nn,1+K/2,psf);
%% Computing weight coefficients of the integral term
 [tt,w, theta]=weight_biss_pr(phi,K,g,alph,g1,g2,mu,nn,ord,psf);
 %% iterative computing of the received radiance
 [radiance]=gauss_rte_biss_p(w,theta,length_y,length_z,step_y,step_z,K,q,b,c);
 intensity=zeros(I,J);
  phi(K+1)=2*pi+0.0015/2;
    for k=1:K
%% Field of view test
if ((phi(k)>=0 && phi(k)<=pi))
intensity(:,:)=intensity(:,:)+radiance(:,:,k)*(phi(k+1)-phi(k));
end
    end
%% power computation
for jj=1:J
powerr(jj)=s*intensity((I-1)/2+1:(M-1)/2+(I-1)/2+1,jj);
end
for ii=1:1:J-start/step_z
powertrc(ii)=powerr(ii+start/step_z);
end
 tm=toc
z=start:step_z:length_z;
semilogy(z,powertrc_7./(pi*10^-6),'r-.')
grid on;
