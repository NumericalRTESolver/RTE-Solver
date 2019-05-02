%% Non-uniform angles discretization
function y=rte_angle_dist_org_pr_biss(g,alph,g1,g2,mu,nn,N,psf)
%% 1-initializing directions uniformly
maxit=550; % max number of iterations: MSE convergence criterion
if (strcmp(psf,'STHG')==1)
%%% STHG fct
% 2D
bta1=@(teta) teta.*(1-g.^2)./(2.*pi.*(1+g.^2-2.*g.*cos(teta)));

bta2=@(teta) (1-g.^2)./(2.*pi.*(1+g.^2-2.*g.*cos(teta)));
elseif (strcmp(psf,'TTHG')==1)
%% 3D
% %% TTHG fct
bta1=@(teta) teta.*(alph*(1-g1.^2)./(2.*pi.*(1+g1.^2-2.*g1.*cos(teta)))+(1-alph)*(1-g2.^2)./(2.*pi.*(1+g2.^2-2.*g2.*cos(teta))));

bta2=@(teta) alph*(1-g1.^2)./(2.*pi.*(1+g1.^2-2.*g1.*cos(teta)))+(1-alph)*(1-g2.^2)./(2.*pi.*(1+g2.^2-2.*g2.*cos(teta)));
%   %% Fournier-Forand fct
elseif (strcmp(psf,'FF')==1)
 v=(3-mu)/2;
dlt=@(teta) 4/(3*(nn-1)^2).*(sin(teta/2)).^2;
dlt_pi=4/(3*(nn-1)^2);
 bta1=@(teta) teta.* 1./(4*pi*(1-dlt(teta)).^2.*dlt(teta).^v).*(v.*(1-dlt(teta))-(1-dlt(teta).^v)+(dlt(teta).*(1-dlt(teta).^v)-v.*(1-dlt(teta))).*(sin(teta./2)).^(-2))+...
     (1-dlt_pi^v)./(16*pi*(dlt_pi-1).*dlt_pi^v).*(3.*(cos(teta)).^2-1);

 bta2=@(teta)  1./(4*pi*(1-dlt(teta)).^2.*dlt(teta).^v).*(v.*(1-dlt(teta))-(1-dlt(teta).^v)+(dlt(teta).*(1-dlt(teta).^v)-v.*(1-dlt(teta))).*(sin(teta./2)).^(-2))+...
     (1-dlt_pi^v)./(16*pi*(dlt_pi-1).*dlt_pi^v).*(3.*(cos(teta)).^2-1);
else
    fprintf('\n Invalid entry for the PSF\n');
end
%% 1-  Uniform distribution of angles
phi=rte_unif_dist_org(N);

%% 2- computing t(k) vector for k=1:N-1
t=zeros(1,2*(N-1)+1);
t(1)=0.0015/2; % t(0)=0


for kk=1:maxit

    for jj=1:2*(N-1)-1
   t(jj+1)=(phi(jj)+phi(jj+1))/2;
    end
    t(2*(N-1)+1)=(phi(2*(N-1))+2*pi)/2;

    for ii=1:2*(N-1)
   phi(ii)=integral (bta1,t(ii),t(ii+1))./integral(bta2,t(ii),t(ii+1));
    end

end
y=phi;
end
