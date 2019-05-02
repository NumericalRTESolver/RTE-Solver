 %% computes the weight coefficients associated to the quadrature method
 function [tt,w,theta]=weight_biss_pr(x,N_angle,g,alph,g1,g2,mu,nn,ord,psf)
tic
w=zeros(N_angle,N_angle);
m=(N_angle+2)/2;
%%% STHG fct
if(strcmp(psf,'STHG'))
f=@(x) (1-g^2)./(2*pi*(1+g^2-2*g*cos(x))); %2D
elseif(strcmp(psf,'TTHG'))
%%% TTHG fct
 f=@(x) alph*(1-g1^2)./(2*pi*(1+g1^2-2*g1*cos(x)))+(1-alph)*(1-g2^2)./(2*pi*(1+g2^2-2*g2*cos(x)));
elseif(strcmp(psf,'FF'))
%% Fournier For. fct
v=(3-mu)/2;
dlt=@(x) 4/(3*(nn-1)^2).*(sin(x/2)).^2;
dlt_pi=4/(3*(nn-1)^2);
f=@(x)   1./(4*pi*(1-dlt(x)).^2.*dlt(x).^v).*(v.*(1-dlt(x))-(1-dlt(x).^v)+(dlt(x).*(1-dlt(x).^v)-v.*(1-dlt(x))).*(sin(x./2)).^(-2))+...
     (1-dlt_pi^v)./(16*pi*(dlt_pi-1).*dlt_pi^v).*(3.*(cos(x)).^2-1);
else
fprintf('\nInvalid entry for the PSF\n');
end

for ll=1:m-1
% for ll=1:m
h1=x(ll);
h2=x(ll+1);
M=7;
h=(h2-h1)/(M-1);
    ss=zeros(1,M);
    tt=zeros(1,M);
    uu=zeros(1,M);
    if (ord==3)
%%%%%%%%%% 3 pts scheme
    %%2 pts
    uu(1)=1/6*(2*f(h1)+f(h1+h))*h;
uu(M)=1/6*(f((M-2)*h)+2*f((M-1)*h))*h;
%% 3 points

% 3 points
for ii=2:M-1
    uu(ii)=1/12*(f(h1-(ii-2)*h)+4*f(h1+(ii-1)*h)+f(h1+ii*h))*2*h;
end
w(1,ll)=sum(uu);


    elseif (ord==5)
 %%%%% 5 pts scheme
 %%2 pts
tt(1)=1/12*(2*f(h1)+f(h1+h))*h;
tt(M)=1/12*(f((M-2)*h)+2*f((M-1)*h))*h;
%% 3 points
tt(2)=1/12*(f(h1)+4*f(h1+h)+f(h1+2*h))*2*h;
tt(M-1)=1/12*(f(h2-2*h)+4*f(h2-h)+f(h2))*2*h;

%% 5 points
for ii=3:M-2
    tt(ii)=1/360*(7*f(h1+(ii-3)*h)+32*f(h1+(ii-2)*h)+12*f(h1+(ii-1)*h)+32*f(h1+(ii)*h)+7*f(h1+(ii+1)*h))*(4*h);
end
w(1,ll)=sum(tt);

    elseif (ord==7)
% %%%%%% 7 pts scheme

 %%2 points
ss(1)=1/18*(2*f(h1)+f(h1+h))*h;
ss(M)=1/18*(f((M-2)*h)+2*f((M-1)*h))*h;
%% 3 points
ss(2)=1/36*(f(h1)+4*f(h1+h)+f(h1+2*h))*2*h;
ss(M-1)=1/36*(f(h2-2*h)+4*f(h2-h)+f(h2))*2*h;
%% 5 points
ss(3)=1/90*1/2*(7*f(h1)+32*f(h1+h)+12*f(h1+2*h)+32*f(h1+3*h)+7*f(h1+4*h))*(4*h);
ss(M-2)=1/90*1/2*(7*f(h2-4*h)+32*f(h2-3*h)+12*f(h2-2*h)+32*f(h2-h)+7*f(h2))*(4*h);
%% 7 points
for ii=4:M-3
    ss(ii)=1/840*1/6*(41*f(h1+(ii-4)*h)+216*f(h1+(ii-3)*h)+27*f(h1+(ii-2)*h)+272*f(h1+(ii-1)*h)+27*f(h1+(ii)*h)+216*f(h1+(ii+1)*h)+41*f(h1+(ii+2)*h))*(6*h);
end

w(1,ll)=sum(ss(1:M-1));

    else
        fprintf('the provided order number is neither 3, nor 5/7 !!! provide one of those integers')
    end
w(1,2*m-1-ll)=w(1,ll);

end

w(1,:)=w(1,:)/sum(w(1,:));
for ii=2:N_angle
for j=1:N_angle
w(ii,j)=w(1,abs(ii-j)+1);
end
end

theta=zeros(N_angle,2);

for ii=1:N_angle

theta(ii,1)=cos(x(ii));
 theta(ii,2)=sin(x(ii));
end
tt=toc
end