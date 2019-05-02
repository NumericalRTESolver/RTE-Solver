%% iterative computing of the received radiance
function [y]=gauss_rte_biss_p(w,theta,length_y,length_z,step_y,step_z,K,q,b,c)
I=(length_y/step_y)+1;
J=(length_z/step_z)+1;
radiance=zeros(I,J,K);
radiance_temp=zeros(I,J,K);

for l=1:320
for k=1:K
denominator1(k)=2*theta(k,2)/(3*step_y)+2*theta(k,1)/(3*step_z)+c;
denominator2(k)=2*theta(k,2)/(3*step_y)-2*theta(k,1)/(3*step_z)+c;
denominator3(k)=-2*theta(k,2)/(3*step_y)-2*theta(k,1)/(3*step_z)+c;
denominator4(k)=-2*theta(k,2)/(3*step_y)+2*theta(k,1)/(3*step_z)+c;

if theta(k,1)>0 && theta(k,2)>0

%%case 1
for i=1:I
for j=1:J
for n=1:K
sum1(n)=radiance(i,j,n)* w(k,n);
end

if (i==1)
    deriv_i=0;
elseif (i==2)
      deriv_i=radiance(i-1,j,k);
else
     deriv_i=radiance(i-1,j,k)+radiance(i-2,j,k);
end

if (j==1)
    deriv_j=0;
elseif (j==2)
      deriv_j=radiance(i,j-1,k);
else
     deriv_j=radiance(i,j-1,k)+radiance(i,j-2,k);
end

numerator(i,j,k)=sum(sum1)*b+(deriv_i)*(theta(k,2)/(3*step_y))+(deriv_j)*(theta(k,1)/(3*step_z))+c*q(i,j,k);
radiance_temp(i,j,k)=numerator(i,j,k)/denominator1(k);
end
end

elseif theta(k,1)<0 && theta(k,2)>0
%%case 2

for i=1:I
for j=1:J

for n=1:K
sum1(n)=radiance(i,j,n)* w(k,n);
end

if (i==1)
    deriv_i=0;
elseif (i==2)
      deriv_i=radiance(i-1,j,k);
else
     deriv_i=radiance(i-1,j,k)+radiance(i-2,j,k);
end

if (j==J)
    deriv_j=0;
elseif (j==J-1)
      deriv_j=radiance(i,j+1,k);
else
     deriv_j=radiance(i,j+1,k)+radiance(i,j+2,k);
end


numerator(i,j,k)=sum(sum1)*b+(deriv_i)*(theta(k,2)/(3*step_y))-(deriv_j)*(theta(k,1)/(3*step_z))+c*q(i,j,k);

radiance_temp(i,j,k)=numerator(i,j,k)/denominator2(k);
end
end
elseif theta(k,1)<0 && theta(k,2)<0
%%case 3
for i=1:I
for j=1:J
for n=1:K
sum1(n)=radiance(i,j,n)* w(k,n);
end

if (i==I)
    deriv_i=0;
elseif (i==I-1)
      deriv_i=radiance(i+1,j,k);
else
     deriv_i=radiance(i+1,j,k)+radiance(i+2,j,k);
end

if (j==J)
    deriv_j=0;
elseif (j==J-1)
      deriv_j=radiance(i,j+1,k);
else
     deriv_j=radiance(i,j+1,k)+radiance(i,j+2,k);
end


numerator(i,j,k)=sum(sum1)*b-(deriv_i)*(theta(k,2)/(3*step_y))-(deriv_j)*(theta(k,1)/(3*step_z))+c*q(i,j,k);

radiance_temp(i,j,k)=numerator(i,j,k)/denominator3(k);
end
end
else

%%case 4

for i=1:I
for j=1:J
for n=1:K
sum1(n)=radiance(i,j,n)* w(k,n);
end

if (i==I)
    deriv_i=0;
elseif (i==I-1)
      deriv_i=radiance(i+1,j,k);
else
     deriv_i=radiance(i+1,j,k)+radiance(i+2,j,k);
end

if (j==1)
    deriv_j=0;
elseif (j==2)
      deriv_j=radiance(i,j-1,k);
else
     deriv_j=radiance(i,j-1,k)+radiance(i,j-2,k);
end

numerator(i,j,k)=sum(sum1)*b-(deriv_i)*(theta(k,2)/(3*step_y))+(deriv_j)*(theta(k,1)/(3*step_z))+c*q(i,j,k);

radiance_temp(i,j,k)=numerator(i,j,k)/denominator4(k);
end
end
end
end
radiance=radiance_temp;
end

y=radiance;