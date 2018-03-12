function es=enerspan(N,center,deltaE)
es=zeros(N,1);
for i=1:N
   es(i)=i*deltaE;
end
es=es-(es(1)+es(N))/2;
es=es+center;
end
