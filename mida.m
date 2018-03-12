clear;
Nk=100;
Nl=100;
N=Nk+Nl+1;
ground_ket=ket(N,1);
ground_ene=10;
deltaE=100;
couple_onek=0.5;
couple_lk_ene=0.0;
ek=enerspan(Nk,ground_ene,deltaE);
el=enerspan(Nl,ground_ene,deltaE);
hamiall=ground_ene*(ground_ket*ground_ket');
%%% k manifold%%%
for i=1:Nk
   hamiall=hamiall+ek(i)*ket(N,i+1)*(ket(N,i+1))';
end
%%% l manifold%%%
for i=1:Nl
   hamiall=hamiall+el(i)*ket(N,i+Nk+1)*(ket(N,i+Nk+1))';
end
%%% coulping over k manifold %%%
couple=zeros(N);
for i=1:Nk
    couple=couple+couple_onek*(ket(N,1)*(ket(N,i+1))'+ket(N,i+1)*(ket(N,1))');
end
%%% coulping between k and l manifold %%%
couplekl=zeros(N);
for l=1+1:Nk+1
    for k=1+1+Nk:Nk+Nl+1
        couplekl=couplekl+couple_lk_ene*sign(rand(1)-0.5)*(ket(N,l)*(ket(N,k))'+ket(N,k)*(ket(N,l))');
    end
end
hamiall=hamiall+couple+couplekl;
t=[0:0.1:10];
for i=1:1:length(t)
    c(i)=abs(ground_ket'*expm(-1*sqrt(-1)*hamiall*t(i))*ground_ket)^2;
end
plot(t,c);
xlabel("time/s")
ylabel("possibility on 1")