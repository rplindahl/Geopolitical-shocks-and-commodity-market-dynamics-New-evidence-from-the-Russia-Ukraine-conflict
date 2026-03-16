function [bhat,badj,bias] = varest_k(y,p);
[bigt,n]=size(y);
w=y(p+1:bigt,:);
q=ones(bigt-p,1); for i=1:p; q=[q y(p+1-i:bigt-i,:)]; end;
bhat=inv(q'*q)*q'*w;
res=w-(q*bhat);

for imc=1:1000;
    sboot=ceil(rand(1,1)*(bigt-p+1));
    jboot=ceil(rand(bigt-p,1)*(bigt-p));
    resboot=res(jboot,:);
    for j=1:bigt-p;
        if j==1;
            q00=1; for i=1:p; q00=[q00 y(sboot+p-i,:)]; end; qboot(j,:)=q00;
        else;
            qboot(j,:)=[1 yboot(j-1,:) qboot(j-1,2:end-n)];
        end;
        yboot(j,:)=(qboot(j,:)*bhat)+resboot(j,:);
    end;
    bboot(:,:,imc)=inv(qboot'*qboot)*qboot'*yboot;
end;

bias=squeeze(mean(bboot,3))-bhat;
badj=adjfunc(bhat,bias,n,p);   %This does the adjustment that prevents the bias-adjustment from pushing the system into the nonstationary region
