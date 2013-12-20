clear all
close all
clc

% run /Users/sebastien/Desktop/cvx/cvx_setup

N = 8;

s = sym('s',[N  ,1]);
u = sym('u',[N-1,1]);

v = [1];
for k = 1:N-1
    v = [v;s(k);u(k)];
end
v = [v;s(N)];

X = v*v.';

x0 = 5;
umin = -21;
umax =  21;

%Stage constraint
f = [s(1) - x0];h = [];
for k = 1:N-1
    f = [f;
         s(k+1) - s(k)^2 - u(k)]; %Dynamics
     
    h = [h;
         u(k) - umax;               %Input bounds
         umin - u(k)];
end

Cost = s.'*s + 0.1*u.'*u;

% Construct matrix state constraint
for k = 1:N
    df(:,k)     = jacobian(f(k),   v(2:end)).';
    d2f(:,:,k)  = jacobian(df(:,k),v(2:end));
    Cstf(k)      = simple(f(k) - df(:,k).'*v(2:end) - 0.5*v(2:end).'*d2f(:,:,k)*v(2:end));
end

for k = 1:N-1
    dh(:,k)     = jacobian(h(k),   v(2:end)).';
    d2h(:,:,k)  = jacobian(dh(:,k),v(2:end));
    Csth(k)      = simple(h(k) - dh(:,k).'*v(2:end) - 0.5*v(2:end).'*d2h(:,:,k)*v(2:end));
end

% Construct Cost matrix
dCost   = jacobian(Cost,  v(2:end)).';
d2Cost  = jacobian(dCost, v(2:end));

for k = 1:N-1
    eval(['s',num2str(k),' = 0;'])
    eval(['u',num2str(k),' = 0;'])
end
eval(['s',num2str(N),' = 0;'])


for k = 1:N-1                     
    Qiconst(:,:,k) = eval([Csth(k)      0.5*dh(:,k).';
                          0.5*dh(:,k) 0.5*d2h(:,:,k)]);
                      
    Check = simple(h(k) - trace(Qiconst(:,:,k)*X))
end

for k = 1:N
    Qconst(:,:,k)  = eval([Cstf(k)      0.5*df(:,k).';
                          0.5*df(:,k) 0.5*d2f(:,:,k)]);
                                           
    %Check = simple(f(k) - v.'*Qstage(:,:,k)*v)
    Check = simple(f(k) - trace(Qconst(:,:,k)*X))
end

Qcost = eval(0.5*[0       dCost.';
                 dCost   d2Cost]);
             
Check = simple(Cost - trace(Qcost*X))

Qunitary = zeros(size(X));Qunitary(1,1) = 1;
Check = simple(1 - trace(Qunitary*X))

%
% SDP Relaxation:
%

Nsol = 10;

%plotting stuff
sIndex = [1:2:length(v)-1];
uIndex = [2:2:length(v)-1];


Xsize = size(X,1);

aTable = logspace(-6,1,Nsol);

for sol = 1:Nsol
    a = aTable(sol);

    cvx_begin

        variable X(Xsize,Xsize);
        minimize( trace(Qcost*X) + a*norm_nuc(X));
        subject to
        for k = 1:size(Qconst,3)
            trace(Qconst(:,:,k)*X ) == 0;
        end
        for k = 1:size(Qiconst,3)
            trace(Qiconst(:,:,k)*X) <= 0;
        end
        trace(Qunitary*X) == 1;

        X == semidefinite(Xsize);
    
    cvx_end
    
    S(:,sol) = svd(X);
    display(['SVD ratio: ',num2str(S(1,sol)/S(2,sol)),' (should be >= 10)'])

    vsol(:,sol) = X(2:end,1);

    figure(1);
    subplot(2,1,1)
    plot([1:1:N],vsol(sIndex,sol));axis tight;hold on
    subplot(2,1,2)
    stairs([1:1:N-1],vsol(uIndex,sol));hold on
    line([1,N-1],[umin umin],'color','k');hold on
    line([1,N-1],[umax umax],'color','k')
    
end

figure(2)
semilogy(S(1,:)./S(2,:))



