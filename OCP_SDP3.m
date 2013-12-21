clear all
close all
clc

% run /Users/sebastien/Desktop/cvx/cvx_setup

N = 10;
Nstate = 2;
Ninput = 1;

s = sym('s',[N  ,Nstate]);
u = sym('u',[N-1,Ninput]);

x0 = [5;0];
umin = -21;
umax =  21;

%Stage constraint
f = [s(1,:).' - x0];h = [];
for k = 1:N-1
    f = [f;
         s(k+1,1) -     s(k,1)^3 - u(k)   ;
         s(k+1,2) - 0.2*s(k,2)   - u(k)^2/400]; %Dynamics
     
    h = [h;
         u(k) - umax;               %Input bounds
         umin - u(k)];
end

%Define lifted variables
v = [1];
for k = 1:N-1
    v = [v;s(k,:).';u(k)];
end
v = [v;s(N,:).'];

X = v*v.';
return
%Cost
Cost = trace(s.'*s) + 0.1*u.'*u;

% Construct matrix state constraint
for k = 1:length(f)
    df(:,k)     = jacobian(f(k),   v(2:end)).';
    d2f(:,:,k)  = jacobian(df(:,k),v(2:end));
    Cstf(k)     = simple(f(k) - df(:,k).'*v(2:end) - 0.5*v(2:end).'*d2f(:,:,k)*v(2:end));
end

for k = 1:length(h)
    dh(:,k)     = jacobian(h(k),   v(2:end)).';
    d2h(:,:,k)  = jacobian(dh(:,k),v(2:end));
    Csth(k)     = simple(h(k) - dh(:,k).'*v(2:end) - 0.5*v(2:end).'*d2h(:,:,k)*v(2:end));
end

% Construct Cost matrix
dCost   = jacobian(Cost,  v(2:end)).';
d2Cost  = jacobian(dCost, v(2:end));

for k = 1:N-1
    for i = 1:Nstate
        eval(['s',num2str(k),'_',num2str(i),' = 0;'])
    end
    eval(['u',num2str(k),' = 0;'])
end
for i = 1:Nstate
    eval(['s',num2str(N),'_',num2str(i),' = 0;'])
end

for k = 1:length(h)                    
    Qiconst(:,:,k) = eval([Csth(k)      0.5*dh(:,k).';
                          0.5*dh(:,k) 0.5*d2h(:,:,k)]);
                      
    Check = simple(h(k) - trace(Qiconst(:,:,k)*X))
end

for k = 1:length(f)
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
sIndex = [];
for i = 1:size(s,2)
    sIndex = [sIndex;
              [i:size(s,2)+1:length(v)-1]];
end
        
uIndex = [size(s,2)+1:size(s,2)+1:length(v)-1];


Xsize = size(X,1);

aTable = [logspace(-4,0,4)]

for sol = 1:length(aTable)
    a = aTable(sol);

    cvx_begin

        variable X(Xsize,Xsize);
        minimize( trace(Qcost*X) + a*norm_nuc(X)/length(X)^2);
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
    for k = 1:size(s,2)
        subplot(size(s,2)+1,1,k)
        plot([1:1:N],vsol(sIndex(k,:),sol));axis tight;hold on
    end
    subplot(size(s,2)+1,1,size(s,2)+1)
    stairs([1:1:N],[vsol(uIndex,sol);vsol(uIndex(end),sol)]);hold on
    line([1,N],[umin umin],'color','k');hold on
    line([1,N],[umax umax],'color','k')
    
end

figure(2)
loglog(aTable,S(1,:)./S(2,:),'linestyle','none','marker','.')
grid

%Check sol 
for sol = 1:size(vsol,2)
    for k = 1:N-1
        for i = 1:size(s,2)
            eval(['s',num2str(k),'_',num2str(i),' = vsol(',num2str(sIndex(i,k)),',',num2str(sol),');'])
        end
        eval(['u',num2str(k),' = vsol(',num2str(uIndex(k)),',',num2str(sol),');'])
    end
    for i = 1:size(s,2)
        eval(['s',num2str(N),'_',num2str(i),' = vsol(',num2str(sIndex(i,N)),',',num2str(sol),');'])
    end
    max(abs(eval(f)))
end
