function [MECS,U_MECS,MECS_times,B] = TargetControl(A,X_t0,X_tf,t0,tf,driver_nodes,time_step);
% X_t0: column vector of state space for initial time point of the control analysis
% A: causal/direct network (dX/dt = AX)
% X_tf: column vector of final/desired state space
% t0: initial time point for control analysis
% tf: final time point for control analysis
% driver_nodes (optional): list of nodes that will be used as drivers
% time_step (optional): integration step
% -------------------------------------------------------------------------%
% MECS: final Minimum Energy Control Strategy, from Klickstein et al.,
% 2016.
% --------------------------------------------------------------------------%
% Yasser Iturria Medina, 03/05/2016.

N_nodes = size(A,1);
if nargin < 6
    driver_nodes = 1:N_nodes;
end
if nargin < 7
    time_step = (tf-t0)/1000;
end
for i = 1:length(driver_nodes)
    disp(['Driver -> ' num2str(driver_nodes(i))]);
    driver_node = driver_nodes(i);
    
    if i == 1,
        Times_state = t0:time_step:tf; % interval of the controlling process, where A
        % is assumed to be constant (e.g. for a given cognitive/clinical state)
        X_in_Times = multi_expv(Times_state-t0, A, X_t0); % Calculating state space trajectories in the time interval of interest.
        B0 = zeros(N_nodes,1);
    end
    
    B  = B0;
    % Input signal for driver node
    B(driver_node,i) = 1;
    
    % Calculating Gramian
    % dt = mean(diff(Times_state));
    Wcontrol = ControlGramian(A,B,t0,tf); % version from Felix Carbonell (faster, exactly the same result).
    % Wcontrol = CtrGram(A,B); % Robust for stable and unestable systems, but not time limited.
    inv_control_gramian = inv(Wcontrol);
    
    % Calculating MECS
    U_MECS = zeros(size(B,2),length(Times_state));
    diff_output = (X_in_Times(:,end) - X_tf);
    U_MECS = -B'*multi_expv(tf - Times_state, A, inv_control_gramian*diff_output);
    MECS_times = sum(U_MECS.^2)*time_step; % control energy associated with the optimal control strategy.
    
    % Total Energy:
    MECS(i) = sum(MECS_times);
end
return;

function W = multi_expv(t, A, v)

[t_sort, ind_t] = sort(t);
nt = length(t);
W = zeros(length(v), nt);
t = [0, t_sort];
vi = v;
for i=2:nt+1
    w = expv(t(i)-t(i-1), A, vi);
    W(:,i-1) = w;
    vi = w;
end
W = W(:,ind_t(ind_t));
return;

function W = ControlGramian(A,B,t0,tf)

h = tf-t0;
n = size(A,1);
C = [-A,B*B';zeros(n),A'];
D = expm(C*h);
W = D(n+1:end,n+1:end)'*D(1:n,n+1:end);
return;

function [P,Tr] = CtrGram(A,B,t0,tf)

%Returns controllability gramian P for an unstable system and the transformation Tr that splits A
%into stable and unstable parts.
%This program uses the results from the following two publications:

%[1] K. Zhou, G. Salomon and E. Wu, 'Balanced realization and model
%reduction for unstable systems", International Journal of Robust and
%Nonlinear Control, vol. 9, pp. 183 - 198, 1999.

%[2] S.K. Nagar and S.K. Singh, 'An algorithmic approach for system
%decomposition and balanced realized model reduction', Journal of the
%Franklin Institute vol 341, pp. 615630, 2004.

%Last updated: 24/09/2012
%Email: chrislbowden@hotmail.com

%Check dimensions
[n n0]   = size(A);
if (n ~= n0)
    error('CtrGram:Dim','The state matrix A must be a square matrix')
end
[n1 p1] = size(B);
if (n1 ~= n)
    error('CtrGram:Dim','A and B must have the same number of rows')
end

%Check the number of positive, zero and negative eigenvalues of A
Evals = eig(A);
Pos = 0;
Neg = 0;
Zer = 0;
for i = 1:n
    if real(Evals(i)) > 0
        Pos = Pos +1;
    end
    if real(Evals(i)) < 0
        Neg = Neg +1;
    end
    if real(Evals(i)) == 0
        Zer = Zer +1;
    end
end

%This method is only valid for systems without poles on the imaginary axis
if (Zer ~= 0)
    error('CtrGram:ImagAxisPole','The state matrix has imaginary axis poles. ')
end
%If A is stable then the function 'gram' from the control system toolbox
%can be used.
if Neg == n
    warning('CtrGram:StabMat','The state matrix is stable. Gramian could also be computed using control system toolbox.')
end

%The transformation T in Theorem 1 in Zhou (1999) can be obtained using the
%algorithm of Nagar and Singh (2004)
%Find the Schur form:
[U,T] = schur(A);

%Reorder the Schur factorization A = U*T*U' of  matrix A so that the
%negative (i.e. stable) cluster of eigenvalues appears in the
%leading (upper left) diagonal blocks of the Schur matrix T
%(Nagar and Singh eqn (2))
[US,TS] = ordschur(U,T,'lhp');

%Solve lyapunov equation - Nagar and Singh eqn(3)]
%This will enable us to decouple the system in order to separate into
%stable and unstable components
A11     =   TS(1:Neg,1:Neg);
A12     =   TS(1:Neg,Neg+1:n);
A22     =   TS(Neg+1:n,Neg+1:n);
S       =   lyap(A11,-A22,A12);

%Use W and it's inverse to find decoupled system - Nagar and Singh eqn(4,5)
W       =   [eye(Neg) S; zeros(n-Neg,Neg) eye(n-Neg)];
Winv    =   [eye(Neg) -S; zeros(n-Neg,Neg) eye(n-Neg)];

%The above can be used to obtain the transformation T (labelled T1 here) in
%Zhou 1999, Theorem 1.
T1      = Winv*US';
T1inv   = US*W;

%To get parts of B corresponding to transformed stable and unstable system:
BB=T1*B;
B1=BB(1:Neg,1:p1);
B2=BB(Neg+1:n,1:p1);

%Gd=[A11 0; 0 A22] where A11 is stable and A22 is unstable
Gd=T1*A*T1inv;

%As and B1 form stable part, Au and Bu form unstable part
As = Gd(1:Neg,1:Neg);
Au = Gd(Neg+1:n,Neg+1:n);

%P1 and P2 are controllability gramians of (As,Bs) and (-Au,Bu)
%- Zhou (1999) p187
P1 = lyapchol(As,B1)'*lyapchol(As,B1);
P2 = lyapchol(-Au,B2)'*lyapchol(-Au,B2);

% P1 = ControlGramian(As,B1,t0,tf);
% P2 = ControlGramian(-Au,B2,t0,tf);

%P is controllability Gramian of the system, the 'larger' P the
%smaller the required control energy
%- Zhou (1999) p187 and p195

P = T1inv*[P1 zeros(Neg,n-Neg); zeros(n-Neg,Neg) P2]*T1inv';
Tr=T1;
return;