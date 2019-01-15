function [MECS,Us_MECS,MECS_times,B] = TargetControl(A,X_t0,X_tf,t0,tf,driver_nodes,time_step);
% X_t0: column vector of state space for initial time point of the control analysis
% A: causal/direct network (dX/dt = AX)
% X_tf: column vector of final/desir    ed state space
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

Us_MECS = zeros(N_nodes, length(t0:time_step:tf));

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
    B(driver_node,1) = 1;
    
    % Calculating Gramian
    % dt = mean(diff(Times_state));
    Wcontrol = ControlGramian(A,B,t0,tf); % version from Felix Carbonell (faster, exactly the same result).
    % Wcontrol = CtrGram(A,B); % Robust for stable and unestable systems, but not time limited.
    inv_control_gramian = inv(Wcontrol);
    
    % Calculating MECS
    U_MECS = zeros(size(B,2),length(Times_state));
    diff_output = (X_in_Times(:,end) - X_tf);
    U_MECS = -B'*multi_expv(tf - Times_state, A, inv_control_gramian*diff_output);
    
    Us_MECS(i,:) = U_MECS;
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

%  [w, err, hump] = expv( t, A, v, tol, m )
%  EXPV computes an approximation of w = exp(t*A)*v for a
%  general matrix A using Krylov subspace  projection techniques.
%  It does not compute the matrix exponential in isolation but instead,
%  it computes directly the action of the exponential operator on the 
%  operand vector. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = expv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = expv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = expv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1. However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and 
%  ||w(t)||/||v|| are of the same order of magnitude (further details in 
%  reference below).
%
%  Example 1:
%  ----------
%    n = 100;
%    A = rand(n);
%    v = eye(n,1);
%    w = expv(1,A,v);
%
%  Example 2:
%  ----------
%    % generate a random sparse matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%        end;
%    end;
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = expv(1,A,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    w_matlab = expm(full(A))*v;
%    toc
%
%    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
%    gap = norm(w-w_matlab)/norm(w_matlab);
%    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
%
%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also MEXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials. 
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function [w, err, hump] = expv( t, A, v, tol, m )

[n,n] = size(A);
if nargin == 3,
  tol = 1.0e-7;
  m = min(n,30);
end;
if nargin == 4,
  m = min(n,30);
end;

anorm = norm(A,'inf'); 
mxrej = 10;  btol  = 1.0e-7; 
gamma = 0.9; delta = 1.2; 
mb    = m; t_out   = abs(t);
nstep = 0; t_new   = 0;
t_now = 0; s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(v); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s; 
sgn = sign(t); nstep = 0;

w = v;
hump = normv;
while t_now < t_out
  nstep = nstep + 1;
  t_step = min( t_out-t_now,t_new );
  V = zeros(n,m+1); 
  H = zeros(m+2,m+2);

  V(:,1) = (1/beta)*w;
  for j = 1:m
     p = A*V(:,j);
     for i = 1:j
        H(i,j) = V(:,i)'*p;
        p = p-H(i,j)*V(:,i);
     end;
     s = norm(p); 
     if s < btol,
        k1 = 0;
        mb = j;
        t_step = t_out-t_now;
        break;
     end;
     H(j+1,j) = s;
     V(:,j+1) = (1/s)*p;
  end; 
  if k1 ~= 0, 
     H(m+2,m+1) = 1;
     avnorm = norm(A*V(:,m+1)); 
  end;
  ireject = 0;
  while ireject <= mxrej,
     mx = mb + k1;
     F = expm(sgn*t_step*H(1:mx,1:mx));
     if k1 == 0,
	err_loc = btol; 
        break;
     else
        phi1 = abs( beta*F(m+1,1) );
        phi2 = abs( beta*F(m+2,1) * avnorm );
        if phi1 > 10*phi2,
           err_loc = phi2;
           xm = 1/m;
        elseif phi1 > phi2,
           err_loc = (phi1*phi2)/(phi1-phi2);
           xm = 1/m;
        else
           err_loc = phi1;
           xm = 1/(m-1);
        end;
     end;
     if err_loc <= delta * t_step*tol,
        break;
     else
        t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
        if ireject == mxrej,
           error('The requested tolerance is too high.');
        end;
        ireject = ireject + 1;
     end;
  end;
  mx = mb + max( 0,k1-1 );
  w = V(:,1:mx)*(beta*F(1:mx,1));
  beta = norm( w );
  hump = max(hump,beta);

  t_now = t_now + t_step;
  t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
  s = 10^(floor(log10(t_new))-1); 
  t_new = ceil(t_new/s) * s;

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end;
err = s_error;
hump = hump / normv;


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