function S = MCM_Rec(sc, vc, fc, A, S0, Z, w, N_steps)
% Multifactorial causal model incorporating neurotransmitter receptors
% Forward simulation
%   sc: structural connectivity [src, tar]
%   vc: vascular connectivity [src, tar]
%   fc: functional connectivity [src, tar]
%   w: factorwise connectivity weighting [fac]
%   A: 3D tensor of factor-receptor interactions [j,m,n] [reg, fac, params]
%   S0: initial factors state [regions, factors]
%       [CBF, A-beta, f-conn, glucose metabolism, GM density, tau]
%   Z0: initial receptors state [regions, receptors]
%   N_steps: number of simulation steps

% TODO: add vascular connectivity
% TODO: intraregional factor self-interaction
N_regs = size(A, 1);
N_facs = size(A, 2);

% Initial states 
S = S0;
dS = zeros(size(S0));

% dS/dt = [sum_(n) alpha_i^n->m S_i^n] + spreading + inputs 
% dS/dt = AS + Bu
% a_i^(n->m)  = f(GE_i, NT-R_i)
% = beta_0^n->m + sum_r beta_r^(n->m) NT-R_(i,r)
% spreading = D^m ( sum_j C_(j->i)S_j^m   - sum_j C_(i->j) S_i^m )
% = D^m *sum_j C_j,i (S_j^m - S_i^m)

% Simplify
% m = n, i=/=j [same factor, different regions]
% A^(m,n)_(i,j) = beta_0^(m,n) + sum_j^N_RE beta_j^(m,n) RE_(j,i) - sum_k^N_roi C_i->k S^m

% m=/=n, i=/=j [different factors, regions]
% A^(m,n)_(i,j) = 0

% m=/n, i=j [different factors, same region]
% A^(m,n)_(i,j) =  b0^(m,n) + sum_j^(N_RE) beta_j^(m,n) RE_(j,i)

% Homogeneous population regression across regions for same factor to get
% fn->m (RE)

% Timestep
h = 0.01;

for i=1:N_steps
   for reg=1:N_regs
      for fac=1:N_facs 
          other_facs = squeeze(S(reg, setdiff(1:N_facs, fac)));
          other_regs = squeeze(S(setdiff(1:N_regs, reg), fac))';
          state = [squeeze(Z(reg,:)) other_facs other_regs];
          dS(reg, fac) = dot(squeeze(A(reg, fac, :)), state);
          
      end
   end    
   
   S = S + (dS * h);
end

% for i=1:N_steps    
%     for reg=1:N_regs
%         % Augment interaction matrix with regional factor interactions
%         alphas = alphas_XZ;
%         alphas(1:N_facs, 1:N_facs) = squeeze(alphas_XX(reg, :, :));
%         % Adjust for double counting local + global interactions
%         alphas(:,1:N_facs) = alphas(:,1:N_facs) ./ 2;
%         
%         % [Intraregional effects of all other factors] 
%         dS = alphas * S(:, reg);
%         
%         for fac=1:(N_facs + N_recs)
%             % [Interregional with the same factor] 
%             if fac <= N_facs
%                 in_rate = (w(fac) * vc(:, reg)) + ((1-w(fac)) * sc(:, reg));
%                 out_rate = (w(fac) * vc(reg, :)) + ((1-w(fac)) * sc(reg, :));
%             else
%                 in_rate = sc(:, reg);
%                 out_rate = sc(reg, :);
%             end
%             
%             inflow = in_rate' * S(fac, :)' * S(fac, reg);
%             outflow = out_rate * S(fac,:) * S(fac, reg);
%             dS(1, fac) = dS(1, fac) + inflow - outflow;
%         end
%                 
%         S = S + (dS + h);
%     end
% end

end

