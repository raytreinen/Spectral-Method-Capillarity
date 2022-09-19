function [] = symmetric_capillary_disc()
% Ray Treinen, April 2022
% 
% Compute capillary surfaces 
% with inclination angles psia at radius r = a and psib at radius r = b
% 
% The output can be set to [R, U, Psi, ell] to return the solution to the
% problem.  The default is set to merely plot the solution.
%
% This function needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.
%
% If this file is used as a function, replace the first line with
% function [R, U, Psi ,ell, n, res_bvp] = symmetric_capillary_disc_damped()
% and delete (or comment) the physical parameters block as well as the figure plotting
% at the end of the program.

%% physical parameters
% this block can be commented out if the program is used inside of a larger
% program and these values are passed as inputs to this program.
kappa = 1;
b = 5;
psib_actual = 7*pi/8;

%% Computational parameters
% computational domain
X = [-1,1];

% maximum loop counts and the tolerances for each loop
max_iter_newton = 100000;
max_iter_bvp = 10000;
tol_newton = 1e-13;
tol_bvp = 1e-12;

% initialize the number of points.  
% this is chosen so that the singularity in the ODE is avoided.
k = 7;
n = 2*k + 1;

%% Initial guesses

% assign to psia and psib the actual values if they are less than pi/2 in
% magnitude, and if not, then replace them with +- pi/2
psib = sign(psib_actual)*min(abs(psib_actual), pi/2);

% compute the Chebyshev points (called tau in the paper)
ss = chebpts(n);

% compute initial guesses for psi, R and U.  These may be replaced below.
Psi0 = @(s) s*psib;

if sqrt(kappa)*b >= 1
    R = abs(b/sin(psib));
    u0 = 2/R;
else
    gma = pi/2 - abs(psib);
    % u0 = 2*cos(gma)/(kappa*b) - b/cos(gma) + (2*b/3)*(1-sin(gma)^3)/cos(gma)^3;
    u0 = 2*cos(gma)/(kappa*b) - b*cos(gma)*(1 + 2*sin(gma))/(3*(1 + sin(gma))^2);
    % average of the two radii from the upper and lower estimates
    R = 0.5*(2/(u0*kappa) + abs(b/sin(psib)));
end
ell0 = abs(R*psib);
v = [sign(psib)*R*sin(Psi0(ss));sign(psib)*( R*(1 - cos(Psi0(ss))) + u0 );Psi0(ss); ell0];

%% plot the initial guesses
%
% figure(3)
% plot(v(1:n),v(n+1:2*n))
% axis equal
% hold on
% % plot(v(1:n),v(2*n+1:3*n))
% figure(4)
% plot(v(1:n),Psi0(s))


%% Solving the problem
% the cases depend on if the solution is a graph over a base domain
if (abs(psib_actual) > pi/2)
    [v, ~] = cheb_engine(v, n);
    
    n = (length(v) - 1)/3;
    psi_vec = linspace(sign(psib_actual)*pi/2,psib_actual,11)';
    psi_vec = psi_vec(2:end);
    
    for i = 1:10
        psib = psi_vec(i);
        [v, n] = cheb_engine(v, n);
    end
else
    [v, n] = cheb_engine(v, n);
end

%% the main computational engine of the file
    function [v, n] = cheb_engine(v, n)
        
        % intialize the residual
        res_bvp = 1;
        
        while res_bvp > tol_bvp
            
            % initialize the differential operator components
            %
            % D0 and D1 are spectral operators corresponding to a
            % downsampling matrix and a differentiation matrix,
            % respectively
            D0 = diffmat([n-1 n],0,X);
            D1 = diffmat([n-1 n],1,X);
            Z0 = sparse(n-1,n);
            D01 = spalloc(n-1, 3*n + 1, n*(n-1));
            D02 = D01;
            D03 = D01;
            D11 = D01;
            D12 = D01;
            D13 = D01;
            D01(1:n-1, 1:n) = D0;
            D11(1:n-1, 1:n) = D1;
            D02(1:n-1, n+1:2*n) = D0;
            D12(1:n-1, n+1:2*n) = D1;
            D03(1:n-1, 2*n+1:3*n) = D0;
            D13(1:n-1, 2*n+1:3*n) = D1;
                        
            % Evaluate the computational vector to check the boundary
            % conditions
            dT1n1 = sparse(1,3*n+1);
            dT1p1 = dT1n1;
            dT3n1 = dT1n1;
            dT3p1 = dT1n1;
            dT1n1(1) = 1;
            dT1p1(n) = 1;
            dT3n1(2*n+1) = 1;
            dT3p1(end -1) = 1;
            
            % building the nonlinear operator N 
            % and the linear operator L
            
            if sqrt(kappa)*b <= 1
                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    (D01*v).*(D13*v) + v(end).*sin(D03*v) - kappa*v(end).*(D02*v).*(D01*v)
                    dT1n1*v + b
                    dT1p1*v - b
                    dT3n1*v + psib
                    dT3p1*v - psib ];
                
                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(D13*v - kappa*v(end).*(D02*v),0,n-1,n-1)*D0, spdiags(-kappa*v(end)*(D01*v),0,n-1,n-1)*D0, spdiags(D01*v,0,n-1,n-1)*D1 + spdiags(v(end)*cos(D03*v),0,n-1,n-1)*D0, sin(D03*v) - kappa*(D02*v).*(D01*v)
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            else
                
                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    D13*v + v(end).*sin(D03*v)./(D01*v) - kappa*v(end).*D02*v
                    dT1n1*v + b
                    dT1p1*v - b
                    dT3n1*v + psib
                    dT3p1*v - psib ];
                
                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(-v(end)*sin(D03*v)./((D01*v).^2),0,n-1,n-1)*D0, -kappa*v(end)*D0, D1 + (spdiags(v(end)*cos(D03*v)./(D01*v),0,n-1,n-1))*D0, sin(D03*v)./(D01*v) - kappa*D02*v
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            end
            
            % initialize a counter and 
            % the residual for the Newton's method loop
            kk = 1;
            res_newton = 1;
                 
            %% Newton's method loop
            while res_newton > tol_newton
                
                lastwarn('', '');
                % the main steps
                dv = -L(v)\N(v);
                
                % warning syntax to catch badly scaled matrices, which
                % happened while developing the code and should no longer
                % happen
                [warnMsg, warnID] = lastwarn();
                if(isempty(warnID))
                else
                    warning('Radii and inclination angles lead to a muilt-scale problem.')
                    return
%                     % plot the current configuration if there is a problem.
%                     R10 = v(1:n)
%                     U10 = v(n+1:2*n)
%                     Psi10 = v(2*n+1:end-1)
%                     ell10 = v(end)
%                     figure(12)
%                     plot(chebfun(R10),chebfun(U10),'.-k')
%                     axis equal
%                     format longe
% 
%                     pause
                end
                
                % Newton's step
                v = v + dv;
                
                % Build a barrier if the inclination angle leaves the
                % region that has a solution

                temp_psi = v(2*n+1:3*n);
                temp_psi(temp_psi > 3.5) = pi;
                temp_psi(temp_psi < -3.5) = -pi;
                v(2*n+1:3*n) = temp_psi;
                
                
                % computing the residual of the Newton's method
                res_newton = norm(dv,'fro')/norm(v,'fro');
                
                % updating the counter and if it grows above the prescribed
                % mazimum, stop this loop and report to the outer loop.
                kk = kk+1;
                
                if kk > max_iter_newton
                    disp('Maximum number of Newton iterations reached')
                    res_newton = res_newton
                    break
                end
                
            end
            
            %% residual and other loop conditions

            % the relative norm of the residual
            res_bvp = norm(N(v),'fro')/norm(v,'fro');
            
            % adaptively add more Chebyshev points and resample the state 
            % of the problem if the tolerance is not met.  
            % Additionally, if there is too much numerical oscillation, add
            % more Chebyshev points and resample the state 
            % of the problem and reset the residual.
            S2 = diffmat([2*n n],0,X);
            if res_bvp > tol_bvp
                nn = n;
                n = n + 4;
                % Sample solutions on the new grid
                S0 = diffmat([n nn],0,X);
                S01 = spalloc(n, 3*nn + 1,n*nn);
                S01(:,1:nn) = S0;
                S02 = spalloc(n, 3*nn + 1,n*nn);
                S02(:,nn+1:2*nn) = S0;
                S03 = spalloc(n, 3*nn + 1,n*nn);
                S03(:,2*nn+1:3*nn) = S0;
                vv = zeros(3*n+1,1);
                vv(1:n) = S01*v;
                vv(n+1:2*n) = S02*v;
                vv(2*n+1:3*n) = S03*v;
                vv(end) = v(end);
                v = vv;
                
            elseif length(find(diff(sign(diff((S2*v(2*n+1:3*n))))))) >= 2
                nn = n;
                n = 2*(n - 1) -1;
                % Sample solutions on the new grid
                S0 = diffmat([n nn],0,X);
                S01 = spalloc(n, 3*nn + 1,n*nn);
                S01(:,1:nn) = S0;
                S02 = spalloc(n, 3*nn + 1,n*nn);
                S02(:,nn+1:2*nn) = S0;
                S03 = spalloc(n, 3*nn + 1,n*nn);
                S03(:,2*nn+1:3*nn) = S0;
                vv = zeros(3*n+1,1);
                vv(1:n) = S01*v;
                vv(n+1:2*n) = S02*v;
                vv(2*n+1:3*n) = S03*v;
                vv(end) = v(end);
                v = vv;
                
                res_bvp = 1;
            else
                break
            end
            
            % if the function exceeds the maximum number of iterations,
            % break with an error statement.
            if n > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                res_bvp = res_bvp
                break
            end
        end
    end
%% Assigning the output into variables that have a physical meaning

R = v(1:n);
U = v(n+1:2*n);
Psi = v(2*n+1:end-1);
ell = v(end);

%% Plotting
% delete or comment this block if the function is to be used inside some
% other algortihm.
figure(1)
plot(chebfun(R),chebfun(U),'.-k')
hold on
axis equal
plot([-b;-b],ylim,'-.k')
plot([b;b],ylim,'-.k')

end