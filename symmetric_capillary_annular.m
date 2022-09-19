function [] = symmetric_capillary_annular()
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
% function [R, U, Psi ,ell, n, res_bvp] = symmetric_capillary_annular(a,b,psia_actual,psib_actual,kappa)
% and delete (or comment) the physical parameters block as well as the figure plotting
% at the end of the program.

%% physical parameters
% this block can be commented out if the program is used inside of a larger
% program and these values are passed as inputs to this program.
kappa = 1;
a = 0.5;
b = 3;
psia_actual = -7.7*pi/8;
psib_actual = -8*pi/8;

%tic
%% Computational parameters
% computational domain
X = [-1,1];

% maximum loop counts and the tolerances for each loop
max_iter_newton = 100000;
max_iter_bvp = 10000;
tol_newton = 1e-13;
tol_bvp = 1e-12;

% initialize the number of points.  
k = 7;
n = 2*k + 1;

%% Initial guesses

% assign to psia and psib the actual values if they are less than pi/2 in
% magnitude, and if not, then replace them with +- pi/2
psia = sign(psia_actual)*min(abs(psia_actual), pi/2);
psib = sign(psib_actual)*min(abs(psib_actual), pi/2);

% compute the Chebyshev points (called tau in the paper)
s = chebpts(n);

% compute initial guesses for psi, R and U.  These may be replaced below.
Psi0 = @(s) (1 + s)*psib/2 + (1 - s)*psia/2;
R0 = @(s) (1 + s)*b/2 + (1 - s)*a/2;
U0 = @(s) s/(b-a);

% depending on the expected geometry, update the initial guess vector v
% this is detailed in the paper.
if (psia*psib <= 0) && (psib >= 0)
    r = (b - a)/(sin(-psia) + sin(psib));
    ell0 = r*(psib - psia);
    v = [(a+b)/2 + r*sin(Psi0(s)); 1/r + r*(1-cos(Psi0(s))); Psi0(s); ell0];
elseif (psia*psib <= 0) && (psib < 0)
    r = (b - a)/(sin(psia) + sin(-psib));
    ell0 = r*(-psib + psia);
    v = [(a+b)/2 + r*sin(-Psi0(s)); -1/r - r*(1-cos(Psi0(s))); Psi0(s); ell0];
elseif (psia*psib >= 0) && (psib >= 0)
    r = (b - a)/(sin(psia) + sin(psib));
    ell0 = r*(psib + psia);
    v = [R0(s); U0(s); Psi0(s); ell0];
else
    r = (b - a)/(-sin(psia) - sin(psib));
    ell0 = abs(r*(psib + psia));
    v = [R0(s); -U0(s); Psi0(s); ell0];
end

%% plot the initial guesses
% 
% figure(3)
% plot(v(1:n),v(n+1:2*n))
% axis equal
% hold on
% % plot(v(1:n),v(2*n+1:3*n))
% figure(4)
% plot(v(1:n),Psi0(s))


%% solving the problem if it is a graph over a base domain and configuring the problem if it is not

% solve the problem if the solution is a graph over a base domain
[v, n] = cheb_engine(v, n);

% solve the problem if the solution is not a graph over a base domain
if ((abs(psib_actual) > pi/2)||(abs(psia_actual) > pi/2))
    n = (length(v) - 1)/3;
    psia_vec = linspace(sign(psia_actual)*pi/2,psia_actual,11)';
    psia_vec = psia_vec(2:end);
    psib_vec = linspace(sign(psib_actual)*pi/2,psib_actual,11)';
    psib_vec = psib_vec(2:end);
    
    for i = 1:10
        psia = psia_vec(i);
        psib = psib_vec(i);
        [v, n] = cheb_engine(v, n);
    end
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
                    dT1n1*v - a
                    dT1p1*v - b
                    dT3n1*v - psia
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
                    dT1n1*v - a
                    dT1p1*v - b
                    dT3n1*v - psia
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
                
                % warning syntax to catch badly scaled matrices.  This
                % happens when b is too small.
                [warnMsg, warnID] = lastwarn();
                if(isempty(warnID))
                else
                    warning('Radii and inclination angles lead to a muilt-scale problem.')
                    return
                    % plot the current configuration if there is a problem.
%                    R10 = v(1:n)
%                     U10 = v(n+1:2*n)
%                     Psi10 = v(2*n+1:end-1)
%                     ell10 = v(end)
%                     figure(12)
%                     plot(chebfun(R10),chebfun(U10),'.-k')
% 
%                     pause
%                     
%                     temp_psi = interp1(X,[psia;psib],chebpts(length(2*n+1:3*n)));
%                     v(2*n+1:3*n) = temp_psi;
                end     
                
                % Newton's step
                v = v + dv;
                
                % the barrier if the solution strays too far from the
                % expected solution.  This forces the inclination angle of
                % the solution to remain within reasonable bounds
                temp_psi = v(2*n+1:3*n);
                temp_psi(temp_psi > 3.5) = pi;
                temp_psi(temp_psi < -3.5) = -pi;
                v(2*n+1:3*n) = temp_psi;
                
                temp_R = v(1:n);
                temp_R(temp_R <= 0) = a/2;
                v(1:n) = temp_R;               
                                
                if v(end) <= (b - a)/2
                    v(end) = b-a;
                elseif v(end) > pi*(b + a)
                    v(end) = ell0;
                end
                
                % computing the residual of the Newton's method
                res_newton = norm(dv,'fro')/norm(v,'fro');
                
                % updating the counter and if it grows above the prescribed
                % mazimum, stop this loop and report to the outer loop.
                kk = kk+1;
                if kk > max_iter_newton
                    disp('Maximum number of Newton iterations reached')
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
                n = 2*(n - 1) - 1;
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
                break
            end
            
            
        end
    end
%% Assigning the output into variables that have a physical meaning

R = v(1:n);
U = v(n+1:2*n);
Psi = v(2*n+1:end-1);
ell = v(end);

%length(R)
% toc

%% Plotting
% delete or comment this block if the function is to be used inside some
% other algortihm.
figure(1)
plot(chebfun(R),chebfun(U),'.-k')
axis equal
hold on
plot([a;a],ylim,'-.k')
plot([b;b],ylim,'-.k')

end