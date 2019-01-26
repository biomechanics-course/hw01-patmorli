% Unit test for finite difference fjacobian function
% using Matlab's unit testing framework.
% To execute tests, execute runtests('testfjacobian')

% define anonymous functions use for jacobian testing
samplefunction0 = @(x) x(1)^2+x(2)*x(3);           % scalar function of 3x1 x
samplefunction1 = @(x) [x(1)^3; x(1)^2+x(2)*x(3)]; % 2x1 function, 3x1 x
samplefunction2 = @(x) [tan(x(1))*x(2); sin(x(1))+x(2); cos(x(2))]; % 3x1 function, 2x1 x
% (anonymous functions allow simple functions to be defined in-line)

% define the analytical jacobians to compare against
analyticaljacobian0 = @(x) [2*x(1) x(3) x(2)];
analyticaljacobian1 = @(x) [3*x(1)^2 0 0; 2*x(1) x(3) x(2)];
analyticaljacobian2 = @(x) [x(2)*sec(x(1))^2 tan(x(1)); cos(x(1)) 1; 0 -sin(x(2))];

% define points to evaluate Jacobian about
x0 = [1 2 3];
x1 = [1 2 3];
x2 = [1 2];

abstol = 1e-5; % a rough absolute tolerance for finite differencing
% All of the variables above are shared between the test cells below.

%% Test 0: Gradient of function 0
% assuming the finite differencing step size is small, the error in
% finite differencing should not be very large.
jacobian0 = fgradient(samplefunction0, x0); % take derivative about x0
assert(isequal(size(jacobian0), [1 3]), 'Jacobian 0 wrong size')
assert(max(max(abs(jacobian0 - analyticaljacobian0(x0)))) < abstol, ...
  'Numerical Jacobian 0 not very accurate')

%% Test 1: Jacobian of function 1
% assuming the finite differencing step size is small, the error in
% finite differencing should not be very large.
jacobian1 = fjacobian(samplefunction1, x1); % take derivative about x1
assert(isequal(size(jacobian1), [2 3]), 'Jacobian 1 wrong size')
assert(max(max(abs(jacobian1 - analyticaljacobian1(x1)))) < abstol, ...
  'Numerical Jacobian 1 not very accurate')

%% Test 2: Jacobian of function 2
jacobian2 = fjacobian(samplefunction2, x2); % take derivative about x1
assert(isequal(size(jacobian2), [3 2]), 'Jacobian 2 wrong size')
assert(max(max(abs(jacobian2 - analyticaljacobian2(x2)))) < abstol, ...
  'Numerical Jacobian 2 not very accurate')

%% Test 3: Change finite difference step size
smalltol = abstol * 1e-2; % try a tighter tolerance
jacobian3 = fjacobian(samplefunction2, x2, 1e-8); % use smaller step
assert(max(max(abs(jacobian3 - analyticaljacobian2(x2)))) < smalltol, ...
  'Numerical Jacobian 1 not very accurate')

%% Test 4: Demonstrate error 
analyticaljacobian = analyticaljacobian2(x2);
exponents = 0:12;
for i = 1:length(exponents)
    jacobian = fjacobian(samplefunction2, x2, 10^(-exponents(i)));
    generror(i) = sqrt(sum(sum((jacobian - analyticaljacobian).^2)) / ...
        prod(size(jacobian)));
    
end
semilogy(exponents, generror)
title('Error analysis of finite difference Jacobian')
xlabel('Finite difference step, 10^{-i}'); ylabel('Root-mean-square error');