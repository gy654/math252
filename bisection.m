f = @(x) x.^3;
  
tol = 1e-5;

tic
[root, nits]= bisection_method(f,a,b,tol)
toc

publish('bisection.m','pdf')

function [root, nits] = bisection_method(f,a,b,tol)
%%%%%%%%
% Inputs: 
% - f: is a function handle for a function we want to find a root of
% - a,b: are the left and right endpoints of the search domain (respectively)
% - tol: is an estimate of how close we would like to be to our root
%
% Outputs:
% -root: is the root of the function
% -nits: the numer of itteration the program took 
%%%%%%%%

if f(a)*f(b) > 0 % check if there is a root in the interval
    %if not, output garbage
    root = NaN; 
    nits = 0;
    disp('you did bad and you should feel bad')
else
    nits = 0; % iteration counter
    a_k = a; % sequence of left endpoints
    b_k = b; % sequence of right endpoints
    m_k = (b_k+a_k)/2; % midpoint of [a_k,b_k]
    while abs(f(m_k)) > tol % loop until midpoint is at MOST 'tol' from the root
        nits = nits + 1; % increment iteration counter
        m_k = (b_k+a_k)/2; % midpoint of [a_k,b_k]
        if sign(f(m_k)) == sign(f(a_k))
            % if f(m_k) has the same sign as 
            % f(a_k), then we can safely replace
            % a_k by m_k and maintain that
            % f(a_k)*f(b_k) < 0 
            a_k = m_k;
        else
            % otherwise replace b_k with m_k
            b_k = m_k;
        end
    end
    root = m_k; %set the output to m_k   
end

end
