function Xs = NewtonRoot(Fun,Funder,Xest,Err,imax)
% NewtonRoot finds the root of Fun = 0 near the point Xest using Newton's method.
% Input variables:
% Fun       Name (string) of a function file that calculates Fun for a given x.
% FunDer    Name (string) of a function file that calculates the derivative of Fun for given x
% Xest      Initial estimate of the solution
% Err       Maxiumum Error.
% imax      Maxiumum number of iterations
% Output variable:
% Xs        Solution

for i=1:imax
    Xi = Xest - feval(Fun,Xest)/feval(FunDer,Xest);
    if abs((Xi-Xest)/Xest) < Err
        Xs = Xi;
        break
    end
    Xest = Xi;
end
if i==imax
    fprintf('Solution was not obtained in %i iterations.\n',imax)
    Xs = ('No answer')
end
