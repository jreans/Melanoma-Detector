function [J, grad] = F15_MLcostFunctionReg(theta, X, y, lambda)
%COSTFUNCTIONREG Compute cost and gradient for logistic regression with regularization
%   J = COSTFUNCTIONREG(theta, X, y, lambda) computes the cost of using
%   theta as the parameter for regularized logistic regression and the
%   gradient of the cost w.r.t. to the parameters. 

% Initialize some useful values
m = length(y); % number of training examples
X = [ones(m, 1), X];
% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta.
%               You should set J to the cost.
%               Compute the partial derivatives and set grad to the partial
%               derivatives of the cost w.r.t. each parameter in theta
h=F24_MLsigmoid(X*theta);
theta(1)=0;
J=(1./m).*(sum(-y.*log(h)-(1-y).*log(1-h)))+lambda./(2.*m).*sum(theta.^2);
grad=(1./m).*(transpose(X)*(h-y))+lambda./m.*theta;


% =============================================================

end
