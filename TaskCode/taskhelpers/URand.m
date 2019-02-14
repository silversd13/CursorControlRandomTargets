function R = URand(X)
% function R = URand(X)
% draws uniformly from columns in X
% returns rounded numbers

N = size(X,2);
R = zeros(1,N);

% gen uniform rand numbers
r = rand(1,N);

% move distr to loc specified in X columns & round
for i=1:N,
    R(i) = round(X(1,i) + diff(X(:,i))*r(i));
end

end % URand