function R = URandDist(X,D)
% function R = URandDist(X,D)
% draws uniformly from columns in X
% returns rounded numbers
% rejects loc if too close (<D px away) to current cursor pos

global Cursor
P = Cursor.State(1:2);

dist = 0;
while dist<D,
    R = URand(X);
    dist = norm(R(:)-P(:));
end

end % URandDist