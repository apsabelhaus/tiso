% Copyright Albert Li 2018

% Gets Hs and Hr if elements are ordered by [s;r]
function [Hs, Hr] = get_H(s,r)
    Hs = [eye(s) zeros(s,r)];
    Hr = [zeros(r,s) eye(r)];
end
