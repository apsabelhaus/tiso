% Andrew P. Sabelhaus 2019

function H = get_H(s, r)
%% get_H 
%   get_H returns the H matrix used to separate cables from bars at various
%   points in the inverse statics optimization process.
%
%   See paper for justification. Applicable to both 2D and 3D cases.

H = [eye(s);
     zeros(r, s)];

end

