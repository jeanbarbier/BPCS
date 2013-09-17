function [dump] = dumping(old,new,m)
% compute the dumping of b and c by a factor m

dump = m .* old + (1 - m) .* new;

end

