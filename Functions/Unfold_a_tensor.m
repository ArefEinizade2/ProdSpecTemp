function [X] = Unfold_a_tensor( T, dim, i )

%% A function for mode-i unfolding a tensor

X = reshape(shiftdim(T,i-1), dim(i), []);

end