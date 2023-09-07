function out = argmin(in,k)
    if nargin == 1
       k=1; 
    end
    [~,out] = mink(in,k);
end