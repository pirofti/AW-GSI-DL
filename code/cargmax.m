function out = argmax(in,k)
    if nargin == 1
       k=1; 
    end
    [~,out] = maxk(in,k);
end