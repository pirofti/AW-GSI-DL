function [l,x_nom_orig,x_leak_orig]=load_leak(leaky_node,Leaktionary,N)
l = leaky_node;
x_nom_orig = Leaktionary{N}.Head;
x_leak_orig = Leaktionary{leaky_node}.Head;
end