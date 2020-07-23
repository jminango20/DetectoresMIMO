function [valor , pos] = second_min( x )
valor = min(x(x>min(x)));
pos = find(x==valor);