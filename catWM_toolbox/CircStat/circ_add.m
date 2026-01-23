function y = circ_add(x,theta)
%y = circ_add(x,theta)
    y = mod(x+theta+pi,2*pi)-pi;
end