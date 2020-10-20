function y = sensf(x)
 y = 2.0*((x - 5).^4 / 600 + exp(-((x - 5)/0.5).^2) + x/10. - 2);
end
