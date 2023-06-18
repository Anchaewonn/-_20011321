function v0 = getnu(E,e)
a = sqrt(1-e^2)*sin(E) / (1-e*cos(E));
b = (cos(E)-e) / (1-e*cos(E));
v0 = atan2(a,b);
end
