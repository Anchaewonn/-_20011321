function E = getE(M,e)
Eold = M;
while(1)
    Enew = M + e*sin(Eold);
    if abs(Enew - Eold) < 0.0001
        break;
    else
        Eold = Enew;
    end
end
E = Enew;
end