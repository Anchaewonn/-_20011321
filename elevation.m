%%항공우주공학과 20011321 안채원
%%우주궤도역학 Week#13 HW -3

function el = elevation(ENU, el_mask)
Rrel = diag(sqrt(diag(ENU*ENU')));
Ru = diag(ENU*[0 0 1]');
el = asind(Ru/Rrel);
    if el <= el_mask
        el = NaN;
    end

end