%%항공우주공학과 20011321 안채원
%%우주궤도역학 Week#13 HW -2

function az = azimuth(ENU)

ENU1 = ENU*[1 0 ; 0 1 ; 0 0];
Rrel = diag(diag(sqrt(ENU1*ENU')));
Rn = diag(ENU1*[0 1]');
az = acosd(diag(Rn/Rrel)');

end