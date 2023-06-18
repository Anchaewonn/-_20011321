%%항공우주공학과 20011321 안채원
%%우주궤도역학 Week#13 HW -(1)

function DCM = ECI2ECEF_DCM(time)

% time ([YYYY,MM,DD,hh,mm,ss] format)
jd = juliandate(time);
thetag = siderealTime(jd);
DCM = [ cosd(thetag) sind(thetag) 0 ; 
        -sind(thetag) cosd(thetag) 0 ;
        0 0 1];
end
