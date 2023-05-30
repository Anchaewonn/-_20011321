%%항공우주공학과 20011321 안채원
%%우주궤도역학 week#12 HW

function rangeInPQW = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomoly)

nu = deg2rad(true_anomoly);
r = semimajor_axis*(1-eccentricity^2) / (1+eccentricity*cos(nu));
rangeInPQW = [r*cos(nu) ; r*sin(nu) ; 0];

end
