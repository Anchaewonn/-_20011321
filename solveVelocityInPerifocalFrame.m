%%항공우주공학과 20011321 안채원
%%우주궤도역학 week#12 HW

function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomoly)

mu = 3.986004418 * 10^5; %km^3/s^-2
nu = deg2rad(true_anomoly);
p = semimajor_axis * (1-eccentricity^2);
velocityInPQW = sqrt(mu/p) * [ -sin(nu) ; eccentricity+cos(nu) ; 0];

end