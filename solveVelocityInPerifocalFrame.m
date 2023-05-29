%%항공우주공학과 20011321 안채원
%%우주궤도역학 week#12 HW

function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomoly)

mu = 3.986004418 * 10^14; %m^3/s^-2
p = semimajor_axis * (1-eccentricity^2);
velocityInPQW = sqrt(mu/p) * [ -sin(true_anomoly) ; eccentricity+cos(true_anomoly) ; 0];

end