%%우주궤도역학 Term project2
%%항공우주공학과 20011321 안채원

clear all
clc

%% 0.data
load('nav.mat');
nav.GPS.a = nav.GPS.a *10^-3; %(km)
nav.QZSS.a = nav.QZSS.a*10^-3; %(km)
nav.BDS.a = nav.BDS.a*10^-3; %(km)

%% 1. v0(true anomaly at t0) (rad)
% 1.(1) E0(eccentricity anomaly at t0) (rad)
GPS.E0 = getE(nav.GPS.M0,nav.GPS.e);
QZSS.E0 = getE(nav.QZSS.M0,nav.QZSS.e);
BDS.E0 = getE(nav.BDS.M0,nav.BDS.e);

% 1.(2) v0(true anomoaly at t0) (rad)
GPS.nu0 = getnu(GPS.E0,nav.GPS.e);
QZSS.nu0 = getnu(QZSS.E0,nav.QZSS.e);
BDS.nu0 = getnu(BDS.E0,nav.BDS.e);

%% 2. v(true anomaly at t)
% 2.(1) n(mean motion) (rad/s)
mu = 3.986004418 * 10^5; %km^3/s^-2
GPS.n = sqrt(mu/nav.GPS.a^3);
QZSS.n = sqrt(mu/nav.QZSS.a^3);
BDS.n = sqrt(mu/nav.BDS.a^3);

% 임의 시간 t
GPS.t = datetime(nav.GPS.toc) + minutes(1:1440);
QZSS.t = datetime(nav.QZSS.toc) + minutes(1:1440);
BDS.t = datetime(nav.BDS.toc) + minutes(1:1440);

for i = 1 : 1440
    % 2.(2) Mean anomaly at t (rad)
    GPS.M(i) = GPS.n*(60*i) + nav.GPS.M0 ;
    QZSS.M(i) = QZSS.n*(60*i) + nav.QZSS.M0;
    BDS.M(i) = BDS.n*(60*i) + nav.BDS.M0;

    %회전수 2pi*k 고려
    GPS.M(i) = rem(GPS.M(i),2*pi);     
    QZSS.M(i) = rem(QZSS.M(i),2*pi);    
    BDS.M(i) = rem(BDS.M(i),2*pi);

    % 2.(3) Eccentricity anomaly at t (rad)
    GPS.E(i) = getE(GPS.M(i),nav.GPS.e);
    QZSS.E(i) = getE(QZSS.M(i),nav.QZSS.e);
    BDS.E(i) = getE(BDS.M(i),nav.BDS.e);

    % 2.(4) true anomaly at t (rad)
    GPS.nu(i) = getnu(GPS.E(i),nav.GPS.e);
    QZSS.nu(i) = getnu(QZSS.E(i),nav.QZSS.e);
    BDS.nu(i) = getnu(BDS.E(i),nav.BDS.e);
end

%% 3. ECI r 구하기

GPS_ECEF=[];
QZSS_ECEF = [];
BDS_ECEF = [];
for i = 1 : 1440 
    % 3.(1) perifocal r
    GPS_PQ = solveRangeInPerifocalFrame(nav.GPS.a, nav.GPS.e, GPS.nu(i));
    QZSS_PQ = solveRangeInPerifocalFrame(nav.QZSS.a, nav.QZSS.e, QZSS.nu(i));
    BDS_PQ = solveRangeInPerifocalFrame(nav.BDS.a, nav.BDS.e, BDS.nu(i));

    % 3.(2) perifocal -> ECI 
    GPS_ECI = PQW2ECI(nav.GPS.omega, nav.GPS.i, nav.GPS.OMEGA) * GPS_PQ;
    QZSS_ECI = PQW2ECI(nav.QZSS.omega, nav.QZSS.i, nav.QZSS.OMEGA) * QZSS_PQ;
    BDS_ECI = PQW2ECI(nav.BDS.omega, nav.BDS.i, nav.BDS.OMEGA) * BDS_PQ;

%% 4. ECI -> ECEF
    GPS_ECEF = [GPS_ECEF,ECI2ECEF_DCM(GPS.t(i))*GPS_ECI];
    QZSS_ECEF = [QZSS_ECEF,ECI2ECEF_DCM(QZSS.t(i))*QZSS_ECI];
    BDS_ECEF = [BDS_ECEF,ECI2ECEF_DCM(BDS.t(i))*BDS_ECI];
end
GPS.ECEF = GPS_ECEF';
QZSS.ECEF = QZSS_ECEF';
BDS.ECEF = BDS_ECEF';

%% 5. ECEF -> GEODETIC,ECI
wgs84 = wgs84Ellipsoid('kilometer');
GPS_geo = [];
QZSS_geo = [];
BDS_geo = [];
GPS_ENU = [];
QZSS_ENU = [];
BDS_ENU = [];

for i = 1 : 1440
    % 5.(1) ECEF->GEODETIC
    [lat1,lon1,h1] = ecef2geodetic(wgs84,GPS.ECEF(i,1),GPS.ECEF(i,2),GPS.ECEF(i,3));
     GPS_geo = [GPS_geo,[lat1,lon1,h1]'];

    [lat2,lon2,h2] = ecef2geodetic(wgs84,QZSS.ECEF(i,1),QZSS.ECEF(i,2),QZSS.ECEF(i,3));
    QZSS_geo = [QZSS_geo,[lat2,lon2,h2]'];

    [lat3,lon3,h3] = ecef2geodetic(wgs84,BDS.ECEF(i,1),BDS.ECEF(i,2),BDS.ECEF(i,3));
    BDS_geo = [BDS_geo, [lat3,lon3,h3]'];

    % 5.(2)  ECEF->ENU
    [east1,north1,up1] = ecef2enu(GPS.ECEF(i,1),GPS.ECEF(i,2),GPS.ECEF(i,3),37,128,0,wgs84);
    GPS_ENU = [GPS_ENU, [east1,north1,up1]'];

    [east2,north2,up2] = ecef2enu(QZSS.ECEF(i,1),QZSS.ECEF(i,2),QZSS.ECEF(i,3),37,128,0,wgs84);
    QZSS_ENU = [QZSS_ENU, [east2,north2,up2]'];

    [east3,north3,up3] = ecef2enu(BDS.ECEF(i,1),BDS.ECEF(i,2),BDS.ECEF(i,3),37,128,0,wgs84);
    BDS_ENU = [BDS_ENU, [east3,north3,up3]'];
end

GPS.GEO = GPS_geo';
QZSS.GEO = QZSS_geo';
BDS.GEO = BDS_geo';
GPS.ENU = GPS_ENU';
QZSS.ENU = QZSS_ENU';
BDS.ENU = BDS_ENU';

%% 6. Ground Track 
% geoplot((GPS.GEO(:,1)),(GPS.GEO(:,2)),'*')
% geoplot((QZSS.GEO(:,1)),(QZSS.GEO(:,2)),'r*')
% geoplot((BDS.GEO(:,1)),(BDS.GEO(:,2)),'g*')

%% 7. Skyplot
GPS.Az = azimuth(GPS.ENU);
GPS.El = elevation(GPS.ENU,10);
figure,skyplot(GPS.Az,GPS.El);

QZSS.Az = azimuth(QZSS.ENU);
QZSS.El = elevation(QZSS.ENU,10);
figure,skyplot(QZSS.Az,QZSS.El);

BDS.Az = azimuth(BDS.ENU);
BDS.El = elevation(BDS.ENU,10);
figure,skyplot(BDS.Az,BDS.El);










