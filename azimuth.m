%%항공우주공학과 20011321 안채원
%%우주궤도역학 Week#13 HW -2
%% input = ENU (nx3) 
function az = azimuth(ENU)
    n = length(ENU);
    az = [];
    for i = 1 : n
        Re = ENU(i,1);
        Rn = ENU(i,2);
        R = sqrt(ENU(i,1)^2+ENU(i,2)^2);
        az_temp = acosd(Rn/R); %deg    

        if Re < 0
            az = [az,360-az_temp];
        else
            az = [az,az_temp];
    end
end