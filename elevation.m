%%항공우주공학과 20011321 안채원
%%우주궤도역학 Week#13 HW -3
%% input = ENU(nx3,km), el_mask(deg) / output = el(rad)
function el = elevation(ENU, el_mask)
   n = length(ENU);
   el = [];
   for i = 1 : n
       Ru = ENU(i,3);
       Rrel = sqrt(ENU(i,1)^2+ENU(i,2)^2+ENU(i,3)^2);
       el_temp = asind(Ru/Rrel);
       if el_temp < el_mask
           el = [el,NaN];
       else
           el = [el,el_temp];
       end
   end   
end