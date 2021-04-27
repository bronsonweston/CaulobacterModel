function [DNAF, DNACtrA,DNACtrAP,DNA2CtrA2,DNA2CtrAP2,DNA2CtrACtrAP] = DNAbinding2(CtrA,CtrAP,kd1,kd2,kd3)
if CtrA < 0
    CtrA=0;
end
if CtrAP < 0
    CtrAP=0;
end

fun= @(DNAf) DNAf + CtrA*DNAf/kd1 + CtrAP*DNAf/kd1 + ((CtrA*DNAf/kd1)^2)/kd2 + ...
    + ((CtrAP*DNAf/kd1)^2)/kd3 + ((CtrA*CtrAP*DNAf^2/kd1^2))/kd2-1;
x0=[0 1];
DNAF=fzero(fun,x0);
DNACtrA=CtrA*DNAF/kd1;
DNACtrAP= CtrAP*DNAF/kd1;
DNA2CtrA2=((CtrA*DNAF/kd1)^2)/kd2;
DNA2CtrAP2=((CtrAP*DNAF/kd1)^2)/kd3;
DNA2CtrACtrAP= ((CtrA*CtrAP*DNAF^2/kd1^2))/kd2;
end

