function [sp,T] = F28_MLgetspwse98(Xalg,Yalg,T)
placese = find(Yalg >= .98);
if isempty(placese)
    placese = length(Yalg);
else
    placese = placese(1);
end
sp = 1-Xalg(placese);
T = T(placese);
end