%function berechnet aus vm, bm und se die Zeiten tb, tv und te
function [tb, tv, te, vm]=calc_t_ramp(vm, bm, se)
if se>=0 
    vm_max = sqrt(bm .* se);
    mask = vm > vm_max;
    %Prüfung ob vm höher ist als vm,max
    vm(mask) = vm_max(mask);

    %Berechne Zeiten
    tb= vm ./ bm;
    te= se ./ vm + tb;
    tv= te - tb;

    
else
    disp('se muss positiv sein!')%Fehlermeldung 
    tb=0;
    te=0;
    tv=0; %Nullsetzen der Ausgabewerte im Fehlerfall
end

 