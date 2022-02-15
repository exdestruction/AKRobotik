function [vm, bm, tb, tv, te] = calc_t_ramp_te_ta(qs, qz, vm_ref, bm_ref)
se = abs(qz-qs);
mask = se > 0;
tb = vm_ref./ bm_ref;
te = se./vm_ref + tb;
tv = te - tb;
[te(:), idx] = max(te);
tb(:) = tb(idx);
tv(:) = tv(idx);
vm = se ./ tv;
bm = vm ./ tb;
end