function test_scattering_potential
FigureInversion('OUTPUT/24-Aug-2017/LOOP-5000-400000-120000-1-false',120000,400000,5000,'tmp');

model=velocity_model();

function calc_dlnv(vref,vnew)
    dlnv=(vnew-vref)/vref;
    fprintf('dlnv = %f\n', dlnv);
end

calc_dlnv(model.vs(2),model.vs(1))
calc_dlnv(model.vs(3),model.vs(2))

calc_dlnv(3300*model.vs(2),2800*model.vs(1));

end

