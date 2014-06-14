function rvar = RelativeVariance( gammabar )

cg = ceil(1./gammabar);
vt1 = -gammabar + 2*gammabar.^2.*cg - gammabar.^3.*cg.^2;
vt2 = 2*(gammabar.*cg - gammabar.^2.*cg.^2);
vt3 = 2/3*gammabar.^3.*cg.^3 + 1/3*gammabar.^3.*cg;

rvar = vt1 + vt2 + vt3;

end
