; Predict the time it will take one healpix to combine
function predictcombtime,glon,glat,nexp
par = [24.847320,  62.468829,  14.928247,  4.2707093,  46.408931,  2.7092888, 51.570365, 2.8514323]
lmcdist = sphdist(280.52,-32.53,glon,glat,/deg)
smcdist = sphdist(302.80,-44.30,glon,glat,/deg)
dt = par[0] * exp(-0.5*(glon^2)/par[1]^2) * exp(-0.5*(glat^2)/par[2]^2) + par[3] + $
       par[4]*exp(-0.5*lmcdist^2/par[5]^2)^0.5 + par[6]*exp(-0.5*smcdist^2/par[7]^2)
dt *= nexp
return,dt
end
