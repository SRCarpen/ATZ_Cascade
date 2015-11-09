# Allochthony model for labeled lakes: Paul Lake 2001
# Stephen R. Carpenter, 2015-09-19

rm(list = ls())
graphics.off()
library(numDeriv)

# Functions for phytoplankton ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Integrated irradiance effect for phytoplankton growth
GAMMA.fun <- function(z,Pbar,dz) {
  eps.total = eps0 + epsDOC*DOC + epsP*Pbar
  Iz <- I0*exp(-z*eps.total)
  rate.at.z <- dz*(1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
  GAMMA = sum(rate.at.z)
  return(GAMMA)
}

# Phytoplankton instantaneous growth rate (losses not included)
Grow.Phyto = function(P0,DOC,Load,Zvec,dz)  {
  Igamma = GAMMA.fun(Zvec,P0,dz)
  Prate = rP*Igamma*P0*Load/(2 + Load) 
  return(Prate)
}

# End of phytoplankton functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Main Program &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# Set up for phytoplankton calculations

I0 <- 600  # Surface irradiance, microEinsteins m-1 s-1; 600 is typical daily median
# P-I curve parameters, median of NTL P-I curves
k_sat <- 0.0194 # per microEinsteins m-1 s-1 (Follows 0.012)
k_inh <- 0.00065 # per microEinsteins m-1 s-1 (Follows mean 0.004, range 0.001-0.007)
# Derived parameter from Follows et al.
Fmax <- ((k_sat + k_inh)/k_sat)*exp(-(k_inh/k_sat)*log(k_inh/(k_sat+k_inh)))
# Light extinction parameters, Carpenter et al. L&O 1998
eps0 <- 0.0213  # Light extinction by pure water
epsDOC = 0.0514  # DOC light extinction coef
epsP <- 0.0177   # Phytoplankton, per m2 (mg phosphorus)-1
rP = 1  # Phytoplankton production per unit P load

# Data for individual whole-lake experiment ++++++++++++++++++++++++++++++++++++++++++++
# Paul Lake 2001  
ZT = 3.5
DOC = 304*12/1000  # umol/L * 12ug/umol * 10^-3 mg/ug
POC = 35.5*12/1000  # umol/L * 12ug/umol * 10^-3 mg/ug
Chl = 4.21
Load = 0.3  # mg P m-2 d-1
ZB = 1.05*0.5 # Estimate converted to g C m-2 from Dry Mass 
Phi.POC = 0.38  # POC allochthony
TPOC = Phi.POC*POC
Phi.Z = 0.36  # Zoop allochthony # Chaob 0.36; Zoop 0.24
GPP = 43.4  # GPP in mmol O2 m-2 d-1 

# End experiment-specific data ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# C:Chl ratio based on autochthonous POC and Chl
APOC = (1-Phi.POC)*POC*1000 # convert mg to ug
CChl = APOC/Chl  # mass:mass

# Areal Phyto Biomass as C 
AC.AR = Chl*ZT*CChl/1000 # Algal g C m-2 using C:Chl ratio

# Algae
I.AC = GPP*12*0.7/1000 # NPP as g C m-2 d-1
Mtot.AC = I.AC/AC.AR  # total mortality of algae due to all losses
s.A = 0.3/ZT  # Phytoplankton sinking (velocity / thermocline depth)
QAZ = (Mtot.AC - s.A)*AC.AR  # grazing flux as g C m-2

print('Results for Paul 2001',quote=F)
print('Phytoplankton',quote=F)
print(c('Phyto C to Chl ratios direct method = ',CChl),quote=F)
print(c('Phyto biomass g C m-2',AC.AR),quote=F)
print(c('NPP g C m-2 d-1',I.AC),quote=F)
print(c('Phyto total mort = ',Mtot.AC),quote=F)
print(c('Sinking component of total Mort',s.A),quote=F)
print(c('Flux to zooplankton ',QAZ),quote=F)

# Scale the Follows et al. production function to observed NPP
# Depth sequence for integration of GPP
nZ <- 20  # Steps for vertical integration
dz <- ZT/nZ
Zvec <- seq(0,ZT,by=dz)
# Compute production for reference conditions
NPP_ref = Grow.Phyto(Chl,DOC,Load,Zvec,dz)
print('',quote=F)
print('Rescaling Follows et al. production function to observed NPP',quote=F)
print('Reference Conditions', quote=F)
print(c('NPP = ',NPP_ref),quote=F)
print('Chl, DOC, Zthermo, Load, rP',quote=F)
print(c(Chl,DOC,ZT,Load,rP),quote=F)
# Rescale rP so that NPP is observed value at reference conditions
rP = I.AC/NPP_ref
NPP_ref = Grow.Phyto(Chl,DOC,Load,Zvec,dz)
print(c('Reference Conditions with rP rescaled to ',rP), quote=F)
print(c('NPP = ',NPP_ref),quote=F)
print('Chl, DOC, ZT, Load, rP',quote=F)
print(c(Chl,DOC,ZT,Load,rP),quote=F)

# Compute attack rate 
# See Chow-Fraser+Sprules_FT_fit.R
handle = 0.005 # handling time in days*animals/algae consumed from Chow-Fraser & Sprules
attack = QAZ/(AC.AR*ZB - handle*AC.AR*QAZ)

print('',quote=F)
print('Grazing',quote=F)
print(c('Attack rate = ',attack),quote=F)
print(c('Handling Time = ',handle),quote=F)

# Compute steady-state algae detritus
pA = 0.3  # egestion coefficient BACK TO D for zoops feeding on algae
pD = 0.3
s.D = 0.5/ZT # Sedimentation loss coefficient = sinking rate/ZT (Reynolds 1984)
Dcoef = c(0,0,0) # vector to hold polynomial coefficients for detritus polynomial
Dcoef[1] = pA*QAZ
Dcoef[2] = pA*QAZ*attack*handle - s.D - (1-pD)*attack*ZB
Dcoef[3] = -1*s.D*attack*handle
Droots = polyroot(Dcoef)
Dstar = max(Re(Droots))

# Flux of algae detritus to zooplankton
QDZ = attack*Dstar*ZB/(1 + attack*handle*Dstar)

print('Detrital algae info',quote=F)
print(c('Detrital algae steady state g C m-2 ',Dstar),quote=F)
print(c('Detrital algae flux to zoopl g C m-2 d-1',QDZ),quote=F)

# Compute TPOC input rate
pT = 0.5 # egestion coefficient for TPOC back to TPOC
TPOCAR = TPOC*ZT # areal TPOC g/m2
QTZ = attack*TPOCAR*ZB/(1 + handle*attack*TPOCAR) # Flux from TPOC to Zoopl
s.T = 0.1/ZT # Sedimentation loss coefficient = sinking rate/ZT (Reynolds 1984)
I.T = s.T*TPOCAR + (1-pT)*QTZ

print('Phyto and TPOC fluxes to Zoopl',quote=F)
print(c(QAZ,QTZ))
print('TPOC fluxes',quote=F)
print(c('TPOC biomass g C m-2',TPOCAR),quote=F)
print(c('TPOC sedimentation loss coefficient = ',s.T),quote=F)
print(c('TPOC input rate g m-2 d-1 = ',I.T),quote=F)

# Compute growth efficiencies on algae and TPOC for zooplankton
gAZ = 0.25 # assumed growth efficiency from algae to zoop 
gDZ = 0.05 # assumed growth efficiency from algal detritus to zoop
gTZ = Phi.Z*(gAZ*QAZ + gDZ*QDZ)/( (1-Phi.Z)*QTZ )

print('Zooplankton',quote=F)
print(c('Efficiencies gAZ, gTZ ',gAZ,gTZ),quote=F)

# Compute Zoop mortality
mZTOT = gAZ*QAZ + gDZ*QDZ + gTZ*QTZ
mZ = mZTOT/ZB

print(c('Total Zoop Mort flux = ',mZTOT),quote=F)
print(c('Zoop Mort Coef = ',mZ),quote=F)
print(c('Zoop biomass, g m-2',ZB),quote=F)

# Parameters for zooplanktivory
mZnp = 0.04 # non-predatory mortality coefficient of zooplankton
QZF = (mZ-mZnp)*ZB # planktivory flux of zooplankton
hF = 1.4*ZT # Based on estimate in the Regime Shift book
cF = QZF*(hF^2 + ZB^2)/(ZB^2) # maximum planktivory rate

# Parameters for zooplankton refuging
D.Z = 0 # Diffusion rate between refuge and foraging arena
Zref = ZB  # Zoop biomass in refuge

print('Zooplankton parameters',quote=F)
print('Planktivory flux, hF, cF',quote=F)
print(c(QZF,hF,cF))

# Analysis of Equilibria $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Function for deviation of A, T and Z from equilibrium
dATZdt.eq = function(lY0) {  
  # unpack state variables (all as g c m-2)
  Y0 = exp(lY0)
  A0 = Y0[1]
  T0 = Y0[2]
  Z0 = Y0[3]
  D0 = Y0[4]
  # Light effects
  vChl = A0*(1/CChl)*1000*(1/ZT) # convert g C/m2 to mg Chl/m3
  NPP = Grow.Phyto(vChl,DOC,Load,Zvec,dz) # NPP
  Anet = NPP - s.A*A0  # Net after sinking
  # Consumption fluxes
  Q.AZ = attack*A0*Z0/(1 - attack*handle*A0)
  Q.TZ = attack*T0*Z0/(1 - attack*handle*T0)
  Q.DZ = attack*D0*Z0/(1 - attack*handle*D0)
  Q.ZF = (cF*Z0^2)/(hF^2 + Z0^2)
  # Dynamics
  dAdt = Anet - Q.AZ
  dTdt = I.T - s.T*T0 -(1-pT)*Q.TZ
  dZdt = gAZ*Q.AZ + gTZ*Q.TZ + gAZ*Q.DZ - mZnp*Z0 - Q.ZF + D.Z*(Zref-Z0)
  dDdt = pA*Q.AZ - s.D*D0 - (1-pD)*Q.DZ
  rates = c(dAdt,dTdt,dZdt,dDdt)
  SSE = sum(rates*rates)  # sum of squared distance from equilibrium
  return(SSE)  # return either rates or SSE 
}

# Function for Jacobian of A, T and Z
dATZdt.jac = function(Y0) {  
  # unpack state variables (all as g c m-2)
  #Y0 = exp(lY0) # no need for transform
  A0 = Y0[1]
  T0 = Y0[2]
  Z0 = Y0[3]
  D0 = Y0[4]
  # Light effects
  vChl = A0*(1/CChl)*1000*(1/ZT) # convert g C/m2 to mg Chl/m3
  NPP = Grow.Phyto(vChl,DOC,Load,Zvec,dz) # NPP
  Anet = NPP - s.A*A0  # Net after sinking
  # Consumption fluxes
  Q.AZ = attack*A0*Z0/(1 - attack*handle*A0)
  Q.TZ = attack*T0*Z0/(1 - attack*handle*T0)
  Q.DZ = attack*D0*Z0/(1 - attack*handle*D0)
  Q.ZF = (cF*Z0^2)/(hF^2 + Z0^2)
  # Dynamics
  dAdt = Anet - Q.AZ
  dTdt = I.T - s.T*T0 -(1-pT)*Q.TZ
  dZdt = gAZ*Q.AZ + gTZ*Q.TZ + gAZ*Q.DZ - mZnp*Z0 - Q.ZF + D.Z*(Zref-Z0)
  dDdt = pA*Q.AZ - s.D*D0 - (1-pD)*Q.DZ
  rates = c(dAdt,dTdt,dZdt,dDdt)
  SSE = sum(rates*rates)  # sum of squared distance from equilibrium
  return(rates)  # return either rates or SSE 
}

# Load regressions to predict ZT from DOC, Chl and P Load
# Best-fitting model predicts ZT from DOC and Chl (ZT_DOC.Chl)
# However prediction from DOC alone is almost as good (ZT_DOC)
# Save line:
# save(ZTvec,DOCvec,ChlVvec,Pvec,ZT_DOC.Chl,ZT_DOC.Load,ZT_DOC,
#      file='ZTmodels.Rdata')
load(file='ZTmodels.Rdata')
ZTb = ZT_DOC$coefficients # intercept and slope for ZT ~ DOC model

# Set up driver gradient ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Scaler = seq(0.4,2.5,length.out=12)
NG = length(Scaler)  # number of gradient steps)
cFgrad = rep(0,NG)  # Vector to hold scaled driver
cFbase = cF  # Save the nominal value

# Vectors to hold results
Avec = rep(0,NG)
Tvec = rep(0,NG)
ZBvec = rep(0,NG)
Dvec = rep(0,NG)
Allovec = rep(0,NG)
Lamvec = rep(0,NG)
Zprod = rep(0,NG)
NPPvec = rep(0,NG)

for(iG in 1:NG)  { # Start gradient over parameter value
  
  # Modify the planktivory parameter
  cFgrad[iG] = cFbase*Scaler[iG]
  cF = cFgrad[iG]
  
  # Find equilibria for Experimental conditions
  Y0 = c(AC.AR,TPOCAR,ZB,Dstar) # guesses
  lY0 = log(Y0)
  ATZeq = optim(lY0,dATZdt.eq,method='Nelder-Mead')
  parest = exp(ATZeq$par)
  Avec[iG] = parest[1]
  Tvec[iG] = parest[2]
  ZBvec[iG] = parest[3]
  Dvec[iG] = parest[4]
  
  # Check stability
  JAC = jacobian(dATZdt.jac,parest)
  JAC.lamda = eigen(JAC,only.values=T)
  Lmods = Mod(JAC.lamda$values)
  iLmax = which.max(Lmods) # which eigenvalue has maximum modulus?
  Lamvec[iG] = JAC.lamda$values[iLmax] # Save the eigenvalue with max modulus
  
  # Compute allochthony for estimates
  gQAZ = gAZ*attack*Avec[iG]*ZBvec[iG]/(1 + attack*handle*Avec[iG])
  gQTZ = gTZ*attack*Tvec[iG]*ZBvec[iG]/(1 + attack*handle*Tvec[iG])
  gQDZ = gDZ*attack*Dvec[iG]*ZBvec[iG]/(1 + attack*handle*Dvec[iG])
  Allovec[iG] = gQTZ/(gQTZ + gQAZ + gQDZ)
  
  # Zooplankton secondary production
  Zprod[iG] = gQTZ + gQAZ + gQDZ - mZnp*ZBvec[iG]
  
  # Compute GPP & NPP
  vChl = Avec[iG]*(1/CChl)*1000*(1/ZT) # convert g C/m2 to mg Chl/m3
  GPPtemp = Grow.Phyto(vChl,DOC,Load,Zvec,dz)
  NPPvec[iG] = GPPtemp - s.A*Avec[iG]
  
}

# Plots

windows()
par(mfrow=c(2,2),cex.axis=1.2,cex.lab=1.2,mar=c(5, 4.2, 4, 2) + 0.1)
plot(cFgrad,Avec,type='l',lwd=2,col='forestgreen',
     xlab = 'Planktivory', ylab = 'Phytos')
plot(cFgrad,Tvec,type='l',lwd=2,col='darkred',
     xlab = 'Planktivory', ylab = 'TPOC')
plot(cFgrad,ZBvec,type='l',lwd=2,col='blue',
     xlab = 'Planktivory', ylab = 'Zoopl')
plot(cFgrad,Allovec,type='l',lwd=2,col='sienna',
     xlab = 'Planktivory', ylab = 'Allochthony')

Lsign = sign(Re(Lamvec))
Lamda = Lsign*Mod(Lamvec)
Lsym = rep(19,NG) # symbol for real vs complex
imLam = Im(Lamvec)
Lsym = ifelse(imLam == 0,19,21)

windows()
par(mfrow=c(1,1),cex.axis=1.5,cex.lab=1.5,mar=c(5, 4.2, 4, 2) + 0.1)
plot(cFgrad,Lamda,type='p',pch=Lsym,col='red',cex=1.5,
     xlab='Planktivory',ylab='Max Eigenvalue',
     main='Solid -> real, Open -> complex')
