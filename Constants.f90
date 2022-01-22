 MODULE Constants
!---------------------------------------------------------------------------------------------------------------------------------!
!===========!
!  GENERAL  !
!===========!
   Integer, parameter :: Nage = 58      ! Modern (Nage = 1), 570 Ma (Nage = 58)
   Integer, parameter :: max  = 2000e6  ! Maximum number for time evolution
   Integer, parameter :: resample = 250 ! Total number of resampling for the MC analysis
   Integer, parameter :: maxL = 1000000 ! Maximum number of times for searching the solution
   double precision, parameter :: pi = dacos(-1d0) ! Circular constant

!========!
!  Flag  !
!========!
   Integer, parameter :: Redfield = 1  ! 1 = Redfield stoichiometry
   Integer, parameter :: Ceq = 0
   Integer, parameter :: output = 0

!============!
!  PHYSICAL  !
!============!
   Integer, parameter :: Nj = 60 ! Number of oceanic layers
   double precision, parameter :: Area0   = 3.62d14      ! Ref. value of ocean area [m2]
   double precision, parameter :: AreaCS0 = 0.2712104d14 ! Ref. value of continental shelf area [m2]
   double precision, parameter :: AreaE   = 5.101d14     ! Earth's surface area [m2]
   double precision, parameter :: rl = 0.9d0         ! Areal fraction of low-mid latitude surface layer (L)
   double precision, parameter :: hm = 100d0         ! Depth of mixed layer [m]
   double precision, parameter :: Sm0 = 35.5d0       ! Ref. value of salinity in the low-mid latitude region [psu]
   double precision, parameter :: Sh0 = 34.5d0       ! Ref. value of salinity in the high latitude region [psu]
   double precision, parameter :: Tcm0 = 15d0        ! Ref. value of surface temperature in the low-mid latitude region [C]
   double precision, parameter :: Tch0 =  0d0        ! Ref. value of surface temperature in the high latitude region [C]
   double precision, parameter :: Ts0  = 288.15d0    ! Ref. value of surface temperature [K]
   double precision, parameter :: Tcs0 = 15d0        ! Ref. value of surface temperature [C]
   double precision, parameter :: rho = 1025d0       ! Density of seawater []
   double precision, parameter :: Ku = 5000d0        ! Vertical diffusion coefficient for the upper ocean [m2/yr]
   double precision, parameter :: Kl = 2500d0        ! Vertical diffusion coefficient for the lower ocean [m2/yr]
   double precision, parameter :: Khor = 1d3         ! Horizontal diffusion coefficient [m2/sec]
   double precision, parameter :: dz = 100d0         ! Spatial resolution [m]
   double precision, parameter :: R = 83.1d0         ! Gas constant []
   double precision, parameter :: Vuh = 100d0        ! Convection rate [m/yr]
   double precision, parameter :: Sdw = 0.25d0       ! Fraction of DeepWater Circulation Box (DCB)

!==================!
!  BIOGEOCHEMICAL  !
!==================!
   double precision, parameter :: pCO2a0 = 280d0     ! Ref. value of atmospheric CO2 level [ppmv]
   double precision, parameter :: pO2a0  = 0.2095d0  ! Ref. value of atmospheric O2 level [atm]
   double precision, parameter :: MaO20  = 38d18     ! Ref. value of mass of O2 in the ocean-atmosphere system [mol]
   double precision, parameter :: CmO20  = 0.25405d0 ! Ref. value of [O2] at L [mM]
   double precision, parameter :: SO40   = 28.9d0    ! Ref. value of [SO42-] at L [mM]
   double precision, parameter :: VpO2  = 1000d0     ! Piston velocity of O2 
   double precision, parameter :: VpH2S = 1072d0     ! Piston velocity of H2S
   double precision, parameter :: VpNH4 = 300d0      ! Piston velocity of NH3
   double precision, parameter :: VpCH4 = 1419d0     ! Piston velocity of CH4
   double precision, parameter :: Sh2s0 = 100d0      ! Solubility of H2S
   double precision, parameter :: Snh40 = 56000d0    ! SOlubility of NH3
   double precision, parameter :: Sch40 = 1.4d0      ! Solubility of CH4
   double precision, parameter :: kHth2s = 2100d0    ! 
   double precision, parameter :: kHtnh4 = 4100d0    ! 
   double precision, parameter :: kHtch4 = 1600d0    ! 

!=================!
!  O2 SATURATION  ! (Garcia and Gordon, 1992)
!=================!
   double precision, parameter :: A0 = 2.00907d0     ! Coefficients for O2 saturation concentration
   double precision, parameter :: A1 = 3.22014d0     ! Coefficients for O2 saturation concentration
   double precision, parameter :: A2 = 4.0501d0      ! Coefficients for O2 saturation concentration
   double precision, parameter :: A3 = 4.94457d0     ! Coefficients for O2 saturation concentration
   double precision, parameter :: A4 =-0.256847d0    ! Coefficients for O2 saturation concentration
   double precision, parameter :: A5 = 3.88767d0     ! Coefficients for O2 saturation concentration
   double precision, parameter :: B0 =-6.24523d-3    ! Coefficients for O2 saturation concentration
   double precision, parameter :: B1 =-7.37614d-3    ! Coefficients for O2 saturation concentration
   double precision, parameter :: B2 =-1.0341d-2     ! Coefficients for O2 saturation concentration
   double precision, parameter :: B3 =-8.17083d-3    ! Coefficients for O2 saturation concentration
   double precision, parameter :: C0 =-4.88682d-7    ! Coefficients for O2 saturation concentration

   double precision, parameter :: alphal = 1d0       ! P uptake efficiency @ low-mid lat. surface layer (L)
   double precision, parameter :: alphah = 0.15d0    ! P uptake efficiency @ high lat. surface layer (H)
   double precision, parameter :: Kpo4 = 1d-6        ! Phosphate half saturation constant [mM]
   double precision, parameter :: MuNF = 87.6d0
   double precision, parameter :: MuOt = 91.25d0
   double precision, parameter :: Mor  = 73d0
   double precision, parameter :: Khpo4 = 3d-5       ! 
   double precision, parameter :: Khno3 = 5d-4       ! 
   double precision, parameter :: Rcp0 = 106d0       ! Canonical Redfield C/P ratio
   double precision, parameter :: Rnp0 =  16d0       ! Canonical Redfield N/P ratio
   double precision, parameter :: RcpMax = 400d0     ! Maximum value of C/P ratio
   double precision, parameter :: RnpMax =  60d0     ! Maximum value of N/P ratio
   double precision, parameter :: m1 = 0.72d0        ! Mass fraction of G1
   double precision, parameter :: m2 = 0.25d0        ! Mass fraction of G2
   double precision, parameter :: m3 = 0.03d0        ! Mass fraction of G3
   double precision, parameter :: Ktox = 0d0         ! Temperature dependency of POM degradation (or = 0.15)
   double precision, parameter :: k1o2  = 0.6d0      ! Decomposition rate of G1 via aerobic respiration [1/day]
   double precision, parameter :: k2o2  = 0.1d0      ! Decomposition rate of G2 via aerobic respiration [1/day]
   double precision, parameter :: k3o2  = 0.0d0      ! Decomposition rate of G3 via aerobic respiration [1/day]
   double precision, parameter :: k1no3 = 0.6d0      ! Decomposition rate of G1 via denitrification [1/day]
   double precision, parameter :: k2no3 = 0.1d0      ! Decomposition rate of G2 via denitrification [1/day]
   double precision, parameter :: k3no3 = 0.0d0      ! Decomposition rate of G3 via denitrification [1/day]
   double precision, parameter :: k1so4 = 0.6d0      ! Decomposition rate of G1 via sulfate reduction [1/day]
   double precision, parameter :: k2so4 = 0.1d0      ! Decomposition rate of G1 via sulfate reduction [1/day]
   double precision, parameter :: k3so4 = 0.0d0      ! Decomposition rate of G1 via sulfate reduction [1/day]
   double precision, parameter :: k1ch4 = 0.6d0      ! Decomposition rate of G1 via methanogenesis [1/day]
   double precision, parameter :: k2ch4 = 0.1d0      ! Decomposition rate of G1 via methanogenesis [1/day]
   double precision, parameter :: k3ch4 = 0.0d0      ! Decomposition rate of G1 via methanogenesis [1/day]
   double precision, parameter :: Ko2  = 0.008d0     ! Half saturation constant for aerobic respiration [mM]
   double precision, parameter :: Ko2d = 0.008d0     ! Half saturation constant for aerobic respiration [mM]
   double precision, parameter :: Kno3 = 0.03d0      ! Half saturation constant for denitrification [mM]
   double precision, parameter :: Kno3d= 0.03d0      ! Half saturation constant for denitrification [mM]
   double precision, parameter :: H2Sk = 3650d0      ! Aerobic oxidation rate of hydrogen sulfide [1/(mMxyr)]
   double precision, parameter :: NH4k = 18250d0     ! Aerobic oxidation rate of ammonium [1/(mMxyr)]
   double precision, parameter :: CH4k = 1d7         ! Aerobic methane oxidation rate [1/(mMxyr)]
   double precision, parameter :: CH4km= 1d0         ! Aerobic methane oxidation rate at mixed layers [1/(mMxyr)]
   double precision, parameter :: kAOM = 0.0003d0    ! Methane oxidation rate via AOM [] // 8d0
   double precision, parameter :: Kso4_aom = 0.093d0 ! Sulfate half saturation constant for AOM []
   double precision, parameter :: oxic = 0.25d0      ! Critical O2 level [mM]
   double precision, parameter :: ccpr = 4240d0      ! 
   double precision, parameter :: CNbur = 10d0       ! C/N ratio of buried sediments (or =37.5 or 1/Rnc
   double precision, parameter :: FsMOR0 = 0.0d12    ! H2S input rate via MOR [mol/yr] (Kagoshima et al., 2015) (or = 0.5d12)
   double precision, parameter :: FsVol0 = 0.7d12    ! SO2, H2S input rate via surface volcanisms [mol/yr] (Kagoshima et al., 2015) (or 0.2d12)
   double precision, parameter :: FsVol_H2S0 = 0.3d12 
   double precision, parameter :: FsVol_SO40 = 0.5d12
   double precision, parameter :: fprefPYR = 1.0d0   ! (or = 0.8)
   double precision, parameter :: Sgyp0 = 200d18     ! Ref. value of the crustal reservoir size of gypsum sulfur [mol]
   double precision, parameter :: Spy0  = 200d18     ! Ref. value of the crustal reservoir size of pyrite sulfur [mol]
   double precision, parameter :: porosity = 0.9d0   ! porosity (or = 0.8)
   double precision, parameter :: msed = 2.7d0       ! 
   double precision, parameter :: density = 2.65d0   ! Dry bulk density [] (or = 2.6)
   double precision, parameter :: nyuSo4 = 0.045d0   ! 
   double precision, parameter :: Dso40 = 149d-4     ! 
   double precision, parameter :: nopd = -0.53d0     ! 
   double precision, parameter :: Aopd = 0.4d0       !
   double precision, parameter :: nsmtz = -0.53d0    !
   double precision, parameter :: Asmtz = 0.4d0      !
   double precision, parameter :: ksmtz = 10d0       !
   double precision, parameter :: be1_oxic = 1d0
   double precision, parameter :: be2_oxic = 50d0
   double precision, parameter :: be1_anox = 5d0
   double precision, parameter :: be2_anox = 70d0
   double precision, parameter :: KorgO2 = 0.334195d0
   double precision, parameter :: KpyrO2 = 0.0175157d0
   double precision, parameter :: nw = 1.5d0

!============================!
!  OXYGEN PENETRATION DEPTH  !
!============================!
   double precision, parameter :: a0opd =-2.24869d0   ! Coefficients for OPD
   double precision, parameter :: a1opd = 0.110645d0  ! Coefficients for OPD
   double precision, parameter :: a2opd = 1.12569d0   ! Coefficients for OPD
   double precision, parameter :: a3opd =-0.281005d0  ! Coefficients for OPD
   double precision, parameter :: a4opd = 0.014827d0  ! Coefficients for OPD
   double precision, parameter :: a5opd =-0.124721d0  ! Coefficients for OPD
   double precision, parameter :: a6opd = 0.0894604d0 ! Coefficients for OPD
   double precision, parameter :: a7opd = 0.00279531d0! Coefficients for OPD
   double precision, parameter :: a8opd =-0.127797d0  ! Coefficients for OPD
   double precision, parameter :: a9opd = 0.0017995d0 ! Coefficients for OPD
   double precision, parameter :: a10opd= 0.0085171d0 ! Coefficients for OPD

   double precision, parameter :: a1ch4 = 0.0030084d0 ! 
   double precision, parameter :: a2ch4 =-0.1655405d0 ! 
   double precision, parameter :: a3ch4 = 3.2305351d0 ! 
   double precision, parameter :: a4ch4 =-25.8343054d0! 
   double precision, parameter :: a5ch4 = 71.5397861d0! 

   double precision, parameter :: sesc = 3.7d-5      ! Coefficients for hydrogen escape
   double precision, parameter :: Btcorg0 = 3d12     ! Ref. value of burial rate of terrigenous organic matter [molC/yr]
   double precision, parameter :: Fmorg0 = 1.25d12   ! Ref. value of organic matter degassing rate [molC/yr]
   double precision, parameter :: Worg0 = 13.0286d12 ! Ref. value of oxidative weathering of organic matter [molC/yr]
   double precision, parameter :: Wgyp0 = 1.6d12     ! Ref. value of gypsum weathering [mol/yr]
   double precision, parameter :: Wpy0  = 1.0d12     ! Ref. value of pyrite weathering [mol/yr]
   double precision, parameter :: WFeO0 = 0d12       ! 
   double precision, parameter :: c0FeO = 0.9713d0
   double precision, parameter :: c1FeO =-0.1095d0
   double precision, parameter :: c2FeO =-0.0737d0
   double precision, parameter :: Pso40 = 2.60d12    ! Ref. value of riverine sulfate input rate [mol/yr]
   double precision, parameter :: Ppo40 = 0.155d12   ! Ref. value of riverine reactive P input rate [mol/yr]
   double precision, parameter :: TBcaso40 = 2.1d12  ! Ref. value of gypsum deposition [mol/yr] // 1.84d12
   double precision, parameter :: kpyrite = 1d-3*1000d0 ! Pyrite precipitation rate
   double precision, parameter :: CjFeII = 0.01d0    ! [FeII]
   double precision, parameter :: epyOxic = 0.1173d0 ! Pyrite burial efficiency for oxic sediments
   double precision, parameter :: Corg0  = 1250d18   ! Crustal organic C reservoir size [mol]
   double precision, parameter :: Ccarb0 = 5000d18   ! Crustal carbonate C reservoir size [mol]
!  S isotope
   double precision, parameter :: d34Sin = 0d0       ! 
   double precision, parameter :: Kmd34S = 0.363d0   ! 
   double precision, parameter :: maxDd34S = 30d0    ! 
   double precision, parameter :: maxDd34S_up = 50d0 ! 
   double precision, parameter :: maxDd34S_low = 15d0! 
   double precision, parameter :: Lam14C = 1d0/8267d0! 
!
   double precision, parameter :: NPPt0   = 5000d12  ! Terrestrial NPP [molC/yr]
   double precision, parameter :: fUV = 0.01d0       ! 
   double precision, parameter :: kfire = 3d0        ! Lenton (2013)
   double precision, parameter :: Jch4t0  = 1d12     ! CH4 flux from terrestrial biosphere
   double precision, parameter :: dr = 13.6d18       ! 
   double precision, parameter :: dd = 0.273d18      ! 
   double precision, parameter :: CPland = 1000d0    ! Corg/Porg ratio of land plants
   double precision, parameter :: k11 = BtCorg0 / (BtCorg0+Ppo40*CPland) ! Fraction of P entering to the terrestrial biosphere
   double precision, parameter :: Wp0 = Ppo40 / (1d0-k11) ! P weathering [molP/yr]
   double precision, parameter :: Pland0 = Wp0 * k11   ! P flux to terrestrial biosphere [molP/yr]

!  pCH4 < 1e-6 bar
   double precision, parameter :: c0l6 =-57265.34429022d0
   double precision, parameter :: c1l6 =-34227.66598214d0
   double precision, parameter :: c2l6 =-8511.63546063d0 
   double precision, parameter :: c3l6 =-1126.69934250d0
   double precision, parameter :: c4l6 =-83.71121714d0
   double precision, parameter :: c5l6 =-3.30919540d0
   double precision, parameter :: c6l6 =-0.05436644d0
   double precision, parameter :: c0m6 =-6.00815080d0
   double precision, parameter :: c1m6 = 18.02266008d0 
   double precision, parameter :: c2m6 = 9.13915748d0
   double precision, parameter :: c3m6 = 2.19756947d0
   double precision, parameter :: c4m6 = 0.28479353d0
   double precision, parameter :: c5m6 = 0.01904814d0
   double precision, parameter :: c6m6 = 0.00051560d0
   double precision, parameter :: c0h6 =-21.50360309d0
   double precision, parameter :: c1h6 =-0.84092722d0
   double precision, parameter :: c2h6 = 0.00101886d0
   double precision, parameter :: c3h6 =-0.02287496d0
!  pCH4 = 1e-5 bar
   double precision, parameter :: c0l5 =-40362.30026998d0
   double precision, parameter :: c1l5 =-23417.68514816d0
   double precision, parameter :: c2l5 =-5652.97793198d0
   double precision, parameter :: c3l5 =-726.55027862d0
   double precision, parameter :: c4l5 =-52.43737190d0
   double precision, parameter :: c5l5 =-2.01499702d0
   double precision, parameter :: c6l5 =-0.03220585d0
   double precision, parameter :: c0m5 =-51.44581130d0
   double precision, parameter :: c1m5 =-35.53559269d0
   double precision, parameter :: c2m5 =-16.38760953d0
   double precision, parameter :: c3m5 =-4.05108269d0
   double precision, parameter :: c4m5 =-0.54278446d0
   double precision, parameter :: c5m5 =-0.03725323d0
   double precision, parameter :: c6m5 =-0.00102639d0
   double precision, parameter :: c0h5 =-21.53920974d0
   double precision, parameter :: c1h5 =-0.77956713d0
   double precision, parameter :: c2h5 =-0.00447583d0
   double precision, parameter :: c3h5 =-0.02413154d0
!  pCH4 = 1e-4 bar
   double precision, parameter :: c0l4 = 28140.72396961d0
   double precision, parameter :: c1l4 = 16376.80570866d0
   double precision, parameter :: c2l4 = 3958.07245924d0
   double precision, parameter :: c3l4 = 508.70002283d0
   double precision, parameter :: c4l4 = 36.66431294d0
   double precision, parameter :: c5l4 = 1.40501860d0
   double precision, parameter :: c6l4 = 0.02236533d0
   double precision, parameter :: c0m4 =-23.30525036d0
   double precision, parameter :: c1m4 =-9.12942967d0
   double precision, parameter :: c2m4 =-7.48270512d0
   double precision, parameter :: c3m4 =-2.84047828d0
   double precision, parameter :: c4m4 =-0.51671076d0
   double precision, parameter :: c5m4 =-0.04486088d0
   double precision, parameter :: c6m4 =-0.00149768d0
   double precision, parameter :: c0h4 =-21.85614193d0
   double precision, parameter :: c1h4 =-0.92312602d0
   double precision, parameter :: c2h4 =-0.05625981d0
   double precision, parameter :: c3h4 =-0.03559932d0
!  pCH4 = 1e-3 bar
   double precision, parameter :: c0l3 =-52833.30554892d0
   double precision, parameter :: c1l3 =-30537.38407390d0
   double precision, parameter :: c2l3 =-7335.27739995d0
   double precision, parameter :: c3l3 =-937.13019049d0
   double precision, parameter :: c4l3 =-67.16841179d0
   double precision, parameter :: c5l3 =-2.56125463d0
   double precision, parameter :: c6l3 =-0.04059850d0
   double precision, parameter :: c0m3 = 60.83212721d0
   double precision, parameter :: c1m3 = 99.29765579d0
   double precision, parameter :: c2m3 = 48.25740585d0
   double precision, parameter :: c3m3 = 11.72742599d0
   double precision, parameter :: c4m3 = 1.52684085d0
   double precision, parameter :: c5m3 = 0.10163283d0
   double precision, parameter :: c6m3 = 0.00271414d0
   double precision, parameter :: c0h3 =-21.87862669d0
   double precision, parameter :: c1h3 =-0.51388277d0
   double precision, parameter :: c2h3 = 0.31136680d0
   double precision, parameter :: c3h3 = 0.03329049d0
!  pCH4 = 2e-3 bar
   double precision, parameter :: c0l23 = 33356.67942747d0
   double precision, parameter :: c1l23 = 18811.76791819d0
   double precision, parameter :: c2l23 = 4411.95131782d0
   double precision, parameter :: c3l23 = 551.11718816d0
   double precision, parameter :: c4l23 = 38.67404407d0
   double precision, parameter :: c5l23 = 1.44556441d0
   double precision, parameter :: c6l23 = 0.02248493d0
   double precision, parameter :: c0m23 = 47.51821659d0
   double precision, parameter :: c1m23 = 87.21966045d0
   double precision, parameter :: c2m23 = 44.52985084d0
   double precision, parameter :: c3m23 = 11.32559623d0
   double precision, parameter :: c4m23 = 1.53947964d0
   double precision, parameter :: c5m23 = 0.10678413d0
   double precision, parameter :: c6m23 = 0.00296709d0
   double precision, parameter :: c0h23 =-22.06530949d0
   double precision, parameter :: c1h23 =-0.97218386d0
   double precision, parameter :: c2h23 = 0.10592109d0
   double precision, parameter :: c3h23 = 0.00207026d0
   
 END MODULE Constants
