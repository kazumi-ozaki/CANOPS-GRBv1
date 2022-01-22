    PROGRAM CANOPS_OPENO2

    USE Constants ! Summarizing constant values used in this model (see Constants.f90).
    IMPLICIT NONE
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!                                                                      !
!   main.f90                                                           !
!   Constants.f90                                                      !
!                                                                      !
!   FINAL UPDATE: 12/30/2021                                           !
!                                                                      !
!   Code credit: Kazumi Ozaki                                          !
!   Department of Environmental Science                                !
!   Toho University, Japan                                             !
!   E-mail: kazumi.ozaki@sci.toho-u.ac.jp                              !
!                                                                      !
!----------------------------------------------------------------------!
!   It is advisable to read the README.txt                             !
!----------------------------------------------------------------------!

!======================!
! VARIABLE DECLARATION !
!======================!
! General variables
    Integer I,J,N,Iw,Isdw,Iredox,Ia,Ie,Ics,Ipo2,Isss,Isst,Ivolc,Isr,Iage,Ix
    Integer river,Idur,Nn,NnS,Imc,seedsize,clock,Nloop,Ierror
    Integer Ic_success,Ic_false
    Integer Idammy,Ipo2init
    Integer(selected_int_kind(r=16)) Ncount
    Integer, allocatable :: seed(:)
    double precision Time,dt,age,Ages(Nage) ! Time
    double precision fR(Nage),fE(Nage),fAD(Nage),fLA(Nage),fD(Nage),fmor(Nage)     ! Forcing parameter
    double precision Geog(Nage),fA(Nage),fL(Nage),fCa(Nage),facs(Nage),fpo2a(Nage) ! Forcing parameter
    double precision dammy,xOut,yOUT,const,const1,const2,const3  ! Others
! Physical variables
    double precision Area,Areal,Areah,AreaCS,Sdwb,Sz(Nj+1),Sdwz(Nj+1) ! Surface area
    double precision Vocean,V,Vm,Vh,Vdw,Vdwz(Nj),Vj(Nj),Vdwj(Nj)      ! Ocean volume
    double precision yz(Nj+1),yz0(Nj+1),Az(Nj+1),Fz(Nj+1),Fdwz(Nj+1),TFz,TFdwz,Ap,ups(Nj),tups,AzDW(Nj+1) ! Geometry
    double precision RI,RIdw,gRI   ! Redox index
    double precision w,wz(Nj),Kz(Nj),Kzh(Nj),dVhor  ! Transport
    double precision Tcir,Sv,Hor(Nj),Wconst,Wcz(Nj) ! Circulation
    double precision Tm(2),Tcm(2),Th(2),Tch(2),Ts(2),Tcs(2),Temp(2,Nj),TempDW(2,Nj),dSST,Tc(Nj),Tk(Nj)  ! Temperature
    double precision Tempc(Nj),TempDWc(Nj),Tempk(Nj),TempDWk(Nj) ! Temperature
    double precision Sm,Sh,dSSS,Salj(2,Nj),SaljDW(2,Nj),Salm(2),Salh(2),S(Nj) ! Salinity
    double precision P(Nj) ! Pressure
    double precision fsr,dlogSR,SR(Nj),SR0(Nj) ! Sedimentation rate
! Gas species in the atmosphere
    double precision pO2a(2),MaO2(2) ! O2
    double precision pH2S,pNH3,pCH4,kch4ox ! Atmospheric H2S, NH3, CH4
    double precision Mch4(2),Mch40  ! CH4
! Gas exchange between the ocean and the atmospehre
    double precision FoaO2m,FoaH2Sm,FoaNH4m,FoaCH4m     ! Gas exchange flux
    double precision FoaO2h,FoaH2Sh,FoaNH4h,FoaCH4h     ! Gas exchange flux
    double precision FoaO2cs,FoaH2Scs,FoaNH4cs,FoaCH4cs ! Gas exchange flux
    double precision Sh2s,Snh4,Sch4 ! Solubility
    double precision Aso4 ! H2S eflux from the ocean to the atmosphere (all H2S will be oxidized to SO4 and goes back to the ocean)
! Biological productivity
    double precision NPPl,NPPh ! Net primary production [molC/yr]
    double precision fexportL,fexportH     ! Export ratio
    double precision Po,Pol,Poh,Fo,Fol,Foh ! Export production [molC/yr]
    double precision PoNH4,PoNH4h ! New production based on NH4
    double precision Rcp,Rnp,Rpc,Rnc,Roc           ! Stoichiometry @ L
    double precision RcpH,RnpH,RncH,RpcH,RocH      ! Stoichiometry @ H
    double precision RcpDW,RnpDW,RncDW,RpcDW,RocDW ! Stoichiometry @ DW
! Fluxes
    double precision Ppo4,Po2,Pno3,Ano3,Pca,Pca0             ! River input flux
    double precision Ppo4M,Pno3M,PcaM,Ppo4H,Pno3H,PcaH,IrivN ! River input flux
    double precision Pso4,Pso4M,Pso4H,Pnh4,Pnh4M,Pnh4H,Anh4  ! River input flux
    double precision BpCa(Nj),BpFe(Nj),BpOrg(Nj)   ! P burial
    double precision Depo(Nj+1),RETURNP(Nj),freeP(Nj+1),BurP(Nj+1) ! P flux
    double precision Worg,Fmorg,Btcorg
    double precision FsMOR,FsVol,FsVol_H2S,FsVol_SO4,Fred ! Volcanic S input rate
    double precision Wgyp                                      ! Gypsum weathering
    double precision WpyBio,WpyAbio,WpyBio0,WpyAbio0,krso4,Wpy ! Pyrite weathering
    double precision fBIFP,BpBIFP ! P scavenging
    double precision Fnfix,Fndef,FnfixH,FndefH
! Distribution of chemical species
    double precision CmO2s,ChO2s,O2js(Nj+1),O2dwjs(Nj+1),SO4Ls(28),SO4Hs(28)  ! Reference value
    double precision CmPO4(2), CjPO4(2,Nj), CdwPO4(2,Nj), ChPO4(2), gPO4j(Nj) ! [PO43-]
    double precision CmNO3(2), CjNO3(2,Nj), CdwNO3(2,Nj), ChNO3(2), gNO3j(Nj) ! [NO3-]
    double precision CmO2(2) ,  CjO2(2,Nj),  CdwO2(2,Nj),  ChO2(2),  gO2j(Nj) ! [O2]
    double precision CmSO4(2), CjSO4(2,Nj), CdwSO4(2,Nj), ChSO4(2), gSO4j(Nj) ! [SO42-]
    double precision CmCH4(2), CjCH4(2,Nj), CdwCH4(2,Nj), ChCH4(2), gCH4j(Nj) ! [CH4]
    double precision CmNH4(2), CjNH4(2,Nj), CdwNH4(2,Nj), ChNH4(2), gNH4j(Nj) ! [NH4+]
    double precision CmH2S(2), CjH2S(2,Nj), CdwH2S(2,Nj), ChH2S(2), gH2Sj(Nj) ! [H2S]
    double precision Cm14C(2), Cj14C(2,Nj), Cdw14C(2,Nj), Ch14C(2), g14Cj(Nj) ! [D14C]
    double precision Dm14C, Dj14C(Nj), Ddw14C(Nj), Dh14C, Dg14Cj(Nj)          ! [D14C]
    double precision PO4j_start(Nj+1,2),NO3j_start(Nj+1,2),NH4j_start(Nj+1,2)
    double precision CH4j_start(Nj+1,2),SO4j_start(Nj+1,2),H2Sj_start(Nj+1,2),O2j_start(Nj+1,2)
    double precision Tempj_start(Nj+1,2),Salj_start(Nj+1,2),C14j_start(Nj+1,2)
! Reservoir size
    double precision mixMo2,mixMso4 ! Mass in the mixed layer
    double precision Mpo4,Mno3,Mo2,Mh2s,Mnh4,Mso4,Mch4_ocn ! Total mass in the ocean
    double precision Corg(2),Sgyp(2),Spy(2) ! Sedimentary reservoirs
! S isotopes
    double precision d34S(2),Dd34S,d34S_up(2),Dd34S_up,d34S_low(2),Dd34S_low,d34S_c(2),Dd34Sc,d34Sgyp(2),d34Spy(2)
    double precision d34S_calc,d34Spy_calc,d34S_calc_up,d34Spy_calc_up,d34S_calc_low,d34Spy_calc_low
! Biological pump in the water column
    double precision Vpom,Ppom1,Ppom2,Ppom3,Bpom ! POM biological pump
    double precision Fpom1(Nj+1),Fpom2(Nj+1),Fpom3(Nj+1),pom,FdwPOM1(Nj+1),FdwPOM2(Nj+1),FdwPOM3(Nj+1) ! POM sinking flux density
    double precision dis1(Nj),dis2(Nj),dis3(Nj),dispom(Nj),dis1DW(Nj),dis2DW(Nj),dis3DW(Nj),dispomDW(Nj)
    double precision disdw1,disdw2,disdw3,disdw ! POM dissolution flux
    double precision   r1o2(Nj),  r1no3(Nj),  r1so4(Nj),  r1ch4(Nj)
    double precision   r2o2(Nj),  r2no3(Nj),  r2so4(Nj),  r2ch4(Nj)
    double precision   r3o2(Nj),  r3no3(Nj),  r3so4(Nj),  r3ch4(Nj)
    double precision r1o2dw(Nj),r1no3dw(Nj),r1so4dw(Nj),r1ch4dw(Nj)
    double precision r2o2dw(Nj),r2no3dw(Nj),r2so4dw(Nj),r2ch4dw(Nj)
    double precision r3o2dw(Nj),r3no3dw(Nj),r3so4dw(Nj),r3ch4dw(Nj)
    double precision     r1(Nj),     r2(Nj),     r3(Nj)
    double precision   r1dw(Nj),   r2dw(Nj),   r3dw(Nj)
    double precision k1o2dw,k1no3dw,k1so4dw,k2o2dw,k2no3dw,k2so4dw,k3o2dw,k3no3dw,k3so4dw,k1ch4dw,k2ch4dw,k3ch4dw !
! Sediments
    double precision dispomSED(Nj),dispomdwSED(Nj)
    double precision Ro2SED(Nj),Rno3SED(Nj),Rso4SED(Nj),Rch4SED(Nj),dRo2SED(Nj),dRno3SED(Nj),dRso4SED(Nj),dRch4SED(Nj)
    double precision Fdeni(Nj),FdeniDW(Nj),Tdeni,TdeniCS ! Benthic denitrification
    double precision Bpyr(Nj),Bdwpyr(Nj),TBpyr,BpyrWC(Nj),BdwpyrWC(Nj),TBpyrWC,fpyrAn,fpyrOx,epyAnox,BpyrWCL,BpyrWCH ! Pyrite burial
    double precision CSj(Nj),CSdwj(Nj) ! C/S ratio
    double precision TBcaso4 ! Gypsum deposition
    double precision cpr(Nj),CPbur,CPburDW ! C/P ratio of buried sediments
    double precision AOMsed(Nj),AOMsedDW(Nj),BSRsed,BSRsedj(Nj),BSRsedjdw(Nj),ksed(Nj)
    double precision topd(Nj),topdDW(Nj),Zopd(Nj),ZopdDW(Nj)
    double precision tsmtz(Nj),tsmtzDW(Nj),Zsmtz(Nj),ZsmtzDW(Nj)
    double precision tb(Nj),tb0(Nj),zbio(Nj)
    double precision beoxic,beanox
    double precision BE(Nj),BEdw(Nj),BE003,BE001,BEdw003,BEdw001,BEs(Nj),BEdws(Nj)
    double precision BEporg(Nj),BEporgDW(Nj),BEporg0(Nj),BEporgDW0(Nj),aveBEc ! Burial efficiencies
    double precision BCorg(Nj),TBCorg,RETURNC(Nj),BCorgDW(Nj),TBCorgDW,RETURNCDW(Nj) ! Benthic C flux
    double precision stafreeP(Nj),TBurP,Trep,TdepoP,Torgp,Tcap,Tfep
    double precision DepoDW(Nj+1),freePdw(Nj+1),TdepoPdw,cprdw(Nj)
    double precision SdepoP,BpCadw(Nj),BpFedw(Nj),BpOrgdw(Nj),RETURNPdw(Nj)
    double precision BurPdw(Nj+1),TburPdw,TrepDw,SBurP
! Reactions
    double precision Ro2(Nj),Rno3(Nj),Rso4(Nj),Rch4(Nj),TRo2,TRno3,TRso4,TRch4
    double precision H2S(Nj),NH4(Nj),O2(Nj),CH4(Nj),H2SDW(Nj),NH4DW(Nj),O2DW(Nj),CH4DW(Nj) ! Reaction
    double precision Kso4 ! 
    double precision KapH2S(Nj),KapNH4(Nj),KapO2(Nj),KapCH4(Nj),KapSO4(Nj),KapNO3(Nj)
    double precision KapNO3DW(Nj),KapH2SDW(Nj),KapNH4DW(Nj),KapO2DW(Nj),KapCH4DW(Nj),KapSO4DW(Nj) !
    double precision KapNH4m,KapH2Sm,KapNH4h,KapH2Sh
    double precision dNH4,dH2S,dO2,dCH4
    double precision dRo2(Nj),dRno3(Nj),dRso4(Nj),dRch4(Nj)
    double precision KapCH4m,KapCH4h,KapSO4m,KapSO4h
    double precision AOMsum
    double precision NH4m,H2Sm,CH4m,AOMm,NH4h,H2Sh,CH4h,AOMh
    double precision tTRo2,tTRno3,tTRso4,tTRch4
    double precision tAOMsed,tAOMwc,tBSRwc,tH2Sox
    double precision KapO2m,KapO2h,KapNO3m,KapNO3h
    double precision Anew,Bnew,Cnew,Dnew,Enew
    double precision anew0,bnew0,cnew0,dnew0,enew0,fnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0,pnew0,qnew0
    double precision AOM(Nj),AOMdw(Nj)
    double precision Ro2det,Rno3det,Rso4det,Rch4det,Ro2detH,Rno3detH,Rso4detH,Rch4detH
! Monte Carlo simulations
    double precision rnd1,rnd2,Xrnd1,Xrnd2,mean,sigma
    double precision Spy_rnd(resample),epy_rnd(resample),Kso4_rnd(resample),W_rnd(resample),Sgyp_rnd(resample)
    double precision pO2a_rnd(resample),Vpom_rnd(resample),Rcp_rnd(resample),fSR_rnd(resample),fBIFP_rnd(resample)
    double precision fexport_rnd(resample),pO2init_rnd(resample),SO4init_rnd(resample),Fred_rnd(resample)
! Other
    double precision RedoxBocn,RedoxBatm,RedoxB ! Global redox balance index
    double precision jSO4,jH2S,tSO4in,tSO4out,tH2Sin,tH2Sout ! Fluxes and residence time of marine S species
    double precision minMo2,maxMh2s,maxT,maxCP ! Min/Max values
    double precision Time_dammy(10000),pO2_dammy(10000),pCH4_dammy(10000),CmSO4_dammy(10000),d34S_dammy(10000),d34Spy_dammy(10000)
    double precision EX_dammy(10000),NPP_dammy(10000),Mp_dammy(10000),Nfix_dammy(10000),Bcorg_dammy(10000),BSRsed_dammy(10000)
    double precision resO2_dammy(10000),resNO3_dammy(10000),resSO4_dammy(10000),resCO2_dammy(10000)
    double precision GRBg_dammy(10000),GRBatm_dammy(10000),GRBocn_dammy(10000)
    double precision NPPt,Vege,Vnpp,VnppO2,VnppCO2,VnppT,Wp,Pland,Jch4t
    double precision Fremt,gamma_t,delta_t,fch4_t,gch4_t,go2_t

!=========!
!  SETUP  !
!=========!
    CALL SETUP(dt,pCH4,Mch40,Mch4,Rcp,Rnp,Rnc,Rpc,Roc,RcpH,RnpH,RncH,RpcH,RocH,RcpDW,RnpDW,RncDW,RpcDW,RocDW &
              ,Sgyp,Spy,d34S,d34S_up,d34S_low,d34S_c,d34Sgyp,FsMOR,FsVol,Salm,Salh,Wgyp,WpyBio0,WpyAbio0,Wpy &
              ,WpyBio,WpyAbio,TBcaso4,Pso4,Ts,Tcs,Tcm,Tm,Tch,Th,Ppo4)

!=======================!
!  READING INPUT FILES  !
!=======================!
    CALL READ(PO4j_start,NO3j_start,NH4j_start,H2Sj_start,CH4j_start,O2j_start,SO4Ls,SO4Hs,SO4j_start  &
             ,Ages,fR,fE,fAD,fLA,fD,fmor,Geog,fA,fL,fCa,facs,fpo2a,tb,Tk,Tc,S,P,BE,Kz,SR0              &
             ,CjPO4,CdwPO4,CjNO3,CdwNO3,CjO2,CdwO2,CjSO4,CdwSO4,CjNH4,CdwNH4,CjH2S,CdwH2S,Temp,Tempdw  &
             ,Dj14C,Ddw14C,Cj14C,Cdw14C,TempC,Tempk,TempDWc,TempDWk,yz0,Kzh,Ku,Kl,Nj,Nage,CjCH4,CdwCH4 &
             ,O2js,O2dwjs,Tempj_start,Salj_start,C14j_start)

!========================!
!  OPENING OUTPUT FILES  !
!========================!
    OPEN( 1, File = './output/random.dat')
    OPEN( 2, File = './output/all.dat')
    OPEN( 3, File = './output/false.dat')
    OPEN( 4, File = './output/success_story.dat')
    OPEN( 7, File = './output/Xj.dat')
    OPEN( 8, File = './output/Xm.dat')
    OPEN( 9, File = './output/O2_H2S.dat')
    OPEN(10, File = './output/pXa.dat')
    OPEN(11, File = './output/NPP.dat')
    OPEN(12, File = './output/mass.dat')
    OPEN(13, File = './output/Tair.dat')
    OPEN(14, File = './output/river.dat')
    OPEN(15, File = './output/crust.dat')
    OPEN(16, File = './output/GRBatm.dat')
    OPEN(17, File = './output/GRBocn.dat')
    OPEN(18, File = './output/S_balance.dat')
    OPEN(19, File = './output/S_flux.dat')
    OPEN(20, File = './output/MSR_AOM.dat')
    OPEN(21, File = './output/MSRsedj.dat')
    OPEN(22, File = './output/AOMj.dat')
    OPEN(23, File = './output/Foa.dat')
    OPEN(24, File = './output/fpom.dat')
    OPEN(25, File = './output/Bcorgj.dat')
    OPEN(26, File = './output/Bcorg.dat')
    OPEN(27, File = './output/Bpj.dat')
    OPEN(28, File = './output/Bp.dat')
    OPEN(29, File = './output/BEjs.dat')
    OPEN(30, File = './output/BEj.dat')
    OPEN(31, File = './output/CorgPorgj.dat')
    OPEN(32, File = './output/CorgPreaj.dat')
    OPEN(33, File = './output/tCorgPrea.dat')
    OPEN(34, File = './output/BenthicFlux.dat')
    OPEN(35, File = './output/OPD.dat')
    OPEN(36, File = './output/respiration.dat')
    OPEN(37, File = './output/Nfix.dat')
    OPEN(38, File = './output/denitrification.dat')
    OPEN(39, File = './output/Tdeni.dat')
    OPEN(40, File = './output/d34S.dat')
    OPEN(41, File = './output/RedoxIndex.dat')
    OPEN(42, File = './output/Ir.dat')
    OPEN(43, File = './output/Ir=0.dat')
    OPEN(44, File = './output/Ir=1.dat')
    OPEN(45, File = './output/Ir=2.dat')
    OPEN(46, File = './output/Ir=3.dat')
    OPEN(47, File = './output/Ir=4.dat')
    OPEN(48, File = './output/Ir=5.dat')
    OPEN(49, File = './output/land.dat')
    OPEN(97, File = './output/check.dat')
    OPEN(98, File = './output/setting.dat')
    OPEN(99, File = './output/file_count.dat')

    Ix = 1
! Age setting
    Do Iage = 1, 1, 1

    Ic_success = 0
    Ic_false   = 0
!=============================== Start Monte Carlo analysis ===============================!
!    Do Imc = 1, resample
    Do Imc = 1, 1
    !____________________________!
    ! Modern condition      !
      W_rnd(Imc)    = 10d0**(-0d0/10d0)
      Fred_rnd(Imc) = 0d0
      Spy_rnd(Imc)  = 200d0
      Sgyp_rnd(Imc) = 200d0
      fsr_rnd(Imc)  = 1d0
      epy_rnd(Imc)  = 1d0
      Kso4_rnd(Imc) = 0.2d0
      Vpom_rnd(Imc) = 100d0
      fBIFP_rnd(Imc) = 0d0
      fexport_rnd(Imc) = 0.2d0
      pO2init_rnd(Imc) = 1d0
      SO4init_rnd(Imc) = SO40
    !____________________________!

      Corg(1)   = Corg0*1d0
      Fmorg = Fmorg0 * Corg(1)/Corg0
      Fred      = Fred_rnd(Imc) * 1d12
      fexportL  = fexport_rnd(Imc)
      fexportH  = fexport_rnd(Imc)
      Spy(1)    = Spy_rnd(Imc) *1d18
      Sgyp(1)   = Sgyp_rnd(Imc)*1d18
      FsVol_H2S = FsVol_H2S0 * Spy(1)/Spy0
      FsVol_SO4 = FsVol_SO40 * Sgyp(1)/Sgyp0
      d34Spy(1) = -d34Sgyp(1) * (Sgyp(1)+Mso4)/Spy(1)
      Wpy       = Wpy0    * Spy(1)/Spy0 * fsr_rnd(Imc)
      WpyBio    = WpyBio0 * Spy(1)/Spy0 * fsr_rnd(Imc) * fprefPYR
      WpyAbio   = WpyAbio0* Spy(1)/Spy0 * fsr_rnd(Imc)
      Wgyp      = Wgyp0 * Sgyp(1)/Sgyp0 * fsr_rnd(Imc)
      epyAnox   = epy_rnd(Imc)
      Kso4 = Kso4_rnd(Imc)
      Vpom = Vpom_rnd(Imc)
      Rcp = Rcp0
      Rnp = Rnp0
      Rnc = Rnp/Rcp
      Rpc = 1d0/Rcp
      Roc = (Rcp+2d0*Rnp)/Rcp
      CPbur = 2d0*Rcp
      CPburDW = 2d0*Rcp
      fBIFP = fBIFP_rnd(Imc)

!________________________Do loop for sedimentation rate__________________________!
!                                                                                !
      Do Isr = 10, 10, -2
        fsr = Isr/10d0 ! Sedimentation rate factor
        Do J = 1, Nj
           SR(J) = fsr_rnd(Imc) * SR0(J) * fR(Iage) ! Sedimentation rate [cm/yr]
           ksed(J) = 0.02d0
           tb(J)  = 10d0/SR(J)  ! Timescale for burial [yr] 
           tb0(J) = 10d0/SR0(J) ! Reference value of burial timescale [yr]
        End Do

        RcpDW = Rcp
        CPbur   = 2d0*Rcp
        CPburDW = 2d0*RcpDW
        const  = 1d0/oxic
        const1 = CPbur*ccpr
        const2 = const*ccpr
        const3 = CPburDW*ccpr
        Do J = 1, Nj
          cpr(J) = const1/(O2js(J+1)*const2+(1d0-O2js(J+1)*const)*CPbur)*(1d0+dexp(-tb0(J)*1d-4))*0.5d0
          If(O2js(J+1) > oxic) then
            cpr(J) = CPbur*(1d0+dexp(-tb0(J)*1d-4))*0.5d0
          End If
          cprDW(J) = const3/(O2dwjs(J+1)*const2+(1d0-O2dwjs(J+1)*const)*CPburDW)*(1d0+dexp(-tb0(J)*1d-4))*0.5d0
          If(O2dwjs(J+1) > oxic) then
            cprDW(J) = CPburDW*(1d0+dexp(-tb0(J)*1d-4))*0.5d0
          End If
        End Do
        const  = (1d0-porosity)*density
        const1 = 1d0/0.05d0
        beoxic = be1_oxic - be2_oxic
        beanox = be1_anox - be2_anox
        Do J = 1, Nj
          If(O2js(J+1) >= 0.2d0) then
            BEs(J) = (SR0(J)**0.4d0)/2.1d0
          Else If((O2js(J+1) < 0.2d0).and.(O2js(J+1) >= 0.03d0)) then
            BE003 = (SR0(J)**0.4d0)/2.1d0
            BE001 = (beanox/(1d0+SR0(J)*const*200d0)+be2_anox)*0.01d0
            BEs(J)= 10d0**((dlog10(BE003)-dlog10(BE001))/(dlog10(0.2d0)-dlog10(0.03d0))  &
                          *(dlog10(O2js(J+1))-dlog10(0.03d0))+dlog10(BE001))
          Else
            BEs(J) = (beanox/(1d0+SR0(J)*const*200d0)+be2_anox)*0.01d0
          End If
          If(O2dwjs(J+1) >= 0.2d0) then
            BEdws(J) = (SR0(J)**0.4d0)/2.1d0
          Else If((O2dwjs(J+1) < 0.2d0).and.(O2dwjs(J+1) >= 0.03d0)) then
            BEdw003 = (SR0(J)**0.4d0)/2.1d0
            BEdw001 = (beanox/(1d0+SR0(J)*const*200d0)+be2_anox)*0.01d0
            BEdws(J)= 10d0**((dlog10(BEdw003)-dlog10(BEdw001))/(dlog10(0.2d0)-dlog10(0.03d0))  &
                            *(dlog10(O2dwjs(J+1))-dlog10(0.03d0))+dlog10(BEdw001))
          Else
            BEdws(J) = (beanox/(1d0+SR0(J)*const*200d0)+be2_anox)*0.01d0
          End If
          BE(J) = (SR0(J)**0.4d0)/2.1d0
          BEdw(J) = (SR0(J)**0.4d0)/2.1d0
          If(BE(J) > 1d0) then
            BE(J) = 1d0
          End If
          BEporg0(J)  = BEs(J)  * Rcp0/cpr(J)
          BEporgDW0(J)= BEdws(J)* Rcp0/cprDW(J)
          Write(29,706) -hm-J*dz,SR0(J),BEs(J),BEdws(J),BEporg0(J),BEporgDW0(J)
          cpr(J)     = CPbur * (1d0+dexp(-tb(J)*0.0001d0))*0.5d0
          cprDW(J)   = CPbur * (1d0+dexp(-tb(J)*0.0001d0))*0.5d0
          BEporg(J)  = BE(J)  *Rcp0/cpr(J)
          BEporgDW(J)= BEdw(J)*Rcp0/cprDW(J)
        End Do

!________________________Do loop for sea surface temperature (SST)_________________________!
!                                                                                          !
      Do Isst = 0,0,5
        dSST   = Isst/10d0
      ! @ low-mid lat.
        Tcm(1) = Tcm0 + dSST       ! [C]
        Tm(1)  = Tcm(1) + 273.15d0 ! [K]
      ! @ high lat.
        Tch(1) = Tch0 + dSST       ! [C]
        Th(1)  = Tch(1) + 273.15d0 ! [K]

!________________________Do loop for sea surface salinity (SSS)____________________________!
!                                                                                          !
      Do Isss = 0, 0, 15
        dSSS = Isss/10d0
        Sm = Sm0 + dSSS  ! @ low-mid lat.
        Sh = Sh0 + dSSS  ! @ high lat.
        Salm(1) = Sm
        Salh(1) = Sh

!________________________Do loop for atmospheric O2 levels_________________________________!
!                                                                                          !
      Do Ipo2 = 0, 0, 1
        Nn = 18
        pO2a(1) = pO2init_rnd(Imc)*pO2a0
        MaO2(1) = MaO20*pO2a(1)/pO2a0
        write(*,*) 'pO2a(1)=',pO2a(1)

!________________________Do loop for areal extent of continental shelves___________________!
!                                                                                          !
      Do Ics = 100, 100, -3
        AreaCS = (Ics/100d0)*AreaCS0      ! Surface area of continental shelves [m^2]
        Area   = Area0 + AreaCS - AreaCS0 ! Ocean surface area [m^2]
        Areal  = Area*rl                  ! Ocean surface area of low-mid lat. region [m^2]
        Areah  = Area - Areal             ! Ocean surface area of high lat. region [m^2]
        Vm = Area*rl*hm                   ! Volume of low-mid lat. surface water [m^3]
        Vh = Area*(1d0-rl)*hm             ! Volume of high lat. surface water [m^3]

!________________________Do loop for riverine input flux___________________________________!
!                                                                                          !
      Do river = -5, -5, 5
      ! Reactive P input rate: Ppo4 [molP/yr]
        Wp = Wp0 * W_rnd(Imc)
        Ppo4  = (1d0-k11) * Wp
        Ppo4M = rl * Ppo4
        Ppo4H = (1d0-rl) * Ppo4
      ! Nitrogen input rate [molN/yr]
        IrivN = 0d0
        Pno3  = IrivN * Ppo4
        Pno3M = rl * Pno3
        Pno3H = (1d0-rl) * Pno3
        Ano3  = 0d0
        Pnh4  = 0d0
        Pnh4M = rl * Pnh4
        Pnh4H = (1d0-rl) * Pnh4

!________________________Do loop for ocean circulation rate________________________________!
!                                                                                          !
      Do Iw = 40, 40, -3
        w = Iw/2d0*3.1536d13/(Sdw*Area) ! [m/yr]
!__________________________________________________________________________________________!
        Write(*,*) '------------------------'
        Write(*,'(A17,f8.4,A9)') 'P input rate   = ',Ppo4/1d12,' TmolP/yr'
        Write(*,'(A17,f8.3,A8)') 'Shelf area     = ',AreaCS/AreaCS0*100d0,' %modern'
        Write(*,'(A17,f8.4,A4)') 'Atmospheric O2 = ',pO2a(1),' atm'
        Write(*,'(A17,f8.4,A6)') 'Spy            = ',Spy(1) /1d18,' EmolS'
        Write(*,'(A17,f8.4,A6)') 'Sgyp           = ',Sgyp(1)/1d18,' EmolS'
        Write(*,*) 'epy  =',epyAnox
        Write(*,*) 'Kso4 =',Kso4
        Write(*,*) 'Vpom =',Vpom
        Write(*,*) 'Rcp  =',Rcp
        Write(*,*) 'fsr  =',fsr_rnd(Imc)
        Write(*,*) 'fBIFP=',fBIFP
        Write(*,*) '------------------------'

!=========================================================!
!  INITIAL VALUE OF DISSOLVED ELEMENTS IN SURFACE WATERS  !
!=========================================================!
!  Concentration, [X]m !
      CmPO4(1) = PO4j_start(1,1)
      CmNO3(1) = NO3j_start(1,1)
      CmNH4(1) = NH4j_start(1,1)
      CmSO4(1) = SO4init_rnd(Imc)
      CmH2S(1) = H2Sj_start(1,1)
      CmCH4(1) = CH4j_start(1,1)
      Cm14C(1) = C14j_start(1,1)
      ChPO4(1) = PO4j_start(1,2)
      ChNO3(1) = NO3j_start(1,2)
      ChNH4(1) = NH4j_start(1,2)
      ChSO4(1) = SO4init_rnd(Imc)
      ChH2S(1) = H2Sj_start(1,2)
      ChCH4(1) = CH4j_start(1,2)
      Ch14C(1) = C14j_start(1,2)
!  [O2]m = (pO2[atm]) x (S_O2[mol/m^3/atm])
      const   = 1d0 / 22.3916d0
      CmO2(1) = const * dexp(lm(Tcm,Sm))*pO2a(1)/pO2a0
      ChO2(1) = const * dexp(lm(Tch,Sh))*pO2a(1)/pO2a0
      Write(*,'(A11,f9.5)') '[O2]m =    ',CmO2(1)
      Write(*,'(A11,f9.5)') '[O2]h =    ',ChO2(1)

!  Sulfur riverine input rate [molS/yr]
      Pso4  = Wpy + Wgyp + FsVol_H2S + FsVol_SO4
      Pso4M = rl*Pso4       ! @ low-mid lat.
      Pso4h = (1d0-rl)*Pso4 ! @ high lat.
      Write(*,*) 'Riverine SO4 input=',Pso4/1d12,' TmolS/yr'

!--------------------!
!  Initial Profiles  !
!--------------------!
      Do J = 1, Nj
         CjPO4(1,J) = PO4j_start(J+1,1)
         CjNO3(1,J) = NO3j_start(J+1,1)
         CjNH4(1,J) = NH4j_start(J+1,1)
         CjSO4(1,J) = SO4j_start(J+1,1)
         CjH2S(1,J) = H2Sj_start(J+1,1)
         CjCH4(1,J) = CH4j_start(J+1,1)
         CjO2(1,J)  =  O2j_start(J+1,1)
         If(CjO2(1,J) <= 1d-5) then
            CjO2(1,J) = 1d-5
         End If
         CdwPO4(1,J) = PO4j_start(J+1,2)
         CdwNO3(1,J) = NO3j_start(J+1,2)
         CdwNH4(1,J) = NH4j_start(J+1,2)
         CdwSO4(1,J) = SO4j_start(J+1,2)
         CdwH2S(1,J) = H2Sj_start(J+1,2)
         CdwCH4(1,J) = CH4j_start(J+1,2)
         CdwO2(1,J)  =  O2j_start(J+1,2)
         If(CdwO2(1,J) <= 1d-5) then
            CdwO2(1,J)=1d-5
         End If
         Temp(1,J) = Tempj_start(J+1,1)
         Tempc(J)  = Temp(1,J)
         Tempk(J)  = Temp(1,J) + 273.15d0
         Cj14C(1,J) = C14j_start(J+1,1)
         Salj(1,J) = Salj_start(J+1,1)
         TempDW(1,J) = Tempj_start(J+1,2)
         TempDWc(J)  = TempDW(1,J)
         TempDWk(J)  = TempDW(1,J) + 273.15d0
         Cdw14C(1,J) = C14j_start(J+1,2)
         SaljDW(1,J) = Salj_start(J+1,2)
         write(97,*) J,CjPO4(1,J),CjNO3(1,J),CjNH4(1,J),CjSO4(1,J),CjH2S(1,J),CjCH4(1,J),CjO2(1,J)
        BpyrWC(J) = 0d0
        BdwpyrWC(J) = 0d0
      End Do

!=====================!
!  MASS IN THE OCEAN  !
!=====================!
      mixMo2  = hm*(CmO2(1)*Areal+ChO2(1)*Areah)
      mixMso4 = hm*(CmSO4(1)*Areal+ChSO4(1)*Areah)
      Mpo4 = hm*(CmPO4(1)*Areal+ChPO4(1)*Areah)
      Mno3 = hm*(CmNO3(1)*Areal+ChNO3(1)*Areah)
      Mo2  = mixMo2
      Mh2s = hm*(CmH2S(1)*Areal+ChH2S(1)*Areah)
      Mnh4 = hm*(CmNH4(1)*Areal+ChNH4(1)*Areah)
      Mso4 = mixMso4
      Do J = 1, Nj
        Mpo4 = Mpo4 + dz*CjPO4(1,J)*Sz(J) + dz*CdwPO4(1,J)*Sdwz(J)
        Mno3 = Mno3 + dz*CjNO3(1,J)*Sz(J) + dz*CdwNO3(1,J)*Sdwz(J)
        Mo2  = Mo2  + dz*CjO2(1,J)*Sz(J)  + dz*CdwO2(1,J)*Sdwz(J)
        Mh2s = Mh2s + dz*CjH2S(1,J)*Sz(J) + dz*CdwH2S(1,J)*Sdwz(J)
        Mnh4 = Mnh4 + dz*CjNH4(1,J)*Sz(J) + dz*CdwNH4(1,J)*Sdwz(J)
        Mso4 = Mso4 + dz*CjSO4(1,J)*Sz(J) + dz*CdwSO4(1,J)*Sdwz(J)
      End Do

!=================================!
!  Ocean Structure & Water Cycle  !
!=================================!
      CALL STRUCTURE(Az,AzDW,yz,yz0,Nj,Sdw,hm,dz,Sz,Sdwz,Area,Areal,Areah,Fz,Fdwz,TFdwz,Area0,rl,TFz,Sdwb  &
                    ,Wconst,w,V,Vdwz,Vdw,Ap,ups,tups,Wcz,Tcir,wz,Sv,Ics,Ppo4,Ppo40,Vj,Vdwj,dVhor,Hor       &
                    ,Kz,Kzh,Ku,Kl,Vocean,Khor)

      Ierror = 0
      Idammy = 0
      Do Idur  = 4, 4, 1
      Do Ivolc = 4, 4, 1
!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!                        START TIME EVOLUTION                                              !
!                                                                                          !
!==========================================================================================!
!==========================================================================================!
      Ncount = 0
      Time = 0d0
      Do
!====================!
!  TIME MEASUREMENT  !
!====================!
      Ncount = Ncount + 1
      dt = 0.1d0
! TIME [yr]
      Time = Ncount*dt
!
!===================!
!  PROGRESS REPORT  !
!===================!
      If(Mod(Ncount*dt,5d4) == 0) then
        Write(*,*) '----------------------------------------------'
        Write(*,'(f12.2,A11,f9.5,A4)') Time,' yr, pO2a= ',pO2a(2),' atm'
        Write(*,'(A7,f8.4,A10,f8.4)') '[O2]m= ',CmO2(1),' : [O2]h= ',ChO2(1)
        Write(*,'(A22,f8.4,A9)')      'Export production    = ',(Po+Foh*Areah)/1d12,  ' TmolC/yr'
        Write(*,'(A22,f8.4,A7)')      'Corg burial rate     = ',(TBCorg+TBCorgDW)/1d12,'TmolC/yr'
        Write(*,'(A22,f8.4,A2)')      'Corg preservation    = ',(TBCorg+TBCorgDW)/(Po+Foh*Areah)*100d0,' %'
        Write(*,'(A22,E10.3e2,A9)')   'SO4 riverine input   =',Pso4/1d12,   ' TmolS/yr'
        Write(*,'(A22,E10.3e2,A9)')   'Pyrite S burial rate =',TBpyr/1d12,  ' TmolS/yr'
        Write(*,'(A22,E10.3e2,A9)')   'Gypsum S deposition  =',TBcaso4/1d12,' TmolS/yr'
        Write(*,'(A34,E10.3e2,A9)')   'Syngenetic pyrite precipitation =',TBpyrWC/1d12,' TmolS/yr'
        Write(*,'(A16,f10.4)')        'S balance      =',1d0-(TBpyr+TBpyrWC+TBcaso4)/(Pso4+FsMOR)
        Write(*,'(A16,f10.4)')        '[SO4]m balance =',1d0-CmSO4(2)/CmSO4(1)
        Write(*,'(A16,f10.4,A3)')     '[SO4]m         =',CmSO4(1),   ' mM'
        Write(*,'(A16,f10.4,A3)')     '[H2S](10)      =',CjH2S(2,10),' mM'
        Write(*,'(A16,f10.4,A3)')     '[CH4](10)      =',CjCH4(2,10),' mM'
        Write(*,'(A10,f10.4,A9)')     'P scav   =',BpBIFP/1d12,        ' TmolP/yr'
        Write(*,'(A10,f10.4,A9)')     'P input  =',Ppo4/1d12,          ' TmolP/yr'
        Write(*,'(A10,f10.4,A9)')     'P output =',(SburP+BpBIFP)/1d12,' TmolP/yr'
        Write(*,'(A10,f10.4)')        'P balance=',dAbs(1d0-(SburP+BpBIFP)/Ppo4)
        Write(*,'(A10,f10.4,A5)')     'pCH4     =',pCH4*1d6,' ppmv'
        Write(*,'(A10,f10.4)')        'MP/MP0   =',Mpo4/1d15/3d0
        Write(*,'(A16,f10.4)')        'Corg/Preact    =',(TBCorg+TBCorgDW)/(SburP+BpBIFP)
        Write(*,'(A16,f10.4)')        'Corg/Preact@sed=',(TBCorg+TBCorgDW)/SburP
        Write(*,'(A16,f10.4)')        'Worg0          =',(Areal*FoaO2m+Areah*FoaO2h-2d0*(Aso4+Ano3)-2d0*Wpy  &
                                                        -2d0*kch4ox*MaO2(1)*Mch4(1)-sesc*Mch4(1)-Fred &
                                                        -2d0*FsVol_H2S-Fmorg+Btcorg)/1d12
        Write(*,'(A18,f10.4,A7)')     'Nitrogen fixation=',(Fnfix+FnfixH)/1d12*14d0,' TgN/yr'
        Write(*,'(A18,f10.4)')        'Deni@surf        =',(Rno3det*(NPPl-Po)*0.8d0+Rno3detH*(NPPh-Poh)*0.8d0)/1d12
        Write(*,*) '----------------------------------------------'
      End If

!===================================!
!  JUDGMENT OF STEADY-STATE OR NOT  !
!===================================!
      If((Ncount*dt >= 10d9).and.(Mod(Ncount*dt,1d4) == 0)           & ! Time > 100 Myr
        .and.(abs(1d0-(TBcaso4+TBpyr+TBpyrWC)/(Pso4+FsMOR)) <= 1d-2) & ! S budget imbalance <1%P
        .and.(RedoxB <= 0.01d0)) then                                  ! Global Redox Budget <0.01 TmolO2 equiv/yr
        Go to 9999
      End If
      If(N == max) then
        Write(*,*) 'REACHING MAXIMUM TIME! ',TIME/1d6,' Myr'
        Write(*,'(f7.3,f7.3,A9)') Ppo4/1d12, (SburP+BpBIFP)/1d12,' TmolP/yr'
        Ierror = 9
        Go to 9998
      End If

!==================================!
!  BURIAL EFFICIENCY OF ORGANIC C  !
!==================================!
      const = (1d0-porosity)*density
            const1= 1d0/0.05d0
            Do J = 1, Nj
              If(CjO2(1,J) >= 0.2d0) then
                BE(J) = (SR(J)**0.4d0)/2.1d0
              Else If((CjO2(1,J) < 0.2d0).and.(CjO2(1,J) >= 0.03d0)) then
                BE003 = (SR(J)**0.4d0)/2.1d0
                BE001 = (beanox/(1d0+SR(J)*const*200d0)+be2_anox)*0.01d0
                BE(J) = 10d0**((dlog10(BE003)-dlog10(BE001))/(dlog10(0.2d0)-dlog10(0.03d0))  &
                              *(dlog10(CjO2(1,J))-dlog10(0.03d0))+dlog10(BE001))
              Else
                BE(J) = (beanox/(1d0+SR(J)*const*200d0)+be2_anox)*0.01d0
              End If
              If(CdwO2(1,J) >= 0.2d0) then
                BEdw(J) = (SR(J)**0.4d0)/2.1d0
              Else If((CdwO2(1,J) < 0.2d0).and.(CdwO2(1,J) >= 0.03d0)) then
                BEdw003 = (SR(J)**0.4d0)/2.1d0
                BEdw001 = (beanox/(1d0+SR(J)*const*200d0)+be2_anox)*0.01d0
                BEdw(J) = 10d0**((dlog10(BEdw003)-dlog10(BEdw001))/(dlog10(0.2d0)-dlog10(0.03d0))  &
                                *(dlog10(CdwO2(1,J))-dlog10(0.03d0))+dlog10(BEdw001))
              Else
                BEdw(J) = (beanox/(1d0+SR(J)*const*200d0)+be2_anox)*0.01d0
              End If
            End do

!==============!
!  DEGGASSING  !
!==============!
!  Metamorphic oxidation of crustal organic C
    Fmorg = Fmorg0 * Corg(1)/Corg0
!  Pyrite degassing
    FsVol_H2S = FsVol_H2S0 * Spy(1)/Spy0
!  Gypsum degassing
    FsVol_SO4 = FsVol_SO40 * Sgyp(1)/Sgyp0

!==============!
!  WEATHERING  !
!==============!
!  Oxidative weathering of organic C [molC/yr] (ref)
    Worg = Worg0 * (1d0+KorgO2) * (pO2a(1)/pO2a0)/(pO2a(1)/pO2a0+KorgO2*fsr_rnd(Imc)**nw) * fsr_rnd(Imc) * Corg(1)/Corg0

!  Pyrite weathering [molS/yr]
    Wpy     = Wpy0     * Spy(1)/Spy0 * fsr_rnd(Imc) * (1d0+KpyrO2) * (pO2a(1)/pO2a0)/(pO2a(1)/pO2a0+KpyrO2*fsr_rnd(Imc)**nw)
    WpyBio  = WpyBio0  * Spy(1)/Spy0 * fprefPYR * fsr_rnd(Imc) ! Biotic
    WpyAbio = WpyAbio0 * Spy(1)/Spy0 * fsr_rnd(Imc)           ! Abiotic

!  Gypsum weathering [molS/yr]
    Wgyp    = Wgyp0 * Sgyp(1)/Sgyp0 * fsr_rnd(Imc)

!  Riverine SO4 flux (weathering + volcanic flux) [molS/yr]
    Pso4    = Wpy + Wgyp + FsVol_H2S + FsVol_SO4 ! Global
    Pso4M   = rl * Pso4       ! Low-mid lat.
    Pso4h   = (1d0-rl) * Pso4 ! High lat.

!  NH4 input rate
    Pnh4  = (TBcorg+TBcorgDW)/CNbur
    Pnh4M = rl*Pnh4
    Pnh4H = (1d0-rl)*Pnh4

!=======================!
!  Corg burial on land  !
!=======================!
  ! O2-effect
    VnppO2 = 1.5d0 - 0.5d0*pO2a(1)/pO2a0
    If(VnppO2 < 0d0) then
      VnppO2 = 0d0
    End If

  ! NPP w/o fire-effect
    Vnpp = VnppO2

  ! Land biomass with fire-effect
    Vege = Vnpp * kfire / (kfire-1d0+ignit(pO2a))

  ! Landbiomass with UV-effect
    Vege = Vege * dtanh(pO2a(1)/pO2a0/fUV)

  ! NPP of terrestrial biosphere (NPPt0 = 5000 TmolC/yr)
    NPPt = NPPt0 * Vege

  ! P flux to terrestrial biosphere
    Pland = k11 * Vege * Wp

  ! Corg burial on land 
    Btcorg = CPland * Pland
    
  ! CH4 flux from the terrestrial biosphere to the atmosphere
    Jch4t = Jch4t0 * (BtCorg/BtCorg0) ! Jch4t0 = 1 Tmol/yr

    Fremt = NPPt  - Btcorg  ! on the land

  ! Respiration pathway
    gamma_t = MaO2(1) / (MaO2(1)+dr)
    delta_t = MaO2(1) / (MaO2(1)+dd)
    fch4_t  = 1d0 - gamma_t          ! 
    If((NPPt == 0d0).or.(BtCorg == 0d0)) then
      gch4_t = 0
    Else
      gch4_t = 2d0 * Jch4t / ((NPPt-BtCorg)*(1d0-delta_t))
    End If
    If(gch4_t > 1d0) then
      gch4_t = 1d0
      Jch4t = ((NPPt-BtCorg)*(1d0-delta_t))*0.5d0
    End If
    If(gch4_t < 0d0) then
      gch4_t = 0d0
    End If
    go2_t = 1d0 - gch4_t

!==================!
! Riverine P flux  !  
!==================! 
    Ppo4  = (1d0-k11*Vege) * Wp
    Ppo4M = rl*Ppo4
    Ppo4H = (1d0-rl)*Ppo4

!------------------------------------!
!  Gypsum deposition flux [molS/yr]  !
!------------------------------------!
      TBcaso4 = TBcaso40 * fCa(Iage) * CmSO4(1)/SO40

!=============!
! TEMPERATURE !
!=============!
      Tcm(2) = Tcm(1)
      Tm(2)  = Tcm(2) + 273.15d0
      Tch(2) = Tch(1)
      Th(2)  = Tch(2) + 273.15d0

!=============!
! TEMPERATURE !
!=============!
      Salm(2) = Salm(1)
      Salh(2) = Salh(1)
      Sm = Salm(2)
      Sh = Salh(2)

      Cm14C(2) = Cm14C(1)
      Ch14C(2) = Ch14C(1)

!==============!
!  ATMOSPHERE  !
!==============!
!------!
!  O2  !
!------!
! Reaction rate (Claire et al., 2006 Geobiology)
    Call PHOTOCHEM(pO2a,pCH4,kch4ox)
!  Mass [mol]
      MaO2(2) = MaO2(1)+dt*(Areal*FoaO2m + Areah*FoaO2h - 2d0*(Aso4+Ano3)  & ! Gas-exchange
                           - Worg - 2d0*Wpy                                & ! Oxidative weathering
                           - 2d0*kch4ox*MaO2(1)*Mch4(1) - sesc*Mch4(1)     & ! CH4 photolysis and H escape
                           - Fred - 2d0*FsVol_H2S - Fmorg                  &
                           + NPPt - go2_t*Fremt -delta_t*gch4_t*Fremt)       ! External reductant input
!  Partial pressure [atm]
      pO2a(2) = pO2a0*MaO2(2)/MaO20
      If(pO2a(2) <= 0d0) then
        Ierror = -1
        Go to 9998
      End If

!-------!
!  CH4  !
!-------!
!  Mass [mol]
      Mch4(2) = Mch4(1) + dt*(Areal*FoaCH4m + Areah*FoaCH4h                & ! Gas-exchange
                              - kch4ox*MaO2(1)*Mch4(1) - sesc*Mch4(1)      &
                              + 0.5d0*gch4_t*Fremt - 0.5d0*delta_t*gch4_t*Fremt) ! CH4 photolysis and H escape
!  Partial pressure [atm]
      pCH4    = Mch4(2)/(143.818d18+pO2a(2)/pO2a0*38d18)

!---------!
!  Fred !
!---------!
      Fred = sesc*Mch4(2)

!----------------------------------------------------!
!  Atmospheric deposition flux to the ocean surface  !
!----------------------------------------------------!
!  SO4
      Aso4 = Areal*FoaH2Sm + Areah*FoaH2Sh
!  NO3
      Ano3 = Areal*FoaNH4m + Areah*FoaNH4h

!===========================================================================!
!===========================================================================!
!                      EVOLUTION OF SURFACE LAYERS                          !
!===========================================================================!
!===========================================================================!
!-------------------------------------------------------!
! O2 saturation concentration at ocean surface [mol/m3] !
!-------------------------------------------------------!
      const   = 1d0 / 22.3916d0
      CmO2s   = const * dexp(lm(Tcm,Sm))*pO2a(1)/pO2a0
      ChO2s   = const * dexp(lm(Tch,Sh))*pO2a(1)/pO2a0

!===========================!
!  BIOLOGICAL PRODUCTIVITY  !
!===========================!
!  Export production flux density [molC/m2/yr] (Yamanaka and Tajika, 1996)
      Fol = alphal * Rcp  * hm * CmPO4(1)*CmPO4(1)/(CmPO4(1)+Kpo4) ! Low-mid latitude
      Foh = alphah * RcpH * hm * ChPO4(1)*ChPO4(1)/(ChPO4(1)+Kpo4) ! High latitude
      If((Fol < 0d0).or.(Foh < 0d0)) then
        Ierror = 2
        Go to 9998
      End If

!  Export production [molC/yr]
      Po  = Areal * Fol
      Poh = Areah * Foh
      Fo  = Fol

!  Net Primary Production (NPP) [molC/yr]
      NPPl = Po  / fexportL
      NPPh = Poh / fexportH

!  Multi-G model [molC/yr]
      Ppom1 = m1 * Po ! Labile
      Ppom2 = m2 * Po ! Semi-labile
      Ppom3 = m3 * Po ! Inert

!========================================!
!  RESPIRATION PATHWAY IN SURFACE LAYERS !
!========================================!
!  Aerobic respiration
      Ro2det   = CmO2(1) / (Ko2+CmO2(1))
      Ro2detH  = ChO2(1) / (Ko2+ChO2(1))

!  Denitrification
      Rno3det  = Ko2d/(Ko2d+CmO2(1)) * CmNO3(1)/(Kno3+CmNO3(1))
      Rno3detH = Ko2d/(Ko2d+ChO2(1)) * ChNO3(1)/(Kno3+ChNO3(1))

!  Sulfate reduction
      Rso4det  = Ko2d/(Ko2d+CmO2(1)) * Kno3d/(Kno3d+CmNO3(1)) * CmSO4(1)/(Kso4+CmSO4(1))
      Rso4detH = Ko2d/(Ko2d+ChO2(1)) * Kno3d/(Kno3d+ChNO3(1)) * ChSO4(1)/(Kso4+ChSO4(1))

!  Methanogenesis
      Rch4det  = 1d0 - Ro2det - Rno3det - Rso4det
      Rch4detH = 1d0 - Ro2detH - Rno3detH - Rso4detH

!  Solubility
      const = 1d0/298.15d0
      const1= 1d0/Tm(1)
      Sh2s = Sh2s0 * dexp(kHth2s*(const1-const)) ! H2S
      Snh4 = Snh40 * dexp(kHtnh4*(const1-const)) ! NH3
      Sch4 = Sch40 * dexp(kHtch4*(const1-const)) ! CH4

!  Gas-exchange flux density [mol/m2/yr]
      FoaO2m  = VpO2  * (CmO2(1)-CmO2s)
      FoaNH4m = VpNH4 * (CmNH4(1)*0.96d0-Snh4*pNH3)
      FoaH2Sm = VpH2S * (CmH2S(1)*0.026785d0-Sh2s*pH2S)
      FoaCH4m = VpCH4 * (CmCH4(1)-Sch4*pCH4)

!-------!
!  NH4  !
!-------!
      const = Area * Wconst
      const1= dt/(Areal*hm)
      KapNH4m = CmNH4(1) + const1*(                              &
                           const*(CjNH4(1,1)-CmNH4(1))           & ! Advection
                         + Sz(1)*Kz(1)*(CjNH4(1,1)-CmNH4(1))/dz  & ! Diffusion
                         - FoaNH4m*Areal                         & ! Gas-exchange
                         - Rnc*NPPl + Rnc*(NPPl-Po)              & ! NPP and respiration
                         + Pnh4M)                                  ! External input
      If(KapNH4m <= 0d0) then
        CmNH4(2) = 0d0
        Fndef = -(const*CjNH4(1,1)+Sz(1)*Kz(1)*CjNH4(1,1)/dz-FoaNH4m*Areal-Rnc*Po+Pnh4M)
      ! NPP_NH4
        PoNH4 = (CmNH4(1)*Areal*hm/dt+const*(CjNH4(1,1)-CmNH4(1))+Sz(1)*Kz(1)*(CjNH4(1,1)-CmNH4(1))/dz  &
                 -FoaNH4m*Areal+Rnc*(NPPl-Po)+Pnh4M)/Rnc
      Else
        Fndef = 0d0
        PoNH4 = NPPl
      End If

!------!
!  O2  !
!------!
      KapO2m  = CmO2(1)  + const1*(                            &
                           const*(CjO2(1,1)-CmO2(1))           &  ! Advection
                         + Sz(1)*Kz(1)*(CjO2(1,1)-CmO2(1))/dz  &  ! Diffuxion
                         - FoaO2m*Areal                        &  ! Gas-exchange
                         + PoNH4 + Roc*(NPPl-PoNH4)            &  ! Biological production
                         - 5d0/4d0*Fnfix - Ro2det*(NPPl-Po))      ! Biological consumption
!-------!
!  H2S  !
!-------!
      BpyrWCL = kpyrite * CjFeII * CmH2S(1) * Areal * dz ! Pyrite precipitation in the surface layer
      KapH2Sm = CmH2S(1) + const1*(                              &
                           const*(CjH2S(1,1)-CmH2S(1))           & ! Advection
                         + Sz(1)*Kz(1)*(CjH2S(1,1)-CmH2S(1))/dz  & ! Diffuxion
                         - FoaH2Sm*Areal                         & ! Gas-exchange
                         + Rso4det*(NPPl-Po)*0.5d0               & ! Biological production
                         - BpyrWCL)                                ! Pyrite precipitation
!-------!
!  CH4  !
!-------!
      KapCH4m = CmCH4(1) + const1*(                              &
                           const*(CjCH4(1,1)-CmCH4(1))           & ! Advection
                         + Sz(1)*Kz(1)*(CjCH4(1,1)-CmCH4(1))/dz  & ! Diffusion
                         - FoaCH4m*Areal                         & ! Gas-exchange
                         + Rch4det*(NPPl-Po)*0.5d0)                ! Biological production
!-------!
!  SO4  !
!-------!
      KapSO4m = CmSO4(1) + const1*(                              &
                           const*(CjSO4(1,1)-CmSO4(1))           & ! Advection
                         + Sz(1)*Kz(1)*(CjSO4(1,1)-CmSO4(1))/dz  & ! Diffusion
                         + Pso4M + Aso4                          & ! External input
                         - TBcaso4                               & ! Gypsum deposition
                         - Rso4det*(NPPl-Po)*0.5d0)                ! Biological consumption

!-------------------------------------!
!  IMPLICIT SCHEME FOR O2-NH4-H2S-CH4 ! (ref)
!-------------------------------------!
      If(KapNH4m > 0d0) then
        anew0 = KapO2m
        bnew0 =-2d0*dt*NH4k
        cnew0 =-2d0*dt*H2Sk
        dnew0 =-2d0*dt*CH4km
        enew0 = KapNH4m
        fnew0 =-dt*NH4k
        gnew0 = KapH2Sm
        hnew0 =-dt*H2Sk
        inew0 = KapCH4m
        jnew0 =-dt*CH4km
        knew0 =-dt*kAOM
        lnew0 = KapSO4m
        mnew0 =-dt*kAOM
        nnew0 = dt*kAOM
        onew0 = dt*H2Sk
        CALL JacobM(anew0,bnew0,cnew0,dnew0,enew0,fnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0  &
                   ,CmO2,CmNH4,CmH2S,CmSO4,CmCH4,Nj,Kso4_aom,Nloop)
        If(Nloop == maxL) then
          Ierror = 3
          Go to 9998
        End If
      Else
        anew0 = KapO2m
        cnew0 =-2d0*dt*H2Sk
        dnew0 =-2d0*dt*CH4km
        gnew0 = KapH2Sm
        hnew0 =-dt*H2Sk
        inew0 = KapCH4m
        jnew0 =-dt*CH4km
        knew0 =-dt*kAOM
        lnew0 = KapSO4m
        mnew0 =-dt*kAOM
        nnew0 = dt*kAOM
        onew0 = dt*H2Sk
        CALL JacobM2(anew0,cnew0,dnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0  &
                    ,CmO2,CmH2S,CmSO4,CmCH4,Nj,Kso4_aom,Nloop)
        If(Nloop == maxL) then
          Ierror = 3
          Go to 9998
        End If
      End If
!      CmSO4(2) = SO40
      NH4m = NH4k*CmNH4(2)*CmO2(2)
      H2Sm = H2Sk*CmH2S(2)*CmO2(2)
      CH4m = CH4km*CmCH4(2)*CmO2(2)
      AOMm = kAOM*CmSO4(2)/(CmSO4(2)+Kso4_aom)*CmCH4(2)

!-------------------------------!
!  P burial rate by scavenging  !
!-------------------------------!
      fBIFP  = fBIFP_rnd(Imc)
      BpBIFP = fBIFP*(const*CjPO4(1,1)+Sz(1)*Kz(1)*(CjPO4(1,1)-CmPO4(1))/dz)

!-----------------!
!  PO4 [molP/m3]  !
!-----------------!
      CmPO4(2) = CmPO4(1) + const1*(                              &
                            const*(CjPO4(1,1)-CmPO4(1))           & ! Advection
                          + Sz(1)*Kz(1)*(CjPO4(1,1)-CmPO4(1))/dz  & ! Diffusion
                          - BpBIFP                                & ! Scavenging
                          - Rpc*Fol*Areal                         & ! Export production
                          + Ppo4M)                                  ! Riverine input

!-------!
!  NO3  !
!-------!
      CmNO3(2) = CmNO3(1) + const1*(                              &
                            const*(CjNO3(1,1)-CmNO3(1))           & ! Advection
                          + Sz(1)*Kz(1)*(CjNO3(1,1)-CmNO3(1))/dz  & ! Diffusion
                          - Fndef - Rno3det*(NPPl-Po)*0.8d0       & ! Biological consumption
                          + Pno3M + Ano3)                         & ! External input
                          + dt*NH4m                                 ! Nitrification
      If(CmNO3(2) <= 0d0) then
        CmNO3(2) = 0d0
        Fnfix = -(const*(CjNO3(1,1)-CmNO3(1))                      &
                 + Sz(1)*Kz(1)*(CjNO3(1,1)-CmNO3(1))/dz            &
                 - Fndef + Pno3M + Ano3 - Rno3det*(NPPl-Po)*0.8d0  &
                 + (CmNO3(1)+dt*NH4m)/const1)
      Else
        Fnfix = 0d0
      End If

!------------------------------!
!  C/N/P stoichiometry of POM  !
!------------------------------!
    If(Redfield == 1) then
    ! Redfield stoichiometry
      Rcp = Rcp0
      Rnp = Rnp0
    Else
    ! Dynamic stoichiometry
      Rcp = Rcp0 + (RcpMax-Rcp0)*0.5d0*(1d0+dtanh((0.1d-3-CmPO4(1))/0.03d-3))
      Rnp = Rnp0 + (RnpMax-Rnp0)*0.5d0*(1d0+dtanh((0.1d-3-CmPO4(1))/0.03d-3))
    End If
    Rnc = Rnp/Rcp
    Rpc = 1d0/Rcp
    Roc = (Rcp+2d0*Rnp)/Rcp

!--------------------------------!
!  High-latitude surface region  !
!--------------------------------!
!  Solubility
      const = 1d0/298.15d0
      const1= 1d0/Th(1)
      Sh2s = Sh2s0 * dexp(kHth2s*(const1-const)) ! H2S
      Snh4 = Snh40 * dexp(kHtnh4*(const1-const)) ! NH3
      Sch4 = Sch40 * dexp(kHtch4*(const1-const)) ! CH4

!  Gas-exchange flux density [mol/m2/yr]
      FoaO2h  = VpO2  * (ChO2(1)-ChO2s)
      FoaNH4h = VpNH4 * (ChNH4(1)*0.96d0-Snh4*pNH3)
      FoaH2Sh = VpH2S * (ChH2S(1)*0.026785d0-Sh2s*pH2S)
      FoaCH4h = VpCH4 * (ChCH4(1)-Sch4*pCH4)

!-------!
!  NH4  !
!-------!
      const = Area * Wconst
      const1= 1d0/Areah
      const2= dt/hm
      KapNH4h = ChNH4(1) + const2*(                                  &
                           Vuh*(CdwNH4(1,1)-ChNH4(1))                & ! Advection
                         + const*const1*(CmNH4(1)-ChNH4(1))          & ! Diffusion
                         - FoaNH4h                                   & ! Gas-exchange
                         - RncH*NPPh*const1 + RncH*(NPPh-Poh)*const1 & ! Biological production/consumption
                         + Pnh4H*const1)                               ! External input
      If(KapNH4h <= 0d0) then
        ChNH4(2) = 0d0
        FndefH = -(Vuh*CdwNH4(1,1)+const*const1*CmNH4(1)-FoaNH4h   &
                  -RncH*NPPh*const1+RncH*(NPPh-Poh)*const1+Pnh4H*const1)*Areah
        PoNH4h = (ChNH4(2)*Areah*hm/dt+Areah*Vuh*(CdwNH4(1,1)-ChNH4(1))  &
                 +const*(CmNH4(1)-ChNH4(1))-FoaNH4h*Areah+Pnh4H          &
                 +RncH*(NPPh-Poh))/RncH
      Else
        FndefH = 0d0
        PoNH4h = NPPh
      End If

!------!
!  O2  !
!------!
      KapO2h  = ChO2(1)  + const2*(                           &
                           Vuh*(CdwO2(1,1)-ChO2(1))           &
                         + const*const1*(CmO2(1)-ChO2(1))     &
                         - FoaO2h                             &
                         + (PoNH4h+RocH*(NPPh-PoNH4h))*const1 &
                         - 5d0/4d0*FnfixH*const1              &
                         - Ro2detH*(NPPh-Poh)*const1)

!-------!
!  H2S  !
!-------!
      BpyrWCH = kpyrite * CjFeII * ChH2S(1) * Areah * dz
      KapH2Sh = ChH2S(1) + const2*(                         &
                           Vuh*(CdwH2S(1,1)-ChH2S(1))       &
                         + const*const1*(CmH2S(1)-ChH2S(1)) &
                         - FoaH2Sh                          &
                         + Rso4detH*(NPPh-Poh)*0.5d0*const1 &
                         - BpyrWCH*const1)

!-------!
!  CH4  !
!-------!
      KapCH4h = ChCH4(1) + const2*(                          &
                           Vuh*(CdwCH4(1,1)-ChCH4(1))        &
                         + const*const1*(CmCH4(1)-ChCH4(1))  &
                         - FoaCH4h                           &
                         + Rch4detH*(NPPh-Poh)*0.5d0*const1)

!-------!
!  SO4  !
!-------!
      KapSO4h = ChSO4(1) + const2*(                          &
                           Vuh*(CdwSO4(1,1)-ChSO4(1))        &
                         + const*const1*(CmSO4(1)-ChSO4(1))  &
                         + Pso4H*const1 - Rso4detH*(NPPh-Poh)*0.5d0*const1)

!=====================================!
!  IMPLICIT SCHEME FOR O2-NH4-H2S-CH4 ! (ref)
!=====================================!
      If(KapNH4h > 0d0) then
        anew0 = KapO2h
        bnew0 =-2d0*dt*NH4k
        cnew0 =-2d0*dt*H2Sk
        dnew0 =-2d0*dt*CH4km
        enew0 = KapNH4h
        fnew0 =-dt*NH4k
        gnew0 = KapH2Sh
        hnew0 =-dt*H2Sk
        inew0 = KapCH4h
        jnew0 =-dt*CH4km
        knew0 =-dt*kAOM
        lnew0 = KapSO4h
        mnew0 =-dt*kAOM
        nnew0 = dt*kAOM
        onew0 = dt*H2Sk
        CALL JacobM(anew0,bnew0,cnew0,dnew0,enew0,fnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0  &
                   ,ChO2,ChNH4,ChH2S,ChSO4,ChCH4,Nj,Kso4_aom,Nloop)
        If(Nloop == maxL) then
          Ierror = 3
          Go to 9998
        End If
      Else
        anew0 = KapO2h
        cnew0 =-2d0*dt*H2Sk
        dnew0 =-2d0*dt*CH4km
        gnew0 = KapH2Sh
        hnew0 =-dt*H2Sk
        inew0 = KapCH4h
        jnew0 =-dt*CH4km
        knew0 =-dt*kAOM
        lnew0 = KapSO4h
        mnew0 =-dt*kAOM
        nnew0 = dt*kAOM
        onew0 = dt*H2Sk
        CALL JacobM2(anew0,cnew0,dnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0  &
                    ,ChO2,ChH2S,ChSO4,ChCH4,Nj,Kso4_aom,Nloop)
        If(Nloop == maxL) then
          Ierror = 3
          Go to 9998
        End If
      End If
!      ChSO4(2) = SO40
      NH4h = NH4k*ChNH4(2)*ChO2(2)
      H2Sh = H2Sk*ChH2S(2)*ChO2(2)
      CH4h = CH4km*ChCH4(2)*ChO2(2)
      AOMh = kAOM*ChSO4(2)/(ChSO4(2)+Kso4_aom)*ChCH4(2)

!-------!
!  PO4  !
!-------!
      ChPO4(2) = ChPO4(1) + const2*(                         &
                            Vuh*(CdwPO4(1,1)-ChPO4(1))       & ! Advection
                          + const*const1*(CmPO4(1)-ChPO4(1)) & ! Diffuxion
                          - RpcH*Foh                         & ! Export production
                          + Ppo4h*const1)                      ! Riverine input

!-------!
!  NO3  !
!-------!
      ChNO3(2) = ChNO3(1) + const2*(                          &
                            Vuh*CdwNO3(1,1)-Vuh*ChNO3(1)      &
                          + const*const1*(CmNO3(1)-ChNO3(1))  &
                          - FndefH*const1                     &
                          + Pno3h*const1                      &
                          - Rno3detH*(NPPh-Poh)*0.8d0*const1) &
                          + dt*NH4h
      If(ChNO3(2) < 0d0) then
        ChNO3(2) = 0d0
        FnfixH   = -(Vuh*CdwNO3(1,1)-Vuh*ChNO3(1)+const*const1*(CmNO3(1)-ChNO3(1))  &
                    -FndefH*const1+Pno3h*const1-Rno3detH*(NPPh-Poh)*0.8d0*const1+(ChNO3(1)+dt*NH4h)/const2)*Areah
      Else
        FnfixH = 0d0
      End If

!------------------------------!
!  C/N/P stoichiometry of POM  !
!------------------------------!
    If(Redfield == 1)then
   ! Redfield ratio
      RcpH = Rcp0
      RnpH = Rnp0
    Else
   ! Dynamic stoichiometry
      RcpH = Rcp0 + (RcpMax-Rcp0)*0.5d0*(1d0+dtanh((0.1d-3-ChPO4(1))/0.03d-3))
      RnpH = Rnp0 + (RnpMax-Rnp0)*0.5d0*(1d0+dtanh((0.1d-3-ChPO4(1))/0.03d-3))
    End If
    RncH = RnpH/RcpH
    RpcH = 1d0/RcpH
    RocH = (RcpH+2d0*RnpH)/RcpH

!====================!
!  SOFT TISSUE PUMP  !
!====================!
      CALL ORG(Ro2,Rno3,Rso4,CjO2,CjNO3,CjSO4,Kso4,r1o2,r1no3                                &
              ,r1so4,r2o2,r2no3,r2so4,r3o2,r3no3                                             &
              ,r3so4,r1,r2,r3,Ppom1,Ppom2,Ppom3,Area,Fpom1,Fpom2,Fpom3,Vpom,pom              &
              ,dis1,dis2,dis3,Az,TRo2,TRno3,TRso4,dispom,Bpom,Fz,Rpc,Po,N,Depo               &
              ,FreeP,cpr,BpCa,CmO2,BpFe,RETURNP,BurP,TBurP,Trep,TdepoP,DepoDW,freePdw        &
              ,TdepoPdw,Fdwz,cprdw,SdepoP,BpCadw,BpFeDW,BpOrgDW,RETURNPdw,CdwO2,BurPdw       &
              ,TBurPdw,TrepDW,SburP,Sdwb,Foh,Fol,FdwPOM1,FdwPOM2,FdwPOM3,BpOrg               &
              ,r1dw,r2dw,r3dw,Areah,Areal,dt,Ppo4,BCorg,TBCorg,RETURNC,BCorgDW,TBCorgDW      &
              ,RETURNCDW,Sz,CPbur,CPburDW,tb,BE,BEdw,dispomDW,Sdwz,dRo2,dRno3,dRso4,CdwNO3   &
              ,CdwSO4,dis1dw,dis2dw,dis3dw,Temp,TempDW,dispomSED,dispomdwSED,Ro2SED,Rno3SED  &
              ,Rso4SED,dRo2SED,dRno3SED,dRso4SED,Fdeni,FdeniDW,Bpyr,Bdwpyr,TBpyr,CSj         &
              ,CSdwj,fpyrAn,Rch4,dRch4,r1ch4,r1ch4dw,r2ch4,r2ch4dw,r3ch4                     &
              ,r3ch4dw,Rch4SED,dRch4SED,TRch4,fpyrOx,O2js,O2dwjs,AOM,AOMdw,tTRo2,tTRno3      &
              ,tTRso4,tTRch4,epyAnox,Zopd,ZopdDW,Zsmtz,ZsmtzDW,topd,topdDW,tsmtz             &
              ,tsmtzDW,AOMsed,AOMsedDW,ksed,BSRsed,BSRsedj,BSRsedjdw                         &
              ,SR,tAOMsed,Rcp,RcpH,Rnc,Rnp,Roc,BEporg,BEporgDW,CjPO4,CdwPO4                  &
              ,aveBEc,RcpDW,RnpDW,RncDW,RpcDW,RocDW,RnpH,RncH,RpcH,RocH,BEporg0,BEporgDW0    &
              ,BEs,BEdws,Ncount,beoxic,beanox)

!===========================================================================!
!===========================================================================!
!                      EVOLUTION OF DEEP OCEANS                             !
!===========================================================================!
!===========================================================================!
!  Temperature  !
!---------------!
      CALL PHYS(Nj,Temp,Tempdw,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,Tcm,Tch,Sdwz,Areah,Vuh,Sdw)
      Do J = 1, Nj
        Tempc(J) = Temp(2,J)
        Tempk(J) = Temp(2,J) + 273.15d0
      End Do
      
!------------!
!  Salinity  !
!------------!
      CALL PHYS(Nj,Salj,Saljdw,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,Salm,Salh,Sdwz,Areah,Vuh,Sdw)

!-------!
!  14C  !
!-------!
      CALL PHYS(Nj,Cj14C,Cdw14C,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,Cm14C,Ch14C,Sdwz,Areah,Vuh,Sdw)
      Do J = 1, Nj
        Cj14C(2,J) = Cj14C(2,J) - dt*Lam14C*Cj14C(1,J)
        Cdw14C(2,J)= Cdw14C(2,J)- dt*Lam14C*Cdw14C(1,J)
      End Do

!--------------------!
!  Chemical species  !
!--------------------!
      CALL PHYS(Nj,CjNH4,CdwNH4,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmNH4,ChNH4,Sdwz,Areah,Vuh,Sdw) ! NH4
      CALL PHYS(Nj,CjO2,CdwO2,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmO2,ChO2,Sdwz,Areah,Vuh,Sdw)     ! O2
      CALL PHYS(Nj,CjH2S,CdwH2S,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmH2S,ChH2S,Sdwz,Areah,Vuh,Sdw) ! H2S
      CALL PHYS(Nj,CjCH4,CdwCH4,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmCH4,ChCH4,Sdwz,Areah,Vuh,Sdw) ! CH4
      CALL PHYS(Nj,CjSO4,CdwSO4,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmSO4,ChSO4,Sdwz,Areah,Vuh,Sdw) ! SO42-
      const = dt/dz
      const1= 1d0/(Rnc*CNbur)
      const2= 1d0/(RncDW*CNbur)
      const3= dt*FsMOR/(V+Vdw)
      Do J = 1, Nj
        KapNH4(J)   = CjNH4(2,J)  + const*(Rnc*dispom(J)+Rnc*(Fpom1(J+1)+Fpom2(J+1)+Fpom3(J+1)) &
                                             *(1d0-BE(J)*const1)*Fz(J+1)/Sz(J))
        KapNH4DW(J) = CdwNH4(2,J) + const*(RncDW*dispomDW(J)+RncDW*(Fdwpom1(J+1)+Fdwpom2(J+1)+Fdwpom3(J+1)) &
                                             *(1d0-BEdw(J)*const2)*Fdwz(J+1)/Sdwz(J))
        KapO2(J)    = CjO2(2,J)   + const*(-Ro2(J)*dispom(J)-Ro2SED(J)*dispomSED(J))
        KapO2DW(J)  = CdwO2(2,J)  + const*(-dRo2(J)*dispomDW(J)-dRo2SED(J)*dispomDWSED(J))
        KapH2S(J)   = CjH2S(2,J)  + const*(Rso4(J)*dispom(J)*0.5d0+Rso4SED(J)*dispomSED(J)*0.5d0 &
                                            +AOMsed(J)*Fz(J+1)/Sz(J)-Bpyr(J)/Sz(J)-BpyrWC(J)/Sz(J))+const3
        KapH2SDW(J) = CdwH2S(2,J) + const*(dRso4(J)*dispomDW(J)*0.5d0+dRso4SED(J)*dispomDWSED(J)*0.5d0 &
                                            +AOMsedDW(J)*Fdwz(J+1)/Sdwz(J) &
                                            -Bdwpyr(J)/Sdwz(J)-BdwpyrWC(J)/Sdwz(J))+const3
        KapCH4(J)   = CjCH4(2,J)  + const*(Rch4(J)*dispom(J)*0.5d0+Rch4SED(J)*dispomSED(J)*0.5d0 &
                                          -AOMsed(J)*Fz(J+1)/Sz(J))
        KapCH4DW(J) = CdwCH4(2,J) + const*(dRch4(J)*dispomDW(J)*0.5d0+dRch4SED(J)*dispomDWSED(J)*0.5d0 &
                                          -AOMsedDW(J)*Fdwz(J+1)/Sdwz(J))
        KapSO4(J)   = CjSO4(2,J)  + const*(-Rso4(J)*dispom(J)*0.5d0-Rso4SED(J)*dispomSED(J)*0.5d0 &
                                          -AOMsed(J)*Fz(J+1)/Sz(J))
        KapSO4DW(J) = CdwSO4(2,J) + const*(-dRso4(J)*dispomDW(J)*0.5d0-dRso4SED(J)*dispomDWSED(J)*0.5d0 &
                                          -AOMsedDW(J)*Fdwz(J+1)/Sdwz(J))
      End Do

!----------------------------!
!  Secondary redox reaction  !
!----------------------------!
      Do J = 1, Nj
        anew0 = KapO2(J)
        bnew0 =-2d0*dt*NH4k
        cnew0 =-2d0*dt*H2Sk
        dnew0 =-2d0*dt*CH4k
        enew0 = KapNH4(J)
        fnew0 =-dt*NH4k
        gnew0 = KapH2S(J)
        hnew0 =-dt*H2Sk
        inew0 = KapCH4(J)
        jnew0 =-dt*CH4k
        knew0 =-dt*kAOM
        lnew0 = KapSO4(J)
        mnew0 =-dt*kAOM
        nnew0 = dt*kAOM
        onew0 = dt*H2Sk
        CALL Jacob(anew0,bnew0,cnew0,dnew0,enew0,fnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0  &
                  ,CjO2,CjNH4,CjH2S,CjSO4,CjCH4,Nj,J,Kso4_aom,Nloop)
        If(Nloop == maxL) then
          Ierror = 3
          Go to 9998
        End If
        NH4(J) = NH4k*CjNH4(2,J)*CjO2(2,J)
        H2S(J) = H2Sk*CjH2S(2,J)*CjO2(2,J)
        CH4(J) = CH4k*CjCH4(2,J)*CjO2(2,J)
        AOM(J) = kAOM*CjSO4(2,J)/(CjSO4(2,J)+Kso4_aom)*CjCH4(2,J)
        BpyrWC(J) = kpyrite*CjFeII*CjH2S(2,J)*Sz(J)*dz
      End Do
      Do J = 1, Nj
        anew0 = KapO2DW(J)
        bnew0 =-2d0*dt*NH4k
        cnew0 =-2d0*dt*H2Sk
        dnew0 =-2d0*dt*CH4k
        enew0 = KapNH4DW(J)
        fnew0 =-dt*NH4k
        gnew0 = KapH2SDW(J)
        hnew0 =-dt*H2Sk
        inew0 = KapCH4DW(J)
        jnew0 =-dt*CH4k
        knew0 =-dt*kAOM
        lnew0 = KapSO4DW(J)
        mnew0 =-dt*kAOM
        nnew0 = dt*kAOM
        onew0 = dt*H2Sk
        CALL Jacob(anew0,bnew0,cnew0,dnew0,enew0,fnew0,gnew0,hnew0,inew0,jnew0,knew0,lnew0,mnew0,nnew0,onew0  &
                  ,CdwO2,CdwNH4,CdwH2S,CdwSO4,CdwCH4,Nj,J,Kso4_aom,Nloop)
        If(Nloop == maxL) then
          Ierror = 3
          Go to 9998
        End If
        NH4DW(J) = NH4k*CdwNH4(2,J)*CdwO2(2,J)
        H2SDW(J) = H2Sk*CdwH2S(2,J)*CdwO2(2,J)
        CH4DW(J) = CH4k*CdwCH4(2,J)*CdwO2(2,J)
        AOMdw(J) = kAOM*CdwSO4(2,J)/(CdwSO4(2,J)+Kso4_aom)*CdwCH4(2,J)
        BdwpyrWC(J) = kpyrite*CjFeII*CdwH2S(2,J)*Sdwz(J)*dz
      End Do
!
!--------------------------------------------------!
!  Total pyrite precipitation in the water column  !
!--------------------------------------------------!
      TBpyrWC = BpyrWCL + BpyrWCH
      Do J = 1, Nj
        TBpyrWC = TBpyrWC + BpyrWC(J) + BdwpyrWC(J)
      End Do
!
!-------------!
!  Total AOM  !
!-------------!
      AOMsum = Vm*kaom*CmCH4(2)*CmSO4(2)/(CmSO4(2)+Kso4_aom) + Vh*kaom*ChCH4(2)*ChSO4(2)/(ChSO4(2)+Kso4_aom)
      Do J = 1, Nj
        AOMsum = AOMsum + AOM(J)*Sz(J)*dz + AOMdw(J)*Sdwz(J)*dz
      End Do

!---------------!
!  [PO4],[NO3]  !
!---------------!
      CALL PHYS(Nj,CjPO4,CdwPO4,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmPO4,ChPO4,Sdwz,Areah,Vuh,Sdw)
      CALL PHYS(Nj,CjNO3,CdwNO3,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,CmNO3,ChNO3,Sdwz,Areah,Vuh,Sdw)
      const = dt/dz
      Do J = 1, Nj
        CjPO4(2,J) = CjPO4(2,J)   + const*(Rpc*dispom(J)+Rpc*((Fpom1(J+1)+Fpom2(J+1)+Fpom3(J+1)))*Fz(J+1)/Sz(J)  &
                                            -BpOrg(J)/Sz(J)-BpCa(J)/Sz(J)-BpFe(J)/Sz(J))
        CdwPO4(2,J) = CdwPO4(2,J) + const*(RpcDW*dispomDW(J)+RpcDW*((Fdwpom1(J+1)+Fdwpom2(J+1)+Fdwpom3(J+1)))*Fdwz(J+1)/Sdwz(J)  &
                                            -BpOrgDW(J)/Sdwz(J)-BpCaDW(J)/Sdwz(J)-BpFeDW(J)/Sdwz(J))
        CjNO3(2,J)  = CjNO3(2,J)  + const*(-0.8d0*Rno3(J)*dispom(J)-0.8d0*Rno3SED(J)*dispomSED(J))  &
                                  + dt*NH4(J)
        CdwNO3(2,J) = CdwNO3(2,J) + const*(-0.8d0*dRno3(J)*dispomDW(J)-0.8d0*dRno3SED(J)*dispomDWSED(J))  &
                                  + dt*NH4DW(J)
        If(RETURNP(J) < 0) then
          CjPO4(2,J) = CjPO4(2,J) - const*RETURNP(J)/Sz(J)
        End If
        If(CjPO4(2,J) < 0) then
          Write(*,*) 'CjPO4(',J,') is negative!',CjPO4(2,J)
          Ierror = 4
          Go to 9998
        End If
        If(RETURNPDW(J) < 0) then
          Write(*,*) 'RETURNPDW(',J,') is negative !',RETURNPDW(J)
          CdwPO4(2,J) = CdwPO4(2,J) - const*RETURNPDW(J)/Sdwz(J)
        End If
        If(CdwPO4(2,J) < 0) then
          Write(*,*) 'CdwPO4(',J,') is negative!',CdwPO4(2,J)
          Ierror = 4
          Go to 9998
        End If
      End Do

!---------!
!  Check  !
!---------!
      Do J = 1, Nj
        If(CjNO3(2,J) < 0) then
          Write(*,*) 'CjNO3(',J,') is negative !'
          Ierror = 5
          Go to 9998
        End If
        If(CjNH4(2,J) < 0) then
          Write(*,*) 'CjNH4(',J,') is negative !'
          Ierror = 6
          Go to 9998
        End If
        If(CjH2S(2,J) < 0) then
          Write(*,*) 'CjH2S(',J,') is negetive !'
          Ierror = 7
          Go to 9998
        End If
        If(CdwNO3(2,J) < 0) then
          Write(*,*) 'CdwNO3(',J,') is negative !'
          Ierror = 5
          Go to 9998
        End If
        If(CdwNH4(2,J) < 0) then
          Write(*,*) 'CdwNH4(',J,') is negative !'
          Ierror = 6
          Go to 9998
        End If
        If(CdwH2S(2,J) < 0) then
          Write(*,*) 'CdwH2S(',J,') is negetive !'
          Ierror = 7
          Go to 9998
        End If
      End Do

!==================!
!  SULFUR ISOTOPE  !
!==================!
      const = Mso4/Vocean
      Dd34S     = maxDd34S*const/(const+Kmd34S)
      Dd34Sc    = maxDd34S*SO40/(SO40+Kmd34S)
      Dd34S_up  = maxDd34S_up*const/(const+Kmd34S)
      Dd34S_low = maxDd34S_low*const/(const+Kmd34S)

      d34S(2) = d34S(1) + dt*((d34Sin-d34S(1))*FsMOR                       &
                             +(d34Spy(1)-d34S(1))*Wpy  &
                             +(d34Sgyp(1)-d34S(1))*Wgyp                    &
                             +(d34Spy(1)-d34S(1))*FsVol_H2S                &
                             +(d34Sgyp(1)-d34S(1))*FsVol_SO4               &
                             +Dd34S*(TBpyr+TBpyrWC))/Mso4
      d34S_up(2) = d34S_up(1) + dt*((d34Sin-d34S_up(1))*FsMOR                       &
                                   +(d34Spy(1)-d34S_up(1))*Wpy  &
                                   +(d34Sgyp(1)-d34S_up(1))*Wgyp                    &
                                   +(d34Spy(1)-d34S(1))*FsVol_H2S                   &
                                   +(d34Sgyp(1)-d34S(1))*FsVol_SO4                  &
                                   +Dd34S_up*(TBpyr+TBpyrWC))/Mso4
      d34S_low(2)=d34S_low(1) + dt*((d34Sin-d34S_low(1))*FsMOR                       &
                                   +(d34Spy(1)-d34S_low(1))*Wpy  &
                                   +(d34Sgyp(1)-d34S_low(1))*Wgyp                    &
                                   +(d34Spy(1)-d34S(1))*FsVol_H2S                    &
                                   +(d34Sgyp(1)-d34S(1))*FsVol_SO4                   &
                                   +Dd34S_low*(TBpyr+TBpyrWC))/Mso4
      d34S_c(2)=d34S_c(1) + dt*((d34Sin-d34S_c(1))*FsMOR                       &
                               +(d34Spy(1)-d34S_c(1))*Wpy  &
                               +(d34Sgyp(1)-d34S_c(1))*Wgyp                    &
                               +(d34Spy(1)-d34S(1))*FsVol_H2S                  &
                               +(d34Sgyp(1)-d34S(1))*FsVol_SO4                 &
                               +Dd34Sc*(TBpyr+TBpyrWC))/Mso4

      d34Spy(2)  = d34Spy(1) + dt*((d34S(1)-Dd34S-d34Spy(1))*(TBpyr+TBpyrWC))/Spy(1)
      d34Sgyp(2) = d34Sgyp(1) + dt*((d34S(1)-d34Sgyp(1))*TBcaso4)/Sgyp(1)

!=========!
!  CRUST  !
!=========!
      Sgyp(2) = Sgyp(1) + dt*(TBcaso4 - Wgyp - FsVol_SO4)                ! Gypsum S
      Spy(2)  = Spy(1)  + dt*(TBpyr + TBpyrWC - Wpy - FsVol_H2S - FsMOR) ! Pyrite S
      Corg(2) = Corg(1) + dt*(TBCorg + TBCorgDW + Btcorg - Worg - Fmorg) ! Organic C
      If(Sgyp(2) < 0d0) then
        Sgyp(2) = 0d0
      End if
      If(Spy(2) < 0d0) then
        Spy(2) = 0d0
      End if
      If(Corg(2) < 0d0) then
        Corg(2) = 0d0
      End If

!=====================!
!  MASS IN THE OCEAN  !
!=====================!
      mixMo2  = hm*(CmO2(2)*Areal+ChO2(2)*Areah)
      mixMso4 = hm*(CmSO4(2)*Areal+ChSO4(2)*Areah)
      Mpo4 = hm*(CmPO4(2)*Areal+ChPO4(2)*Areah)
      Mno3 = hm*(CmNO3(2)*Areal+ChNO3(2)*Areah)
      Mo2  = mixMo2
      Mh2s = hm*(CmH2S(2)*Areal+ChH2S(2)*Areah)
      Mnh4 = hm*(CmNH4(2)*Areal+ChNH4(2)*Areah)
      Mso4 = mixMso4
      Mch4_ocn = hm*(CmCH4(2)*Areal+ChCH4(2)*Areah)
      Do J = 1, Nj
        Mpo4 = Mpo4 + dz*CjPO4(2,J)*Sz(J) + dz*CdwPO4(2,J)*Sdwz(J)
        Mno3 = Mno3 + dz*CjNO3(2,J)*Sz(J) + dz*CdwNO3(2,J)*Sdwz(J)
        Mo2  = Mo2  + dz*CjO2(2,J)*Sz(J)  + dz*CdwO2(2,J)*Sdwz(J)
        Mh2s = Mh2s + dz*CjH2S(2,J)*Sz(J) + dz*CdwH2S(2,J)*Sdwz(J)
        Mnh4 = Mnh4 + dz*CjNH4(2,J)*Sz(J) + dz*CdwNH4(2,J)*Sdwz(J)
        Mso4 = Mso4 + dz*CjSO4(2,J)*Sz(J) + dz*CdwSO4(2,J)*Sdwz(J)
        Mch4_ocn = Mch4_ocn + dz*CjCH4(2,J)*Sz(J) + dz*CdwCH4(2,J)*Sdwz(J)
      End Do

!--------------------------------------------------------------------!
!                      Output (transition)                           !
!--------------------------------------------------------------------!
       If(((output == 0).and.(Ncount*dt >= 1d3).and.(Ncount*dt < 1d4).and.(Mod(Ncount*dt,1d2) == 0))  &
      .or.((output == 0).and.(Ncount*dt >= 1d4).and.(Ncount*dt < 1d5).and.(Mod(Ncount*dt,1d3)) == 0)  &
      .or.((output == 0).and.(Ncount*dt >= 1d5).and.(Ncount*dt < 1d6).and.(Mod(Ncount*dt,1d4)) == 0)  &
      .or.((output == 0).and.(Ncount*dt >= 1d6).and.(Ncount*dt < 1d7).and.(Mod(Ncount*dt,1d5)) == 0)  &
      .or.((output == 0).and.(Ncount*dt >= 1d7).and.(Ncount*dt < 1d8).and.(Mod(Ncount*dt,1d5)) == 0)  &
      .or.((output == 0).and.(Ncount*dt >= 1d8).and.(Ncount*dt < 1d9).and.(Mod(Ncount*dt,1d6)) == 0)  &
      .or.((output == 0).and.(Ncount*dt >= 1d9).and.(Ncount*dt < 1d10).and.(Mod(Ncount*dt,1d6)) == 0)) then

    ! RECORDING
      Idammy = Idammy + 1
      Time_dammy(Idammy)  = Time
      pO2_dammy(Idammy)   = pO2a(2)
      pCH4_dammy(Idammy)  = pCH4
      CmSO4_dammy(Idammy) = CmSO4(2)
      d34S_dammy(Idammy)  = d34S(2)
      d34Spy_dammy(Idammy)= d34S(2)-Dd34S
      EX_dammy(Idammy)    = (Po+Foh*Areah)*12d0/1d15
      NPP_dammy(Idammy)   = (NPPl+NPPh)*12d0/1d15
      Mp_dammy(Idammy)    = Mpo4/1d15
      Nfix_dammy(Idammy)  = (Fnfix+FnfixH)/1d12*14d0
      Bcorg_dammy(Idammy) = (TBCorg+TBCorgDW)/1d12
      BSRsed_dammy(Idammy)= BSRsed/1d12
      resO2_dammy(Idammy) = tTRo2
      resNO3_dammy(Idammy)= tTRno3
      resSO4_dammy(Idammy)= tTRso4
      resCO2_dammy(Idammy)= tTRch4

      xOUT = Time
      yOUT = Time
!  tBSRwc, tBSRsed, H2Sox, tAOMwc, tAOMsed
      tBSRwc = Rso4det*(NPPl-Po)*0.5d0 + Rso4detH*(NPPh-Poh)*0.5d0
      tH2Sox = Vm*H2Sk*CmO2(2)*CmH2S(2) + Vh*H2Sk*ChO2(2)*ChH2S(2)
      tAOMwc = Vm*kaom*CmCH4(2)*CmSO4(2)/(CmSO4(2)+Kso4_aom) + Vh*kaom*ChCH4(2)*ChSO4(2)/(ChSO4(2)+Kso4_aom)
      tAOMsed = 0d0
      Do J = 1, Nj
        tBSRwc  = tBSRwc  + dispom(J)*Rso4(J)*Sz(J)*0.5d0 + dispomDW(J)*dRso4(J)*Sdwz(J)*0.5d0
        tH2Sox  = tH2Sox  + Sz(J)*dz*H2Sk*CjO2(2,J)*CjH2S(2,J) + Sdwz(J)*dz*H2Sk*CdwO2(2,J)*CdwH2S(2,J)
        tAOMwc  = tAOMwc  + Sz(J)*dz*kaom*CjCH4(2,J)*CjSO4(2,J)/(CjSO4(2,J)+Kso4_aom)  &
                          + Sdwz(J)*dz*kaom*CdwCH4(2,J)*CdwSO4(2,J)/(CdwSO4(2,J)+Kso4_aom)
        tAOMsed = tAOMsed + AOMsed(J)*Fz(J+1) + AOMsedDW(J)*Fdwz(J+1)
      End Do
      jH2S = FsMOR-TBpyr-TBpyrWC-tH2Sox+tAOMwc+tAOMsed-Aso4+BSRsed+tBSRwc
      jSO4 = Pso4-TBcaso4+tH2Sox-tAOMwc-tAOMsed+Aso4-BSRsed+tBSRwc
      tSO4in  = Mso4/(Pso4+Aso4+tH2Sox)
      tSO4out = Mso4/(TBcaso4+BSRsed+tBSRwc+tAOMwc+tAOMsed)
      tH2Sin  = Mh2s/(FsMOR+tAOMwc+tAOMsed+BSRsed+tBSRwc)
      tH2Sout = Mh2s/(TBpyr+TBpyrWC+tH2Sox+Aso4)
    ! MSR_AOM.dat
      Write(20,706) xOUT,yOUT,BSRsed/1d12,tAOMsed/1d12,tBSRwc/1d12,tAOMwc/1d12,AOMsum/1d12
    ! S_balance.dat
      Write(18,710) xOUT,yOUT,jH2S/1d12,jSO4/1d12,tH2Sox/1d12,Aso4/1d12,tSO4in/1d12,tSO4out/1d12,tH2Sin/1d12,tH2Sout/1d12
  
      Do J = 1, Nj
    ! OPD.dat
        Write(35,709) -hm-dz*J,Zopd(J),ZopdDW(J),Zsmtz(J),ZsmtzDW(J),topd(J),topddw(J),tsmtz(J),tsmtzdw(J)
    ! MSRsedj.dat
        Write(21,704) -hm-dz*J,BSRsedj(J),BSRsedjdw(J)
      End Do
      write(35,*) ' '
      write(21,*) ' '
    ! Foa.dat
      Write(23,705) xOUT,yOUT,(Areal*FoaCH4m+Areah*FoaCH4h)/1d12 &
                             ,(Areal*FoaH2Sm+Areah*FoaH2Sh)/1d12 &
                             ,(Areal*FoaNH4m+Areah*FoaNH4h)/1d12
    ! crust.dat
      Write(15,706) xOUT,yOUT,Corg(2)/1d18,Sgyp(2)/1d18,Spy(2)/1d18,(Sgyp(2)+Spy(2))/1d18
    ! Nfix.dat [TgN/yr]
      Write(37,705) xOUT,yOUT,Fnfix/1d12*14d0,FnfixH/1d12*14d0,(Fnfix+FnfixH)/1d12*14d0
! Redox balance
      RedoxBatm = (Areal*FoaO2m+Areah*FoaO2h+sesc*Mch4(1)-Worg-2d0*Wpy  &
                   -2d0*(Areal*FoaCH4m+Areah*FoaCH4h)-2d0*(Areal*FoaH2Sm+Areah*FoaH2Sh)  &
                   -3d0/4d0*(Areal*FoaNH4m+Areah*FoaNH4h)-5d0/4d0*Ano3 &
                   -Fred-2d0*FsVol_H2S-Fmorg+Btcorg)/1d12
      RedoxBocn = (TBcorg+TBcorgDW+2d0*(TBpyr+TBpyrWC)+2d0*(Areal*FoaCH4m+Areah*FoaCH4h)+2d0*(Areal*FoaH2Sm+Areah*FoaH2Sh) &
                   +3d0/4d0*(Areal*FoaNH4m+Areah*FoaNH4h)+5d0/4d0*Ano3-(Areal*FoaO2m+Areah*FoaO2h)-2d0*FsMOR)/1d12
      RedoxB    = (TBcorg+TBCorgDW+2d0*(TBpyr+TBpyrWC)+sesc*Mch4(1)-Worg-2d0*Wpy &
                  -Fred-2d0*FsMOR-2d0*FsVol_H2S-Fmorg+Btcorg)/1d12
!  @atmosphere
    ! RedoxBatm.dat
      Write(16,713) xOUT,yOUT,(Areal*FoaO2m+Areah*FoaO2h)/1d12,sesc*Mch4(1)/1d12,Worg/1d12,2d0*Wpy/1d12  &
                   ,2d0*(Areal*FoaCH4m+Areah*FoaCH4h+Areal*FoaH2Sm+Areah*FoaH2Sh+Areal*FoaNH4m+Areah*FoaNH4h)/1d12 &
                   ,Fred/1d12,2d0*FsVol_H2S/1d12,0d0,Fmorg/1d12,Btcorg/1d12,RedoxBatm
!  @ocean
    ! RedoxBocn.dat
      Write(17,708) xOUT,yOUT,(TBcorg+TBcorgDW)/1d12,2d0*(TBpyr+TBpyrWC)/1d12  &
                   ,2d0*(Areal*FoaCH4m+Areah*FoaCH4h+Areal*FoaH2Sm+Areah*FoaH2Sh+Areal*FoaNH4m+Areah*FoaNH4h)/1d12  &
                   ,(Areal*FoaO2m+Areah*FoaO2h)/1d12,RedoxBocn,RedoxB

      d34S_calc      = d34Sin + Dd34S*(TBpyr+TBpyrWC)/(TBpyr+TBpyrWC+TBcaso4)
      d34Spy_calc    = d34S_calc - Dd34S
      d34S_calc_up   = d34Sin + Dd34S_up*(TBpyr+TBpyrWC)/(TBpyr+TBpyrWC+TBcaso4)
      d34Spy_calc_up = d34S_calc_up - Dd34S_up
      d34S_calc_low  = d34Sin + Dd34S_low*(TBpyr+TBpyrWC)/(TBpyr+TBpyrWC+TBcaso4)
      d34Spy_calc_low= d34S_calc_low - Dd34S_low
    ! d34S.dat
      Write(40,718) xOUT,yOUT,d34S(2),d34Spy(2),d34Sgyp(2),d34S_up(2),d34S_low(2),d34S_c(2)         &
                   ,d34S(2)-Dd34S,d34S_up(2)-Dd34S_up,d34S_low(2)-Dd34S_low,d34S_c(2)-Dd34Sc        &
                   ,d34S_calc,d34Spy_calc,d34S_calc_up,d34Spy_calc_up,d34S_calc_low,d34Spy_calc_low
    ! all.dat
      Write(2,729) xOUT,pO2a(1)/pO2a0,epyAnox,Kso4,Spy(1)/1d18,Sgyp(1)/1d18,Ppo4/Ppo40,Vpom,Rcp,fsr_rnd(Imc),fBIFP,BpBIFP/1d12 &
                  ,abs(1d0-(TBcaso4+TBpyr+TBpyrWC)/(Pso4+FsMOR)),1d0-(SburP+BpBIFP)/Ppo4,pCH4*1d6 &
                  ,Wpy/1d12,Wgyp/1d12,fexportL &
                  ,CmSO4(2),ChSO4(2),Mpo4/1d15,(Po+Foh*Areah)*12d0/1d15,(NPPl+NPPh)*12d0/1d15,(Fnfix+FnfixH)/1d12*14d0 &
                  ,(TBCorg+TBCorgDW)/1d12,Corg(1)/1d18,RedoxBatm,RedoxBocn,RedoxB
    ! BEj.dat
      Do J = 1, Nj
         write(30,707) xOUT,yOUT,-(J+1d0)*dz,BE(J),BEdw(J),BEporg(J),BEporgDW(J)
      End do
    ! land.dat
      write(49,713) xOUT,yOUT,Vege,NPPt/1d12,Wp/1d12,Pland/1d12,Ppo4/1d12,k11,Jch4t/1d12,gamma_t,delta_t,go2_t,gch4_t

      CALL OUTPT(xOUT,yOUT,Worg  &
                ,Tcm,Tch,Sv,Ppo4,Ppo40,AreaCS,AreaCS0,pO2a,pO2a0,dSST,dSSS,Vpom,Pno3,Ano3  &
                ,Po,Foh,Areah,ChPO4,ChNO3,ChO2,ChSO4  &
                ,Mh2s,Mo2,Mpo4,CjO2,CjH2S,hm,Temp,TempDW,Sm,Sh,CmPO4,CmNO3,CmO2,CmSO4,CmNH4,CmH2S  &
                ,Dm14C,Cm14C,Dh14C,Ch14C,dz,CjPO4,CdwPO4,CjNO3,CdwNO3  &
                ,CdwO2,CjSO4,CdwSO4,CjNH4,CdwNH4,CdwH2S,Dj14C,Cj14C,Ddw14C,Cdw14C  &
                ,gO2j,gPO4j,gNO3j,gH2Sj,gNH4j,gSO4j,Dg14Cj,Sz,Sdwz,Fpom1,Fpom2,Fpom3,Fdwpom1,Fdwpom2,Fdwpom3  &
                ,Fz,Fdwz,Mno3,Mso4  &
                ,Mnh4,TBCorg,TBCorgDW,BCorg,BCorgDW,RETURNC,RETURNCdw,RETURNP,BurP,SburP,cpr,BpOrg,BpCa,BpFe  &
                ,BpOrgDW,BpCaDW,BpFeDW,Torgp,Tcap,Tfep,RETURNPdw,Ia,Ie,Iredox,Nj,Vj,Vdwj,RI,RIdw,gRI,Pno3M  &
                ,ChNH4,ChH2S,dispomSED,dispomdwSED,Ro2SED,Rno3SED,Rso4SED  &
                ,dRo2SED,dRno3SED,dRso4SED,Fdeni,FdeniDW,Tdeni,TdeniCS,Salj,SaljDW  &
                ,TBcaso4,TBpyr,CmCH4,ChCH4,CjCH4,CdwCH4,AOM,AOMdw,TRo2,TRno3,TRso4,TRch4  &
                ,tTRo2,tTRno3,tTRso4,tTRch4,pCH4,TBpyrWC,NPPl,NPPh,Vm,Vh,Salm,Salh,Mch4_ocn)
      End If
!--------------------------------------------------------------------!
!                  End of output (transition)                        !
!--------------------------------------------------------------------!

!-----------!
!  Recycle  !
!-----------!
      MaO2(1)    = MaO2(2)
      pO2a(1)    = pO2a(2)
      Mch4(1)    = Mch4(2)
      Ts(1)      = Ts(2)
      Tm(1)      = Tm(2)
      Th(1)      = Th(2)
      Tcs(1)     = Tcs(2)
      Tcm(1)     = Tcm(2)
      Tch(1)     = Tch(2)
      Salm(1)    = Salm(2)
      Salh(1)    = Salh(2)
      d34S(1)    = d34S(2)
      d34S_up(1) = d34S_up(2)
      d34S_low(1)= d34S_low(2)
      d34S_c(1)  = d34S_c(2)
      Corg(1)    = Corg(2)
      Spy(1)     = Spy(2)
      Sgyp(1)    = Sgyp(2)
      d34Sgyp(1) = d34Sgyp(2)
      d34Spy(1)  = d34Spy(2)
! [X]m
      CmPO4(1)  = CmPO4(2)
      CmNO3(1)  = CmNO3(2)
      CmSO4(1)  = CmSO4(2)
      CmO2(1)   = CmO2(2)
      CmNH4(1)  = CmNH4(2)
      CmCH4(1)  = CmCH4(2)
      CmH2S(1)  = CmH2S(2)
! [X]h
      ChPO4(1)  = ChPO4(2)
      ChNO3(1)  = ChNO3(2)
      ChNH4(1)  = ChNH4(2)
      ChSO4(1)  = ChSO4(2)
      ChH2S(1)  = ChH2S(2)
      ChCH4(1)  = ChCH4(2)
      ChO2(1)   = ChO2(2)
! [X]j
      Do J = 1, Nj
        CjPO4(1,J)  = CjPO4(2,J)
        CjNO3(1,J)  = CjNO3(2,J)
        CjNH4(1,J)  = CjNH4(2,J)
        CjSO4(1,J)  = CjSO4(2,J)
        CjH2S(1,J)  = CjH2S(2,J)
        CjCH4(1,J)  = CjCH4(2,J)
        CjO2(1,J)   = CjO2(2,J)
        CdwPO4(1,J) = CdwPO4(2,J)
        CdwNO3(1,J) = CdwNO3(2,J)
        CdwNH4(1,J) = CdwNH4(2,J)
        CdwSO4(1,J) = CdwSO4(2,J)
        CdwH2S(1,J) = CdwH2S(2,J)
        CdwCH4(1,J) = CdwCH4(2,J)
        CdwO2(1,J)  = CdwO2(2,J)
         Temp(1,J)  =  Temp(2,J)
        Cj14C(1,J)  = Cj14C(2,J)
         Salj(1,J)  = Salj(2,J)
         TempDW(1,J)= TempDW(2,J)
         Cdw14C(1,J)= Cdw14C(2,J)
         SaljDW(1,J)= SaljDW(2,J)
      End Do

!     END TIME EVOLUTION
      End Do

!==================================================================================================================================!
 9999 Continue
      Ic_success = Ic_success + 1
      Write(*,'(f10.2,A17,f8.4,f8.4,A8)') Ncount*dt,'yr. Steady State ',Ppo4/1d12,SburP/1d12,' TmolP/yr'
      Write(*,'(A24,f9.3,A3)') 'Residense time of DIP = ',Mpo4/Ppo4,' yr'
      Write(*,*) 'C/S=',(TBcorg+TBcorgDW)/(TBpyr+TBpyrWC)
    ! success_story.dat
      Do I = 1, Idammy
         Write(4,726) Time_dammy(I),pO2_dammy(I),epyAnox,Kso4,Spy(1)/1d18,Sgyp(1)/1d18,Ppo4/Ppo40,Vpom,Rcp   &
                      ,fsr_rnd(Imc),fBIFP,fexportL,pCH4_dammy(I),CmSO4_dammy(I),GRBg_dammy(I),GRBatm_dammy(I) &
                      ,GRBocn_dammy(I),d34S_dammy(I),d34Spy_dammy(I),EX_dammy(I),NPP_dammy(I),Mp_dammy(I)     &
                      ,Nfix_dammy(I),Bcorg_dammy(I),BSRsed_dammy(I),Corg(1)/1d18
      End Do
      Write(4,*) ' '
!--------------------------------------------------------------------!
!                      Output (Steady State)                         !
!--------------------------------------------------------------------!
      If(output == 1) then
      xOUT = pO2a(1)
      yOUT = Ppo4/Ppo40
! Redox balance
      RedoxBatm = (Areal*FoaO2m+Areah*FoaO2h+sesc*Mch4(1)-Worg-2d0*Wpy  &
                   -2d0*(Areal*FoaCH4m+Areah*FoaCH4h)-2d0*(Areal*FoaH2Sm+Areah*FoaH2Sh)  &
                   -3d0/4d0*(Areal*FoaNH4m+Areah*FoaNH4h)-5d0/4d0*Ano3 &
                   -Fred-2d0*FsVol_H2S-Fmorg+Btcorg)/1d12
      RedoxBocn = (TBcorg+TBcorgDW+2d0*(TBpyr+TBpyrWC)+2d0*(Areal*FoaCH4m+Areah*FoaCH4h)+2d0*(Areal*FoaH2Sm+Areah*FoaH2Sh) &
                   +3d0/4d0*(Areal*FoaNH4m+Areah*FoaNH4h)+5d0/4d0*Ano3-(Areal*FoaO2m+Areah*FoaO2h)-2d0*FsMOR)/1d12
      RedoxB    = (TBcorg+TBCorgDW+2d0*(TBpyr+TBpyrWC)+sesc*Mch4(1)-Worg-2d0*Wpy &
                  -Fred-2d0*FsMOR-2d0*FsVol_H2S-Fmorg+Btcorg)/1d12
!  @atmosphere
    ! RedoxBatm.dat
      Write(16,713) xOUT,yOUT,(Areal*FoaO2m+Areah*FoaO2h)/1d12,sesc*Mch4(1)/1d12,Worg/1d12,2d0*Wpy/1d12  &
                   ,2d0*(Areal*FoaCH4m+Areah*FoaCH4h+Areal*FoaH2Sm+Areah*FoaH2Sh+Areal*FoaNH4m+Areah*FoaNH4h)/1d12 &
                   ,Fred/1d12,2d0*FsVol_H2S/1d12,0d0,Fmorg/1d12,Btcorg/1d12,RedoxBatm
!  @ocean
    ! RedoxBocn.dat
      Write(17,709) xOUT,yOUT,(TBcorg+TBcorgDW)/1d12,2d0*TBpyr/1d12,2d0*TBpyrWC/1d12  &
                   ,2d0*(Areal*FoaCH4m+Areah*FoaCH4h+Areal*FoaH2Sm+Areah*FoaH2Sh+Areal*FoaNH4m+Areah*FoaNH4h)/1d12  &
                   ,(Areal*FoaO2m+Areah*FoaO2h)/1d12,RedoxBocn,RedoxB
      Do J = 1, Nj
      ! OPD.dat
        Write(35,709) -hm-dz*J,Zopd(J),ZopdDW(J),Zsmtz(J),ZsmtzDW(J),topd(J),topddw(J),tsmtz(J),tsmtzdw(J)
      ! MSRsedj.dat
        Write(21,704) -hm-dz*J,BSRsedj(J),BSRsedjdw(J)
      End Do
    ! Foa.dat
      Write(23,705) xOUT,yOUT,(Areal*FoaCH4m+Areah*FoaCH4h)/1d12  &
                             ,(Areal*FoaH2Sm+Areah*FoaH2Sh)/1d12  &
                             ,(Areal*FoaNH4m+Areah*FoaNH4h)/1d12
    ! crust.dat
      Write(15,706) xOUT,yOUT,Corg(2)/1d18,Sgyp(2)/1d18,Spy(2)/1d18,(Sgyp(2)+Spy(2))/1d18
    ! Nfix.dat [TgN/yr]
      Write(37,*) xOUT,yOUT,Fnfix/1d12*14d0,FnfixH/1d12*14d0,(Fnfix+FnfixH)/1d12*14d0
! S balance
!  tBSRwc, tBSRsed, H2Sox, tAOMwc, tAOMsed
      tBSRwc = Rso4det*(NPPl-Po)*0.5d0 + Rso4detH*(NPPh-Poh)*0.5d0
      tH2Sox = Vm*H2Sk*CmO2(2)*CmH2S(2)+Vh*H2Sk*ChO2(2)*ChH2S(2)
      tAOMwc = Vm*kaom*CmCH4(2)*CmSO4(2)/(CmSO4(2)+Kso4_aom)  &
              + Vh*kaom*ChCH4(2)*ChSO4(2)/(ChSO4(2)+Kso4_aom)
      tAOMsed = 0d0
      Do J = 1, Nj
        tBSRwc  = tBSRwc  + dispom(J)*Rso4(J)*Sz(J)*0.5d0 + dispomDW(J)*dRso4(J)*Sdwz(J)*0.5d0
        tH2Sox  = tH2Sox  + Sz(J)*dz*H2Sk*CjO2(2,J)*CjH2S(2,J)  &
                          + Sdwz(J)*dz*H2Sk*CdwO2(2,J)*CdwH2S(2,J)
        tAOMwc  = tAOMwc  + Sz(J)*dz*kaom*CjCH4(2,J)*CjSO4(2,J)/(CjSO4(2,J)+Kso4_aom)  &
                          + Sdwz(J)*dz*kaom*CdwCH4(2,J)*CdwSO4(2,J)/(CdwSO4(2,J)+Kso4_aom)
        tAOMsed = tAOMsed + AOMsed(J)*Fz(J+1) + AOMsedDW(J)*Fdwz(J+1)
      End Do
      jH2S    = FsMOR-TBpyr-TBpyrWC-tH2Sox+tAOMwc+tAOMsed-Aso4+BSRsed+tBSRwc
      jSO4    = Pso4-TBcaso4+tH2Sox-tAOMwc-tAOMsed+Aso4-BSRsed+tBSRwc
      tSO4in  = Mso4/(Pso4+Aso4+tH2Sox)
      tSO4out = Mso4/(TBcaso4+BSRsed+tBSRwc+tAOMwc+tAOMsed)
      tH2Sin  = Mh2s/(FsMOR+tAOMwc+tAOMsed+BSRsed+tBSRwc)
      tH2Sout = Mh2s/(TBpyr+TBpyrWC+tH2Sox+Aso4)
    ! MSR_AOM.dat
      Write(20,707) xOUT,yOUT,BSRsed/1d12,tAOMsed/1d12,tBSRwc/1d12,tAOMwc/1d12,AOMsum/1d12
    ! S_balance.dat
      Write(18,710) xOUT,yOUT,jH2S/1d12,jSO4/1d12,tH2Sox/1d12,Aso4/1d12,tSO4in/1d12,tSO4out/1d12,tH2Sin/1d12,tH2Sout/1d12

      d34S_calc      = d34Sin + Dd34S*(TBpyr+TBpyrWC)/(TBpyr+TBpyrWC+TBcaso4)
      d34Spy_calc    = d34S_calc - Dd34S
      d34S_calc_up   = d34Sin + Dd34S_up*(TBpyr+TBpyrWC)/(TBpyr+TBpyrWC+TBcaso4)
      d34Spy_calc_up = d34S_calc_up - Dd34S_up
      d34S_calc_low  = d34Sin + Dd34S_low*(TBpyr+TBpyrWC)/(TBpyr+TBpyrWC+TBcaso4)
      d34Spy_calc_low= d34S_calc_low - Dd34S_low
    ! d34S.dat
      Write(40,716) xOUT,yOUT,d34S(2),d34S_up(2),d34S_low(2),d34S_c(2)                               &
                   ,d34S(2)-Dd34S,d34S_up(2)-Dd34S_up,d34S_low(2)-Dd34S_low,d34S_c(2)-Dd34Sc         &
                   ,d34S_calc,d34Spy_calc,d34S_calc_up,d34Spy_calc_up,d34S_calc_low,d34Spy_calc_low
    ! all.dat
      Write(2,729) pO2a(1)/pO2a0,epyAnox,Kso4,Spy(1)/1d18,Sgyp(1)/1d18,Ppo4/Ppo40,Vpom,Rcp,fsr_rnd(Imc),fBIFP,BpBIFP/1d12 &
                  ,abs(1d0-(TBcaso4+TBpyr+TBpyrWC)/(Pso4+FsMOR)),1d0-(SburP+BpBIFP)/Ppo4,pCH4*1d6 &
                  ,Wpy/1d12,Wgyp/1d12,fexportL &
                  ,CmSO4(2),ChSO4(2),Mpo4/1d15,(Po+Foh*Areah)*12d0/1d15,(NPPl+NPPh)*12d0/1d15,(Fnfix+FnfixH)/1d12*14d0 &
                  ,(TBCorg+TBCorgDW)/1d12,TBcaso4/1d12,0d0,RedoxBatm,RedoxBocn,RedoxB
    ! BEj.dat
      Do J = 1, Nj
         write(30,707) xOUT,yOUT,-(J+1d0)*dz,BE(J),BEdw(J),BEporg(J),BEporgDW(J)
      End do

      CALL OUTPT(xOUT,yOUT,Worg  &
                ,Tcm,Tch,Sv,Ppo4,Ppo40,AreaCS,AreaCS0,pO2a,pO2a0,dSST,dSSS,Vpom,Pno3,Ano3  &
                ,Po,Foh,Areah,ChPO4,ChNO3,ChO2,ChSO4  &
                ,Mh2s,Mo2,Mpo4,CjO2,CjH2S,hm,Temp,TempDW,Sm,Sh,CmPO4,CmNO3,CmO2,CmSO4,CmNH4,CmH2S  &
                ,Dm14C,Cm14C,Dh14C,Ch14C,dz,CjPO4,CdwPO4,CjNO3,CdwNO3  &
                ,CdwO2,CjSO4,CdwSO4,CjNH4,CdwNH4,CdwH2S,Dj14C,Cj14C,Ddw14C,Cdw14C  &
                ,gO2j,gPO4j,gNO3j,gH2Sj,gNH4j,gSO4j,Dg14Cj,Sz,Sdwz,Fpom1,Fpom2,Fpom3,Fdwpom1,Fdwpom2,Fdwpom3  &
                ,Fz,Fdwz,Mno3,Mso4  &
                ,Mnh4,TBCorg,TBCorgDW,BCorg,BCorgDW,RETURNC,RETURNCdw,RETURNP,BurP,SburP,cpr,BpOrg,BpCa,BpFe  &
                ,BpOrgDW,BpCaDW,BpFeDW,Torgp,Tcap,Tfep,RETURNPdw,Ia,Ie,Iredox,Nj,Vj,Vdwj,RI,RIdw,gRI,Pno3M  &
                ,ChNH4,ChH2S,dispomSED,dispomdwSED,Ro2SED,Rno3SED,Rso4SED  &
                ,dRo2SED,dRno3SED,dRso4SED,Fdeni,FdeniDW,Tdeni,TdeniCS,Salj,SaljDW  &
                ,TBcaso4,TBpyr,CmCH4,ChCH4,CjCH4,CdwCH4,AOM,AOMdw,TRo2,TRno3,TRso4,TRch4  &
                ,tTRo2,tTRno3,tTRso4,tTRch4,pCH4,TBpyrWC,NPPl,NPPh,Vm,Vh,Salm,Salh,Mch4_ocn)
      End If
!-------------------------------------------------------------------!
!                 End of OutPut (Steady State)                      !
!-------------------------------------------------------------------!

9998  Continue

      Ic_false = Ic_false + 1
      If(Imc == resample) then
      ! file_count.dat
        write(99,*) Ic_success,Ic_false
      End If
      ! false.dat
        Write(3,917) Ierror,Time,pO2a(2)/pO2a0,epyAnox,Kso4,Spy(1)/1d18,Sgyp(1)/1d18,Ppo4/Ppo40,Vpom,Rcp,fsr_rnd(Imc),fBIFP,BpBIFP &
                    ,abs(1d0-(TBcaso4+TBpyr+TBpyrWC)/(Pso4+FsMOR)),1d0-(SburP+BpBIFP)/Ppo4,pCH4,fexportL
      Do I = 7, 41
         Write(I,*) ' '
      End Do

      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do

!=================!
!  CLOSING FILES  !
!=================!
      Do J = 1, 99
        CLOSE(J)
      End Do

!====================!
!  FORMAT STATEMENT  !
!====================!
      703   Format(3E13.4E2)
      704   Format(4E13.4E2)
      705   Format(5E13.4E2)
      706   Format(6E13.4E2)
      707   Format(7E13.4E2)
      708   Format(8E13.4E2)
      709   Format(9E13.4E2)
      710   Format(10E13.4E2)
      711   Format(11E13.4E2)
      712   Format(12E13.4E2)
      713   Format(13E13.4E2)
      714   Format(14E13.4E2)
      715   Format(15E13.4E2)
      716   Format(16E13.4E2)
      717   Format(17E13.4E2)
      718   Format(18E13.4E2)
      719   Format(19E13.4E2)
      720   Format(20E13.4E2)
      721   Format(21E13.4E2)
      722   Format(22E13.4E2)
      723   Format(23E13.4E2)
      724   Format(24E13.4E2)
      726   Format(26E13.4E2)           
      729   Format(29E13.4E2)

      916   Format(I3,16E13.4E2)
      917   Format(I3,17E13.4E2)

!=============!
!  FUNCTIONS  !
!=============!
      CONTAINS
!  lm
      FUNCTION lm(Temp,Sal)
        USE Constants
        Implicit none
        double precision:: lm,lTm
        double precision, Intent(IN) :: Temp(2),Sal
        lTm = dlog((298.15d0-Temp(1))/(273.15d0+Temp(1)))
        lm  = A0 + A1*lTm + A2*lTm*lTm + A3*lTm**3d0 + A4*lTm**4d0 + A5*lTm**5d0  &
             + Sal*(B0+B1*lTm+B2*lTm**2d0+B3*lTm**3d0) + C0*Sal**2d0
      END FUNCTION lm
      
!  O2-dependency of abiotic pyrite weathering
      FUNCTION fO2Wpy(pO2)
        USE Constants
        Implicit none
        double precision:: fO2Wpy
        double precision, Intent(IN) :: pO2(2)
        fO2Wpy = (1d0+KpyrO2*fsr_rnd(Imc)) * (pO2a(1)/pO2a0)/(pO2a(1)/pO2a0+KpyrO2*fsr_rnd(Imc))
      END FUNCTION fO2Wpy

!  Reaction rate of CH4 photolysis
      FUNCTION PSI(pO2)
        USE Constants
        Implicit none
        double precision:: PSI,psiO2
        double precision, Intent(IN):: pO2(2)
        psiO2   = dlog10(pO2(1)*(143.818d0+38d0)*1d18)
        PSI     = 10d0**(a1ch4*psiO2**4d0+a2ch4*psiO2**3d0+a3ch4*psiO2**2d0+a4ch4*psiO2+a5ch4)
      END FUNCTION

! ignition factor
      FUNCTION ignit(mO2)
        USE Constants
        double precision:: ignit
        double precision, Intent(IN) :: mO2(2)
        ignit = 48d0*mO2(1)-9.08d0             ! Ignition factor (Lenton, 2013)
        If(ignit < 0d0) then
          ignit = 0d0
        End If
        If(ignit > 5d0) then
          ignit = 5d0
        End If
      End FUNCTION ignit

  End Program CANOPS_OPENO2

!_____________________________________________________________________________________________________!
!_____________________________________________________________________________________________________!
!                                                                                                     !
!   SUBROUTINES                                                                                       !
!_____________________________________________________________________________________________________!
!_____________________________________________________________________________________________________!

    SUBROUTINE SETUP(dt,pCH4,Mch40,Mch4,Rcp,Rnp,Rnc,Rpc,Roc,RcpH,RnpH,RncH,RpcH,RocH &
                    ,RcpDW,RnpDW,RncDW,RpcDW,RocDW,Sgyp,Spy,d34S,d34S_up,d34S_low,d34S_c,d34Sgyp     &
                    ,FsMOR,FsVol,Salm,Salh,Wgyp,WpyBio0,WpyAbio0,Wpy,WpyBio,WpyAbio,TBcaso4,Pso4     &
                    ,Ts,Tcs,Tcm,Tm,Tch,Th,Ppo4)
    USE Constants
    double precision dt,pCH4,Mch40,Mch4(2),Rcp,Rnp,Rnc,Rpc,Roc,RcpH,RnpH,RncH,RpcH,RocH &
                    ,RcpDW,RnpDW,RncDW,RpcDW,RocDW,Sgyp(2),Spy(2),d34S(2),d34S_up(2),d34S_low(2)   &
                    ,d34S_c(2),d34Sgyp(2),FsMOR,FsVol,Wgyp,WpyBio0,WpyAbio0,Wpy,WpyBio,WpyAbio     &
                    ,TBcaso4,Pso4,Ts(2),Tcs(2),Tcm(2),Tm(2),Tch(2),Th(2),Salm(2),Salh(2),Ppo4
                    
! Time step [yr]
    dt = 0.02d0

! Atmospheric composition (CH4,NH3,H2S)
    pCH4    = 10d-6
    Mch40   = pCH4*(143.818d0+38d0)*1d18
    Mch4(1) = 0d0
      
! Redfield ratio [mol/mol]
    Rcp   = Rcp0
    Rnp   = Rnp0
    Rnc   = Rnp/Rcp
    Rpc   = 1d0/Rcp
    Roc   = (Rcp+2d0*Rnp)/Rcp
    RcpH  = Rcp0
    RnpH  = Rnp0
    RncH  = RnpH/RcpH
    RpcH  = 1d0/RcpH
    RocH  = (RcpH+2d0*RnpH)/RcpH
    RcpDW = Rcp0
    RnpDW = Rnp0
    RncDW = RnpH/RcpH
    RpcDW = 1d0/RcpH
    RocDW = (RcpH+2d0*RnpH)/RcpH
      
! Crustal S reservoir [molS]
    Sgyp(1) = 200d18
    Spy(1)  = 200d18
      
! Sulfur isotope [permil]
    d34S(1)     = 10d0
    d34S_up(1)  = 23.2d0
    d34S_low(1) = 6.96d0
    d34S_c(1)   = d34S(1)
    d34Sgyp(1)  = d34S(1)
      
! Outgassing flux of S [molS/yr]
    FsMOR = FsMOR0 * 1d0
    FsVol = FsVol0 * 1d0
      
! Gypsum weathering [molS/yr]
    Wgyp = Wgyp0 * Sgyp(1)/Sgyp0
      
! Pyrite weathering [molS/yr]
    WpyBio0  = Wpy0 * 0.4d0
    WpyAbio0 = Wpy0 * 0.6d0
    Wpy      = Wpy0     * Spy(1)/Spy0
    WpyBio   = WpyBio0  * Spy(1)/Spy0
    WpyAbio  = WpyAbio0 * Spy(1)/Spy0
      
! Gypsum deposition rate [molS/yr]
    TBcaso4 = TBcaso40
      
! Pso4
    Pso4 = Pso40
      
! Temprature
    Ts(1)   = Ts0
    Tcs(1)  = Tcs0
    Tcm(1)  = Tcm0
    Tm(1)   = Tcm(1) + 273.15d0
    Tch(1)  = Tch0
    Th(1)   = Tch(1) + 273.15d0
      
! Salinity
    Salm(1) = Sm0
    Salm(2) = Sm0
    Salh(1) = Sh0
    Salh(2) = Sh0
      
! Riverine P input rate [molP/yr]
    Ppo4 = Ppo40

    RETURN
    End SUBROUTINE SETUP

!_____________________________________________________________________________________________________!
!_____________________________________________________________________________________________________!

    SUBROUTINE READ(PO4j_start,NO3j_start,NH4j_start,H2Sj_start,CH4j_start,O2j_start,SO4Ls,SO4Hs   &
                   ,SO4j_start,Ages,fR,fE,fAD,fLA,fD,fmor,Geog,fA,fL,fCa,facs,fpo2a,tb,Tk,Tc,S,P   &
                   ,BE,Kz,SR0,CjPO4,CdwPO4,CjNO3,CdwNO3,CjO2,CdwO2,CjSO4,CdwSO4,CjNH4,CdwNH4       &
                   ,CjH2S,CdwH2S,Temp,Tempdw,Dj14C,Ddw14C,Cj14C,Cdw14C,TempC,Tempk,TempDWc,TempDWk &
                   ,yz0,Kzh,Ku,Kl,Nj,Nage,CjCH4,CdwCH4,O2js,O2dwjs,Tempj_start,Salj_start,C14j_start)
    Integer Nj,Nage,J
    double precision PO4j_start(Nj+1,2),NO3j_start(Nj+1,2),NH4j_start(Nj+1,2),H2Sj_start(Nj+1,2) &
                    ,CH4j_start(Nj+1,2),O2j_start(Nj+1,2),SO4Ls(28),SO4Hs(28),SO4j_start(Nj+1,2)    &
                    ,Ages(Nage),fR(Nage),fE(Nage),fAD(Nage),fLA(Nage),fD(Nage),fmor(Nage),Geog(Nage)         &
                    ,fA(Nage),fL(Nage),fCa(Nage),facs(Nage),fpo2a(Nage),tb(Nj),Tk(Nj),Tc(Nj),S(Nj),P(Nj)     &
                    ,BE(Nj),Kz(Nj),SR0(Nj)                                                                   &
                    ,CjPO4(2,Nj),CdwPO4(2,Nj),CjNO3(2,Nj),CdwNO3(2,Nj),CjO2(2,Nj),CdwO2(2,Nj)                &
                    ,CjSO4(2,Nj),CdwSO4(2,Nj),CjNH4(2,Nj),CdwNH4(2,Nj),CjH2S(2,Nj),CdwH2S(2,Nj)              &
                    ,Temp(2,Nj),Tempdw(2,Nj),Dj14C(Nj),Ddw14C(Nj),Cj14C(2,Nj),Cdw14C(2,Nj)                   &
                    ,TempC(Nj),Tempk(Nj),TempDWc(Nj),TempDWk(Nj),yz0(Nj+1),Kzh(Nj),Ku,Kl                     &
                    ,dammy,CjCH4(2,Nj),CdwCh4(2,Nj),O2js(Nj+1),O2dwjs(Nj+1) &
                    ,Tempj_start(Nj+1,2),Salj_start(Nj+1,2),C14j_start(Nj+1,2)

    OPEN(1, File = './input/Xj_reference.dat')
    Do J = 1, Nj+1
      Read(1,*)  dammy,PO4j_start(J,1),PO4j_start(J,2),NO3j_start(J,1),NO3j_start(J,2),O2j_start(J,1) ,O2j_start(J,2) &
                ,SO4j_start(J,1),SO4j_start(J,2),NH4j_start(J,1),NH4j_start(J,2),H2Sj_start(J,1),H2Sj_start(J,2) &
                ,CH4j_start(J,1),CH4j_start(J,2),Tempj_start(J,1),Tempj_start(J,2),Salj_start(J,1),Salj_start(J,2) &
                ,C14j_start(J,1),C14j_start(J,2)
      O2js(J)   = O2j_start(J,1)
      O2dwjs(J) = O2j_start(J,2)
    End Do
    Close(1)

    OPEN(1, File = './input/f_GEOCARB1.dat')
    OPEN(2, File = './input/f_GEOCARB2.dat')
    OPEN(3, File = './input/fCa.dat')
    OPEN(4, File = './input/Acs_pO2.dat')
    Do J = 1, Nage
      Read(1,*) Ages(J),fR(J),fE(J),fAD(J),fLA(J)
      Read(2,*) fD(J),fmor(J),Geog(J),fA(J),fL(J)
      Read(3,*) dammy,fCa(J)
      Read(4,*) dammy,facs(J),fpo2a(J)
    End Do
    Do I = 1, 4
      Close(I)
    End Do

    OPEN( 1, File = './input/tb.dat')
    OPEN( 2, File = './input/tk.dat')
    OPEN( 3, File = './input/s.dat')
    OPEN( 4, File = './input/p.dat')
    OPEN( 7, File = './input/y.dat')
    OPEN( 8, File = './input/be.dat')
    OPEN( 9, File = './input/kz.dat')
    OPEN(10, File = './input/sr.dat')
    Do J = 1, Nj
      Read( 1,*) tb(J)
      Read( 2,*) Tk(J)
                 Tc(J)=Tk(J)-273.15d0
      Read( 3,*)  S(J)
      Read( 4,*)  P(J)
      Read( 8,*) BE(J)
      Read( 9,*) Kz(J)
      Read(10,*) SR0(J)
    End Do
! Seabed topography
    Do J = 1, Nj+1
      Read(7,*) yz0(J)
    End Do

    Do I = 1, 4
      Close(I)
    End Do
    Do I = 7, 10
      Close(I)
    End Do

    RETURN
    End SUBROUTINE READ

!_____________________________________________________________________________________________________!
!_____________________________________________________________________________________________________!

    SUBROUTINE STRUCTURE(Az,AzDW,yz,yz0,Nj,Sdw,hm,dz,Sz,Sdwz,Area,Areal,Areah,Fz,Fdwz,TFdwz,Area0  &
                        ,rl,TFz,Sdwb,Wconst,w,V,Vdwz,Vdw,Ap,ups,tups,Wcz,Tcir,wz,Sv,Ics,Ppo4,Ppo40 &
                        ,Vj,Vdwj,dVhor,Hor,Kz,Kzh,Ku,Kl,Vocean,Khor)
    Integer J,Nj,Ics
    double precision Az(Nj+1),AzDW(Nj+1),yz(Nj+1),yz0(Nj+1),Sdw,hm,dz,Sz(Nj+1),Sdwz(Nj+1)          &
                    ,Area,Area0,Areal,Areah,Fz(Nj+1),Fdwz(Nj+1),TFdwz,TFz,Sdwb,Wconst,w,V,Vdwz(Nj) &
                    ,Vdw,Ap,rl,ups(Nj),tups,Wcz(Nj),Tcir,wz(Nj),Sv,Ppo4,Ppo40,Vj(Nj),Vdwj(Nj)      &
                    ,const,dVhor,Hor(Nj),Kz(Nj),Kzh(Nj),Ku,Kl,Vocean,Khor

!=======================!
!  SEAFLOOR TOPOGRAPHY  !
!=======================!
!  Areal fraction of each water layer
    const = 1d0 - Sdw
    Do J = 1, Nj
      Az(J) = (1d0-yz0(J)*0.01d0)*const
    End Do
    Az(Nj+1) = 0d0

!  Upper surface area of each water layer, Sz, Sdwz [m2]
    Sz(1)   = (1d0-Sdw)*Area ! low-mid lat.
    Sdwz(1) = Sdw*Area       ! high lat.
    const = Sdw/(1d0-Sdw)*Area0
    Do J = 2, Nj+1
      Sz(J)   = Az(J)*Area0 ! low-mid lat.
      Sdwz(J) = Az(J)*const ! high lat.
    End Do

!  Seabed area, Fz, Fdwz [m2]
    Fz(1)   = 0d0  ! 100m water depth @ low-mid lat.
    Fdwz(1) = 0d0  ! 100m water depth @ high lat.
    Do J = 2, Nj+1
      Fz(J)   = Sz(J-1)   - Sz(J)   ! low-mid lat.
      Fdwz(J) = Sdwz(J-1) - Sdwz(J) ! high lat.
    End Do

!  Total seabed area
    TFz   = 0d0
    TFdwz = 0d0
    Do J = 1, Nj+1
      TFz   = TFz + Fz(J) + Fdwz(J) ! global
      TFdwz = TFdwz + Fdwz(J)       ! high lat.
    End Do

!  Upper surface area of high lat. deep water
    Sdwb = Area * Sdw

!  Volume [m3]
    V   = 0d0
    Vdw = 0d0
    Do J = 1, Nj
      V   = V   + Sz(J)*dz   ! Total volume @ low-mid lat.
      Vdw = Vdw + Sdwz(J)*dz ! Total volume @ high lat.
      Vj(J)   = Sz(J)*dz     ! low-mid lat.
      Vdwj(J) = Sdwz(J)*dz   ! high lat.
    End Do
    ! Total volume
     Vocean = V + Vdw + Area*dz ! [m^3]

!=====================!
!  OCEAN CIRCULATION  !
!=====================!
!  Overturning rate
    Wconst = w * Sdw
    Sv = Area*Wconst/3.1536d13
    Write(*,'(A10,f6.2,A2)') 'Sv      = ',Area*Wconst/3.1536d13,'Sv'

!  Ocean circulation time [yr]
    Tcir = (V+Vdw+hm*Area)/(Wconst*Area)
    Write(*,'(A13,f8.3,A3)') 'Tcir       = ', Tcir,' yr'
    Write(*,*) '----------------------------------------------'

!  Isopicnal flow from high lat. to low-mid lat.
!  Fraction
    Ap = 0d0
    Do J = 11, Nj
      Ap = Ap + Sz(J)*dz
    End Do
    Do J = 1, 10
      ups(J) = 0d0
    End Do
    Do J = 11, Nj
      ups(J) = (Sz(J)-Sz(J+1))/Sz(11)
    End Do
    tups = 0d0
    Do J = 1, Nj
      tups = tups + ups(J) ! tups should be 1
    End Do

!  Upwelling/downwelling rate [m/yr]
    Do J = 1, 10
      Wcz(J) = Wconst
    End Do
    Wcz(Nj) = ups(Nj)*Wconst
    Do J = 2, 50
      Wcz(Nj+1-J) = Wcz(Nj+2-J) + ups(Nj+1-J)*Wconst
    End Do
    Do J = 1, Nj
      wz(J) = Area*Wcz(J)/Sz(J)
    End Do

!  Isopicnal mixing
    dVhor = Khor*3.1536d7*dz
    Do J = 1, 4
       Hor(J) = 0.001d0
    End Do
    Do J = 5, Nj
       Hor(J) = 3d-5
    End Do

!  Vertical eddy diffusion
    Do J = 1, 4
      Kz(J) = Ku
    End Do
    Do J = 5, 15
      Kz(J) = Kl
    End Do
    Do J = 1, Nj
      Kzh(J) = Kz(J)*2d0
    End Do

    RETURN
    End SUBROUTINE STRUCTURE

!_____________________________________________________________________________________________________!
!_____________________________________________________________________________________________________!
      
    SUBROUTINE PHYS(Nj,X,Xdw,dt,dz,Wcz,Area,Sz,Kz,Kzh,dVhor,ups,Wconst,Xm,Xh,Sdwz,Areah,Vuh,Sdw)
    Integer J,Nj
    double precision X(2,Nj),Xdw(2,Nj),dt,dz,Wcz(Nj),Area,Sz(Nj+1),Sdwz(Nj+1)  &
                    ,Kz(Nj),Kzh(Nj),dVhor,ups(Nj),Wconst,Xm(2),Xh(2),Areah,Vuh,Sdw &
                    ,const,const1,const2,const3

!========================================================!
!  CONCENTRATION CHANGE DUE TO OCEAN CIRCULATION/MIXING  !
!========================================================!
    const = dt/dz
    const1= 1d0/dz
    const2= Area*Wconst
    const3= const*dVhor
!  Surface layers
    X(2,1)   = X(1,1)   + const*((- Wcz(1)*X(1,1) + Wcz(2)*X(1,2))*Area/Sz(1)   &
                                  + ups(1)*const2*Xdw(1,1)/Sz(1)                &
                                  + (Kz(1)*const1)*(Xm(1)-X(1,1))               &
                                  -(Kz(2)*const1)*(X(1,1)-X(1,2))*Sz(2)/Sz(1))  &
                        + const3*(Xdw(1,1)-X(1,1))/Sz(1)
    Xdw(2,1) = Xdw(1,1) + const*(const2*(Xh(1)-Xdw(1,1))/Sdwz(1)                         &
                                 + Areah*Vuh*(Xh(1)-Xdw(1,1))/(Area*Sdw)                 &
                                 - (Kzh(2)*const1)*(Xdw(1,1)-Xdw(1,2))*Sdwz(2)/Sdwz(1))  &
                        - const3*(Xdw(1,1)-X(1,1))/Sdwz(1)
!  Abyssal layers
    Do J = 2, Nj-1
      X(2,J)   =  X(1,J)  + const*((-Wcz(J)*X(1,J)+Wcz(J+1)*X(1,J+1))*Area/Sz(J)         &
                                    + ups(J)*const2*Xdw(1,J)/Sz(J)                       &
                                    + (Kz(J)*const1)*(X(1,J-1)-X(1,J))                   &
                                    - (Kz(J+1)*const1)*(X(1,J)-X(1,J+1))*Sz(J+1)/Sz(J))  &
                          + const3*(Xdw(1,J)-X(1,J))/Sz(J)
      Xdw(2,J) = Xdw(1,J) + const*(Area*(Wcz(J)*Xdw(1,J-1)-Wcz(J+1)*Xdw(1,J))/Sdwz(J)           &
                                  - ups(J)*const2*Xdw(1,J)/Sdwz(J)                              &
                                  + (Kzh(J)*const1)*(Xdw(1,J-1)-Xdw(1,J))                       &
                                  - (Kzh(J+1)*const1)*(Xdw(1,J)-Xdw(1,J+1))*Sdwz(J+1)/Sdwz(J))  &
                          - const3*(Xdw(1,J)-X(1,J))/Sdwz(J)
    End Do
!  Bottom layers
    X(2,Nj)   = X(1,Nj)   + const*(- Wcz(Nj)*X(1,Nj)*Area/Sz(Nj)           &
                                   + ups(Nj)*const2*Xdw(1,Nj)/Sz(Nj)       &
                                   + (Kz(Nj)*const1)*(X(1,Nj-1)-X(1,Nj)))  &
                          + const3*(Xdw(1,Nj)-X(1,Nj))/Sz(Nj)
    Xdw(2,Nj) = Xdw(1,Nj) + const*(Area*Wcz(Nj)*Xdw(1,Nj-1)/Sdwz(Nj)            &
                                   - ups(Nj)*const2*Xdw(1,Nj)/Sdwz(Nj)          &
                                   + (Kzh(Nj)*const1)*(Xdw(1,Nj-1)-Xdw(1,Nj)))  &
                          - const3*(Xdw(1,Nj)-X(1,Nj))/Sdwz(Nj)

    RETURN
    End SUBROUTINE PHYS

!_____________________________________________________________________________________________________!
!_____________________________________________________________________________________________________!
    SUBROUTINE ORG(Ro2,Rno3,Rso4,CjO2,CjNO3,CjSO4,Kso4,r1o2,r1no3                                &
                  ,r1so4,r2o2,r2no3,r2so4,r3o2,r3no3                                             &
                  ,r3so4,r1,r2,r3,Ppom1,Ppom2,Ppom3,Area,Fpom1,Fpom2,Fpom3,Vpom,pom              &
                  ,dis1,dis2,dis3,Az,TRo2,TRno3,TRso4,dispom,Bpom,Fz,Rpc,Po,N,Depo               &
                  ,FreeP,cpr,BpCa,CmO2,BpFe,RETURNP,BurP,TBurP,Trep,TdepoP,DepoDW,freePdw        &
                  ,TdepoPdw,Fdwz,cprdw,SdepoP,BpCadw,BpFeDW,BpOrgDW,RETURNPdw,CdwO2,BurPdw       &
                  ,TBurPdw,TrepDW,SburP,Sdwb,Foh,Fol,FdwPOM1,FdwPOM2,FdwPOM3,BpOrg               &
                  ,r1dw,r2dw,r3dw,Areah,Areal,dt,Ppo4,BCorg,TBCorg,RETURNC,BCorgDW,TBCorgDW      &
                  ,RETURNCDW,Sz,CPbur,CPburDW,tb,BE,BEdw,dispomDW,Sdwz,dRo2,dRno3,dRso4,CdwNO3   &
                  ,CdwSO4,dis1dw,dis2dw,dis3dw,Temp,TempDW,dispomSED,dispomdwSED,Ro2SED,Rno3SED  &
                  ,Rso4SED,dRo2SED,dRno3SED,dRso4SED,Fdeni,FdeniDW,Bpyr,Bdwpyr,TBpyr,CSj         &
                  ,CSdwj,fpyrAn,Rch4,dRch4,r1ch4,r1ch4dw,r2ch4,r2ch4dw,r3ch4                     &
                  ,r3ch4dw,Rch4SED,dRch4SED,TRch4,fpyrOx,O2js,O2dwjs,AOM,AOMdw,tTRo2,tTRno3      &
                  ,tTRso4,tTRch4,epyAnox,Zopd,ZopdDW,Zsmtz,ZsmtzDW,topd,topdDW,tsmtz             &
                  ,tsmtzDW,AOMsed,AOMsedDW,ksed,BSRsed,BSRsedj,BSRsedjdw                         &
                  ,SR,tAOMsed,Rcp,RcpH,Rnc,Rnp,Roc,BEporg,BEporgDW,CjPO4,CdwPO4                  &
                  ,aveBEc,RcpDW,RnpDW,RncDW,RpcDW,RocDW,RnpH,RncH,RpcH,RocH,BEporg0,BEporgDW0    &
                  ,BEs,BEdws,Ncount,beoxic,beanox)
    USE Constants
    Implicit none
    Integer I,J,N
    Integer(selected_int_kind(r=16)) Ncount
    double precision Ro2(Nj),Rno3(Nj),Rso4(Nj),CjO2(2,Nj),CjNO3(2,Nj)                                 &
                    ,CjSO4(2,Nj),Kso4,r1o2(Nj),r1no3(Nj),r1so4(Nj),r2o2(Nj)                           &
                    ,r2no3(Nj),r2so4(Nj),r3o2(Nj),r3no3(Nj),r3so4(Nj)                                 &
                    ,r1(Nj),r2(Nj),r3(Nj),Ppom1,Ppom2,Ppom3,Area,Fpom1(Nj+1),Fpom2(Nj+1),Fpom3(Nj+1)  &
                    ,Vpom,pom,pomdw,dis1(Nj),dis2(Nj),dis3(Nj),Az(Nj+1),TRo2,TRno3,TRso4              &
                    ,dispom(Nj),Bpom,Fz(Nj+1),Rpc,Po,Depo(Nj+1),FreeP(Nj+1),cpr(Nj),BpCa(Nj)          &
                    ,CmO2(2),BpFe(Nj),BpOrg(Nj),RETURNP(Nj),BurP(Nj+1),TBurP,Trep,TdepoP              &
                    ,DepoDW(Nj+1),freePdw(Nj+1),TdepoPdw,Fdwz(Nj+1),cprdw(Nj),SdepoP,BpCadw(Nj)       &
                    ,BpFeDw(Nj),RETURNPdw(Nj),CdwO2(2,Nj),BurPdw(Nj+1),TBurPdw,TrepDW,SBurP,Sdwb      &
                    ,Foh,Fol,FdwPOM1(Nj+1),FdwPOM2(Nj+1),FdwPOM3(Nj+1),r1dw(Nj),r2dw(Nj),r3dw(Nj)     &
                    ,Areah,Areal,dt,Ppo4,BCorg(Nj),TBCorg,RETURNC(Nj),BpOrgDW(Nj),BCorgDW(Nj)         &
                    ,TBCorgDW,RETURNCDW(Nj),Sz(Nj+1),CPbur,CPburDW,tb(Nj),BE(Nj),BEdw(Nj),dispomDW(Nj)&
                    ,Sdwz(Nj+1),dRo2(Nj),dRno3(Nj),dRso4(Nj),r1o2dw(Nj),r1no3dw(Nj),r1so4dw(Nj)       &
                    ,r2o2dw(Nj),r2no3dw(Nj),r2so4dw(Nj),r3o2dw(Nj),r3no3dw(Nj),r3so4dw(Nj)            &
                    ,CdwNO3(2,Nj),CdwSO4(2,Nj),dis1dw(Nj),dis2dw(Nj),dis3dw(Nj),Temp(2,Nj)            &
                    ,TempDW(2,Nj),dispomSED(Nj),dispomdwSED(Nj),Ro2SED(Nj),Rno3SED(Nj)                &
                    ,Rso4SED(Nj),dRo2SED(Nj),dRno3SED(Nj),dRso4SED(Nj),Fdeni(Nj),FdeniDW(Nj)          &
                    ,Tfpom(Nj+1),TfpomDW(Nj+1),Bpyr(Nj),Bdwpyr(Nj),TBpyr,CSj(Nj),CSdwj(Nj),fpyrAn     &
                    ,Rch4(Nj),dRch4(Nj),r1ch4(Nj),r1ch4dw(Nj),r2ch4(Nj),r2ch4dw(Nj)                   &
                    ,r3ch4(Nj),r3ch4dw(Nj),Rch4SED(Nj),dRch4SED(Nj),TRch4,fpyrOx,fpyr,O2js(Nj+1)      &
                    ,O2dwjs(Nj+1),AOM(Nj),AOMdw(Nj),tTRo2,tTRno3,tTRso4,tTRch4,pomSED,pomdwSED        &
                    ,epyAnox,fpyrDW,Zopd(Nj),ZopdDW(Nj),Zsmtz(Nj),ZsmtzDW(Nj),topd(Nj),topdDW(Nj)     &
                    ,tsmtz(Nj),tsmtzDW(Nj),AOMsed(Nj),AOMsedDW(Nj),BSRsed                             &
                    ,BSRsedj(Nj),BSRsedjdw(Nj),ksed(Nj)                                               &
                    ,SR(Nj),x,y,z,tAOMsed,tt,ttdw                                                     &
                    ,Rcp,RcpH,Rnc,Rnp,Roc,BEporg(Nj),BEporgDW(Nj),CjPO4(2,Nj),CdwPO4(2,Nj),aveBEc     &
                    ,TdepoC,RcpDW,RnpDW,RncDW,RpcDW,RocDW,RnpH,RncH,RpcH,RocH,BEporg0(Nj)             &
                    ,BEporgDW0(Nj),BEs(Nj),BEdws(Nj),const,const1,const2,const3                       &
                    ,BEporg_oxic(Nj),BEporgDW_oxic(Nj),beoxic,beanox,BE_oxic(Nj),BEdw_oxic(Nj)        &
                    ,Tdeni,TorgP,TfeP,Tcap,TreturnP,Tdeniwc

!==================================!
!  BIOLOGICAL PUMP & SEDIMENTATION !
!==================================!
!--------------------------!
!  POM decomposition rate  !
!--------------------------!
      Do J = 1, Nj
!  Respiration pathway
      ! low-mid lat.
        Ro2(J)   = CjO2(1,J)/(Ko2+CjO2(1,J))
        Rno3(J)  = Ko2d/(Ko2d+CjO2(1,J))*CjNO3(1,J)/(Kno3+CjNO3(1,J))
        Rso4(J)  = Ko2d/(Ko2d+CjO2(1,J))*Kno3d/(Kno3d+CjNO3(1,J))*CjSO4(1,J)/(Kso4+CjSO4(1,J))
        Rch4(J)  = 1d0-Ro2(J)-Rno3(J)-Rso4(J)
      ! high lat.
        dRo2(J)  = CdwO2(1,J)/(Ko2+CdwO2(1,J))
        dRno3(J) = Ko2d/(Ko2d+CdwO2(1,J))*CdwNO3(1,J)/(Kno3+CdwNO3(1,J))
        dRso4(J) = Ko2d/(Ko2d+CdwO2(1,J))*Kno3d/(Kno3d+CdwNO3(1,J))*CdwSO4(1,J)/(Kso4+CdwSO4(1,J))
        dRch4(J) = 1d0-dRo2(J)-dRno3(J)-dRso4(J)

!  Remineralization rate of labile organic matter (G1)
        r1o2(J)    = k1o2  * Ro2(J)   * dexp(Ktox*(Temp(1,J)-15d0))
        r1no3(J)   = k1no3 * Rno3(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r1so4(J)   = k1so4 * Rso4(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r1ch4(J)   = k1ch4 * Rch4(J)  * dexp(Ktox*(Temp(1,J)-15d0))
         r1o2dw(J) = k1o2  * dRo2(J)  * dexp(Ktox*(TempDW(1,J)-15d0))
        r1no3dw(J) = k1no3 * dRno3(J) * dexp(Ktox*(TempDW(1,J)-15d0))
        r1so4dw(J) = k1so4 * dRso4(J) * dexp(Ktox*(TempDW(1,J)-15d0))
        r1ch4dw(J) = k1ch4 * dRch4(J) * dexp(Ktox*(TempDW(1,J)-15d0))
! Remineralization rate of semi-labile organic matter (G2)
        r2o2(J)    = k2o2  * Ro2(J)   * dexp(Ktox*(Temp(1,J)-15d0))
        r2no3(J)   = k2no3 * Rno3(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r2so4(J)   = k2so4 * Rso4(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r2ch4(J)   = k2ch4 * Rch4(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r2o2dw(J)  = k2o2  * dRo2(J)  * dexp(Ktox*(TempDW(1,J)-15d0))
        r2no3dw(J) = k2no3 * dRno3(J) * dexp(Ktox*(TempDW(1,J)-15d0))
        r2so4dw(J) = k2so4 * dRso4(J) * dexp(Ktox*(TempDW(1,J)-15d0))
        r2ch4dw(J) = k2ch4 * dRch4(J) * dexp(Ktox*(TempDW(1,J)-15d0))
! Remineralization rate of inert organic matter (G3)
        r3o2(J)    = k3o2  * Ro2(J)   * dexp(Ktox*(Temp(1,J)-15d0))
        r3no3(J)   = k3no3 * Rno3(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r3so4(J)   = k3so4 * Rso4(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r3ch4(J)   = k3ch4 * Rch4(J)  * dexp(Ktox*(Temp(1,J)-15d0))
        r3o2dw(J)  = k3o2  * dRo2(J)  * dexp(Ktox*(TempDW(1,J)-15d0))
        r3no3dw(J) = k3no3 * dRno3(J) * dexp(Ktox*(TempDW(1,J)-15d0))
        r3so4dw(J) = k3so4 * dRso4(J) * dexp(Ktox*(TempDW(1,J)-15d0))
        r3ch4dw(J) = k3ch4 * dRch4(J) * dexp(Ktox*(TempDW(1,J)-15d0))

        If(Temp(1,J) < 15d0) then
          r1o2(J)  = k1o2  * Ro2(J)
          r1no3(J) = k1no3 * Rno3(J)
          r1so4(J) = k1so4 * Rso4(J)
          r1ch4(J) = k1ch4 * Rch4(J)
          r2o2(J)  = k2o2  * Ro2(J)
          r2no3(J) = k2no3 * Rno3(J)
          r2so4(J) = k2so4 * Rso4(J)
          r2ch4(J) = k2ch4 * Rch4(J)
          r3o2(J)  = k3o2  * Ro2(J)
          r3no3(J) = k3no3 * Rno3(J)
          r3so4(J) = k3so4 * Rso4(J)
          r3ch4(J) = k3ch4 * Rch4(J)
        End If
        If(TempDW(1,J) < 15d0) then
          r1o2dw(J)  = k1o2  * dRo2(J)
          r1no3dw(J) = k1no3 * dRno3(J)
          r1so4dw(J) = k1so4 * dRso4(J)
          r1ch4dw(J) = k1ch4 * dRch4(J)
          r2o2dw(J)  = k2o2  * dRo2(J)
          r2no3dw(J) = k2no3 * dRno3(J)
          r2so4dw(J) = k2so4 * dRso4(J)
          r2ch4dw(J) = k2ch4 * dRch4(J)
          r3o2dw(J)  = k3o2  * dRo2(J)
          r3no3dw(J) = k3no3 * dRno3(J)
          r3so4dw(J) = k3so4 * dRso4(J)
          r3ch4dw(J) = k3ch4 * dRch4(J)
        End If
!  Remineralization rate
        r1(J)   = r1o2(J)   + r1no3(J)   + r1so4(J)   + r1ch4(J)
        r2(J)   = r2o2(J)   + r2no3(J)   + r2so4(J)   + r2ch4(J)
        r3(J)   = r3o2(J)   + r3no3(J)   + r3so4(J)   + r3ch4(J)
        r1dw(J) = r1o2dw(J) + r1no3dw(J) + r1so4dw(J) + r1ch4dw(J)
        r2dw(J) = r2o2dw(J) + r2no3dw(J) + r2so4dw(J) + r2ch4dw(J)
        r3dw(J) = r3o2dw(J) + r3no3dw(J) + r3so4dw(J) + r3ch4dw(J)
      End Do

!----------------------------!
!  C/P ratio of sinking POM  !
!----------------------------!
       RcpDW = (Fol*(Sdwb-Areah)+Foh*Areah)/(Fol*Rpc*(Sdwb-Areah)+Foh*RpcH*Areah)
       RpcDW = 1d0/RcpDW
       RncDW = (Rnc*Fol*(Sdwb-Areah)+RncH*Foh*Areah)/(Fol*(Sdwb-Areah)+Foh*Areah)
       RnpDW = RncDW*RcpDW

!------------------------------------------!
!  POM settling flux density [molC/m2/yr]  !
!------------------------------------------!
       Fpom1(1)   = Ppom1/Areal
       Fpom2(1)   = Ppom2/Areal
       Fpom3(1)   = Ppom3/Areal
       FdwPOM1(1) = m1*(Areah*Foh+(Sdwb-Areah)*Fol)/Sdwb
       FdwPOM2(1) = m2*(Areah*Foh+(Sdwb-Areah)*Fol)/Sdwb
       FdwPOM3(1) = m3*(Areah*Foh+(Sdwb-Areah)*Fol)/Sdwb
       const = dz/Vpom
       Do J = 2, Nj+1
          Fpom1(J)   = Fpom1(J-1) * dexp(-r1(J-1)*const) ! G1
          Fpom2(J)   = Fpom2(J-1) * dexp(-r2(J-1)*const) ! G2
          Fpom3(J)   = Fpom3(J-1) * dexp(-r3(J-1)*const) ! G3
          FdwPOM1(J) = FdwPOM1(J-1) * dexp(-r1dw(J-1)*const) ! G1
          FdwPOM2(J) = FdwPOM2(J-1) * dexp(-r2dw(J-1)*const) ! G2
          FdwPOM3(J) = FdwPOM3(J-1) * dexp(-r3dw(J-1)*const) ! G3
       End Do

!  POM settling flux density
       Do J = 1, Nj+1
          Tfpom(J)   = Fpom1(J)   + Fpom2(J)   + Fpom3(J)   ! low-mid lat.
          TfpomDW(J) = FdwPOM1(J) + FdwPOM2(J) + FdwPOM3(J) ! high lat.
       End Do

!==============!
!  OPD & SMTZ  !
!==============!
      const = (1d0-porosity)*density
      Do J = 1, Nj
       ! OPD @ low-mid lat.
         x = dlog10(SR(J))
         y = dlog10(CjO2(1,J)*1d3)
         z = dlog10(Tfpom(J+1)*0.1d0)
         tt= Temp(1,J)
         Zopd(J) = 10d0**(a0opd + a1opd*x + a2opd*y + a3opd*z + a4opd*x**2d0 + a5opd*y**2d0 &
                        + a6opd*z**2d0 + a7opd*x*y + a8opd*y*z + a9opd*x*z + a10opd*tt)
       ! OPD @ high lat.
         y = dlog10(CdwO2(1,J)*1d3)
         z = dlog10(TfpomDW(J+1)*0.1d0)
         ttdw = TempDW(1,J)
         ZopdDW(J) = 10d0**(a0opd + a1opd*x + a2opd*y + a3opd*z + a4opd*x**2d0 + a5opd*y**2d0 &
                          + a6opd*z**2d0 + a7opd*x*y + a8opd*y*z + a9opd*x*z + a10opd*ttdw)
         If((CjO2(1,J) <= 1d-3).or.(Zopd(J) <= 0d0)) then
           Zopd(J) = 0d0
         End If
         If((CdwO2(1,J) <= 1d-3).or.(ZopdDW(J) <= 0d0)) then
           ZopdDW(J) = 0d0
         End If
       ! SMTZ
         Zsmtz(J)   = ksmtz * CjSO4(1,J)**Asmtz  * (const*SR(J))**nopd
         ZsmtzDW(J) = ksmtz * CdwSO4(1,J)**Asmtz * (const*SR(J))**nopd
         topd(J)    = Zopd(J)/SR(J)
         topdDW(J)  = ZopdDW(J)/SR(J)
         tsmtz(J)   = (Zsmtz(J)-Zopd(J))/SR(J)
         tsmtzDW(J) = (ZsmtzDW(J)-ZopdDW(J))/SR(J)
      End Do

!------------------------!
!  POM dissolution flux  !
!------------------------!
      pom   = 0d0
      pomdw = 0d0
      pomSED   = 0d0
      pomdwSED = 0d0
      Do J = 1, Nj
      ! Flux density
        dis1(J)    = Fpom1(J)   - Fpom1(J+1)
        dis2(J)    = Fpom2(J)   - Fpom2(J+1)
        dis3(J)    = Fpom3(J)   - Fpom3(J+1)
        dis1DW(J)  = Fdwpom1(J) - Fdwpom1(J+1)
        dis2DW(J)  = Fdwpom2(J) - Fdwpom2(J+1)
        dis3DW(J)  = Fdwpom3(J) - Fdwpom3(J+1)
        dispom(J)  = dis1(J)   + dis2(J)   + dis3(J)
        dispomDW(J)= dis1DW(J) + dis2DW(J) + dis3DW(J)
      ! Total decomposition flux
        pom   = pom   + dispom(J)   * Sz(J)
        pomdw = pomdw + dispomDW(J) * Sdwz(J)
      ! Benthic remineralization flux density
        dispomSED(J)   = Tfpom(J+1)   * (1d0-BE(J))   * Fz(J+1)/Sz(J)
        dispomdwSED(J) = TfpomDW(J+1) * (1d0-BEdw(J)) * Fdwz(J+1)/Sdwz(J)
      ! Total Benthic remineralization flux
        pomSED   = pomSED   + Tfpom(J+1)   * (1d0-BE(J))   * Fz(J+1)
        pomdwSED = pomdwSED + TfpomDW(J+1) * (1d0-BEdw(J)) * Fdwz(J+1)
      End Do

!--------------------------------!
!  Denitrification [molC/m2/yr]  ! ()
!--------------------------------!
      const = 1d0/3.65d0
      Do J = 1, Nj
        Fdeni(J) = 3.65d0*10d0**(                                                 &
                     - 2.2567d0 - 0.185d0*dlog10(Tfpom(J+1)*const)                &
                     - 0.221d0*dlog10(Tfpom(J+1)*const)*dlog10(Tfpom(J+1)*const)  &
                     - 0.3995d0*dlog10(CjNO3(1,J)*rho)*dlog10(CjO2(1,J)*rho)      &
                     + 1.25d0*dlog10(CjNO3(1,J)*rho)                              &
                     + 0.4721d0*dlog10(CjO2(1,J)*rho)                             &
                     - 0.0996d0*dlog10(J*dz+hm)                                   &
                     + 0.4256d0*dlog10(Tfpom(J+1)*const)*dlog10(CjO2(1,J)*rho))
        FdeniDW(J) = 3.65d0*10d0**(                                                   &
                     - 2.2567d0 - 0.185d0*dlog10(TfpomDW(J+1)*const)                  &
                     - 0.221d0*dlog10(TfpomDW(J+1)*const)*dlog10(TfpomDW(J+1)*const)  &
                     - 0.3995d0*dlog10(CdwNO3(1,J)*rho)*dlog10(CdwO2(1,J)*rho)        &
                     + 1.25d0*dlog10(CdwNO3(1,J)*rho)                                 &
                     + 0.4721d0*dlog10(CdwO2(1,J)*rho)                                &
                     - 0.0996d0*dlog10(J*dz+hm)                                       &
                     + 0.4256d0*dlog10(TfpomDW(J+1)*const)*dlog10(CdwO2(1,J)*rho))
      End Do

!-------------------------------!
!  P deposition flux [molP/yr]  !
!-------------------------------!
      Depo(1)   = 0d0
      DepoDW(1) = 0d0
      TdepoP    = 0d0
      TdepoPdw  = 0d0
      Do J = 2, Nj+1
       ! Deposition flux
         Depo(J)   = Rpc   * (Fpom1(J)+Fpom2(J)+Fpom3(J))*Fz(J)
         DepoDW(J) = RpcDW * (Fdwpom1(J)+Fdwpom2(J)+Fdwpom3(J))*Fdwz(J)
       ! Total deposition flux
         TdepoP    = TdepoP   + Depo(J)   ! low-mid lat.
         TdepoPdw  = TdepoPdw + DepoDW(J) ! high lat.
      End Do
      SdepoP = TdepoP + TdepoPdw ! global

!----------------------------------------!
!  Corg deposition/burial flux [molC/yr] !
!----------------------------------------!
      TBCorg   = 0d0
      TBCorgDW = 0d0
      TdepoC   = 0d0
      Do J = 1, Nj
       ! Total deposition flux
         TdepoC       = TdepoC + (Depo(J+1)+DepoDW(J+1))*Rcp
       ! Burial flux
         BCorg(J)     = BE(J)*Depo(J+1)*Rcp
         TBCorg       = TBCorg + BCorg(J)
         BCorgDW(J)   = BEdw(J)*DepoDW(J+1)*RcpDW
         TBCorgDW     = TBCorgDW + BCorgDW(J)
       ! Benthic remineralization rate
         RETURNC(J)   = (Fpom1(J+1)+Fpom2(J+1)+Fpom3(J+1))*Fz(J+1) - BCorg(J)
         RETURNCDW(J) = (Fdwpom1(J+1)+Fdwpom2(J+1)+Fdwpom3(J+1))*Fdwz(J+1) - BCorgDW(J)
      End Do

!---------!
!  %deni  !
!---------!
      Do J = 1, Nj
       Rno3SED(J)  = Fdeni(J)/(Tfpom(J+1)*(1d0-BE(J)))
       dRno3SED(J) = FdeniDW(J)/(TfpomDW(J+1)*(1d0-BEdw(J)))
         If(Rno3SED(J) >= 0.9d0) then
           Rno3SED(J) = 0.9d0
           Fdeni(J) = Rno3SED(J)*(Tfpom(J+1)*(1d0-BE(J)))
         End If
         If(dRno3SED(J) >= 0.9d0) then
           dRno3SED(J) = 0.9d0
           FdeniDW(J) = dRno3SED(J)*(TfpomDW(J+1)*(1d0-BEdw(J)))
         End If
         Ro2SED(J) = (1d0-Rno3SED(J))*(1d0-dexp(-ksed(J)*topd(J)))
         If((1d0-dexp(-ksed(J)*topd(J))) >= 1d0) then
           Rso4SED(J) = 0d0
           Rch4SED(J) = 0d0
         Else
           Rso4SED(J) = (1d0-Ro2SED(J)-Rno3SED(J))*CjSO4(1,J)/(Kso4+CjSO4(1,J))
           Rch4SED(J) = 1d0 - Rno3SED(J) - Ro2SED(J) - Rso4SED(J)
         End If
         dRo2SED(J)=(1d0-dRno3SED(J))*(1d0-dexp(-ksed(J)*topdDW(J)))
         If((1d0-dexp(-ksed(J)*topdDW(J))) >= 1d0) then
           dRso4SED(J) = 0d0
           dRch4SED(J) = 0d0
         Else
           dRso4SED(J) = (1d0-dRo2SED(J)-dRno3SED(J))*CdwSO4(1,J)/(Kso4+CdwSO4(1,J))
           dRch4SED(J) = 1d0 - dRno3SED(J) - dRo2SED(J) - dRso4SED(J)
         End If
      End Do
! AOMsed
      tAOMsed = 0d0
      Do J = 1, Nj
        AOMsed(J) = Rch4SED(J)*Tfpom(J+1)*(1d0-BE(J))*0.5d0*CjSO4(1,J)/(Kso4+CjSO4(1,J))
        If(AOMsed(J) <= 0d0) then
          AOMsed(J) = 0d0
        End If
        AOMsedDW(J) = dRch4SED(J)*TfpomDW(J+1)*(1d0-BEdw(J))*0.5d0*CdwSO4(1,J)/(Kso4+CdwSO4(1,J))
        If(AOMsedDW(J) <= 0d0) then
          AOMsedDW(J) = 0d0
        End If
        tAOMsed = tAOMsed + AOMsed(J)*Fz(J+1) + AOMsedDW(J)*Fdwz(J+1)
      End Do

!=================!
!  Pyrite burial  !
!=================!
      BSRsed = 0d0
      Do J = 1, Nj
       ! Pyrite burial efficiency
         fpyr   = epyAnox - (epyAnox-epyOxic) * dtanh(CjO2(1,J)*1d3)
         fpyrDW = epyAnox - (epyAnox-epyOxic) * dtanh(CdwO2(1,J)*1d3)
       ! C/S ratio
         CSj(J)   = Bcorg(J)/(fpyr*(Rso4SED(J)*Tfpom(J+1)*(1d0-BE(J))*Fz(J+1)*0.5d0+AOMsed(J)*Fz(J+1)))
         CSdwj(J) = BcorgDW(J)/(fpyrDW*(dRso4SED(J)*TfpomDW(J+1)*(1d0-BEdw(J))*Fdwz(J+1)*0.5d0+AOMsedDW(J)*Fdwz(J+1)))
       ! Benthic sulfate reduction
         BSRsedj(J)   = Rso4SED(J)  * Tfpom(J+1)   * (1d0-BE(J))   * Fz(J+1)   * 0.5d0
         BSRsedjdw(J) = dRso4SED(J) * TfpomDW(J+1) * (1d0-BEdw(J)) * Fdwz(J+1) * 0.5d0
         BSRsed = BSRsed + BSRsedj(J) + BSRsedjdw(J)
      End Do
    ! Pyrite burial
      TBpyr = 0d0
      Do J = 1, Nj
        Bpyr(J)   = BCorg(J)/CSj(J)         ! low-mid lat.
        Bdwpyr(J) = BCorgDW(J)/CSdwj(J)     ! high lat.
        TBpyr = TBpyr + Bpyr(J) + Bdwpyr(J) ! global
      End Do

      const  = (1d0-porosity)*density
      const1 = 1d0/0.07d0
      Do J = 1, Nj
        BE_oxic(J)   = (beoxic/(1d0+SR(J)*const*const1)+50d0)*0.01d0
        BEdw_oxic(J) = (beoxic/(1d0+SR(J)*const*const1)+50d0)*0.01d0
      End Do

!----------------!
!  (Corg/Porg)b  !
!----------------!
      CPbur   = 2d0*Rcp
      CPburDW = 2d0*RcpDW
      const = 1d0/oxic
      const1= CPbur*ccpr
      const2= const*ccpr
      const3= CPburDW*ccpr
      Do J = 1, Nj
        cpr(J) = const1/(CjO2(1,J)*const2+(1d0-CjO2(1,J)*const)*CPbur)*(1d0+dexp(-tb(J)*1d-4))*0.5d0
          If(CjO2(1,J) > oxic) then
            cpr(J) = CPbur*(1d0+dexp(-tb(J)*1d-4))*0.5d0
          End If
        cprDW(J) = const3/(CdwO2(1,J)*const2+(1d0-CdwO2(1,J)*const)*CPburDW)*(1d0+dexp(-tb(J)*1d-4))*0.5d0
          If(CdwO2(1,J) > oxic) then
            cprDW(J) = CPburDW*(1d0+dexp(-tb(J)*1d-4))*0.5d0
          End If
      End Do

!============!
!  P burial  !
!============!
      const = 1d0/oxic
      Do J = 1, Nj
!  Burial efficiency of organic P
        BEporg(J)   = BE(J)   * Rcp/cpr(J)
        BEporgDW(J) = BEdw(J) * RcpDW/cprDW(J)
        BEporg_oxic(J)   = BE_oxic(J)   * Rcp/(CPbur*(1d0+dexp(-tb(J)*1d-4))*0.5d0)
        BEporgDW_oxic(J) = BEdw_oxic(J) * RcpDW/(CPburDW*(1d0+dexp(-tb(J)*1d-4))*0.5d0)
!  Burial rate @ low-mid lat.
        If(CjO2(1,J) < oxic) then
          BpOrg(J)= BEporg(J)      * Depo(J+1)  ! Organic P burial
          BpCa(J) = BEporg_oxic(J) * Depo(J+1) * 2d0*(1d0+(1d0-1d0)*CjO2(1,J)*const)   ! Authigenic P burial
          BpFe(J) = BEporg_oxic(J) * Depo(J+1) * CjO2(1,J)*const                       ! Fe-bound P burial
        Else
          BpOrg(J) = BEporg(J)      * Depo(J+1)     ! Organic P burial
          BpCa(J)  = BEporg_oxic(J) * Depo(J+1)*2d0 ! Authigenic P burial
          BpFe(J)  = BEporg_oxic(J) * Depo(J+1)     ! Fe-bound P burial
        End If
!  Benthic efflux
        freeP(J+1) = Depo(J+1) - BpOrg(J) - BpCa(J) - BpFe(J)
        RETURNP(J) = freeP(J+1)
        If(RETURNP(J) <= 0d0) then
          write(*,*) 'RETURNP<0!!',J
          BpCa(J) = freeP(J+1) - BpFe(J)
          RETURNP(J) = 0d0
        End If
!  Burial rate @ high lat.
        If(CdwO2(1,J) < oxic) then
          BpOrgDW(J)= BEporgDW(J)      * DepoDW(J+1)  ! Organic P burial
          BpCaDW(J) = BEporgDW_oxic(J) * DepoDW(J+1) * 2d0*(1d0+(1d0-1d0)*CdwO2(1,J)*const)   ! Authigenic P burial
          BpFeDW(J) = BEporgDW_oxic(J) * DepoDW(J+1) * CdwO2(1,J)*const                       ! Fe-bound P burial
        Else
          BpOrgDW(J) = BEporgDW(J)      * DepoDW(J+1)       ! Organic P burial
          BpCaDW(J)  = BEporgDW_oxic(J) * DepoDW(J+1) * 2d0 ! Authigenic P burial
          BpFeDW(J)  = BEporgDW_oxic(J) * DepoDW(J+1)       ! Fe-bound P burial
        End If
!  Benthic efflux
        freePDW(J+1) = DepoDW(J+1) - BpOrgDW(J) - BpCaDW(J) - BpFeDW(J)
        RETURNPDW(J) = freePDW(J+1)
        If(RETURNPDW(J) <= 0d0) then
          write(*,*) 'RETURNP<0!!'
          BpCaDW(J) = freePDW(J+1) - BpFeDW(J)
          RETURNPDW(J) = 0d0
        End If
      End Do
!  Total burial rate
      BurP(1)   = 0d0
      BurPDW(1) = 0d0
      Do J = 2, Nj+1
        BurP(J)   = BpOrg(J-1)   + BpFe(J-1)   + BpCa(J-1)   ! low-mid lat.
        BurPDW(J) = BpOrgDW(J-1) + BpFeDW(J-1) + BpCaDW(J-1) ! high lat.
      End Do
      TBurP   = 0d0
      TBurPdw = 0d0
      Trep    = 0d0
      TrepDW  = 0d0
      Do J = 1, Nj+1
        TBurP   = TBurP   + BurP(J)
        TBurPdw = TBurPdw + BurPDW(J)
      End Do
      Do J = 1, Nj
        If(RETURNP(J) < 0d0) then
          Write(*,*) '!!'
          Trep = Trep + RETURNP(J)
        End If
        If(RETURNPDW(J) < 0d0) then
          Write(*,*) '!!'
          TrepDW = TrepDW + RETURNPDW(J)
        End If
      End Do
      TBurP   = TBurP   + Trep   ! low-mid lat.
      TBurPdw = TBurPdw + TrepDW ! high lat.
      SBurP = TBurP + TBurPdw    ! global

      Torgp = 0d0
      Tcap  = 0d0
      Tfep  = 0d0
      Do J = 1, Nj
        Tcap  = Tcap  + BpCa(J)  + BpCadw(J)
        Tfep  = Tfep  + BpFe(J)  + BpFedw(J)
        Torgp = Torgp + BpOrg(J) + BpOrgDW(J)
        If(RETURNP(J) < 0) then
          Tcap = Tcap + RETURNP(J)
        End If
        If(RETURNPdw(J) < 0) then
          Tcap = Tcap + RETURNPdw(J)
        End If
      End Do
      TreturnP = 0d0
      Do J = 1, Nj
        TreturnP = TreturnP + returnP(J) + returnPdw(J)
      End Do
      Tdeni = 0d0
      Tdeniwc = 0d0
      Do J = 1, Nj
        Tdeni = Tdeni + 4d0/5d0*(Fdeni(J)*Fz(J+1)+FdeniDW(J)*Fdwz(J+1))
        Tdeniwc = Tdeniwc + 4d0/5d0*(dispom(J)*Rno3(J)*Sz(J)+dispomDW(J)*dRno3(J)*Sdwz(J))
      End Do

!-------------------!
!  Progress report  !
!-------------------!
      If(Mod(Ncount*dt,5d4) == 0) then
        TRo2   = 0d0
        TRno3  = 0d0
        TRso4  = 0d0
        TRch4  = 0d0
        tTRo2  = 0d0
        tTRno3 = 0d0
        tTRso4 = 0d0
        tTRch4 = 0d0
        Do J = 1, Nj
          TRo2  = TRo2  + (dispom(J)*Ro2(J)*Sz(J)+dispomDW(J)*dRo2(J)*Sdwz(J))/(pom+pomdw)
          TRno3 = TRno3 + (dispom(J)*Rno3(J)*Sz(J)+dispomDW(J)*dRno3(J)*Sdwz(J))/(pom+pomdw)
          TRso4 = TRso4 + (dispom(J)*Rso4(J)*Sz(J)+dispomDW(J)*dRso4(J)*Sdwz(J))/(pom+pomdw)
          TRch4 = TRch4 + (dispom(J)*Rch4(J)*Sz(J)+dispomDW(J)*dRch4(J)*Sdwz(J))/(pom+pomdw)
        End Do
        Do J = 1, Nj
          tTRo2 = tTRo2 + (dispom(J)*Ro2(J)*Sz(J) + dispomDW(J)*dRo2(J)*Sdwz(J)  &
                          + Tfpom(J+1)*(1d0-BE(J))*Ro2SED(J)*Fz(J+1)             &
                          + TfpomDW(J+1)*(1d0-BEdw(J))*dRo2SED(J)*Fdwz(J+1))
          tTRno3 = tTRno3 + (dispom(J)*Rno3(J)*Sz(J) + dispomDW(J)*dRno3(J)*Sdwz(J)  &
                            + Tfpom(J+1)*(1d0-BE(J))*Rno3SED(J)*Fz(J+1)              &
                            + TfpomDW(J+1)*(1d0-BEdw(J))*dRno3SED(J)*Fdwz(J+1))
          tTRso4 = tTRso4 + (dispom(J)*Rso4(J)*Sz(J) + dispomDW(J)*dRso4(J)*Sdwz(J)  &
                            + Tfpom(J+1)*(1d0-BE(J))*Rso4SED(J)*Fz(J+1)              &
                            + TfpomDW(J+1)*(1d0-BEdw(J))*dRso4SED(J)*Fdwz(J+1))
          tTRch4 = tTRch4 + (dispom(J)*Rch4(J)*Sz(J) + dispomDW(J)*dRch4(J)*Sdwz(J)  &
                            + Tfpom(J+1)*(1d0-BE(J))*Rch4SED(J)*Fz(J+1)              &
                            + TfpomDW(J+1)*(1d0-BEdw(J))*dRch4SED(J)*Fdwz(J+1))
        End Do
        tTRo2  = tTRo2/(pom+pomdw+pomSED+pomdwSED)
        tTRno3 = tTRno3/(pom+pomdw+pomSED+pomdwSED)
        tTRso4 = tTRso4/(pom+pomdw+pomSED+pomdwSED)
        tTRch4 = tTRch4/(pom+pomdw+pomSED+pomdwSED)
        Write(*,*) '------------------------------'
        Write(*,'(A18,f7.3)') 'aerobic oxidation:',TRo2
        Write(*,'(A18,f7.3)') 'denitrification  :',TRno3
        Write(*,'(A18,f7.3)') 'sulfate reduction:',TRso4
        Write(*,'(A18,f7.3)') 'methanogenesis   :',TRch4
        Write(*,*) '----sediment----'
        Write(*,'(A18,f7.3)') 'aerobic oxidation:',tTRo2
        Write(*,'(A18,f7.3)') 'denitrification  :',tTRno3
        Write(*,'(A18,f7.3)') 'sulfate reduction:',tTRso4
        Write(*,'(A18,f7.3)') 'methanogenesis   :',tTRch4
        Write(*,*) 'Benthic MSR =',BSRsed/1d12,' Tmol/yr'
        Write(*,*) 'P deposiiton=',SdepoP/1d12,'TmolP/yr'
        Write(*,*) 'POM deposition=',(pom+pomdw)/1d12,'Tmol/yr'
        Write(*,*) 'POM decomposition@sed=',(pomSED+pomdwSED)/1d12,'Tmol/yr'
        Write(*,*) 'Porg burial=',TorgP/1d12
        Write(*,*) 'Fe-bound P burial=',Tfep/1d12
        Write(*,*) 'Ca-bound P burial=',Tcap/1d12
        Write(*,*) 'P efflux=',TreturnP/1d12
        Write(*,*) 'Deni@wc =',Tdeniwc/1d12
        Write(*,*) 'Deni@sed=',Tdeni/1d12
        Write(*,*) '------------------------------'
      End If

      RETURN
      End SUBROUTINE ORG

!=================================================================================================!
!=================================================================================================!

      SUBROUTINE Jacob(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o  &
                      ,CjO2,CjNH4,CjH2S,CjSO4,CjCH4,Nj,Ji,Kso4_aom,Nloop)
      Integer Ii,Nj,max,Ji,coln,rown,MAXmtx,III,JJJ,Ir,Ic,Nloop
      Parameter(max=1e6,rown=5,coln=6,MAXmtx=50)
      double precision a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q  &
                      ,CjO2(2,Nj),CjNH4(2,Nj),CjH2S(2,Nj),CjSO4(2,Nj),CjCH4(2,Nj),f1,f2,f3,f4,f5,f6  &
                      ,O2,NH4,H2S,CH4,SO4,NO3,dO2,dNH4,dH2S,dCH4,dSO4,dNO3  &
                      ,dmax,am(MAXmtx,MAXmtx),Kso4_aom,Nmax,Nmin &
                      ,const,const1,const2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      O2  = CjO2(1,Ji)
      CH4 = CjCH4(1,Ji)
      NH4 = CjNH4(1,Ji)
      H2S = CjH2S(1,Ji)
      SO4 = CjSO4(1,Ji)

      const  = 2d0*f
      const1 = 2d0*h
      const2 = 2d0*j

      Nloop = 0
      Do Ii = 1, max
        Nloop = Nloop + 1
        If(Ii == max) then
          Write(*,*) 'Ii = max !!: ',Ji,O2,CH4,NH4,H2S,SO4
          Go to 1
        End If
! f
        f1 = -O2  + a + const*O2*NH4 + const1*O2*H2S + const2*O2*CH4
        f2 = -NH4 + e + f*O2*NH4
        f3 = -H2S + g + h*O2*H2S - m*CH4*SO4/(SO4+Kso4_aom)
        f4 = -CH4 + i + j*O2*CH4 + m*CH4*SO4/(SO4+Kso4_aom)
        f5 = -SO4 + l + m*CH4*SO4/(SO4+Kso4_aom) - h*O2*H2S
! aij
        am(1,1) = -1d0 + const*NH4 + const1*H2S + const2*CH4
        am(1,2) = const*O2
        am(1,3) = const1*O2
        am(1,4) = const2*O2
        am(1,5) = 0d0
        am(1,6) = -f1
        am(2,1) = f*NH4
        am(2,2) = -1d0 + f*O2
        am(2,3) = 0d0
        am(2,4) = 0d0
        am(2,5) = 0d0
        am(2,6) = -f2
        am(3,1) = h*H2S
        am(3,2) = 0d0
        am(3,3) = -1d0 + h*O2
        am(3,4) = -m*SO4/(SO4+Kso4_aom)
        am(3,5) = -m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(3,6) = -f3
        am(4,1) = j*CH4
        am(4,2) = 0d0
        am(4,3) = 0d0
        am(4,4) = -1d0 + j*O2 + m*SO4/(SO4+Kso4_aom)
        am(4,5) = m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(4,6) = -f4
        am(5,1) = -h*H2S
        am(5,2) = 0d0
        am(5,3) = -h*O2
        am(5,4) = m*SO4/(SO4+Kso4_aom)
        am(5,5) = -1d0 + m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(5,6) = -f5
! Matrix solver
        CALL MATRIX(coln,rown,am,MAXmtx)
        dO2  = am(1,coln)
        dNH4 = am(2,coln)
        dH2S = am(3,coln)
        dCH4 = am(4,coln)
        dSO4 = am(5,coln)
905     Format(E13.4E3,E13.4E3,E13.4E3,E13.4E3,E13.4E3)
906     Format(E13.4E3,E13.4E3,E13.4E3,E13.4E3,E13.4E3,E13.4E3)

        dmax = f1*f1 + f2*f2 + f3*f3 + f4*f4 + f5*f5

        If(dmax <= 1d-20) then
          CjO2(2,Ji)  = O2
          CjCH4(2,Ji) = CH4
          CjNH4(2,Ji) = NH4
          CjH2S(2,Ji) = H2S
          CjSO4(2,Ji) = SO4
          Go to 1
        Else
          O2  = O2  + dO2
          NH4 = NH4 + dNH4
          H2S = H2S + dH2S
          CH4 = CH4 + dCH4
          SO4 = SO4 + dSO4
        End If
      End Do
1     Continue

      RETURN
      End SUBROUTINE Jacob

!=================================================================================================!
!=================================================================================================!

      SUBROUTINE JacobM(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o  &
                       ,CmO2,CmNH4,CmH2S,CmSO4,CmCH4,Nj,Kso4_aom,Nloop)
      Integer Ii,Nj,max,coln,rown,MAXmtx,III,JJJ,Ir,Ic,Nloop
      Parameter(max=1e6,rown=5,coln=6,MAXmtx=50)
      double precision a,b,c,d,e,f,g,h,i,j,k,l,m,n,o  &
                      ,CmO2(2),CmNH4(2),CmH2S(2),CmSO4(2),CmCH4(2),f1,f2,f3,f4,f5,f6  &
                      ,O2,NH4,H2S,CH4,SO4,NO3,dO2,dNH4,dH2S,dCH4,dSO4,dNO3  &
                      ,dmax,am(MAXmtx,MAXmtx),Kso4_aom,Nmax,Nmin &
                      ,const,const1,const2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      O2  = CmO2(1)
      CH4 = CmCH4(1)
      NH4 = CmNH4(1)
      H2S = CmH2S(1)
      SO4 = CmSO4(1)

      const  = 2d0*f
      const1 = 2d0*h
      const2 = 2d0*j

      Nloop = 0d0
      Do Ii = 1, max
        Nloop = Nloop + 1
        If(Ii == max) then
          Go to 111
        End If
        f1 = -O2  + a + const*O2*NH4 + const1*O2*H2S + const2*O2*CH4
        f2 = -NH4 + e + f*O2*NH4
        f3 = -H2S + g + h*O2*H2S - m*CH4*SO4/(SO4+Kso4_aom)
        f4 = -CH4 + i + j*O2*CH4 + m*CH4*SO4/(SO4+Kso4_aom)
        f5 = -SO4 + l + m*CH4*SO4/(SO4+Kso4_aom) - h*O2*H2S
! aij
        am(1,1) = -1d0 + const*NH4 + const1*H2S + const2*CH4
        am(1,2) = const*O2
        am(1,3) = const1*O2
        am(1,4) = const2*O2
        am(1,5) = 0d0
        am(1,6) = -f1
        am(2,1) = f*NH4
        am(2,2) = -1d0 + f*O2
        am(2,3) = 0d0
        am(2,4) = 0d0
        am(2,5) = 0d0
        am(2,6) = -f2
        am(3,1) = h*H2S
        am(3,2) = 0d0
        am(3,3) = -1d0 + h*O2
        am(3,4) = -m*SO4/(SO4+Kso4_aom)
        am(3,5) = -m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(3,6) = -f3
        am(4,1) = j*CH4
        am(4,2) = 0d0
        am(4,3) = 0d0
        am(4,4) = -1d0 + j*O2+m*SO4/(SO4+Kso4_aom)
        am(4,5) = m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(4,6) = -f4
        am(5,1) = -h*H2S
        am(5,2) = 0d0
        am(5,3) = -h*O2
        am(5,4) = m*SO4/(SO4+Kso4_aom)
        am(5,5) = -1d0 + m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(5,6) = -f5
! Matrix solver
        CALL MATRIX(coln,rown,am,MAXmtx)
        dO2  = am(1,coln)
        dNH4 = am(2,coln)
        dH2S = am(3,coln)
        dCH4 = am(4,coln)
        dSO4 = am(5,coln)
        dNO3 = am(6,coln)
905     Format(E13.4E3,E13.4E3,E13.4E3,E13.4E3,E13.4E3)
906     Format(E13.4E3,E13.4E3,E13.4E3,E13.4E3,E13.4E3,E13.4E3)

        dmax = f1*f1 + f2*f2 + f3*f3 + f4*f4 + f5*f5

        If(dmax <= 1d-20) then
          CmO2(2)  = O2
          CmCH4(2) = CH4
          CmNH4(2) = NH4
          CmH2S(2) = H2S
          CmSO4(2) = SO4
          Go to 111
        Else
          O2  = O2  + dO2
          NH4 = NH4 + dNH4
          H2S = H2S + dH2S
          CH4 = CH4 + dCH4
          SO4 = SO4 + dSO4
        End If
      End Do
111   Continue

      RETURN
      End SUBROUTINE JacobM

!=================================================================================================!
!=================================================================================================!

      SUBROUTINE JacobM2(a,c,d,g,h,i,j,k,l,m,n,o  &
                        ,CmO2,CmH2S,CmSO4,CmCH4,Nj,Kso4_aom,Nloop)
      Integer Ii,Nj,max,coln,rown,MAXmtx,III,JJJ,Ir,Ic,Nloop
      Parameter(max=1e6,rown=4,coln=5,MAXmtx=50)
      double precision a,c,d,g,h,i,j,k,l,m,n,o  &
                      ,CmO2(2),CmH2S(2),CmSO4(2),CmCH4(2),f1,f2,f3,f4  &
                      ,O2,H2S,CH4,SO4,NO3,dO2,dH2S,dCH4,dSO4  &
                      ,dmax,am(MAXmtx,MAXmtx),Kso4_aom,Nmax,Nmin &
                      ,const,const1,const2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      O2  = CmO2(1)
      CH4 = CmCH4(1)
      H2S = CmH2S(1)
      SO4 = CmSO4(1)

      const1 = 2d0*h
      const2 = 2d0*j

      Nloop = 0
      Do Ii = 1, max
        Nloop = Nloop + 1
        If(Ii == max) then
          Go to 121
        End If
        f1 = -O2  + a + const1*O2*H2S + const2*O2*CH4
        f2 = -H2S + g + h*O2*H2S - m*CH4*SO4/(SO4+Kso4_aom)
        f3 = -CH4 + i + j*O2*CH4 + m*CH4*SO4/(SO4+Kso4_aom)
        f4 = -SO4 + l + m*CH4*SO4/(SO4+Kso4_aom) - h*O2*H2S
! aij
        am(1,1) = -1d0 + const1*H2S + const2*CH4
        am(1,2) = const1*O2
        am(1,3) = const2*O2
        am(1,4) = 0d0
        am(1,5) = -f1
        am(2,1) = h*H2S
        am(2,2) = -1d0 + h*O2
        am(2,3) = -m*SO4/(SO4+Kso4_aom)
        am(2,4) = -m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(2,5) = -f2
        am(3,1) = j*CH4
        am(3,2) = 0d0
        am(3,3) = -1d0 + j*O2 + m*SO4/(SO4+Kso4_aom)
        am(3,4) = m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(3,5) = -f3
        am(4,1) = -h*H2S
        am(4,2) = -h*O2
        am(4,3) = m*SO4/(SO4+Kso4_aom)
        am(4,4) = -1d0 + m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(4,5) = -f4
! Matrix solver
        CALL MATRIX(coln,rown,am,MAXmtx)
        dO2  = am(1,coln)
        dH2S = am(2,coln)
        dCH4 = am(3,coln)
        dSO4 = am(4,coln)

        dmax = f1*f1 + f2*f2 + f3*f3 + f4*f4

        If(dmax <= 1d-20) then
          CmO2(2)  = O2
          CmCH4(2) = CH4
          CmH2S(2) = H2S
          CmSO4(2) = SO4
          Go to 121
        Else
          O2  = O2  + dO2
          H2S = H2S + dH2S
          CH4 = CH4 + dCH4
          SO4 = SO4 + dSO4
        End If
      End Do
121   Continue

      RETURN
      End SUBROUTINE JacobM2

!=================================================================================================!
!=================================================================================================!

      SUBROUTINE JacobMIX(e,f,g,h,i,j,k,l,m,n,o,CmO2,CmH2S,CmSO4,CmCH4,Kso4_aom,Nloop)
      Integer max,coln,rown,MAXmtx,III,JJJ,Nloop
      Parameter(max=1e6,rown=3,coln=4,MAXmtx=50)
      double precision e,f,g,h,i,j,k,l,m,n,o  &
                      ,CmO2(2),CmH2S(2),CmSO4(2),CmCH4(2),f1,f2,f3,O2,H2S,CH4,SO4  &
                      ,dH2S,dCH4,dSO4,dmax,am(MAXmtx,MAXmtx),Kso4_aom
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      O2  = CmO2(2)
      CH4 = CmCH4(1)
      H2S = CmH2S(1)
      SO4 = CmSO4(1)

      Nloop = 0
      Do Ii = 1, max
        Nloop = Nloop + 1
        If(Ii == max) then
          Go to 1928
        End If
        f1 = -H2S + g + h*O2*H2S - m*CH4*SO4/(SO4+Kso4_aom)
        f2 = -CH4 + i + j*O2*CH4 + m*CH4*SO4/(SO4+Kso4_aom)
        f3 = -SO4 + l + m*CH4*SO4/(SO4+Kso4_aom) - h*O2*H2S
        am(1,1) = -1d0 + h*O2
        am(1,2) = -m*SO4/(SO4+Kso4_aom)
        am(1,3) = -m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(1,4) = -f1
        am(2,1) = 0d0
        am(2,2) = -1d0 + j*O2+m*SO4/(SO4+Kso4_aom)
        am(2,3) = m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(2,4) = -f2
        am(3,1) = -h*O2
        am(3,2) = m*SO4/(SO4+Kso4_aom)
        am(3,3) = -1d0 + m*CH4*Kso4_aom/(SO4+Kso4_aom)**2d0
        am(3,4) = -f3

        CALL MATRIX(coln,rown,am,MAXmtx)
        dH2S = am(1,coln)
        dCH4 = am(2,coln)
        dSO4 = am(3,coln)

        dmax = f1*f1 + f2*f2 + f3*f3

        If(dmax <= 1d-20) then
          CmCH4(2) = CH4
          CmH2S(2) = H2S
          CmSO4(2) = SO4
          Go to 1928
        Else
          H2S = H2S + dH2S
          CH4 = CH4 + dCH4
          SO4 = SO4 + dSO4
        End If
      End Do
1928  Continue

      RETURN
      End SUBROUTINE JacobMIX

!=================================================================================================!
!=================================================================================================!

      SUBROUTINE MATRIX(coln,rown,am,MAXmtx)
      Integer I,J,K,MAXmtx,row2,II,JJ,coln,rown,colnall,minn
      double precision am(MAXmtx,MAXmtx),s
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      colnall = coln
! Upper triangular matrix transformation
! CALL gaussjordan(MAXmtx,am,rown,colnall)

      minn = Min(rown,colnall)
      Do I = 1, minn
! Row exchange
         CALL mtxrxcgmax(MAXmtx,am,rown,colnall,I,row2)
         Do J = 1, rown
           If (J /= I) then
             s = am(J,I)/am(I,I)
             Do K = 1, colnall
               am(J,K) = am(J,K) - s*am(I,K)
             End Do
           End If
         End Do
         s = am(I,I)
         Do K = 1, colnall
           am(I,K) = am(I,K)/s
         End Do
      End Do

      RETURN
      End SUBROUTINE MATRIX

!=================================================================================================!
!=================================================================================================!

      SUBROUTINE mtxrxcgmax(MAXmtx,am,rn,cn,c1,maxp)
      Integer cn,rn,c1,maxp,MAXmtx,MAXC,I
      double precision am(MAXmtx, MAXmtx),tr,max,ael
      Parameter(MAXC=20)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      max  = dabs(am(c1,c1))
      maxp = c1
      Do I = c1+1, rn
         ael = dabs(am(i,c1))
         If (max < ael) then
            max  = ael
            maxp = I
         End if
      End Do
      If (maxp > c1) then
         Do I = 1, cn
            tr = am(c1,I)
            am(c1,I)   = am(maxp,I)
            am(maxp,I) = tr
         End Do
      End If

      RETURN
      End SUBROUTINE mtxrxcgmax

!----------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------!
      Subroutine PHOTOCHEM(pO2,pCH4,kch4ox)
        USE Constants
        Implicit none
        double precision:: pO2(2),pCH4,phi_o2,phi_ch4
        double precision:: kch4ox_6,kch4ox_5,kch4ox_4,kch4ox_3,kch4ox_23,kch4ox
    
        phi_o2 = dlog10(pO2(1)/0.987d0)   ! 1 bar = 0.987 atm
        If(phi_o2 <= -9d0) then
          kch4ox_6  = c0l6 + c1l6*phi_o2 + c2l6*phi_o2*phi_o2 + c3l6*phi_o2**3d0     &
                           + c4l6*phi_o2**4d0 + c5l6*phi_o2**5d0 + c6l6*phi_o2**6d0
          kch4ox_5  = c0l5 + c1l5*phi_o2 + c2l5*phi_o2*phi_o2 + c3l5*phi_o2**3d0     &
                           + c4l5*phi_o2**4d0 + c5l5*phi_o2**5d0 + c6l5*phi_o2**6d0
          kch4ox_4  = c0l4 + c1l4*phi_o2 + c2l4*phi_o2*phi_o2 + c3l4*phi_o2**3d0     &
                           + c4l4*phi_o2**4d0 + c5l4*phi_o2**5d0 + c6l4*phi_o2**6d0
          kch4ox_3  = c0l3 + c1l3*phi_o2 + c2l3*phi_o2*phi_o2 + c3l3*phi_o2**3d0     &
                           + c4l3*phi_o2**4d0 + c5l3*phi_o2**5d0 + c6l3*phi_o2**6d0
          kch4ox_23 = c0l23 + c1l23*phi_o2 + c2l23*phi_o2*phi_o2 + c3l23*phi_o2**3d0 &
                            + c4l23*phi_o2**4d0 + c5l23*phi_o2**5d0 + c6l23*phi_o2**6d0
        Else If((phi_o2 > -9d0).and.(phi_o2 <= -3d0)) then
          kch4ox_6  = c0m6 + c1m6*phi_o2 + c2m6*phi_o2*phi_o2 + c3m6*phi_o2**3d0     &
                           + c4m6*phi_o2**4d0 + c5m6*phi_o2**5d0 + c6m6*phi_o2**6d0
          kch4ox_5  = c0m5 + c1m5*phi_o2 + c2m5*phi_o2*phi_o2 + c3m5*phi_o2**3d0     &
                           + c4m5*phi_o2**4d0 + c5m5*phi_o2**5d0 + c6m5*phi_o2**6d0
          kch4ox_4  = c0m4 + c1m4*phi_o2 + c2m4*phi_o2*phi_o2 + c3m4*phi_o2**3d0     &
                           + c4m4*phi_o2**4d0 + c5m4*phi_o2**5d0 + c6m4*phi_o2**6d0
          kch4ox_3  = c0m3 + c1m3*phi_o2 + c2m3*phi_o2*phi_o2 + c3m3*phi_o2**3d0     &
                           + c4m3*phi_o2**4d0 + c5m3*phi_o2**5d0 + c6m3*phi_o2**6d0
          kch4ox_23 = c0m23 + c1m23*phi_o2 + c2m23*phi_o2*phi_o2 + c3m23*phi_o2**3d0 &
                            + c4m23*phi_o2**4d0 + c5m23*phi_o2**5d0 + c6m23*phi_o2**6d0
        Else
          kch4ox_6  = c0h6 + c1h6*phi_o2 + c2h6*phi_o2*phi_o2 + c3h6*phi_o2**3d0
          kch4ox_5  = c0h5 + c1h5*phi_o2 + c2h5*phi_o2*phi_o2 + c3h5*phi_o2**3d0
          kch4ox_4  = c0h4 + c1h4*phi_o2 + c2h4*phi_o2*phi_o2 + c3h4*phi_o2**3d0
          kch4ox_3  = c0h3 + c1h3*phi_o2 + c2h3*phi_o2*phi_o2 + c3h3*phi_o2**3d0
          kch4ox_23 = c0h23 + c1h23*phi_o2 + c2h23*phi_o2*phi_o2 + c3h23*phi_o2**3d0
        End If
        phi_ch4 = dlog10(pCH4/0.987d0)
        If(phi_ch4 <= -6d0) then
          kch4ox = kch4ox_6
        Else If((phi_ch4 > -6d0).and.(phi_ch4 <= -5d0)) then
          kch4ox = kch4ox_6 + (kch4ox_5 - kch4ox_6) * (phi_ch4 + 6d0)
        Else If((phi_ch4 > -5d0).and.(phi_ch4 <= -4d0)) then
          kch4ox = kch4ox_5 + (kch4ox_4 - kch4ox_5) * (phi_ch4 + 5d0)
        Else If((phi_ch4 > -4d0).and.(phi_ch4 <= -3d0)) then
          kch4ox = kch4ox_4 + (kch4ox_3 - kch4ox_4) * (phi_ch4 + 4d0)
        Else If((phi_ch4 > -3d0).and.(phi_ch4 <= dlog10(2d-3))) then
          kch4ox = kch4ox_3 + (kch4ox_23 - kch4ox_3) * (phi_ch4 + 3d0) / (dlog10(2d-3)+3d0)
        Else
          kch4ox = kch4ox_23
        End If
        kch4ox = 10d0**kch4ox

        Return
      End Subroutine PHOTOCHEM
!=================================================================================================!
!=================================================================================================!

      SUBROUTINE OUTPT(xOUT,yOUT,Worg,Tcm,Tch,Sv,Ppo4,Ppo40,AreaCS,AreaCS0,pO2a,pO2a0,dSST,dSSS       &
                      ,Vpom,Pno3,Ano3,Po,Foh,Areah,ChPO4,ChNO3,ChO2,ChSO4                             &
                      ,Mh2s,Mo2,Mpo4,CjO2,CjH2S,hm,Temp,TempDW,Sm,Sh                                  &
                      ,CmPO4,CmNO3,CmO2,CmSO4,CmNH4,CmH2S                                             &
                      ,Dm14C,Cm14C,Dh14C,Ch14C,dz,CjPO4,CdwPO4,CjNO3,CdwNO3                           &
                      ,CdwO2,CjSO4,CdwSO4,CjNH4,CdwNH4,CdwH2S,Dj14C,Cj14C,Ddw14C,Cdw14C               &
                      ,gO2j,gPO4j,gNO3j,gH2Sj,gNH4j,gSO4j,Dg14Cj,Sz,Sdwz                              &
                      ,Fpom1,Fpom2,Fpom3,Fdwpom1,Fdwpom2,Fdwpom3,Fz,Fdwz                              &
                      ,Mno3,Mso4,Mnh4,TBCorg,TBCorgDW,BCorg,BCorgDW,RETURNC,RETURNCdw,RETURNP         &
                      ,BurP,SburP,cpr,BpOrg,BpCa,BpFe,BpOrgDW,BpCaDW,BpFeDW,Torgp,Tcap,Tfep,RETURNPdw &
                      ,Ia,Ie,Iredox,Nj,Vj,Vdwj,RI,RIdw,gRI,Pno3M,ChNH4,ChH2S,dispomSED,dispomdwSED    &
                      ,Ro2SED,Rno3SED,Rso4SED,dRo2SED,dRno3SED,dRso4SED                               &
                      ,Fdeni,FdeniDW,Tdeni,TdeniCS,Salj,SaljDW,TBcaso4,TBpyr                          &
                      ,CmCH4,ChCH4,CjCH4,CdwCH4,AOM,AOMdw                                             &
                      ,TRo2,TRno3,TRso4,TRch4,tTRo2,tTRno3,tTRso4,tTRch4,pCH4,TBpyrWC                 &
                      ,NPPl,NPPh,Vm,Vh,Salm,Salh,Mch4_ocn)
      Integer J,Nj,Ia,Ie,Iredox,gIa,gIe,gIredox,Iadw,Iedw,Idwredox
      double precision  xOUT,yOUT,dz,hm,Worg,Sv,Ppo4,Ppo40,AreaCS,AreaCS0,pO2a(2),pO2a0,pCH4            &
                       ,Vpom,Pno3,Ano3,Po,Foh,Areah,NPPl,NPPh,Vm,Vh                                     &
                       ,Mh2s,Mo2,Mpo4,Mno3,Mnh4,Mso4,Mch4_ocn                                           &
                       ,Tcm(2),Tch(2),Temp(2,Nj),TempDW(2,Nj)                                           &
                       ,Sm,Sh,Salj(2,Nj),SaljDW(2,Nj),Salm(2),Salh(2),dSST,dSSS                         &
                       ,CmPO4(2),CmNO3(2),CmO2(2),CmSO4(2),CmNH4(2),CmH2S(2),CmCH4(2)                   &
                       ,ChPO4(2),ChNO3(2),ChO2(2),ChSO4(2),ChNH4(2),ChH2S(2),ChCH4(2)                   &
                       ,Dm14C,Cm14C(2),Dh14C,Ch14C(2)                                                   &
                       ,CjPO4(2,Nj),CdwPO4(2,Nj),CjNO3(2,Nj),CdwNO3(2,Nj)                               &
                       ,CjO2(2,Nj) , CdwO2(2,Nj),CjSO4(2,Nj),CdwSO4(2,Nj)                               &
                       ,CjNH4(2,Nj),CdwNH4(2,Nj),CjH2S(2,Nj),CdwH2S(2,Nj),CjCH4(2,Nj),CdwCH4(2,Nj)      &
                       ,Dj14C(Nj),Cj14C(2,Nj),Ddw14C(Nj),Cdw14C(2,Nj)                                   &
                       ,gO2j(Nj),gPO4j(Nj),gNO3j(Nj),gH2Sj(Nj),gNH4j(Nj),gSO4j(Nj),gCH4j(Nj),Dg14Cj(Nj) &
                       ,Sz(Nj),Sdwz(Nj),Fpom1(Nj+1),Fpom2(Nj+1),Fpom3(Nj+1),Fdwpom1(Nj+1)               &
                       ,Fdwpom2(Nj+1),Fdwpom3(Nj+1),Fz(Nj+1),Fdwz(Nj+1)                                 &
                       ,TBCorg,TBCorgDW,BCorg(Nj),BCorgDW(Nj)                                           &
                       ,RETURNC(Nj),RETURNCdw(Nj),RETURNP(Nj),RETURNPdw(Nj),BurP(Nj+1),SburP,cpr(Nj)    &
                       ,BpOrg(Nj),BpCa(Nj),BpFe(Nj),BpOrgDW(Nj),BpCaDW(Nj),BpFeDW(Nj),Torgp,Tcap,Tfep   &
                       ,Vj(Nj),Vdwj(Nj),RI,RIdw,gRI,RIj,RIdwj,gRIj,sRI,sRIdw,gsRI,TVj,TVdwj,TV,Pno3M    &
                       ,dispomSED(Nj),dispomdwSED(Nj),Ro2SED(Nj),Rno3SED(Nj),Rso4SED(Nj)                &
                       ,dRo2SED(Nj),dRno3SED(Nj),dRso4SED(Nj)                                           &
                       ,Fdeni(Nj),FdeniDW(Nj),Tdeni,TdeniCS                                             &
                       ,TBcaso4,TBpyr,TBpyrWC,AOM(Nj),AOMdw(Nj)                                         &
                       ,TRo2,TRno3,TRso4,TRch4,tTRo2,tTRno3,tTRso4,tTRch4
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Respiration ! respiration.dat
      Write(36,810) xOUT,yOUT,TRo2,TRno3,TRso4,TRch4,tTRo2,tTRno3,tTRso4,tTRch4

! Anoxygenic Oxidation of Methane
      Do J = 1, Nj
      ! AOMj.dat
        Write(22,805) xOUT,yOUT,-(hm+dz/2d0)-(J-1d0)*dz,AOM(J),AOMdw(J)
      End Do

! S outputs, S_flux.dat
      Write(19,806) xOUT,yOUT,TBpyr/1d12,TBpyrWC/1d12,(TBpyr+TBpyrWC)/1d12,TBcaso4/1d12

!==========================!
! BENTHIC DENITRIFICATION  !
!==========================!
      Tdeni = 0d0
      Do J = 1, Nj
        Tdeni = Tdeni + 4d0/5d0*(Fdeni(J)*Fz(J+1)+FdeniDW(J)*Fdwz(J+1))
      End Do
      TdeniCS = 4d0/5d0*(Fdeni(1)*Fz(2)+FdeniDW(1)*Fdwz(2)) ! @ continental shelves
    ! Tdeni.dat
      Write(39,806) xOUT,yOUT,TdeniCS,Tdeni,TdeniCS*14d0/1d12,Tdeni*14d0/1d12
    ! denitrification.dat
      Do J = 1, Nj
        Write(38,809) xOUT,yOUT,-(J+1d0)*dz,Rno3SED(J),dRno3SED(J)          &
                     ,4d0/5d0*Fdeni(J),4d0/5d0*FdeniDW(J)                   &
                     ,4d0/5d0*Fdeni(J)*Fz(J+1),4d0/5d0*FdeniDW(J)*Fdwz(J+1)
      End Do
      Write(38,*) ' '

! TEMPERATURE: Tair.dat
      Write(13,804) xOUT,yOUT,Tcm(2),Tch(2)

! RIVERINE INPUT RATE: river.dat
      Write(14,805) xOUT,yOUT,Ppo4,Pno3,Ano3

! SETTINGS: setting.dat
      Write(98,807) Sv,Ppo4/Ppo40,AreaCS/AreaCS0,pO2a(2)/pO2a0,dSST,dSSS,Vpom

! ATMOSPHERE: pXa.dat
      Write(10,804) xOUT,yOUT,pO2a(2),pCH4

! Export production & NPP [GtC/yr], NPP.dat
      Write(11,808) xOUT,yOUT,Po*12d0/1d15,Foh*Areah*12d0/1d15,NPPl*12d0/1d15,NPPh*12d0/1d15 &
                   ,(Po+Foh*Areah)*12d0/1d15,(NPPl+NPPh)*12d0/1d15

!============================!
! OCEANIC VERTICAL PROFILES  !
!============================!
    ! Xj.dat
      Dm14C = (Cm14C(2)-1d0)*1000d0
      Dh14C = (Ch14C(2)-1d0)*1000d0
      Write(7,830) xOUT,yOUT,-hm/2d0,CmPO4(2),ChPO4(2),CmNO3(2),ChNO3(2),CmO2(2),ChO2(2)   &
                  ,CmSO4(2),ChSO4(2),CmNH4(2),ChNH4(2),CmH2S(2),ChH2S(2),CmCH4(2),ChCH4(2) &
                  ,Tcm(2),Tch(2),Salm(2),Salh(2),Dm14C,Dh14C                               &
                  ,(CmPO4(2)*Vm+ChPO4(2)*Vh)/(Vm+Vh),(CmNO3(2)*Vm+ChNO3(2)*Vh)/(Vm+Vh)     &
                  ,(CmO2(2)*Vm+ChO2(2)*Vh)/(Vm+Vh),(CmSO4(2)*Vm+ChSO4(2)*Vh)/(Vm+Vh)       &
                  ,(CmNH4(2)*Vm+ChNH4(2)*Vh)/(Vm+Vh),(CmH2S(2)*Vm+ChH2S(2)*Vh)/(Vm+Vh)     &
                  ,(CmCH4(2)*Vm+ChCH4(2)*Vh)/(Vm+Vh)
      Do J = 1, Nj
        Dj14C(J)  = (Cj14C(2,J)-1d0)*1000d0
        Ddw14C(J) = (Cdw14C(2,J)-1d0)*1000d0
      ! Global average
        gO2j(J)  = (CjO2(2,J) *Sz(J)+CdwO2(2,J) *Sdwz(J)) / (Sz(J)+Sdwz(J))
        gPO4j(J) = (CjPO4(2,J)*Sz(J)+CdwPO4(2,J)*Sdwz(J)) / (Sz(J)+Sdwz(J))
        gNO3j(J) = (CjNO3(2,J)*Sz(J)+CdwNO3(2,J)*Sdwz(J)) / (Sz(J)+Sdwz(J))
        gH2Sj(J) = (CjH2S(2,J)*Sz(J)+CdwH2S(2,J)*Sdwz(J)) / (Sz(J)+Sdwz(J))
        gNH4j(J) = (CjNH4(2,J)*Sz(J)+CdwNH4(2,J)*Sdwz(J)) / (Sz(J)+Sdwz(J))
        gSO4j(J) = (CjSO4(2,J)*Sz(J)+CdwSO4(2,J)*Sdwz(J)) / (Sz(J)+Sdwz(J))
        gCH4j(J) = (CjCH4(2,J)*Sz(J)+CdwCH4(2,J)*Sdwz(J)) / (Sz(J)+Sdwz(J))
        Dg14Cj(J)= (Dj14C(J)  *Sz(J)+Ddw14C(J)  *Sdwz(J)) / (Sz(J)+Sdwz(J))
        Write(7,830) xOUT,yOUT,-(hm+dz/2d0)-(J-1d0)*dz,CjPO4(2,J),CdwPO4(2,J),CjNO3(2,J),CdwNO3(2,J) &
                    ,CjO2(2,J),CdwO2(2,J),CjSO4(2,J),CdwSO4(2,J),CjNH4(2,J),CdwNH4(2,J),CjH2S(2,J)   &
                    ,CdwH2S(2,J),CjCH4(2,J),CdwCH4(2,J)                                              &
                    ,Temp(2,J),TempDW(2,J),Salj(2,J),SaljDW(2,J),Dj14C(J),Ddw14C(J)                  &
                    ,gPO4j(J),gNO3j(J),gO2j(J),gSO4j(J),gNH4j(J),gH2Sj(J),gCH4j(J)
      End Do
      Write(7,*) ' '
! @ surface, Xm.dat
      Write(8,816) xOUT,yOUT,CmPO4(2),ChPO4(2),CmNO3(2),ChNO3(2),CmO2(2),ChO2(2),CmSO4(2),ChSO4(2) &
                            ,CmNH4(2),ChNH4(2),CmH2S(2),ChH2S(2),CmCH4(2),ChCH4(2)
! @ intermediate depths, O2_H2S_IW.dat 
      Write(9,808) xOUT,yOUT,CjO2(2,1),CdwO2(2,1),CjO2(2,5),CjO2(2,10),CjO2(2,25),CjO2(2,35),CjO2(2,45),CjO2(2,55) &
                  ,CjH2S(2,1),CdwH2S(2,1),CjH2S(2,5),CjH2S(2,10),CjH2S(2,25),CjH2S(2,35),CjH2S(2,45),CjH2S(2,55)

! Total volume
      TV    = 0d0
      TVj   = 0d0
      TVdwj = 0d0
      Do J = 1, Nj
        TVj   = TVj   + Vj(J)
        TVdwj = TVdwj + Vdwj(J)
        TV    = TV + Vj(J) + Vdwj(J)
      End Do

! Redox Index
      sRI   = 0d0
      sRIdw = 0d0
      gsRI  = 0d0
      Do J = 1, Nj
        RIj   = 0d0
        RIdwj = 0d0
        If(CjO2(2,J) <= 1d-3) then
          RIj = 1d0
        End If
        If(CjH2S(2,J) >= 0.1d0) then
          RIj = 2d0
        End If
        If(CdwO2(2,J) <= 1d-3) then
          RIdwj = 1d0
        End If
        If(CdwH2S(2,J) >= 0.1d0) then
          RIdwj = 2d0
        End If
        sRI   = sRI   + RIj*Vj(J)
        sRIdw = sRIdw + RIdwj*Vdwj(J)
      End Do
      RI   = sRI/(2d0*TVj)
      RIdw = sRIdw/(2d0*TVdwj)
      gRI  = (sRI+sRIdw)/(2d0*TV)
    ! RI.dat
      Write(41,805) xOUT,yOUT,RI,RIdw,gRI

!---------!
!  Fpom   !
!---------!
    ! fpom.dat
      Do J = 1, Nj+1
        Write(24,804) -hm-(J-1d0)*dz, Fpom1(J)+Fpom2(J)+Fpom3(J), Fdwpom1(J)+Fdwpom2(J)+Fdwpom3(J)  &
                                    , (Fpom1(J)+Fpom2(J)+Fpom3(J))/(Fpom1(1)+Fpom2(1)+Fpom3(1))
      End Do
      Write(24,*) ' '
!
!========!
!  MASS  ! mass.dat
!========!
      Write(12,809) xOUT,yOUT,Mpo4,Mno3,Mo2,Mso4,Mnh4,Mh2s,Mch4_ocn

!====================!
!  ORGANIC C BURIAL  !
!====================!
!  Total burial rate, Bcorg.dat
      Write(26,805) xOUT,yOUT,TBCorg/1d12,TBCorgDW/1d12,(TBCorg+TBCorgDW)/1d12

!  Burial flux density, Bcorgj.dat
      Do J = 1, Nj
        Write(25,807) xOUT,yOUT,-hm-dz*J,BCorg(J),BCorgDW(J),BCorg(J)/Fz(J+1),BCorgDW(J)/Fdwz(J+1)
      End Do
      Write(25,*) ' '

!=================!
!  P BURIAL RATE  !
!=================!
!  Burial flux and flux density, Bpj.dat
      Do J = 2, Nj+1
        Write(27,807) xOUT,yOUT,-hm-dz*(J-1),(BpOrg(J-1)+BpCa(J-1)+BpFe(J-1))/1d12          &
                                            ,(BpOrgDW(J-1)+BpCaDW(J-1)+BpFeDW(J-1))/1d12    &
                                            ,(BpOrg(J-1)+BpCa(J-1)+BpFe(J-1))/Fz(J)         &
                                            ,(BpOrgDW(J-1)+BpCaDW(J-1)+BpFeDW(J-1))/Fdwz(J)
      End Do
      Write(27,*) ' '

!  Total burial rate
      Torgp = 0d0
      Tcap  = 0d0
      Tfep  = 0d0
      Do J = 1, Nj
        Tcap  = Tcap  + BpCa(J)  + BpCadw(J)
        Tfep  = Tfep  + BpFe(J)  + BpFedw(J)
        Torgp = Torgp + BpOrg(J) + BpOrgDW(J)
        If(RETURNP(J) < 0) then
          Tcap = Tcap + RETURNP(J)
        End If
        If(RETURNPdw(J) < 0) then
          Tcap = Tcap + RETURNPdw(J)
        End If
      End Do
    ! Bp.dat
      Write(28,806) xOUT,yOUT,Torgp/1d12,Tcap/1d12,Tfep/1d12,(Torgp+Tcap+Tfep)/1d12

!====================!
!  C/P BURIAL RATIO  !
!====================!
!  Corg/Porg  ! CorgPorgj.dat
      Do J = 1, Nj
        Write(31,805) xOUT,yOUT,-hm-dz*J,BCorg(J)/BpOrg(J),BCorgDW(J)/BpOrgDW(J)
      End Do
      Write(31,*) ' '

!  Corg/Preac  ! CorgPreaj.dat
      Do J = 2, Nj+1
        If(ReturnP(J-1) < 0)then
          Write(32,804) xOUT,yOUT,-hm-dz*(J-1),BCorg(J-1)/(BurP(J)-ReturnP(J-1))
        Else
          Write(32,804) xOUT,yOUT,-hm-dz*(J-1),BCorg(J-1)/BurP(J)
        End If
      End Do
      Write(32,*) ' '
    ! Global average, tCorgPrea.dat
      Write(33,803) xOUT,yOUT,(TBCorg+TBCorgDW)/SburP

!================!
!  BENTHIC FLUX  ! BenthicFlux.dat
!================!
      Do J = 1, Nj
        Write(34,811) xOUT,yOUT,-hm-dz*J,ReturnP(J),ReturnPdw(J),ReturnC(J),ReturnCdw(J) &
                     ,ReturnP(J)/Fz(J+1),ReturnPDW(J)/Fdwz(J+1),ReturnC(J)/Fz(J+1),ReturnCdw(J)/Fdwz(J+1)
      End Do
      Write(34,*) ' '

!==========!
!  Iredox  !
!==========!
      Ia = 0
      Ie = 0
      Do J = 1, Nj
        If(Ia == 0) then
          If(CjO2(2,J) < 0.001d0) then
            Ia = 1
          End If
        End If
        If(Ia == 1) then
          If(CjO2(2,J) >= 0.001d0) then
            Ia = 2
          End If
        End If
        If(Ie == 0) then
          If(CjH2S(2,J) > 0.1d0) then
            Ie = 1
          End If
        End If
        If(Ie == 1) then
          If(CjH2S(2,J) <= 0.1d0) then
            Ie = 2
          End If
        End If
      End Do
      Iadw = 0
      Iedw = 0
      Do J = 1, Nj
        If(Iadw == 0) then
          If(CdwO2(2,J) < 0.001d0) then
            Iadw = 1
          End If
        End If
        If(Iadw == 1) then
          If(CdwO2(2,J) >= 0.001d0) then
            Iadw = 2
          End If
        End If
        If(Iedw == 0) then
          If(CdwH2S(2,J) > 0.1d0) then
            Iedw = 1
          End If
        End If
        If(Iedw == 1) then
          If(CdwH2S(2,J) <= 0.1d0) then
            Iedw = 2
          End If
        End If
      End Do
      gIa = 0
      gIe = 0
      Do J = 1, Nj
        If(gIa == 0) then
          If(gO2j(J) < 0.001d0) then
            gIa = 1
          End If
        End If
        If(gIa == 1) then
          If(gO2j(J) >= 0.001d0) then
            gIa = 2
          End If
        End If
        If(gIe == 0) then
          If(gH2Sj(J) > 0.1d0) then
            gIe = 1
          End If
        End If
        If(gIe == 1) then
          If(gH2Sj(J) <= 0.1d0) then
            gIe = 2
          End If
        End If
      End Do
      If((Ia == 0).and.(Ie == 0)) then
        Iredox = 0
      ! Ir=0.dat
        Write(43,'(f7.4,f7.4,i5)') xOUT,yOUT,Iredox
      Else If((Ia == 2).and.(Ie == 0)) then
        Iredox = 1
      ! Ir=1.dat
        Write(44,'(f7.4,f7.4,i5)') xOUT,yOUT,Iredox
      Else If((Ia == 1).and.(Ie == 0)) then
        Iredox = 2
      ! Ir=2.dat
        Write(45,'(f7.4,f7.4,i5)') xOUT,yOUT,Iredox
      Else If((Ia == 1).and.(Ie == 2)) then
        Iredox = 3
      ! Ir=3.dat
        Write(46,'(f7.4,f7.4,i5)') xOUT,yOUT,Iredox
      Else If((Ia == 1).and.(Ie == 1)) then
        Iredox = 4
      ! Ir=4.dat
        Write(47,'(f7.4,f7.4,i5)') xOUT,yOUT,Iredox
      Else If((Ia == 2).and.(Ie == 2)) then
        Iredox = 5
      ! Ir=4.dat
        Write(48,'(f7.4,f7.4,i5)') xOUT,yOUT,Iredox
      Else
        Iredox = 100
      End If
      If((Iadw == 0).and.(Iedw == 0)) then
        Idwredox = 0
      Else If((Iadw == 2).and.(Iedw == 0)) then
        Idwredox = 1
      Else If((Iadw == 1).and.(Iedw == 0)) then
        Idwredox = 2
      Else If((Iadw == 1).and.(Iedw == 2)) then
        Idwredox = 3
      Else If((Iadw == 1).and.(Iedw == 1)) then
        Idwredox = 4
      Else If((Iadw == 2).and.(Iedw == 2)) then
        Idwredox = 5
      Else
        Idwredox = 100
      End If
      If((gIa == 0).and.(gIe == 0)) then
        gIredox = 0
      Else If((gIa == 2).and.(gIe == 0)) then
        gIredox = 1
      Else If((gIa == 1).and.(gIe == 0)) then
        gIredox = 2
      Else If((gIa == 1).and.(gIe == 2)) then
        gIredox = 3
      Else If((gIa == 1).and.(gIe == 1)) then
        gIredox = 4
      Else If((gIa == 2).and.(gIe == 2)) then
        gIredox = 5
      Else
        gIredox = 100
      End If
    ! Ir.dat
      Write(42,'(f7.4,f7.4,i5,i5,i5)') xOUT,yOUT,Iredox,Idwredox,gIredox

      803   Format(3E13.4E3)
      804   Format(4E13.4E3)
      805   Format(5E13.4E3)
      806   Format(6E13.4E3)
      807   Format(7E13.4E3)
      808   Format(8E13.4E3)
      809   Format(9E13.4E3)
      810   Format(10E13.4E3)
      811   Format(11E13.4E3)
      812   Format(12E13.4E3)
      813   Format(13E13.4E3)
      814   Format(14E13.4E3)
      815   Format(15E13.4E3)
      816   Format(16E13.4E3)
      817   Format(17E13.4E3)
      818   Format(18E13.4E3)
      819   Format(19E13.4E3)
      820   Format(20E13.4E3)
      821   Format(21E13.4E3)
      822   Format(22E13.4E3)
      823   Format(23E13.4E3)
      830   Format(30E13.4E3)

      RETURN
      End SUBROUTINE OUTPT
