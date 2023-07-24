!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine POTGEN(LNPT,NPP,VLIM,XO,RM2,VV,NCN,CNN,PPAR,QPAR,NSR,NLR,DSCM,REQ,RREF,PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)

use LEVEL_COMMON, only: bDS, bTT, cDS
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: LNPT, NCN
integer(kind=iwp), intent(in) :: NPP, PPAR, QPAR, NSR, NLR, MMLR(3), NCMM, IVSR, IDSTT
real(kind=wp), intent(in) :: VLIM, XO(NPP), DSCM, REQ, RREF, PARM(4), CMM(3), RHOAB
real(kind=wp), intent(inout) :: RM2(NPP)
real(kind=wp), intent(out) :: VV(NPP), CNN
integer(kind=iwp) :: I, IORD, IORDD, J, LVSR
real(kind=wp) :: BETA, DM(3), PVSR, ULR, ZP, ZQ, ZZ
real(kind=wp), save :: BINF, ULRe

! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!if (NCMM > 4) then
!  write(u6,*) 'potgen has the following at the beginning:'
!  write(u6,*) 'IAN1 = ',IAN1
!  write(u6,*) 'IMN1 = ',IMN1
!  write(u6,*) 'IAN2 = ',IAN2
!  write(u6,*) 'IMN2 = ',IMN2
!  write(u6,*) 'CHARGE = ',CHARGE
!  write(u6,*) 'NUMPOT = ',NUMPOT
!  write(u6,*) 'RH = ',RH
!  write(u6,*) 'RMIN = ',RMIN
!  write(u6,*) 'PRV = ',PRV
!  write(u6,*) 'ARV = ',ARV
!  write(u6,*) 'EPS = ',EPS
!  write(u6,*) 'NTP = ',NTP
!  write(u6,*) 'LPPOT = ',LPPOT
!  write(u6,*) 'IOMEG1(now OMEGA) = ',OMEGA
!  write(u6,*) 'VLIM = ',VLIM
!  write(u6,*) 'IPOTL = ',IPOTL
!  write(u6,*) 'PPAR = ',PPAR
!  write(u6,*) 'QPAR = ',QPAR
!  write(u6,*) 'NSR = ',NSR
!  write(u6,*) 'NLR = ',NLR
!  write(u6,*) 'IBOB = ',IBOB
!end if
!write(u6,*) 'DSCM = ',DSCM
!write(u6,*) 'REQ = ',REQ
!write(u6,*) 'RREF = ',RREF
!write(u6,*) 'NCMM = ',NCMM
!write(u6,*) 'IVSR = ',IVSR
!write(u6,*) 'IDSTT = ',IDSTT
!write(u6,*) 'RHOAB = ',RHOAB
!write(u6,*) 'MMLR = ',MMLR
!write(u6,*) 'CMM = ',CMM
!write(u6,*) 'PARM = ',PARM
!write(u6,*) 'NLEV1 = ',NLEV1
!write(u6,*) 'AUTO1 = ',AUTO1
!write(u6,*) 'LCDC = ',LCDC
!write(u6,*) 'LXPCT = ',LXPCT
!write(u6,*) 'NJM = ',NJM
!write(u6,*) 'JDJR = ',JDJR
!write(u6,*) 'IWF = ',IWF
!write(u6,*) 'LPRWF = ',LPRWF
! Use the RM2 dummy variable:
if (RM2(1) > 0) RM2(1) = RM2(2)
LNPT = 1
IORD = NLR
!=======================================================================
! Generate an MLR potential as per Dattani & Le Roy J.Mol.Spec. 2011
!=======================================================================
write(u6,*) 'Beginning to process MLR potential!'
write(u6,*) ''
!if (IPOTL == 4) then
if (LNPT > 0) then
  NCN = MMLR(1)
  CNN = CMM(1)
  call dampF(REQ,RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)
  ULRe = sum(DM(1:NCMM)*CMM(1:NCMM)/REQ**MMLR(1:NCMM))
  write(u6,*) 'Finished calculating damping functions'
  !end if
  BINF = log(Two*DSCM/ULRe)
  write(u6,602) NCN,PPAR,QPAR,DSCM,REQ
  write(u6,607) PPAR,PPAR,QPAR,NSR,NLR,IORD+1,PARM(1:IORD+1)
  write(u6,613) RREF
  if (RHOAB > Zero) then
    PVSR = Half*IVSR
    if (IDSTT > 0) then
      PVSR = Half*IVSR
      write(u6,664) RHOAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
    else
      LVSR = IVSR/2
      write(u6,666) RHOAB,LVSR,bTT(LVSR)
    end if
  end if
  write(u6,617) BINF,MMLR(1),CMM(1),MMLR(1)
  if (NCMM > 1) then
    do I=2,NCMM !Removed IF stmnt that prints C10 nicely
      write(u6,619) MMLR(I),CMM(I),MMLR(I)
    end do
  end if
end if
!  Loop over distance array XO(I)
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING:
!write(u6,*) 'PPAR=',PPAR
!write(u6,*) 'REQ=',REQ
!do I=1,3
!  write(u6,*) 'XO=',XO(I)
!end do
do I=1,NPP
  !write(u6,*) 'Calculating radial variables.'
  ! (r^n - rx^n)/(r^n + rx^n) for n={p,q},x={eq,ref}:
  ZZ = (XO(i)**PPAR-REQ**PPAR)/(XO(i)**PPAR+REQ**PPAR)
  ZP = (XO(i)**PPAR-RREF**PPAR)/(XO(i)**PPAR+RREF**PPAR)
  ZQ = (XO(i)**QPAR-RREF**QPAR)/(XO(i)**QPAR+RREF**QPAR)
  BETA = Zero
  if (ZZ > 0) then
    IORDD = NLR
  else
    IORDD = NSR
  end if
  do J=IORDD,0,-1
    BETA = BETA*ZQ+PARM(J+1)
  end do
  ! Calculate MLR exponent coefficient
  !write(u6,*) 'Calculating MLR exponent coefficient.'
  BETA = BINF*ZP+(One-ZP)*BETA
  ULR = Zero
  ! Calculate local value of uLR(r)
  if ((NCMM < 3) .or. (MMLR(2) > 0)) then
    ! IVSR gets corrupted so make sure it's -2.
    !IVSR = -2
    !write(u6,*) 'IVSR=',IVSR
    if (RHOAB > Zero) call dampF(XO(I),RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)
    ULR = sum(DM(1:NCMM)*CMM(1:NCMM)/XO(I)**MMLR(1:NCMM))
  end if
  BETA = (ULR/ULRe)*exp(-BETA*ZZ)
  VV(I) = DSCM*(One-BETA)**2-DSCM+VLIM
  !write(u6,*) 'I=',I,'/',NPP
end do
! OPTIONALLY PRINT THESE VARIABLE WHEN DEBUGGING
!write(u6,*) 'NPP=',NPP
!write(u6,*) 'VLIM=',VLIM
!write(u6,*) 'DSCM=',DSCM
!write(u6,*) 'ZZ=',ZZ
!write(u6,*) 'ULRe=',ULRe
!write(u6,*) 'ULR=',ULR
!write(u6,*) 'BETA=',BETA
!end if
! OPTIONALLY PRINT SOME V(R) VALUES WHEN DEBUGGING:
!write(u6,*) 'Finished MLR generation. First/last V(R):'
!do I=1,3
! write(u6,*) 'V(',I,')=',VV(I)
!end do
!write(u6,*) 'V(                 20000)=',VV(20000)
!write(u6,*) 'V(',NPP,')=',VV(NPP)

return

602 format(/' MLR(n=',i1,'; p=',I1,', q=',I1,') Potential with:   De=',F10.3,'[cm-1]    Re=',F12.8,'[A]')
607 format('   with exponent coefficient   beta(r)= beta{INF}*y',I1,' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x, &
           'exponent coefft. power series orders',I4,' for  R < Re  and',I4,' for  R > Re'/6x,'and',i3,' coefficients:',1PD16.8, &
           2d16.8:/(10x,4d16.8:))
613 format(6x,'with radial variables  y_p & y_q  defined w.r.t.  RREF=',F10.7)
!615 format(6x,'radial variables  y_p & y_q  defined w.r.t.  RREF= Re=' F10.7)
!614 format(/' Potential is an SPF expansion in  (r-Re)/(',F5.2,'* r)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',1PD16.9, &
!           '   and',i3,'  b_i  coefficients:'/(5D16.8))
!616 format(/' Potential is an O-T expansion in  (r-Re)/[',f5.2,'*(r+Re)]  with   Re=',f12.9/5x,'De=',G18.10,'   c0=',1PD16.9, &
!           '   and',i3,'  c_i coefficients:'/(5D16.8))
617 format('   while betaINF=',f12.8,'  & uLR defined by  C',i1,' =',1PD13.6,'[cm-1 Ang','^',0P,I1,']')
!618 format(/' Potential is a general GPEF expansion in  (r**',i1,' - Re**',i1,')/(',SP,F5.2,'*r**',SS,i1,SP,F6.2,'*Re**',SS,i1, &
!           ')'/5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
619 format(50x,'C',I1,' =',1PD13.6,'[cm-1 Ang','^',0P,I1,']')
!621 format(50x,'C',I2,'=',1PD13.6,'[cm-1 Ang','^',0P,I2,']')
!620 format(/' Potential is a power series in  r  of  order',i3,' with   V(r=0)=',f11.4/3x,'& coefficients (from linear term):', &
!           1P2d16.8:/(5x,4D16.8:))
!626 format('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]   and'/'     Damping function  D(r)= exp[ -',0Pf6.4,'*(',f7.4, &
!           '/X -1.0)**',f5.2,']')
!636 format(3x,'where   fsw(r) = 1/[1 - exp{',f7.4,'*(r -',f7.4,')}]')
!642 format(' where for  r < Rinn=',F7.4,'   V=',SP,F12.4,1x,1PD13.6,'/R**12' )
!644 format('  and  for  r > Rout=',F7.3,'   V= VLIM ',(SP,1PD14.6,'/r**',SS,I2):/(39x,SP,1PD14.6,'/r**',SS,I2))
!652 format(6x,'where the radial variable   y_',I1,'= (r**',I1,' - RREF**',i1,')/(r**',I1,' + RREF**',i1, ')')
!654 format(10x,'is defined w.r.t.   RREF=',F11.8)
!656 format(10x,'is defined w.r.t.   RREF= Re= ',F11.8)
664 format(4x,'uLR inverse-power terms incorporate DS-type damping with   RHOAB=',f10.7/8x, &
           'defined to give very short-range  Dm(r)*Cm/r^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2, &
           '(RHOAB*r)/m -',f6.3,'(RHOAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
666 format(4x,'uLR inverse-power terms incorporate TT-type damping with   RHOAB=',f10.7/8x, &
           'defined to give very short-range  Dm(r)*Cm/r^m  behaviour   r^{',I2,'}'/8x, &
           'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r)^k/k!}]   where   bTT=',f6.3,'*RHOAB')

end subroutine POTGEN
