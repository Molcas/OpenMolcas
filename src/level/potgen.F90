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
subroutine POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,NCN,CNN,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM,MMLR,CMM, &
                  NCMM,IVSR,IDSTT,RHOAB)

implicit none
integer NBOB
parameter(NBOB=20)
integer I, J, IBOB, IAN1, IAN2, IMN1, IMN2, IORD, IORDD, IPOTL, PPAR, QPAR, NCN, NSR, NLR, NPP, LNPT, NCMM, IVSR, LVSR, IDSTT, &
        KDER, MMLR(3)
!character*2 NAME1,NAME2
real*8 BETA, BINF, CNN, DSCM, REQ, RREF, VLIM, ZZ, ZP, ZQ, ZME, ULR, ULRe, RHOAB, DM(3), CMM(3), RM3, PVSR, PARM(4), XO(NPP), &
       VV(NPP), RM2(NPP), bTT(-1:2), cDS(-2:0), bDS(-2:0)
save IORD, IORDD
save BINF, ZME, ULR, ULRe
! Damping function parameters for printout .....
data bTT/2.44d0,2.78d0,3.126d0,3.471d0/
data bDS/3.3d0,3.69d0,3.95d0/
data cDS/0.423d0,0.40d0,0.39d0/
save bTT, bDS, cDS
! Electron mass, as per 2006 physical constants
data ZME/5.4857990943d-4/

write(6,*) 'ZME=',ZME,'bTT=',bTT ! Make them "referenced"
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
! Also make sure some of these variables are "used" if NCMM>4
if (NCMM > 4) then
  !write(6,*) 'potgen.f has the following at the beginning:'
  write(6,*) 'IAN1 = ',IAN1
  write(6,*) 'IMN1 = ',IMN1
  write(6,*) 'IAN2 = ',IAN2
  write(6,*) 'IMN2 = ',IMN2
  !write(6,*) 'CHARGE = ',CHARGE
  !write(6,*) 'NUMPOT = ',NUMPOT
  !write(6,*) 'RH = ',RH
  !write(6,*) 'RMIN = ',RMIN
  !write(6,*) 'PRV = ',PRV
  !write(6,*) 'ARV = ',ARV
  !write(6,*) 'EPS = ',EPS
  !write(6,*) 'NTP = ',NTP
  !write(6,*) 'LPPOT = ',LPPOT
  !write(6,*) 'IOMEG1(now OMEGA) = ',OMEGA
  !write(6,*) 'VLIM = ',VLIM
  write(6,*) 'IPOTL = ',IPOTL
  !write(6,*) 'PPAR = ',PPAR
  !write(6,*) 'QPAR = ',QPAR
  !write(6,*) 'NSR = ',NSR
  !write(6,*) 'NLR = ',NLR
  write(6,*) 'IBOB = ',IBOB
end if
!write(6,*) 'DSCM = ',DSCM
!write(6,*) 'REQ = ',REQ
!write(6,*) 'RREF = ',RREF
!write(6,*) 'NCMM = ',NCMM
!write(6,*) 'IVSR = ',IVSR
!write(6,*) 'IDSTT = ',IDSTT
!write(6,*) 'RHOAB = ',RHOAB
!write(6,*) 'MMLR = ',MMLR
!write(6,*) 'CMM = ',CMM
!write(6,*) 'PARM = ',PARM
!write(6,*) 'NLEV1 = ',NLEV1
!write(6,*) 'AUTO1 = ',AUTO1
!write(6,*) 'LCDC = ',LCDC
!write(6,*) 'LXPCT = ',LXPCT
!write(6,*) 'NJM = ',NJM
!write(6,*) 'JDJR = ',JDJR
!write(6,*) 'IWF = ',IWF
!write(6,*) 'LPRWF = ',LPRWF
! Use the RM2 dummy variable:
if (RM2(1) > 0) RM2(1) = RM2(2)
LNPT = 1
IORD = NLR
!=======================================================================
! Generate an MLR potential as per Dattani & Le Roy J.Mol.Spec. 2011
!=======================================================================
write(6,*) 'Beginning to process MLR potential!'
write(6,*) ''
!if(IPOTL == 4) then
if (LNPT > 0) then
  NCN = MMLR(1)
  CNN = CMM(1)
  ULRe = 0.d0
  KDER = 0
  if (KDER /= 0) write(6,*) KDER ! Make sure it's "referenced")
  call dampF(REQ,RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)
  do J=1,NCMM
    ULRe = ULRe+DM(J)*CMM(J)/REQ**MMLR(J)
  end do
  write(6,*) 'Finished calculating damping functions'
  !end if
  BINF = dlog(2.d0*DSCM/ULRe)
  write(6,602) NCN,PPAR,QPAR,DSCM,REQ
  write(6,607) PPAR,PPAR,QPAR,NSR,NLR,IORD+1,(PARM(J),J=1,IORD+1)
  write(6,613) RREF
  if (RHOAB > 0.d0) then
    PVSR = 0.5d0*IVSR
    if (IDSTT > 0) then
      PVSR = 0.5d0*IVSR
      write(6,664) RHOAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
    else
      LVSR = IVSR/2
      write(6,666) RHOAB,LVSR,bTT(LVSR)
    end if
  end if
  write(6,617) BINF,MMLR(1),CMM(1),MMLR(1)
  if (NCMM > 1) then
    do I=2,NCMM !Removed IF stmnt that prints C10 nicely
      write(6,619) MMLR(I),CMM(I),MMLR(I)
    end do
  end if
end if
!  Loop over distance array XO(I)
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING:
!write(6,*) 'PPAR=',PPAR
!write(6,*) 'REQ=',REQ
!do I=1,3
!  write(6,*) 'XO=',XO(I)
!end do
do I=1,NPP
  !write(6,*) 'Calculating radial variables.'
  ! (r^n - rx^n)/(r^n + rx^n) for n={p,q},x={eq,ref}:
  ZZ = (XO(i)**PPAR-REQ**PPAR)/(XO(i)**PPAR+REQ**PPAR)
  ZP = (XO(i)**PPAR-RREF**PPAR)/(XO(i)**PPAR+RREF**PPAR)
  ZQ = (XO(i)**QPAR-RREF**QPAR)/(XO(i)**QPAR+RREF**QPAR)
  BETA = 0.d0
  if (ZZ > 0) IORDD = NLR
  if (ZZ <= 0) IORDD = NSR
  do J=IORDD,0,-1
    BETA = BETA*ZQ+PARM(J+1)
  end do
  ! Calculate MLR exponent coefficient
  !write(6,*) 'Calculating MLR exponent coefficient.'
  BETA = BINF*ZP+(1.d0-ZP)*BETA
  ULR = 0.d0
  ! Calculate local value of uLR(r)
  if ((NCMM >= 3) .and. (MMLR(2) <= 0)) then
    RM3 = 1.d0/XO(I)**3
    write(6,*) RM3 !Make it "referenced"
  else
    ! IVSR gets corrupted so make sure it's -2.
    !IVSR = -2
    !write(6,*) 'IVSR=',IVSR
    if (RHOAB > 0.d0) call dampF(XO(I),RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)
    do J=1,NCMM
      ULR = ULR+DM(J)*CMM(J)/XO(I)**MMLR(J)
    end do
  end if
  BETA = (ULR/ULRe)*dexp(-BETA*ZZ)
  VV(I) = DSCM*(1.d0-BETA)**2-DSCM+VLIM
  !write(6,*) 'I=',I,'/',NPP
end do
! OPTIONALLY PRINT THESE VARIABLE WHEN DEBUGGING
!write(6,*) 'NPP=',NPP
!write(6,*) 'VLIM=',VLIM
!write(6,*) 'DSCM=',DSCM
!write(6,*) 'ZZ=',ZZ
!write(6,*) 'ULRe=',ULRe
!write(6,*) 'ULR=',ULR
!write(6,*) 'BETA=',BETA
!end if
! OPTIONALLY PRINT SOME V(R) VALUES WHEN DEBUGGING:
!write(6,*) 'Finished MLR generation. First/last V(R):'
!do I=1,3
! write(6,*) 'V(',I,')=',VV(I)
!end do
!write(6,*) 'V(                 20000)=',VV(20000)
!write(6,*) 'V(',NPP,')=',VV(NPP)

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
