!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine SIZES()

use Index_Functions, only: nTri_Elem, nTri3_Elem
use Symmetry_Info, only: Mul
use PrintLevel, only: USUAL
use SUPERINDEX, only: SUPINI
use Cholesky, only: NumCho
use EQSOLV, only: IFCOUP, NLIST, NLSTOT
use caspt2_global, only: do_csf, do_grad, do_nac, if_invar, if_invaria, if_SSDM, ipea_shift, iPrGlb, NCMO
use general_data, only: nActel, nAsh
use caspt2_module, only: IfChol, IfDW, IfMSCoup, IfProp, IfRMS, IfXMS, MxCI, nAshT, nASup, nBasT, nBSqT, nBTri, &
                         nCases, nConf, nFroT, nG1, nG2, nG3Tot, nIMx, nInDep, nISh, nIshT, nISup, nOrb, nOSqT, nOTri, nSMx, nSsh, &
                         nState, nSym, nTGTU, nTU, Zeta
#ifdef _DEBUGPRINT_
use caspt2_module, only: NAMX
#endif
use stdalloc, only: mma_MaxDBLE
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "warnings.h"
integer(kind=iwp) :: ICASE, ICASE1, ICASE2, ISL1, ISL2, ISL3, ISYM, ISYM1, ISYM2, M, M11, M12, M21, M22, M31, MAXAIS, MC1S1DER, &
                     MC2DER, MEMBASE, MINBUF, MMX, MXLeft, MXLFT, N, NAS, NAS1, NAS2, NBOTTOM, NCH, nCLag, NCX, NDD, NEED, NEED0, &
                     NFIA, NFIT, NFTA, NG02, NG10, NG12, NG20, NG30, NGARR, ngrad, ngrad1, ngrad10, ngrad11, ngrad11_1, ngrad11_2, &
                     ngrad2, ngrad3, ngrad4, ngrad5, ngrad6, ngrad6_1, ngrad6_2, ngrad7, ngrad7_1, ngrad7_2, ngrad8, ngrad9, NIN, &
                     NIN1, NIS, NIS1, NIS2, NMKRHS, NMX, nOLag, nOMax, NPLBUF, NPoly, nPrp, nPrp1, nPrp2, nRHSP, nSgm1, nSgm2, &
                     nSigma, nSigma_inner, nSigma_outer, NSLag, nTG1, nTG2, nTG3, NumChT, nV, nVCUtil, nWLag, NX
real(kind=wp) :: XX, YY
integer(kind=iwp), external :: iParDiv
integer(kind=iwp), parameter :: Magic = 119000

! Available workspace right now:
call mma_MaxDBLE(MXLEFT)
! 250000 words margin, for various purposes.
NBOTTOM = 250000
! POLY package:
! MINBUF=Min acceptable nr of buffers in POLY3 for efficiency.
MINBUF = 1
NPOLY = NBOTTOM
! Sizes of G10,G20,G02,G12,G30: Faked, since no more used.
NG10 = 1
NG20 = 1
NG02 = 1
NG30 = 1
NG12 = 1

if (NACTEL > 0) then
  NGARR = 2*(NG10+NG20+NG02+NG30+NG12)
  NPOLY = NPOLY+NGARR+5*MXCI
  MXLFT = MXLEFT-NPOLY
  XX = Half*real(3*MXCI+3,kind=wp)
  YY = real(MXLFT,kind=wp)
  ! Allowed by memory:
  NPLBUF = int(sqrt(YY+XX**2)-XX)
  ! Preferred size for efficiency:
  NPLBUF = max(MINBUF,NPLBUF)
  ! Actual max needed in-core:
  NPLBUF = min(NPLBUF,nTri_Elem(NASHT))
  NPOLY = NPOLY+NPLBUF*(NPLBUF+3+3*MXCI)
  !write(u6,*) ' Memory requirements for POLY3 (Above SGUGA):'
  !write(u6,*) '   CI vector                      :',MXCI
  !write(u6,*) '   Extra margin for small scratch :',NBOTTOM
  !write(u6,*) '   Ten arrays G10,G20..F12        :',NGARR
  !write(u6,*) '   H0,SGM0,C1,C2                  :',4*MXCI
  !write(u6,*) '   VXYZB,TUBUF,DTUB               :',3*NPLBUF
  !write(u6,*) '   A,B1,B2                        :',3*NPLBUF*MXCI
  !write(u6,*) '   SCR                            :',NPLBUF**2
  !write(u6,*)
  !write(u6,*) ' Total, for POLY3                 :',NPOLY
  !write(u6,*)
end if

! Precompute sizes, offsets etc.
call SUPINI()

! SBMAT need:
!NG3C = 0
!do ISYM=1,NSYM
!  N = NTUV(ISYM)
!  NG3C = NG3C+nTri_Elem(N)
!end do
!NG3C = iPARDIV(NG3TOT,NG2)

! Sizes and addresses to lists:
do ISL1=1,NSYM
  do ISL3=1,NSYM
    ISL2 = Mul(ISL1,ISL3)
    NLIST(ISL1,ISL3,1) = NASH(ISL2)*NTU(ISL3)*2
    NLIST(ISL1,ISL3,2) = NLIST(ISL1,ISL3,1)
    NLIST(ISL1,ISL3,3) = NASH(ISL2)*NTU(ISL3)
    NLIST(ISL1,ISL3,4) = NASH(ISL2)*NTGTU(ISL3)*2
    NLIST(ISL1,ISL3,5) = NLIST(ISL1,ISL3,3)
    NLIST(ISL1,ISL3,6) = NLIST(ISL1,ISL3,4)
    NLIST(ISL1,ISL3,7) = NASH(ISL2)*NASH(ISL3)*2
    NLIST(ISL1,ISL3,8) = NLIST(ISL1,ISL3,7)
    NLIST(ISL1,ISL3,9) = NASH(ISL2)*NASH(ISL3)
    NLIST(ISL1,ISL3,10) = NLIST(ISL1,ISL3,9)
    if (ISL1 == 1) NLIST(ISL1,ISL3,10) = NLIST(ISL1,ISL3,10)-NASH(ISL2)
    NLIST(ISL1,ISL3,12) = NASH(ISL1)*NASH(ISL2)
    NLIST(ISL1,ISL3,13) = NLIST(ISL1,ISL3,12)
    if (ISL3 == 1) NLIST(ISL1,ISL3,13) = NLIST(ISL1,ISL3,13)-NASH(ISL1)
    NLIST(ISL1,ISL3,11) = NLIST(ISL1,ISL3,12)
    if (ISL3 == 1) NLIST(ISL1,ISL3,11) = NLIST(ISL1,ISL3,11)+NASH(ISL1)*NASHT
    NLIST(ISL1,ISL3,14) = NISH(ISL1)*NISH(ISL2)
    NLIST(ISL1,ISL3,15) = NLIST(ISL1,ISL3,14)
    if (ISL3 == 1) NLIST(ISL1,ISL3,15) = NLIST(ISL1,ISL3,15)-NISH(ISL1)
    NLIST(ISL1,ISL3,16) = NSSH(ISL1)*NSSH(ISL2)
    NLIST(ISL1,ISL3,17) = NLIST(ISL1,ISL3,16)
    if (ISL3 == 1) NLIST(ISL1,ISL3,17) = NLIST(ISL1,ISL3,17)-NSSH(ISL1)
  end do
end do
NLSTOT = 4*sum(NLIST(1:NSYM,1:NSYM,:))

! maximum orbitals in an irrep
NOMAX = maxval(NORB(1:NSYM))
if (.not. IFCHOL) then
  ! MKRHS needs:
  NMKRHS = 2*NOMAX**2
  NMX = 0
  do ICASE=1,13
    do ISYM=1,NSYM
      NIN = NINDEP(ISYM,ICASE)
      if (NIN == 0) cycle
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      NV = NAS*NIS
      if (NV == 0) cycle
      if (ICASE > 11) then
        N = NV
      else
        N = 2*NV+nTri_Elem(NAS)
      end if
      NMX = max(N,NMX)
    end do
  end do
  NMKRHS = NMKRHS+NMX
  !if (MAXIT == 0) then
  !  NSIGMA = 0
  !else
  !  NSIGMA = NBOTTOM+NLSTOT+MMX
  !end if
else
  ! RHSALL2 and ADDRHS needs:
  NMKRHS = maxval(NASUP(1:NSYM,:)*NISUP(1:NSYM,:))
  NMKRHS = 2*iPARDIV(NMKRHS,2*NOMAX**2)
end if
NMKRHS = NMKRHS+NBOTTOM

! PCG/(new)SIGMA routine needs: twice a global array RHS size for
! transformations (overrides NSIGMA computed above)
NSIGMA_OUTER = 0
NSIGMA_INNER = 0
NVCUTIL = 0
do ICASE=1,13
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NRHSP = iPARDIV(NAS*NIS,0)
    if (ICASE > 11) then
      NSIGMA_INNER = max(NSIGMA_INNER,NRHSP)
    else
      NSIGMA_INNER = max(NSIGMA_INNER,NAS*NIS+NRHSP)
      NSIGMA_OUTER = max(NSIGMA_OUTER,NAS*NIS+NRHSP)
    end if
    NVCUTIL = max(NVCUTIL,2*NRHSP)
  end do
end do
NSIGMA = max(NSIGMA_INNER+NSIGMA_OUTER,NVCUTIL)+NBOTTOM

! PRPCTL needs:
! In DIADNS alone, NDD words are needed:
#ifdef _DEBUGPRINT_
write(u6,*) ' Memory requirements for PRPCTL (Above SGUGA).'
write(u6,*)
write(u6,*) ' PRP1) First phase of PRPCTL.'
write(u6,*)
write(u6,*) '    A) DIADNS alone, broken up in case/symm:'
#endif
MMX = 0
do ICASE=1,13
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NV = NAS*NIS
    if (NV == 0) cycle
    if (ICASE > 11) then
      M = 2*NV
    else
      M = 3*NV+nTri_Elem(NAS)
    end if
    MMX = max(M,MMX)
  end do
end do

NDD = 0
do ICASE=1,13
  do ISYM=1,NSYM
    NX = 0
    NIN = NINDEP(ISYM,ICASE)
    if (NIN > 0) then
      NIS = NISUP(ISYM,ICASE)
      if (NIS > 0) then
        NAS = NASUP(ISYM,ICASE)
#       ifdef _DEBUGPRINT_
        write(u6,*) ' Case, Symm:',ICASE,ISYM
        write(u6,*) ' NIN,NAS,NIS:',NIN,NAS,NIS
        write(u6,*) ' NIMX,NSMX:',NIMX,NAMX
#       endif
        if ((ICASE == 2) .or. (ICASE == 3)) NX = NIN*NIMX**2
        if ((ICASE == 6) .or. (ICASE == 7)) NX = NIN*NSMX*NIMX**2
        if ((ICASE == 8) .or. (ICASE == 9)) NX = NIN*NSMX**2
        if ((ICASE == 10) .or. (ICASE == 11)) NX = NIN*NIMX*NSMX**2
        if ((ICASE == 12) .or. (ICASE == 13)) NX = max(NAS*NIMX**2,NIS*NSMX**2)
      end if
    end if
#   ifdef _DEBUGPRINT_
    write(u6,'(1x,a,2i4,5x,i8)') '      Case, symm:',ICASE,ISYM,2*NX
#   endif
    NDD = max(2*NX,NDD)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) '       DIADNS needs the maximum, or NDD=',NDD
#endif
NCMO = NBSQT
NPRP1 = NBOTTOM+NCMO+notri+NLSTOT+2*NOSQT+MMX+NDD

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) '    B) Also needed, for 1st phase of PRPCTL:'
write(u6,'(1x,a,i8)') '       NBOTTOM:',NBOTTOM
write(u6,'(1x,a,i8)') '       NCMO   :',NCMO
write(u6,'(1x,a,i8)') '       notri :',notri
write(u6,'(1x,a,i8)') '       NLSTOT :',NLSTOT
write(u6,'(1x,a,i8)') '     2*NOSQT  :',2*NOSQT
write(u6,'(1x,a,i8)') '       MMX    :',MMX
#endif
NPRP2 = NBOTTOM+2*NCMO+notri+NBAST*(NBAST+1)

#ifdef _DEBUGPRINT_
write(u6,*) ' PRP2) Second phase of PRPCTL.'
write(u6,*)
write(u6,'(1x,a,i8)') '       NBOTTOM:',NBOTTOM
write(u6,'(1x,a,i8)') '     2*NCMO   :',2*NCMO
write(u6,'(1x,a,i8)') '       notri :',notri
write(u6,'(1x,a,i8)') 'NBAST*(NBAST+1)',NBAST*(NBAST+1)
#endif

! (rough) memory estimation for gradients
! Cholesky vectors try to use all the available memory, so the
! memory used for them are not counted (as in other steps)
ngrad = 0
if (do_grad) then
  ! memory in caspt2_grad
  nCLag = nConf*nState
  nOLag = NBSQT
  NSLag = nState*nState
  nWLag = NBTRI
  ngrad1 = NBSQT*4+nCLag*2+nOLag*2+nSLag*2+nWLag+NBSQT*2+nAshT**2
  if (IFXMS .or. IFRMS) ngrad1 = ngrad1+NBSQT
  if (IFDW .and. (zeta >= Zero)) ngrad1 = ngrad1+nState
  if (do_nac) ngrad1 = ngrad1+NBSQT
  if (nFroT /= 0) ngrad1 = ngrad1+nFroT**2

  ! "base" memory in dens
  nch = 0
  !if (IfChol) nch = nvloc_chobatch(1)
  ngrad2 = NOSQT*2+NBSQT*9+2*max(nbast**2,nch)+2*nAshT**2
  if ((nFroT /= 0) .or. (.not. if_invaria)) ngrad2 = ngrad2+2*NBSQT
  if (do_csf) ngrad2 = ngrad2+NBSQT
  ngrad2 = ngrad2+nAshT**2
  if (.not. if_invaria) ngrad2 = ngrad2+nAshT**2

  membase = MXLEFT-ngrad1-ngrad2-NBOTTOM

  ! gradient: residual
  ngrad3 = NSIGMA-NBOTTOM

  ! gradient: density (trdns2o)
  ngrad4 = NDD*3/2

  ! gradient: off-diagonal derivatives (sigder)
  NFIT = sum(NASH(1:NSYM)*NISH(1:NSYM))+1
  NFIA = sum(NSSH(1:NSYM)*NISH(1:NSYM))+1
  NFTA = sum(NSSH(1:NSYM)*NASH(1:NSYM))+1

  MMX = 0
  do ICASE1=1,11
    do ISYM1=1,NSYM
      NIS1 = NISUP(ISYM1,ICASE1)
      NAS1 = NASUP(ISYM1,ICASE1)
      NIN1 = NINDEP(ISYM1,ICASE1)
      NSGM2 = NIS1*NAS1
      NSGM1 = 1
      if (ICASE1 == 1) then
        NSGM1 = NASH(ISYM1)*NISH(ISYM1)
      else if (ICASE1 == 4) then
        NSGM1 = NASH(ISYM1)*NSSH(ISYM1)
      else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
        NSGM1 = NIS1
      end if
      M11 = NSGM2+NSGM1
      M12 = 2*iPARDIV(NSGM2,0)+NSGM2
      if (ICASE1 <= 11) then
        M12 = M12+iPARDIV(max(NSGM2,nTri_Elem(NAS1)),0)
      else
        M12 = M12+iPARDIV(NSGM2,0)
      end if
      do ICASE2=ICASE1+1,NCASES
        if (IFCOUP(ICASE2,ICASE1) == 0) cycle
        do ISYM2=1,NSYM
          NIS2 = NISUP(ISYM2,ICASE2)
          NAS2 = NASUP(ISYM2,ICASE2)
          NCX = NIS2*NAS2
          ! first loop
          M21 = M11+2*iPARDIV(NCX,0)
          M31 = M11-NSGM1+iPARDIV(NCX,0)
          if (iCASE1 <= 11) then
            ! C1S1DER
            MC1S1DER = iPARDIV(NAS1*NAS1+NIN1*NIS1+NAS1*NIN1,0)
            M31 = M31+iPARDIV(NSGM2,0)+max(MC1S1DER+NAS1*NAS1,iPARDIV(NSGM2,0))
          end if

          ! second loop
          M22 = NSGM1+iPARDIV(NCX,0)
          if (ICASE2 <= 11) then
            ! C2DER
            MC2DER = iPARDIV(NAS2*(2*NIS2+max(NAS2,NIS2)),0)
            M22 = M22+NAS2*NAS2+MC2DER
          end if
          M = max(M21,M31,M12,M22)
          MMX = max(M,MMX)
        end do
      end do
    end do
  end do
  ngrad5 = NFIT*2+NFIA*2+NFTA*2+MMX

  ! gradient: CI derivatives (clagx)
  ngrad6 = NG1*4+NG2*4+NG3TOT*3
  ! CLagD (hmm, not for the A and C subspaces in parallel)
  MMX = 0
  do ICASE=1,11
    do ISYM=1,NSYM
      NIS = NISUP(ISYM,ICASE)
      NAS = NASUP(ISYM,ICASE)
      NIN = NINDEP(ISYM,ICASE)
      M = 2*NAS**2+3*NIN*NIS+NAS*NIS
      if (IFMSCOUP) M = M+NIN*NIS
      ! inside CLagDX
      M = M+2*NAS**2+NAS*min(NAS,NIS)+NAS*NIN+NIN
      select case (ICASE)
        case (1)
          M = M+nTri_Elem(NAS)+NG3TOT
        case (2,3)
          M = M+2*NASHT**4
          if (IPEA_SHIFT /= Zero) M = M+nTri_Elem(NAS)
        case (4)
          M = M+nTri_Elem(NAS)+NG3TOT
        case (5)
          if (IPEA_SHIFT /= Zero) M = M+nTri_Elem(NAS)
        case (6,7)
          if (IPEA_SHIFT /= Zero) M = M+nTri_Elem(NAS)
        case (8,9)
          if (IPEA_SHIFT /= Zero) M = M+nTri_Elem(NAS)
        case (10,11)
          if (IPEA_SHIFT /= Zero) M = M+nTri_Elem(NAS)
      end select
      MMX = max(M,MMX)
    end do
  end do
  ngrad6_1 = MMX
  ! CnstCLag (derfg3)
  ngrad6_2 = NG3TOT+NCONF
  MXLFT = membase-ngrad6-ngrad6_2
  MXLFT = min(MXCI*(3*NASHT**2+NASHT),MXLFT)
  ngrad6_2 = ngrad6_2+MXLFT
  ngrad6 = ngrad6+max(ngrad6_1,ngrad6_2)

  ! gradient: effective Hamiltonian (DerHEff)
  NTG1 = NASHT**2
  NTG2 = NASHT**4
  NTG3 = nTri3_Elem(NTG1)
  ngrad7 = NTG1+NTG2+NTG3
  ! DerHeffX
  MAXAIS = 0
  do ICASE=1,13
    do ISYM=1,NSYM
      NIS = NISUP(ISYM,ICASE)
      NAS = NASUP(ISYM,ICASE)
      NIN = NINDEP(ISYM,ICASE)
      M = 2*iPARDIV(NAS*NIS,0)
      MAXAIS = max(M,MAXAIS)
    end do
  end do
  ngrad7_1 = MAXAIS
  ! DERTG3
  ngrad7_2 = MXCI*3+2*NASHT**2
  MXLFT = membase-ngrad7-ngrad7_2
  MXLFT = min(MXCI*(2*NASHT**2+4),MXLFT)
  ngrad7_2 = ngrad7_2+MXLFT
  ngrad7 = ngrad7+max(ngrad7_1,ngrad7_2)

  ! gradient: iterative CI derivatives (DEPSAOffC)
  ngrad8 = 0
  if (.not. if_invar) ngrad8 = nConf*nState*5+nConf+nState**3+nAshT**2+nAshT**4

  ! gradient: kappa (orbital) derivatives
  !           (OLagNS_RI, OLagNS2, OLagVVVO)
  NumChT = 0
  if (IFCHOL) then
    call Get_iArray('NumCho',NumCho,nSym)
    NumChT = sum(NumCho(1:nSym))
    ngrad9 = NumChT**2+1
  else
    ngrad9 = (NISHT+NASHT)**2*NBAST**2+1
  end if

  ! gradient: state-specific density matrix
  ngrad10 = 0
  ! IFCHOL is always true
  if (if_SSDM) ngrad10 = NCONF+NumChT**2+2*nBasT**2+2*NumChT+NBSQT

  ! memory in XMS (XMS_Grad)
  ngrad11 = 0
  if (IFMSCOUP) then
    ngrad11_1 = NCONF*4+2*NASHT**2+2*NASHT**4+NASHT**6
    ngrad11_2 = 3*NASHT**2+NBSQT*4+NCONF
    ngrad11 = max(ngrad11_1,ngrad11_2)
  end if

  !! NPRP1 corresponds to trdns2d
  ngrad = NBOTTOM+ngrad1+ngrad2+max(ngrad3,ngrad4,ngrad5,ngrad6,ngrad7,ngrad8,ngrad9,ngrad10,ngrad11,NPRP1-NBOTTOM)

# ifdef _DEBUGPRINT_
  write(u6,*) ' PRP3) Gradient.'
  write(u6,*)
  write(u6,'(1x,a,i12,a)') 'NBOTTOM        :',NBOTTOM,' (essential)'
  write(u6,'(1x,a,i12,a)') 'caspt2_grad    :',ngrad1,' (essential)'
  write(u6,'(1x,a,i12,a)') 'dens           :',ngrad2,' (essential)'
  write(u6,'(1x,a,i12)') '(Part of) PRP1 :',NPRP1-NBOTTOM
  write(u6,'(1x,a,i12)') 'caspt2_res     :',ngrad3
  write(u6,'(1x,a,i12)') 'trnds2o        :',ngrad4
  write(u6,'(1x,a,i12)') 'sigder         :',ngrad5
  write(u6,'(1x,a,i12)') 'clagx          :',ngrad6
  write(u6,'(1x,a,i12)') 'derheff        :',ngrad7
  write(u6,'(1x,a,i12)') 'DEPSAOffC      :',ngrad8
  write(u6,'(1x,a,i12)') 'OLagNS/OLagVVVO:',ngrad9
  write(u6,'(1x,a,i12)') 'CnstAB_SSDM    :',ngrad10
  write(u6,'(1x,a,i12)') 'XMS_Grad       :',ngrad11
  write(u6,'(1x,a,a )') '----------------','------------'
  write(u6,'(1x,a,i12)') 'Total Estimate :',ngrad
# endif
end if

NPRP = 0
if (IFPROP .or. do_grad) NPRP = max(NPRP1,NPRP2,ngrad)
if (IPRGLB >= USUAL) then
  write(u6,'(A)') repeat('-',80)
  write(u6,*) 'Estimated memory requirements:'
  write(u6,'(a,i12)') '  POLY3 :             ',NPOLY
  write(u6,'(a,i12)') '  RHS:                ',NMKRHS
  write(u6,'(a,i12)') '  SIGMA :             ',NSIGMA
  write(u6,'(a,i12)') '  PRPCTL:             ',NPRP
  !SVC: NPRP includes NDD if it is needed, so this is confusing
  !write(u6,'(a,i12)') '  DIADNS:             ',NDD
  write(u6,'(a,i12)') ' Available workspace: ',MXLEFT
  write(u6,*)
end if

NEED0 = max(NPOLY,NMKRHS,NSIGMA)
NEED = max(NEED0,NPRP)
if (NEED > MXLEFT) then
  if (NEED0 <= MXLEFT) then
    write(u6,'(5X,A)') repeat('*',26)
    write(u6,'(5X,A)') ' Memory problem!! The memory is insufficient '
    write(u6,'(5X,A)') ' for the property section.'
    write(u6,'(5X,A,I5,A)') '* Need at least ',2+NEED/Magic,' MB *'
    if (do_grad) then
      write(u6,'(5X,A)') ' The property (including gradient) section will be skipped.'
      write(u6,'(5X,A)') ' If you cannot increase the memory, consider using numerical gradients'
      if (IFPROP) write(u6,'(5X,A)') ' (without PROPerty keyword)'
    else
      write(u6,'(5X,A)') ' The property section will be skipped.'
    end if
    write(u6,'(5X,A)') repeat('*',26)
    if (do_grad) call Quit(_RC_MEMORY_ERROR_)
    NEED = NEED0
    IFPROP = .false.
  else if (.not. IFCHOL) then
    ! not a Cholesky calculation, keep old memory requirements
    write(u6,'(5X,A)') repeat('*',26)
    write(u6,'(5X,A)') '* Insufficient memory !! *'
    write(u6,'(5X,A,I5,A)') '* Need at least ',2+NEED0/Magic,' MB *'
    write(u6,'(5X,A)') repeat('*',26)
    write(u6,*)
    write(u6,*) ' Program execution stops -- sorry!'
    write(u6,*)
    call Quit(_RC_MEMORY_ERROR_)
  else if (max(NPOLY,NSIGMA) <= MXLEFT) then
    ! This is a Cholesky calculation, only give recommended amount
    write(u6,'(5X,A)') repeat('*',40)
    write(u6,'(5X,A,I5,A)') '* Below comfortable memory of ',2+NEED0/Magic,' MB *'
    write(u6,'(5X,A)') repeat('*',40)
    write(u6,*)
    write(u6,*) ' Program will try to adapt'
    write(u6,*) ' If it fails, please raise the memory to at least',2+max(NPOLY,int(0.55_wp*NMKRHS),NSIGMA)/Magic,' MB'
    write(u6,*) ' (Maybe more, this aint rocket science)'
    write(u6,*)
    write(u6,*) ' With print level DEBUG you will get some more'
    write(u6,*) ' informative memory estimates during computation'
    write(u6,*)
  else
    write(u6,'(5X,A)') repeat('*',26)
    write(u6,'(5X,A)') '* Insufficient memory... *'
    write(u6,'(5X,A)') repeat('*',26)
    write(u6,*)
    write(u6,'(A,I6,A)') ' If possible, I would like to have   ',NEED0/Magic,' MB'
    write(u6,'(A,I6,A)') ' Please raise the memory to at least ',max(NPOLY,int(0.55_wp*NMKRHS),NSIGMA)/Magic,' MB'
    write(u6,*) ' (Maybe more, this aint rocket science)'
    write(u6,*)
    call Quit(_RC_MEMORY_ERROR_)
  end if
end if

end subroutine SIZES
