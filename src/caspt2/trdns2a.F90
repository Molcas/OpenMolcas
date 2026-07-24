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

subroutine TRDNS2A(IVEC,JVEC,DPT2,NDPT2)
! Add to the diagonal blocks of transition density matrix,
!    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>,
! where p,q are active indices. Compare TRDNS2D.
! The present solution gives just a reasonable approximation,
! with correct trace.

use Index_Functions, only: iTri
use PrintLevel, only: VERBOSE
use caspt2_global, only: DREF, iPrGlb
use general_data, only: nActel, nAsh
use caspt2_module, only: nAES, nAshT, nInDep, nIsh, nISup, nOrb, nSym
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NDPT2
real(kind=wp), intent(inout) :: DPT2(NDPT2)
integer(kind=iwp) :: ICASE, IOFDPT, ISYM, IT, ITABS, ITQ, ITU, IU, IUABS, IUQ, IUT, lVec1, lVec2, NA, NADIFF, NAHOLE, NI, NIN, &
                     NIS, NO, nVec
real(kind=wp) :: COEF1, COEF2, D, DR, OVL
integer(kind=iwp), parameter :: NACTD(13) = [1,2,2,-1,0,1,1,-2,-2,-1,-1,0,0]
real(kind=wp), external :: RHS_DDOT

if (IPRGLB >= VERBOSE) then
  call WarningMessage(1,'Computing approximated density.')
  write(u6,*) ' The active/active submatrices of the density'
  write(u6,*) ' matrix is roughly approximated only.'
end if

COEF1 = Zero
COEF2 = Zero
NAHOLE = 2*NASHT-NACTEL
do ICASE=1,13
  NADIFF = NACTD(ICASE)
  if (NACTEL+NADIFF < 0) cycle
  if (NAHOLE-NADIFF < 0) cycle
  OVL = Zero
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    NIS = NISUP(ISYM,ICASE)
    NVEC = NIN*NIS
    if (NVEC == 0) cycle
    call RHS_ALLO(NIN,NIS,LVEC1)
    call RHS_ALLO(NIN,NIS,LVEC2)
    call RHS_READ_SR(LVEC1,iCASE,iSYM,IVEC)
    call RHS_READ_SR(LVEC2,iCASE,iSYM,JVEC)
    OVL = OVL+RHS_DDOT(NIN,NIS,LVEC1,LVEC2)
    call RHS_FREE(LVEC1)
    call RHS_FREE(LVEC2)
  end do
  if (NADIFF > 0) then
    COEF1 = COEF1+OVL*real(NADIFF,kind=wp)/real(max(1,NAHOLE),kind=wp)
    COEF2 = COEF2+OVL*real(NAHOLE-NADIFF,kind=wp)/real(max(1,NAHOLE),kind=wp)
  else
    COEF2 = COEF2+OVL*real(NACTEL+NADIFF,kind=wp)/real(max(1,NACTEL),kind=wp)
  end if
end do

IOFDPT = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NO = NORB(ISYM)
  do IT=1,NA
    ITQ = NI+IT
    ITABS = NAES(ISYM)+IT
    do IU=1,IT
      IUQ = NI+IU
      IUABS = NAES(ISYM)+IU
      DR = DREF(iTri(ITABS,IUABS))
      D = COEF2*DR
      if (IT == IU) D = D+Two*COEF1
      ITU = ITQ+NO*(IUQ-1)
      IUT = IUQ+NO*(ITQ-1)
      DPT2(IOFDPT+ITU) = DPT2(IOFDPT+ITU)+D
      DPT2(IOFDPT+IUT) = DPT2(IOFDPT+ITU)
    end do
  end do
  IOFDPT = IOFDPT+NO**2
end do

end subroutine TRDNS2A
