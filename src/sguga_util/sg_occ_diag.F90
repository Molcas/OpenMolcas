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

!> DIAG(I) = sum_lev n_lev(I) ETA(lev) for the CSFs I of symmetry ISYCI of
!> state iState (n the occupation numbers of the levels); from H0DIAG and
!> DIELMV of CASPT2: the upper and lower half walks of each (midvertex,
!> upper symmetry) block separately.
subroutine SG_Occ_Diag(iState,ISYCI,ETA,DIAG,nDiag)

use Symmetry_Info, only: Mul
use sguga, only: CIS, SGS
use CSF_Graph, only: nPack
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iState, ISYCI, nDiag
real(kind=wp), intent(in) :: ETA(*)
real(kind=wp), intent(out) :: DIAG(nDiag)
integer(kind=iwp) :: ICS, IEMU, ISYDWN, ISYUP, JCS, MV, NDWN, NUP

DIAG(:) = Zero
IEMU = 1
do MV=1,CIS(iState)%nMidV
  do ISYUP=1,SGS(iState)%nSym
    NUP = CIS(iState)%NOW(1,ISYUP,MV)
    if (NUP == 0) cycle
    ISYDWN = Mul(ISYUP,ISYCI)
    NDWN = CIS(iState)%NOW(2,ISYDWN,MV)
    if (NDWN == 0) cycle
    ICS = 1+CIS(iState)%IOW(1,ISYUP,MV)
    JCS = 1+CIS(iState)%IOW(2,ISYDWN,MV)
    call Block_Diag(ICS,JCS,NUP,NDWN,DIAG(IEMU))
    IEMU = IEMU+NUP*NDWN
  end do
end do

contains

  subroutine Block_Diag(ICS,JCS,NUP,NDWN,EMU)
    integer(kind=iwp), intent(in) :: ICS, JCS, NUP, NDWN
    real(kind=wp), intent(inout) :: EMU(NUP,NDWN)
    integer(kind=iwp) :: I, IC, IC1, II, IOC, ISTEP, LEV, LV1, nIpWlk, nLev, MidLev
    real(kind=wp) :: rSum

    nIpWlk = CIS(iState)%nIpWlk
    nLev = SGS(iState)%nLev
    MidLev = SGS(iState)%MidLev
    ! the upper half walks
    do I=1,NUP
      II = ICS-1+NIPWLK*(I-1)
      rSum = Zero
      do LV1=MidLev+1,nLev,nPack
        II = II+1
        IC = CIS(iState)%ICase(II)
        do LEV=LV1,min(LV1+(nPack-1),nLev)
          IC1 = IC/4
          ISTEP = IC-4*IC1
          IOC = (ISTEP+1)/2
          rSum = rSum+real(IOC,kind=wp)*ETA(LEV)
          IC = IC1
        end do
      end do
      EMU(I,:) = EMU(I,:)+rSum
    end do
    ! the lower half walks
    do I=1,NDWN
      II = JCS-1+NIPWLK*(I-1)
      rSum = Zero
      do LV1=1,MidLev,nPack
        II = II+1
        IC = CIS(iState)%ICase(II)
        do LEV=LV1,min(LV1+(nPack-1),MidLev)
          IC1 = IC/4
          ISTEP = IC-4*IC1
          IOC = (ISTEP+1)/2
          rSum = rSum+real(IOC,kind=wp)*ETA(LEV)
          IC = IC1
        end do
      end do
      EMU(:,I) = EMU(:,I)+rSum
    end do
  end subroutine Block_Diag

end subroutine SG_Occ_Diag
