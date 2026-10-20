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

subroutine MK_TWDM(mSym,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,iOFF,NBASF,ISY12)
! CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
! AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: MUL
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mSYM, nTDMZZ, nSCR, IOFF(mSYM), NBASF(mSym), ISY12
real(kind=wp), intent(in) :: TDMZZ(NTDMZZ), WDMZZ(NTDMZZ)
real(kind=wp), intent(out) :: SCR(nSCR,4)
integer(kind=iwp) :: I, IJ, IOF, ISY, ISY1, ISY2, ITD, J, NB, NB1, NB2
real(kind=wp) :: TDM, WDM

SCR(:,:) = Zero
if (ISY12 == 1) then
  ! SPECIAL CASE: DIAGONAL SYMMETRY BLOCKS.
  IOF = 0
  ITD = 0
  do ISY=1,mSym
    NB = NBASF(ISY)
    if (NB == 0) cycle
    do J=1,NB
      do I=1,NB
        ITD = ITD+1
        TDM = TDMZZ(ITD)
        WDM = WDMZZ(ITD)
        IJ = IOF+iTri(I,J)
        if (I > J) then
          SCR(IJ,2) = SCR(IJ,2)+TDM
          SCR(IJ,4) = SCR(IJ,4)+WDM
        else if (I < J) then
          SCR(IJ,2) = SCR(IJ,2)-TDM
          SCR(IJ,4) = SCR(IJ,4)-WDM
        end if
        SCR(IJ,1) = SCR(IJ,1)+TDM
        SCR(IJ,3) = SCR(IJ,3)+WDM
      end do
    end do
    IOF = IOF+nTri_Elem(NB)
  end do
else
  ! GENERAL CASE, NON-DIAGONAL SYMMETRY BLOCKS
  ! THEN LOOP OVER ELEMENTS OF TDMZZ
  ITD = 0
  do ISY1=1,mSym
    NB1 = NBASF(ISY1)
    if (NB1 == 0) cycle
    ISY2 = MUL(ISY1,ISY12)
    NB2 = NBASF(ISY2)
    if (NB2 == 0) cycle
    if (ISY1 > ISY2) then
      do J=1,NB2
        do I=1,NB1
          ITD = ITD+1
          TDM = TDMZZ(ITD)
          WDM = WDMZZ(ITD)
          IJ = IOFF(ISY1)+I+NB1*(J-1)
          SCR(IJ,1) = SCR(IJ,1)+TDM
          SCR(IJ,2) = SCR(IJ,2)+TDM
          SCR(IJ,3) = SCR(IJ,3)+WDM
          SCR(IJ,4) = SCR(IJ,4)+WDM
        end do
      end do
    else
      do J=1,NB2
        do I=1,NB1
          ITD = ITD+1
          TDM = TDMZZ(ITD)
          WDM = WDMZZ(ITD)
          IJ = IOFF(ISY2)+J+NB2*(I-1)
          SCR(IJ,1) = SCR(IJ,1)+TDM
          SCR(IJ,2) = SCR(IJ,2)-TDM
          SCR(IJ,3) = SCR(IJ,3)+WDM
          SCR(IJ,4) = SCR(IJ,4)-WDM
        end do
      end do
    end if
  end do
end if

end subroutine MK_TWDM
