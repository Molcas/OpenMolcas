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

subroutine Prepare(nGrdPt,Grid,B,GrdI)
! Some stuff for preparing the gradient integral computation

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iChTbl
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrdPt
real(kind=wp), intent(in) :: Grid(3,nGrdPt), B(nGrdPt)
real(kind=wp), intent(out) :: GrdI(4,nGrdPt)
#include "Molcas.fh"
#include "disp.fh"
integer(kind=iwp) :: iCar, iCnt, iCnttp, iComp, iIrrep, jOper, LuWr, mc, mdc, mDisp, nCnttp_Valence, nDiff, nDisp
logical(kind=iwp) :: DoRys, TstFnc
character, parameter :: xyz(0:2) = ['x','y','z']
integer(kind=iwp), external :: iPrmt

LuWr = u6
DoRys = .true.
nDiff = 3
call IniSew(DoRys,nDiff)

! Copy the grid coordinates and weights in ONE array
! This is the only solution I found to pass info trough oneel_g !

GrdI(1:3,:) = Grid
GrdI(4,:) = B

nCnttp_Valence = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux) exit
  nCnttp_Valence = nCnttp_Valence+1
end do

! Compute number of centers and displacements. Ignore pseudo centers.
! If any pseudo centers disable use of translational and rotational
! invariance.

mDisp = 0
mdc = 0
do iCnttp=1,nCnttp_Valence
  if (dbsc(iCnttp)%pChrg) then
    mdc = mdc+dbsc(iCnttp)%nCntr
  else
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      mDisp = mDisp+3*(nIrrep/dc(mdc)%nStab)
    end do
  end if
end do

! Initialize the Direct array. Why? I don't know.

Direct(:) = .true.

! Generate symmetry adapted cartesian displacements

IndDsp(:,:) = 0
InxDsp(:,:) = 0
Disp_Fac(:,:,:) = One
mult_Disp(:) = 1
nDisp = 0
do iIrrep=0,nIrrep-1
  lDisp(iIrrep) = 0
  ! Loop over basis function definitions
  mdc = 0
  mc = 1
  do iCnttp=1,nCnttp_Valence
    ! Loop over unique centers associated with this basis set.
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      IndDsp(mdc,iIrrep) = nDisp
      ! Loop over the cartesian components
      do iCar=0,2
        iComp = 2**iCar
        if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab) .and. (.not. dbsc(iCnttp)%pChrg)) then
          nDisp = nDisp+1
          if (iIrrep == 0) InxDsp(mdc,iCar+1) = nDisp
          lDisp(iIrrep) = lDisp(iIrrep)+1
          mult_Disp(nDisp) = nIrrep/dc(mdc)%nStab
          if (iIrrep == 0) then
            do jOper=0,nIrrep-1
              Disp_Fac(iCar+1,jOper,mdc) = iPrmt(jOper,iComp)*iChTbl(iIrrep,jOper)
            end do
          end if
          write(ChDisp(nDisp),'(A,1X,A1)') dc(mdc)%LblCnt,xyz(iCar)

        end if
      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do
  end do

end do

if (nDisp /= mDisp) then
  call WarningMessage(2,'Error in espf/prepare')
  write(LuWr,*) ' Wrong number of symmetry adapted displacements',nDisp,'=/=',mDisp
  call Abend()
end if

return

end subroutine Prepare
