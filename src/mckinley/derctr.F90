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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine DerCtr(ldot,JfGrd,IndGrd,JfHss,IndHss,JfG,nSD,iSD4)

use McKinley_global, only: sIrrep
use Index_Functions, only: iTri
use Disp, only: IndDsp
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: ldot
logical(kind=iwp), intent(out) :: JfGrd(3,4), JfHss(4,3,4,3), JfG(4)
integer(kind=iwp), intent(out) :: IndGrd(3,4,0:7), IndHss(4,3,4,3,0:7)
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4)
integer(kind=iwp) :: iAtom, ic1, ic2, iCar, iComp, ii, iIrrep, ij, istop, jAtom, jCar, JndGrd(3,4,0:7), mdci, mdcj, mdck, mdcl, &
                     nDisp, nnIrrep
logical(kind=iwp) :: IfG(4), IfGrd(3,4), IfHss(4,3,4,3)
logical(kind=iwp), external :: TF

mdci = iSD4(10,1)
mdcj = iSD4(10,2)
mdck = iSD4(10,3)
mdcl = iSD4(10,4)

nnIrrep = nIrrep
IfGrd(:,:) = .false.
if (sIrrep) nnIrrep = 1

do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdci,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdci,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,1,iIrrep) = nDisp
      IfGrd(iCar+1,1) = .true.
    else
      IndGrd(iCar+1,1,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdcj,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdcj,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,2,iIrrep) = nDisp
      IfGrd(iCar+1,2) = .true.
    else
      IndGrd(iCar+1,2,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdck,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdck,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,3,iIrrep) = nDisp
      IfGrd(iCar+1,3) = .true.
    else
      IndGrd(iCar+1,3,iIrrep) = 0
    end if
  end do
end do
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(mdcl,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(mdcl,iIrrep,iComp)) then
      nDisp = nDisp+1
      IndGrd(iCar+1,4,iIrrep) = nDisp
      IfGrd(iCar+1,4) = .true.
    else
      IndGrd(iCar+1,4,iIrrep) = 0
    end if
  end do
end do
JndGrd(:,:,0:nnIrrep-1) = IndGrd(:,:,0:nnIrrep-1)
JfGrd(:,:) = IfGrd(:,:)
IndHss(:,:,:,:,0:nirrep-1) = 0
IfHss(:,:,:,:) = .false.
JfHss(:,:,:,:) = .false.
if (.not. ldot) return

do iAtom=1,4
  do jAtom=1,iAtom

    ! Calculate the index for the derivative

    do iIrrep=0,nnIrrep-1
      do iCar=1,3
        if (iAtom == jAtom) then
          iStop = iCar
        else
          iStop = 3
        end if
        do jCar=1,iStop
          IfHss(iAtom,iCar,jAtom,jCar) = .true.
          if ((jndGrd(iCar,iAtom,iIrrep) > 0) .and. (jndGrd(jCar,jAtom,iIrrep) > 0)) then
            IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = iTri(JndGrd(iCar,iAtom,iIrrep),JndGrd(jCar,jAtom,iIrrep))
          else
            IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = 0
          end if
        end do
      end do
    end do
  end do
end do

! Scramble the control array for the hessian

IfG(:) = .true.
do iAtom=1,4
  JfG(iAtom) = IfG(iAtom)
  do jAtom=1,iAtom
    do iCar=1,3
      if (iAtom == jAtom) then
        iStop = iCar
      else
        iStop = 3
      end if
      if (iAtom >= jAtom) then
        JfHss(iAtom,iCar,jAtom,1:iStop) = IfHss(iAtom,iCar,jAtom,1:iStop)
      else if (iAtom < jAtom) then
        JfHss(iAtom,iCar,jAtom,1:iStop) = IfHss(jAtom,1:iStop,iAtom,iCar)
      end if
    end do
  end do
end do
if (sIrrep) then
  do ii=1,4
    do ic1=1,3
      if (IndGrd(ic1,ii,0) == 0) JfGrd(ic1,ii) = .false.
      do ij=1,4
        do ic2=1,3
          if (IndHss(ii,ic1,ij,ic2,0) == 0) JfHss(ii,ic1,ij,ic2) = .false.
        end do
      end do
    end do
  end do
end if
!----------------------------------------------------------------------*
!
!  End Hess
!
!----------------------------------------------------------------------*

end subroutine DerCtr
