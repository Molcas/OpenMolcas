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

subroutine DerCtr(mdci,mdcj,mdck,mdcl,ldot,JfGrd,IndGrd,JfHss,IndHss,JfG)

!#define _OLD_CODE_
use McKinley_global, only: sIrrep
use Index_Functions, only: iTri
use Symmetry_Info, only: nIrrep
#ifdef _OLD_CODE_
use Center_Info, only: dc
use Symmetry_Info, only: iOper
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mdci, mdcj, mdck, mdcl
logical(kind=iwp), intent(in) :: ldot
logical(kind=iwp), intent(out) :: JfGrd(3,4), JfHss(4,3,4,3), JfG(4)
integer(kind=iwp), intent(out) :: IndGrd(3,4,0:7), IndHss(4,3,4,3,0:7)
#include "Molcas.fh"
#include "disp.fh"
integer(kind=iwp) :: iAtom, ic1, ic2, iCar, iComp, ii, iIrrep, ij, istop, jAtom, jCar, JndGrd(3,4,0:7), nDisp, nnIrrep
logical(kind=iwp) :: IfG(4), IfGrd(3,4), IfHss(4,3,4,3)
logical(kind=iwp), external :: TF
#ifdef _OLD_CODE_
integer(kind=iwp) :: i, iCo(4), iCom(0:7,0:7), idcrr(0:7), ielem, iStabM(0:7), iTmp(0:7), j, jOper, LmbdR, nDCRR, nCoM, nMax, nStabM
logical(kind=iwp) :: chck
logical(kind=iwp), external :: TstFnc
#endif

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
#ifdef _OLD_CODE_
iCo(1) = mdci
iCo(2) = mdcj
iCo(3) = mdck
iCo(4) = mdcl
#endif
IndHss(:,:,:,:,0:nirrep-1) = 0
IfHss(:,:,:,:) = .false.
JfHss(:,:,:,:) = .false.
if (.not. ldot) return

do iAtom=1,4
  do jAtom=1,iAtom

    ! This segment of the code is not really needed.
    ! If turned on it should not do much of a difference.

#   ifdef _OLD_CODE_
    call DCR(LmbdR,dc(iCo(iAtom))%iStab,dc(iCo(iAtom))%nStab,dc(iCo(jAtom))%iStab,dc(iCo(jAtom))%nStab,iDCRR,nDCRR)

    ! Find the stabilizer for A and B

    call Inter(dc(iCo(iAtom))%iStab,dc(iCo(iAtom))%nStab,dc(iCo(jAtom))%iStab,dc(iCo(jAtom))%nStab,iStabM,nStabM)

    ! Generate all possible (left) CoSet
    ! To the stabilizer of A and B

    do jOper=0,nStabM-1
      iCoM(0:nIrrep-1,jOper) = ieor(iOper(0:nIrrep-1),iStabM(jOper))
    end do

    ! Order the Coset so we will have the unique ones first
    ! Check uniqueness

    nMax = 1
    outer: do j=1,nIrrep-1
      do i=0,nMax-1
        do ielem=0,nStabM-1
          if (iCoM(i,1) == iCoM(j,ielem)) cycle outer
        end do
      end do

      ! Move unique CoSet

      nMax = nMax+1
      iTmp(0:nStabM-1) = iCoM(nMax-1,0:nStabM-1)
      iCoM(nMax-1,0:nStabM-1) = iCoM(j,0:nStabM-1)
      iCoM(j,0:nStabM-1) = iTmp(0:nStabM-1)
      if (nMax == nIrrep/nStabM) exit outer
    end do outer

    ! Check if the derivative is needed in the present symmetry

    nCoM = nIrrep/nStabM

    do iCar=1,3
      if (iAtom == jAtom) then
        istop = iCar
      else
        iStop = 3
      end if
      do jCar=1,istop
        iComp = ieor(2**(iCar-1),2**(jCar-1))
        Chck = TstFnc(iCoM,0,iComp,nStabM)
        if (Chck) IfHss(iAtom,iCar,jAtom,jCar) = .true.
      end do
    end do
#   endif

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

return

end subroutine DerCtr
