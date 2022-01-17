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
! Copyright (C) 1991,1992,2008, Roland Lindh                           *
!***********************************************************************

subroutine Mk_ChDisp()
!***********************************************************************
!                                                                      *
! Object: To generate the displacement labels
!                                                                      *
! Called from: Gateway
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September '91                                            *
!                                                                      *
!             Modified to complement GetInf, January '92.              *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iCar, iCnt, iCnttp, iComp, iDisp, iIrrep, mc, mdc, mDisp, nCnttp_Valence, nDisp(0:7)
logical(kind=iwp) :: TstFnc
integer(kind=iwp), allocatable :: DegDisp(:)
character(len=LenIn6), allocatable :: ChDisp(:)
character, parameter :: xyz(0:2) = ['x','y','z']

!                                                                      *
!***********************************************************************
!                                                                      *
nCnttp_Valence = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux) exit
  nCnttp_Valence = nCnttp_Valence+1
end do

mDisp = 0
mdc = 0
do iCnttp=1,nCnttp_Valence
  if (dbsc(iCnttp)%pChrg) then
    mdc = mdc+dbsc(iCnttp)%nCntr
    cycle
  end if
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    mDisp = mDisp+3*(nIrrep/dc(mdc)%nStab)
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(ChDisp,mDisp,label='ChDisp')
call mma_allocate(DegDisp,mDisp,label='DegDisp')
iDisp = 0
do iIrrep=0,nIrrep-1
  ! Loop over basis function definitions
  mdc = 0
  mc = 1
  nDisp(iIrrep) = 0
  do iCnttp=1,nCnttp_Valence
    ! Loop over unique centers associated with this basis set.
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      ! Loop over the cartesian components
      do iCar=0,2
        iComp = 2**iCar
        if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab) .and. (.not. dbsc(iCnttp)%pChrg)) then
          iDisp = iDisp+1
          ChDisp(iDisp) = ' '
          write(ChDisp(iDisp)(1:(LenIn6)),'(A,1X,A1)') dc(mdc)%LblCnt(1:LenIn4),xyz(iCar)
          DegDisp(iDisp) = nIrrep/dc(mdc)%nStab
          nDisp(iIrrep) = nDisp(iIrrep)+1
        end if

      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do
  end do

end do

if (iDisp /= mDisp) then
  call WarningMessage(2,' Wrong number of symmetry adapted displacements')
  write(u6,*) iDisp,'=/=',mDisp
  call Abend()
end if

call Put_iScalar('nChDisp',iDisp)
call Put_cArray('ChDisp',ChDisp(1),LenIn6*iDisp)
call Put_iArray('nDisp',nDisp,nIrrep)
call Put_iArray('DegDisp',DegDisp,iDisp)
call mma_deallocate(ChDisp)
call mma_deallocate(DegDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_ChDisp
