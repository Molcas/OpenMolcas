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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2016,2017, Roland Lindh                                *
!***********************************************************************

subroutine R1IntB()
!***********************************************************************
!                                                                      *
!     purpose: Read basis set informations and one-electron integrals  *
!              were not needed so far.                                 *
!                                                                      *
!***********************************************************************

use InfSCF, only: Darwin, KntE, lRel, MssVlc, nBT
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iComp, iOpt, iRC, iSyLbl
character(len=8) :: Label

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for kinetic energy, mass velocity and Darwin integrals

call mma_allocate(KntE,nBT+4,Label='KntE')
call mma_allocate(MssVlc,nBT+4,Label='MssVlc')
call mma_allocate(Darwin,nBT+4,Label='Darwin')

! Read kinetic energy integrals
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'Kinetic '
call RdOne(iRc,iOpt,Label,iComp,KntE,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'R1Intb: Error readin ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if

! Read mass velocity integrals
lRel = .false.
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'MassVel '
call RdOne(iRc,iOpt,Label,iComp,MssVlc,iSyLbl)

if (iRc == 0) then
  ! Read Darwin integrals
  iRc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  iSyLbl = 1
  Label = 'Darwin  '
  call RdOne(iRc,iOpt,Label,iComp,Darwin,iSyLbl)
  if (iRc == 0) lRel = .true.
end if

if (.not. lRel) then
  call mma_deallocate(MssVlc)
  call mma_deallocate(Darwin)
  call mma_allocate(MssVlc,0,Label='MssVlc')
  call mma_allocate(Darwin,0,Label='Darwin')
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine R1IntB
