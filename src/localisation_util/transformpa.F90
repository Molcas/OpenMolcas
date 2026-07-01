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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine transformPA(PA,nOrb2Loc,Umat,forward)
! this subroutine transforms the <s|PA|t> matrices according to PA <- U^T * PA * U

use Localisation_globals, only: BName, Debug, nAtoms, nBas_Start
use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOrb2Loc
real(kind=wp), intent(inout) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(in) :: Umat(nOrb2Loc,nOrb2Loc)
logical(kind=iwp), intent(in) :: forward
integer(kind=iwp) :: iAtom
character(len=LenIn+8) :: PALbl
real(kind=wp), allocatable :: QU(:,:)

call mma_allocate(QU,nOrb2Loc,nOrb2Loc,Label='QU')
do iAtom=1,nAtoms
  QU(:,:) = Zero

  if (forward) then
    ! transform as U^T*PA*U to reset the step
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
                One,PA(:,:,iAtom),nOrb2Loc, &
                    Umat,nOrb2Loc, &
                Zero,QU,nOrb2Loc)
    call dgemm_('T','N',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
                One,Umat,nOrb2Loc, &
                    QU,nOrb2Loc, &
                Zero,PA(:,:,iAtom),nOrb2Loc)

  else
    ! transform as U*PA*U^T
    call dgemm_('N','T',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
                One,PA(:,:,iAtom),nOrb2Loc, &
                    Umat,nOrb2Loc, &
                Zero,QU,nOrb2Loc)
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc, &
                One,Umat,nOrb2Loc, &
                    QU,nOrb2Loc, &
                Zero,PA(:,:,iAtom),nOrb2Loc)
  end if
end do
call mma_deallocate(QU)

if (Debug) then
  write(u6,*) 'In transformPA'
  write(u6,*) '--------------'
  do iAtom=1,nAtoms
    PALbl = 'PA__'//BName(nBas_Start(iAtom))(1:LenIn)
    call RecPrt(PALbl,' ',PA(:,:,iAtom),nOrb2Loc,nOrb2Loc)
  end do
end if

end subroutine transformPA
