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

subroutine transformPA(PA,nOrb2Loc,Umat,Umat_inv)
! this subroutine transforms the <s|PA|t> matrices according to PA <- U^T * PA * U

use stdalloc, only: mma_allocate,mma_deallocate
use definitions, only: wp,iwp,u6
use constants, only: Zero,One
use Molcas, only: LenIn
use Localisation_globals, only: BName,nBas_Start,nAtoms, Debug

implicit none

integer(kind=iwp), intent(in) :: nOrb2Loc
real(kind=wp), intent(in) :: Umat(nOrb2Loc,nOrb2Loc),Umat_inv(nOrb2Loc,nOrb2Loc)
real(kind=wp), intent(inout) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
integer(kind=iwp) :: iAtom
real(kind=wp),allocatable :: QU(:,:)
character(len=LenIn+8) :: PALbl

call mma_allocate(QU,nOrb2Loc,nOrb2Loc,Label="QU")
do iAtom = 1, nAtoms
    QU(:,:) = Zero
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,&
                        One,PA(:,:,iAtom),nOrb2Loc,&
                            Umat,nOrb2Loc,&
                        Zero,QU,nOrb2Loc)
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,&
                        One,Umat_inv,nOrb2Loc,&
                            QU,nOrb2Loc,&
                        Zero,PA(:,:,iAtom),nOrb2Loc)
end do
call mma_deallocate(QU)

if (Debug) then
    write(u6,*) 'In transformPA'
    write(u6,*) '--------------'
    do iAtom=1,nAtoms
        PALbl = 'PA__'//BName(nBas_Start(iAtom))(1:LenIn)
        call RecPrt(PALbl,' ',PA(:,:,iAtom),nOrb2Loc,nOrb2Loc)
    end do
endif

end subroutine transformPA

