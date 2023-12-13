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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine GenerateB(CMO,nBas,nOrb2Loc,Lbl_AO,Lbl,nComp,Debug)
! Author: T.B. Pedersen
!
! Purpose: generate the dipole matrices for Boys localisation, i.e.
!          transform from AO to MO basis.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas, nOrb2Loc, nComp
real(kind=wp), intent(in) :: CMO(*), Lbl_AO(nBas,nBas,nComp)
real(kind=wp), intent(out) :: Lbl(nOrb2Loc,nOrb2Loc,nComp)
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: i, iComp, iMO, j
real(kind=wp) :: Cmp, Tst
real(kind=wp), allocatable :: Dbar(:,:)

if ((nBas < 1) .or. (nOrb2Loc < 1)) return

call mma_allocate(Dbar,nBas,nOrb2Loc,label='Dbar')
do iComp=1,nComp
  call DGEMM_('N','N',nBas,nOrb2Loc,nBas,One,Lbl_AO(:,:,iComp),nBas,CMO,nBas,Zero,Dbar,nBas)
  call DGEMM_('T','N',nOrb2Loc,nOrb2Loc,nBas,One,CMO,nBas,Dbar,nBas,Zero,Lbl(:,:,iComp),nOrb2Loc)
end do
call mma_deallocate(Dbar)

if (Debug) then
  write(u6,*)
  write(u6,*) 'In GenerateB'
  write(u6,*) '------------'
  write(u6,*) '[Assuming doubly occupied orbitals]'
  do iComp=1,nComp
    Cmp = Zero
    do iMO=1,nOrb2Loc
      Cmp = Cmp+Lbl(iMO,iMO,iComp)
    end do
    Cmp = Two*Cmp
    write(u6,'(A,I5,1X,F15.8)') 'Component, Exp. Val.:',iComp,Cmp
    do j=1,nOrb2Loc-1
      do i=j+1,nOrb2Loc
        Tst = Lbl(i,j,iComp)-Lbl(j,i,iComp)
        if (abs(Tst) > 1.0e-14_wp) then
          write(u6,*) 'GenerateB: broken symmetry!'
          write(u6,*) '  Component: ',iComp
          write(u6,*) '  i and j  : ',i,j
          write(u6,*) '  Dij      : ',Lbl(i,j,iComp)
          write(u6,*) '  Dji      : ',Lbl(j,i,iComp)
          write(u6,*) '  Diff.    : ',Tst
          call SysAbendMsg('GenerateB','Broken symmetry!',' ')
        end if
      end do
    end do
  end do
end if

end subroutine GenerateB
