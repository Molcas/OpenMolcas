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

subroutine ORTHO_MOTRA(nSym,nBas,nDel,Ovlp,CMO)
! Objective: Orthonormalize input vectors
!           (Gram-Schmidt orthogonaliztion)

use motra_global, only: N2MAX
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nDel(nSym)
real(kind=wp), intent(in) :: Ovlp(*)
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp) :: IJ, IM, ISYM, NORBI
real(kind=wp), allocatable :: W1(:), W2(:), W3(:)

! Allocate work space

call mma_allocate(W1,N2MAX,label='SCR1')
call mma_allocate(W2,N2MAX,label='SCR2')
call mma_allocate(W3,N2MAX,label='SCR3')

! Loop over symmetries and orthogonalize

IJ = 1
IM = 1
do ISYM=1,NSYM
  NORBI = NBAS(ISYM)-NDEL(ISYM)
  if (NORBI > 0) then
    call SQUARE(OVLP(IJ),W3,1,NBAS(ISYM),NBAS(ISYM))
    call DGEMM_('N','N',NBAS(ISYM),NORBI,NBAS(ISYM),One,W3,NBAS(ISYM),CMO(IM),NBAS(ISYM),Zero,W2,NBAS(ISYM))
    call DGEMM_('T','N',NORBI,NORBI,NBAS(ISYM),One,CMO(IM),NBAS(ISYM),W2,NBAS(ISYM),Zero,W1,NORBI)
    call ORTHOX_MOTRA(W1,CMO(IM),NORBI,NBAS(ISYM))
  end if
  IJ = IJ+NBAS(ISYM)*(NBAS(ISYM)+1)/2
  IM = IM+NBAS(ISYM)*NBAS(ISYM)
end do

! Deallocate work space and exit

call mma_deallocate(W1)
call mma_deallocate(W2)
call mma_deallocate(W3)

return

end subroutine ORTHO_MOTRA
