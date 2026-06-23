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

subroutine DIAFCK(NO,FOCK,IOSTA,IOEND,TSCT,NB,CMO1SCT,CMO2SCT)
! Use orbital rotations on a subblock of active orbitals
! in order to diagonalize that block of a Fock matrix, and
! in addition apply a counterrotation to the CI array that
! keeps the wave function invariant to the rotation.

! FOCK is the Fock matrix of this symmetry, in square form.
! to be tranformed. TSCT is a corresponding section of the
! transformation matrix, CMO1SCT is a section of the CMO array (in)
! and CMO2SCT the same, transformed.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NO, IOSTA, IOEND, NB
real(kind=wp), intent(inout) :: FOCK(NO,NO)
real(kind=wp), intent(out) :: TSCT(IOSTA:IOEND,IOSTA:IOEND), CMO2SCT(NB,*)
real(kind=wp), intent(in) :: CMO1SCT(NB,*)
integer(kind=iwp) :: I, IJ, J, JMX, K, NSCT
real(kind=wp) :: SWAP, V, VMX
real(kind=wp), allocatable :: TMP(:)

! Size of orbital section to process:
NSCT = IOEND+1-IOSTA
! Local scratch area:
call mma_allocate(TMP,NO*NSCT,Label='TMP')
! Put part of FOCK to be diagonalized into triangular scratch area:
IJ = 0
do I=IOSTA,IOEND
  TMP(IJ+1:IJ+I-IOSTA+1) = FOCK(I,IOSTA:I)
  IJ = IJ+I-IOSTA+1
end do
! Put unit matrix into TSCT:
call unitmat(TSCT,NSCT)
! Diagonalize, and order for best submatrix condition:
call Jacob(TMP,TSCT,NSCT,NSCT)

do I=IOSTA,IOEND
  JMX = I
  VMX = abs(TSCT(I,JMX))

  do J=I+1,IOEND
    V = abs(TSCT(I,J))
    if (V > VMX) then
      JMX = J
      VMX = V
    end if
  end do

  if (JMX > I) then
    do K=IOSTA,IOEND
      SWAP = TSCT(K,I)
      TSCT(K,I) = TSCT(K,JMX)
      TSCT(K,JMX) = SWAP
    end do
  end if

  if (TSCT(I,I) < Zero) TSCT(IOSTA:IOEND,I) = -TSCT(IOSTA:IOEND,I)
end do

! Transform the Fock matrix:
call DGEMM_('N','N',NO,NSCT,NSCT,One,FOCK(1,IOSTA),NO,TSCT,NSCT,Zero,TMP,NO)
FOCK(:,IOSTA:IOEND) = reshape(TMP(:),[NO,NSCT])

do I=IOSTA,IOEND
  FOCK(I,1:NO) = TMP(NO*(I-IOSTA)+1:NO*(I-IOSTA+1))
end do

call DGEMM_('T','N',NSCT,NSCT,NSCT,One,TSCT,NSCT,TMP(IOSTA),NO,Zero,FOCK(IOSTA,IOSTA),NO)

! Then transform the CMO coeffs:
call DGEMM_('N','N',NB,NSCT,NSCT,One,CMO1SCT,NB,TSCT,NSCT,Zero,CMO2SCT,NB)

call mma_deallocate(TMP)

end subroutine DIAFCK
