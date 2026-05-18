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

subroutine TRANSFOCK(TORB,NTORB,F,NF,IDIR)
! Purpose: given an orbital transformation array and some
! one-electron matrix in storage format as e.g. FIFA and FIMO
! transform the matrix to use the new orbital basis.

use caspt2_module, only: nIsh, nOMx, nRas1, nRas2, nRas3, nSsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NTORB, NF, IDIR
real(kind=wp), intent(in) :: TORB(NTORB)
real(kind=wp), intent(inout) :: F(NF)
integer(kind=iwp) :: NT, NI, NR1, NR2, NR3, NS, NO, IJOFF, ITOFF, I, J, II, JJ, IJ, ISYM, IOFF
real(kind=wp), allocatable :: FSQ(:), TMP(:), TSQ(:)

NT = 0
NOMX = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NR1 = NRAS1(ISYM)
  NR2 = NRAS2(ISYM)
  NR3 = NRAS3(ISYM)
  NS = NSSH(ISYM)
  NO = NI+NR1+NR2+NR3+NS
  NOMX = max(NOMX,NO)
  NT = NT+NI**2+NR1**2+NR2**2+NR3**2+NS**2
end do

call mma_allocate(FSQ,NOMX**2,LABEL='FSQ')
call mma_allocate(TSQ,NOMX**2,LABEL='TSQ')
call mma_allocate(TMP,NOMX**2,LABEL='TMP')

IJOFF = 0
ITOFF = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NR1 = NRAS1(ISYM)
  NR2 = NRAS2(ISYM)
  NR3 = NRAS3(ISYM)
  NS = NSSH(ISYM)
  NO = NI+NR1+NR2+NR3+NS
  if (NO == 0) cycle

  ! Copy the matrices to square storage: first fill with zeroes.
  TSQ(1:NO**2) = Zero
  ! Copy inactive TORB block to TSQ
  IOFF = 0
  do I=1,NI
    do J=1,NI
      TSQ(I+NO*(J-1)) = TORB(ITOFF+I+NI*(J-1))
    end do
  end do
  ! Copy ras1 T block to TSQ, and so on..
  ITOFF = ITOFF+NI**2
  IOFF = IOFF+NI
  do I=1,NR1
    II = IOFF+I
    do J=1,NR1
      JJ = IOFF+J
      TSQ(II+NO*(JJ-1)) = TORB(ITOFF+I+NR1*(J-1))
    end do
  end do

  ITOFF = ITOFF+NR1**2
  IOFF = IOFF+NR1
  do I=1,NR2
    II = IOFF+I
    do J=1,NR2
      JJ = IOFF+J
      TSQ(II+NO*(JJ-1)) = TORB(ITOFF+I+NR2*(J-1))
    end do
  end do

  ITOFF = ITOFF+NR2**2
  IOFF = IOFF+NR2
  do I=1,NR3
    II = IOFF+I
    do J=1,NR3
      JJ = IOFF+J
      TSQ(II+NO*(JJ-1)) = TORB(ITOFF+I+NR3*(J-1))
    end do
  end do
  ! Finally, the secondary orbitals (non-deleted, virtual).
  ITOFF = ITOFF+NR3**2
  IOFF = IOFF+NR3
  do I=1,NS
    II = IOFF+I
    do J=1,NS
      JJ = IOFF+J
      TSQ(II+NO*(JJ-1)) = TORB(ITOFF+I+NS*(J-1))
    end do
  end do
  ITOFF = ITOFF+NS**2
  ! Now transfer the Fock matrix block to square storage:
  IJ = 0
  do I=1,NO
    do J=1,I
      IJ = IJ+1
      FSQ(J+NO*(I-1)) = F(IJOFF+IJ)
      FSQ(I+NO*(J-1)) = F(IJOFF+IJ)
    end do
  end do
  if (IDIR >= 0) then
    ! T^T F T
    ! Transform, first do FSQ*TSQ -> TMP...
    call DGEMM_('N','N',NO,NO,NO,One,FSQ,NO,TSQ,NO,Zero,TMP,NO)
    ! ... and then do TSQ(transpose)*TMP -> FSQ.
    call DGEMM_('T','N',NO,NO,NO,One,TSQ,NO,TMP,NO,Zero,FSQ,NO)
  else
    ! T F T^T
    ! Or inverse transformation
    call DGEMM_('N','T',NO,NO,NO,One,FSQ,NO,TSQ,NO,Zero,TMP,NO)
    call DGEMM_('N','N',NO,NO,NO,One,TSQ,NO,TMP,NO,Zero,FSQ,NO)
  end if
  ! Transfer FSQ values back to F, in triangular storage.
  IJ = 0
  do I=1,NO
    do J=1,I
      IJ = IJ+1
      F(IJOFF+IJ) = FSQ(I+NO*(J-1))
    end do
  end do
  IJOFF = IJOFF+(NO*(NO+1))/2
  ! and repeat, using next symmetry block.
end do
call mma_deallocate(FSQ)
call mma_deallocate(TSQ)
call mma_deallocate(TMP)

end subroutine TRANSFOCK
