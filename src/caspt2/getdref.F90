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

subroutine GETDREF(DREF,NDREF)

use PrintLevel, only: DEBUG
use caspt2_global, only: iPrGlb
use caspt2_module, only: NASHT,  NG1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NDREF
real(kind=wp), intent(out) :: DREF(NDREF)
integer(kind=iwp) :: I, IJ, J
real(kind=wp), allocatable :: G1(:)

! Get active 1-el density matrix GAMMA1 and
! construct DREF in a tringular storage.

! Remember: NDREF=1 if NASHT=0.
DREF(:) = Zero

if (NASHT == 0) return
! Active 1-el density matrix:
call mma_allocate(G1,NG1,Label='G1')
call PT2_GET(NG1,'GAMMA1',G1)
do I=1,NASHT
  do J=1,I
    IJ = (I*(I-1))/2+J
    DREF(IJ) = G1(I+NASHT*(J-1))
  end do
end do
call mma_deallocate(G1)

if (IPRGLB >= DEBUG) write(u6,*) ' GETDREF has constructed DREF.'

end subroutine GETDREF

