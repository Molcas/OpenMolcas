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
! Copyright (C) 2000, Jonna Stalring                                   *
!***********************************************************************

subroutine rddj(G1r,G1Q,G2r,iestate)
! Jonna 000411
!
! Reads the one and two electron densities for estate
! and returns them in rectangular and single triangular storage

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: LuJob
use input_mclr, only: iTOC, ntAsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: G1r(*), G1Q(*), G2r(*)
integer(kind=iwp), intent(in) :: iestate
integer(kind=iwp) :: i, iB, iDij, iDkl, iIJKL, iR, iRij, iRijkl, iRkl, jB, jDisk, kB, lB, nG1, nG2
real(kind=wp) :: Fact, rdum(1)
real(kind=wp), allocatable :: G2Q(:)

ng1 = nTri_Elem(ntash)
ng2 = nTri_Elem(ng1)

call mma_allocate(G2Q,ng2,Label='G2Q')

! Read one and two el dens for state iestate

iR = iestate
jdisk = itoc(3)
do i=1,iR-1
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
end do
call dDaFile(LUJOB,2,G1q,ng1,jDisk)
call dDaFile(LUJOB,0,rdum,ng1,jDisk)
call dDaFile(LUJOB,2,G2Q,Ng2,jDisk)
call dDaFile(LUJOB,0,rdum,Ng2,jDisk)

! Make one el rectangular and two el singel triang.

do iB=1,ntash
  do jB=1,ntash
    iDij = iTri(ib,jB)
    iRij = jb+(ib-1)*ntash
    do kB=1,ntash
      do lB=1,ntash
        iDkl = iTri(kB,lB)
        iRkl = lb+(kb-1)*ntash
        fact = One
        if ((iDij >= iDkl) .and. (kB == lB)) fact = Two
        if ((iDij < iDkl) .and. (iB == jB)) fact = Two
        iijkl = iTri(iDij,iDkl)
        iRijkl = iTri(iRij,iRkl)
        G2R(iRijkl) = Fact*G2Q(iijkl)
      end do
    end do
  end do
end do
do iB=1,ntash
  do jB=1,ntash
    G1R(iB+(jB-1)*ntash) = G1Q(iTri(iB,jB))
  end do
end do

call mma_deallocate(G2Q)

end subroutine rddj
