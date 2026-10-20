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

subroutine GETDPREF(DREF,NDREF,PREF,NPREF)
! Get active 1-density and 2-density matrices GAMMA1 and
! GAMMA2, and construct DREF and PREF which are in a tringular
! storage.

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG
use caspt2_module, only: NASHT, NG1, NG2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(out) :: DREF(NDREF), PREF(NPREF)
integer(kind=iwp) :: I, IJ, IJKL, IJKLT, IJLK, IJT, J, JI, JIKL, JILK, K, KL, KLT, L, LK, N2
real(kind=wp) :: P1, P2
real(kind=wp), allocatable :: G1(:), G2(:)

! Remember: NDREF=1 if NASHT=0. Similar NPREF.
DREF(1) = Zero
PREF(1) = Zero
if (NASHT == 0) return
! Active density 1-matrix:
call mma_allocate(G1,NG1,LABEL='G1')
call PT2_GET(NG1,'GAMMA1',G1)
do I=1,NASHT
  IJ = nTri_Elem(I-1)
  DREF(IJ+1:IJ+I) = G1(I:I+NASHT*(I-1):NASHT)
end do
call mma_deallocate(G1)

! CONSTRUCT PREF, 2-ELECTRON DENSITY MATRIX:
call mma_allocate(G2,NG2,LABEL='G2')
call PT2_GET(NG2,'GAMMA2',G2)
IJT = 0
IJKLT = 0
N2 = NASHT**2
do I=1,NASHT
  do J=1,I
    IJT = IJT+1
    IJ = I+NASHT*(J-1)
    JI = J+NASHT*(I-1)
    KLT = 0
    do K=1,NASHT
      do L=1,K
        KLT = KLT+1
        if (KLT > IJT) cycle
        IJKLT = IJKLT+1
        KL = K+NASHT*(L-1)
        LK = L+NASHT*(K-1)

        P1 = Half*G2(IJ+N2*(KL-1))
        P2 = Half*G2(IJ+N2*(LK-1))
        IJKL = iTri(IJ,KL)
        IJLK = iTri(IJ,LK)
        JIKL = iTri(JI,KL)
        JILK = iTri(JI,LK)
        PREF(IJKL) = P1
        PREF(IJLK) = P2
        PREF(JIKL) = P2
        PREF(JILK) = P1
      end do
    end do
  end do
end do

call mma_deallocate(G2)

if (IPRGLB >= DEBUG) then
  write(u6,*) ' GETDPREF has constructed DREF and PREF.'
end if

end subroutine GETDPREF
