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

subroutine HMAT(C,HC,HH,HD,NDIM,NDIMH,NTRIAL)
! RASSCF program: version IBM-3090: SX section
!
! Purpose: Calculation of the Davidson Hamiltonian.
! The sigma vector HC is computed in SIGVEC
!
! ********** IBM-3090 Release 88 09 08 **********

use Index_Functions, only: nTri_Elem
use wadr, only: BM, DIA, F1, F2, NLX, SXG, SXH, SXN
use rasscf_global, only: NSXS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NDIM, NDIMH, NTRIAL
real(kind=wp) :: C(*), HC(*), HH(*), HD(NDIM)
integer(kind=iwp) :: I, IJ, IST, J, JST, L1
real(kind=wp), allocatable :: XC(:), XX(:)
real(kind=wp), external :: DDot_

! COMPUTE SIGMA VECTORS

IST = 1+NDIMH*NDIM

call mma_allocate(XX,NLX,Label='XX')
call mma_allocate(XC,NSXS,Label='XC')

call SIGVEC(C(IST),HC(IST),HD,BM,SXN,SXG,SXH,DIA,F1,F2,XX,XC,NTRIAL)

call mma_deallocate(XX)
call mma_deallocate(XC)

! Compute a new block of the HH matrix

IJ = nTri_Elem(NDIMH)
L1 = NDIMH+1
NDIMH = NDIMH+NTRIAL
do I=L1,NDIMH
  JST = 1
  do J=1,I
    IJ = IJ+1
    HH(IJ) = DDOT_(NDIM,HC(IST),1,C(JST),1)
    JST = JST+NDIM
  end do
  IST = IST+NDIM
end do

end subroutine HMAT
