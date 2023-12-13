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

subroutine Print_qEVec(EVec,nH,EVal,nq,rK,qEVec,LuTmp)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nH, nq, LuTmp
real(kind=wp), intent(in) :: EVec(nH,nH), EVal(nTri_Elem(nH))
real(kind=wp), intent(out) :: rK(nq,nH), qEVec(nq,nH)
integer(kind=iwp) :: iiQQ, IncQQ, iq, iQQ, Lu, mQQ
real(kind=wp) :: temp
character(len=14), allocatable :: qLbl(:)
real(kind=wp), parameter :: Thr = 1.0e-4_wp
real(kind=wp), external :: DDot_

call mma_allocate(qLbl,nq,Label='qLbl')

do iq=1,nq
  read(LuTmp) qLbl(iq),rK(iq,:)
end do

call DGEMM_('N','N',nq,nH,nH,One,rK,nq,EVec,nH,Zero,qEVec,nq)

Lu = u6
IncQQ = 5
do iiQQ=1,nH,IncQQ
  mQQ = min(nH,iiQQ+IncQQ-1)
  write(Lu,*)
  write(Lu,'(14X,5I10)') (iQQ,iQQ=iiQQ,mQQ)
  write(Lu,'(1X,A,5F10.6)') 'Eigenvalues   ',(EVal(nTri_Elem(iQQ)),iQQ=iiQQ,mQQ)
  write(Lu,*)
  do iq=1,nq
    temp = sqrt(DDot_(nH,qEVec(iq,1),nq,qEVec(iq,1),nq))
    if (temp > Thr) write(Lu,'(1X,A,5F10.6)') qLbl(iq),qEVec(iq,iiQQ:mQQ)
  end do
  write(Lu,*)
end do

call mma_deallocate(qLbl)

return

end subroutine Print_qEVec
