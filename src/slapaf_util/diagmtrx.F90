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

subroutine DiagMtrx(H,nH,iNeg)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nH
real(kind=wp), intent(in) :: H(nH,nH)
integer(kind=iwp), intent(out) :: iNeg
#include "print.fh"
integer(kind=iwp) :: i, ij, iPrint, iRout, j, Lu, LuTmp, nq, nQQ
real(kind=wp) :: SumHii
logical(kind=iwp) :: Exists
character(len=16) :: filnam
real(kind=wp), allocatable :: EVal(:), EVec(:), qEVec(:), rK(:)

Lu = u6
iRout = 21
iPrint = nPrint(iRout)

call mma_allocate(EVal,nTri_Elem(nH),Label='EVal')
call mma_allocate(EVec,nH**2,Label='EVec')

! Copy elements for H

SumHii = Zero
do i=1,nH
  do j=1,i
    ij = nTri_Elem(i-1)+j
    EVal(ij) = H(i,j)
  end do
  SumHii = SumHii+H(i,i)
end do
!write(Lu,*) ' SumHii=',SumHii

! Set up a unit matrix

call unitmat(EVec,nH)

! Compute eigenvalues and eigenvectors

call NIDiag_new(EVal,EVec,nH,nH)
call Jacord(EVal,EVec,nH,nH)

! Print out the result

iNeg = 0
do i=1,nH
  if (EVal(nTri_Elem(i)) < Zero) iNeg = iNeg+1
end do
if (iprint > 5) then
  write(Lu,*)
  write(Lu,*) '*****************************************************************'
  write(Lu,*) '* Eigenvalues and Eigenvectors of the Hessian                   *'
  write(Lu,*) '*****************************************************************'
end if

filnam = 'SPCINX'
call f_Inquire(filnam,Exists)

if (Exists .and. (iprint > 5)) then

  ! Read linear combinations from disc

  LuTmp = 11
  call molcas_binaryopen_Vanilla(luTmp,filnam)
  !open(luTmp,File=filnam,Form='unformatted',Status='unknown')
  rewind(LuTmp)

  read(LuTmp) nq,nQQ

  if (nQQ == nH) then

    call mma_allocate(rK,nq*nQQ,Label='rK')
    call mma_allocate(qEVec,nq*nH,Label='qEVec')

    call Print_qEVec(EVec,nH,EVal,nq,rK,qEVec,LuTmp)

    call mma_deallocate(qEVec)
    call mma_deallocate(rK)

  else

    write(Lu,*)
    write(Lu,*) 'Eigenvalues of the Hessian'
    write(Lu,*)
    write(Lu,'(1X,10F10.5)') (EVal(nTri_Elem(i)),i=1,nH)
    write(Lu,*)
    write(Lu,*) 'Eigenvectors of the Hessian'
    write(Lu,*)
    do i=1,nH
      write(Lu,'(1X,10F10.5)') (EVec((j-1)*nH+i),j=1,nH)
    end do
  end if

  close(LuTmp)

else if (iprint > 5) then

  call Print_qEVec2(nH,EVal,EVec)

end if

call mma_deallocate(EVec)
call mma_deallocate(EVal)

return

end subroutine DiagMtrx
