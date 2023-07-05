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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nH, iNeg
real(kind=wp) :: H(nH,nH)
#include "print.fh"
integer(kind=iwp) :: i, ij, iPrint, iRout, j, Lu, LuTmp, nq, nQQ
real(kind=wp) :: SumHii
logical(kind=iwp) :: Exists
character(len=16) :: filnam
real(kind=wp), allocatable :: EVal(:), EVec(:), qEVec(:), rK(:)

Lu = u6
iRout = 21
iPrint = nPrint(iRout)

call mma_allocate(EVal,nH*(nH+1)/2,Label='EVal')
call mma_allocate(EVec,nH*nH,Label='EVec')

! Copy elements for H

SumHii = Zero
do i=1,nH
  do j=1,i
    ij = i*(i-1)/2+j
    EVal(ij) = H(i,j)
  end do
  SumHii = SumHii+H(i,i)
end do
!write(Lu,*) ' SumHii=',SumHii

! Set up a unit matrix

call dcopy_(nH*nH,[Zero],0,EVec,1)
call dcopy_(nH,[One],0,EVec,nH+1)

! Compute eigenvalues and eigenvectors

call NIDiag_new(EVal,EVec,nH,nH)
call Jacord(EVal,EVec,nH,nH)

! Print out the result

iNeg = 0
do i=1,nH
  if (EVal(i*(i+1)/2) < Zero) iNeg = iNeg+1
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
    write(Lu,'(1X,10F10.5)') (EVal(i*(i+1)/2),i=1,nH)
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
