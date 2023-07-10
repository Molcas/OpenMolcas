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

#include "compiler_features.h"
#ifdef _DEBUGPRINT_
subroutine DiagMtrx_T(H,nH,iNeg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nH, iNeg
real(kind=wp) :: H(*)
#include "print.fh"
integer(kind=iwp) :: i, iprint, iRout, j, Lu, LuTmp, nq, nQQ
logical(kind=iwp) :: Exists
character(len=16) :: filnam
real(kind=wp), allocatable :: EVal(:), EVec(:), qEVec(:), rK(:)

Lu = u6
iRout = 22
iprint = nPrint(iRout)

call mma_allocate(EVal,nH*(nH+1)/2,Label='EVal')
call mma_allocate(EVec,nH*nH,Label='EVec')

! Copy elements for H

EVal(:) = H(1:nH*(nH+1)/2)

! Set up a unit matrix

call unitmat(EVec,nH)

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
  write(Lu,*) 'Eigenvalues of the Hessian'
  write(Lu,*)
  write(Lu,'(5G20.6)') (EVal(i*(i+1)/2),i=1,nH)
end if

call f_Inquire('SPCINX',Exists)

if (Exists .and. (iprint > 5)) then

  ! Read linear combinations from disc

  LuTmp = 11
  filnam = 'SPCINX'
  !open(luTmp,File=filnam,Form='unformatted',Status='unknown')
  call molcas_binaryopen_vanilla(luTmp,filnam)
  rewind(LuTmp)

  read(LuTmp) nq,nQQ

  if (nQQ == nH) then

    call mma_allocate(rK,nq*nQQ,Label='rK')
    call mma_allocate(qEVec,nq*nH,Label='qEVec')

    call Print_qEVec(EVec,nH,EVal,nq,rK,qEVec,LuTmp)

    call mma_deallocate(qEVec)
    call mma_deallocate(rk)

  else

    write(Lu,*)
    write(Lu,*) 'Eigenvectors of the Hessian'
    write(Lu,*)
    do i=1,nH
      write(Lu,'(10F10.5)') (EVec((j-1)*nH+i),j=1,nH)
    end do
  end if

  close(LuTmp)

else if (iprint > 5) then
  write(Lu,*)
  write(Lu,*) 'Eigenvectors of the Hessian'
  write(Lu,*)
  do i=1,nH
    write(Lu,'(10F10.5)') (EVec((j-1)*nH+i),j=1,nH)
  end do

end if

call mma_deallocate(EVec)
call mma_deallocate(EVal)

return

end subroutine DiagMtrx_T
#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(DiagMtrx_T)

#endif
