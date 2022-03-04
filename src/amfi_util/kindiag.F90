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

subroutine kindiag(TKIN,TKINTRIA,ndim,evec,eval,breit)
!bs determines eigenvectors and -values of TKIN

implicit real*8(a-h,o-z)
dimension tkin(ndim,ndim), TKINTRIA((ndim*ndim+ndim)/2), eval(ndim), evec(ndim,ndim)
logical breit

!bs move symmetric matrix to triangular matrix
itria = 1
do irun2=1,ndim
  do irun1=1,irun2
    TKINTRIA(itria) = TKIN(irun1,irun2)
    itria = itria+1
  end do
end do
do irun2=1,ndim
  do irun1=1,ndim
    evec(irun1,irun2) = 0d0
  end do
end do
do irun1=1,ndim
  evec(irun1,irun1) = 1d0
end do
!bs now diagonalize
call Jacob(TKINTRIA,evec,ndim,ndim)
!bs get the eigenvalues
do irun=1,ndim
  eval(irun) = TKINTRIA((irun*irun+irun)/2)
end do
if (breit) then
  do irun=1,ndim
    eval(irun) = 0d0
  end do
end if
!bs ensure normalization of the vectors.
do IRUN=1,ndim
  fact = 0d0
  do JRUN=1,ndim
    fact = fact+evec(JRUN,IRUN)*evec(JRUN,IRUN)
  end do
  fact = 1d0/sqrt(fact)
  do JRUN=1,ndim
    evec(JRUN,IRUN) = fact*evec(JRUN,IRUN)
  end do
end do

return

end subroutine kindiag
