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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine NmHess(dq,nInter,g,nIter,Hess,Delta,q,FEq,Cubic,DipM,dDipM)
!***********************************************************************
!                                                                      *
! Object: to numerically evaluate the molecular Hessian.               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             May '92                                                  *
!***********************************************************************

use Constants, only: Two, Six
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nInter, nIter
real(kind=wp) :: dq(nInter,nIter), g(nInter,nIter), Hess(nInter,nInter), Delta, q(nInter,nIter+1), FEq(nInter,nInter,nInter), &
                 DipM(3,nIter), dDipM(3,nInter)
logical(kind=iwp) :: Cubic
#include "print.fh"
integer(kind=iwp) :: iCount, iInter, iPrint, iRout, jInter, kInter, kIter, kIter1, kIter2, kIter3, kIter4

iRout = 181
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt('NmHess:  g',' ',g,nInter,nIter)
  call RecPrt('NmHess:  q',' ',q,nInter,nIter)
  call RecPrt('NmHess: dq',' ',dq,nInter,nIter)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Form the derivative of the dipole moment

do iInter=1,nInter
  kIter = 1+(iInter-1)*2
  kIter1 = kIter+1
  kIter2 = kIter+2
  !write(u6,*) kIter1,kIter2
  dDipM(1,iInter) = (DipM(1,kIter1)-DipM(1,kIter2))/(Two*Delta)
  dDipM(2,iInter) = (DipM(2,kIter1)-DipM(2,kIter2))/(Two*Delta)
  dDipM(3,iInter) = (DipM(3,kIter1)-DipM(3,kIter2))/(Two*Delta)
end do
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('DipM',' ',DipM,3,nIter)
call RecPrt('dDipM',' ',dDipM,3,nInter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Evaluate the Hessian

do iInter=1,nInter
  do jInter=1,nInter
    kIter = 1+(jInter-1)*2
    kIter1 = kIter+1
    kIter2 = kIter+2
    !write(u6,*) iInter,jInter,kIter1,kIter2
    ! Observe the sign convention due to the use of forces
    ! rather than gradients!!!
    Hess(iInter,jInter) = -(g(iInter,kIter1)-g(iInter,kIter2))/(Two*Delta)
  end do
end do
if (iPrint >= 99) call RecPrt(' Numerical Hessian',' ',Hess,nInter,nInter)

! Symmetrize

do iInter=1,nInter
  do jInter=1,iInter-1
    Hess(iInter,jInter) = (Hess(iInter,jInter)+Hess(jInter,iInter))/Two
    Hess(jInter,iInter) = Hess(iInter,jInter)
  end do
end do
if (iPrint >= 99) call RecPrt(' Symmetrized Hessian',' ',Hess,nInter,nInter)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the numerical cubic constants if data is available.

if (Cubic) then

  ! F(i,j,j)

  do iInter=1,nInter
    do jInter=1,nInter
      kIter = 1+(jInter-1)*2
      kIter1 = kIter+1
      kIter2 = kIter+2
      ! Observe the sign convention due to the use of forces
      ! rather than gradients!!!
      FEq(iInter,jInter,jInter) = -(g(iInter,kIter1)+g(iInter,kIter2))/(Delta**2)
    end do
  end do

  ! F(i,j,k); j>k

  do iInter=1,nInter
    iCount = 0
    do jInter=1,nInter
      do kInter=1,jInter-1
        iCount = iCount+1
        kIter = 1+2*nInter+(iCount-1)*4
        kIter1 = kIter+1
        kIter2 = kIter+2
        kIter3 = kIter+3
        kIter4 = kIter+4
        ! Observe the sign convention due to the use of forces
        ! rather than gradients!!!
        FEq(iInter,jInter,kInter) = -(g(iInter,kIter1)-g(iInter,kIter2)-g(iInter,kIter3)+g(iInter,kIter4))/(Two*Delta)**2
      end do
    end do
  end do

  ! Symmetrize

  do iInter=1,nInter
    do jInter=1,iInter
      do kInter=1,jInter
        FEq(iInter,jInter,kInter) = (FEq(iInter,jInter,kInter)+FEq(iInter,kInter,jInter)+FEq(jInter,iInter,kInter)+ &
                                     FEq(jInter,kInter,iInter)+FEq(kInter,jInter,iInter)+FEq(kInter,iInter,jInter))/Six
        FEq(iInter,jInter,kInter) = FEq(iInter,jInter,kInter)
        FEq(iInter,kInter,jInter) = FEq(iInter,jInter,kInter)
        FEq(jInter,iInter,kInter) = FEq(iInter,jInter,kInter)
        FEq(jInter,kInter,iInter) = FEq(iInter,jInter,kInter)
        FEq(kInter,iInter,jInter) = FEq(iInter,jInter,kInter)
      end do
    end do
  end do

end if

return

end subroutine NmHess
