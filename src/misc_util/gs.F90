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

subroutine GS(drdq,nLambda,T,nInter,Swap,RD)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8 drdq(nInter,nLambda), T(nInter,nInter)
logical Swap, RD
real*8, allocatable :: Temp(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_

! Be careful here so that noise is not converted to a basis!

Thr = 1.0D-12
#ifdef _DEBUGPRINT_
call RecPrt('GS: dRdQ',' ',drdq,nInter,nLambda)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Initial check to see if the dRdQ vectors are independent.

call mma_Allocate(Temp,nInter,nLambda,Label='Temp')
Temp(:,:) = drdq(:,:)
call GS_(drdq,nInter,nLambda,Thr)
jLambda = 0
do i=1,nLambda
  XX = sqrt(DDot_(nInter,drdq(1,i),1,drdq(1,i),1))
  !write(6,*) 'i,XX=',i,XX
  if (XX > Thr) then
    jLambda = jLambda+1
    ! RD = remove degeneracies
    if (RD) then
      if (jLambda /= i) then
        call dCopy_(nInter,drdq(1,i),1,drdq(1,jLambda),1)
      end if
    end if
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('GS: dRdQ(orth)',' ',drdq,nInter,nLambda)
#endif
if ((.not. RD) .and. (jLambda /= nLambda)) then
  write(6,*) ' Constraints are linear dependent!'
  call abend()
end if
nLambda = jLambda
!                                                                      *
!***********************************************************************
!                                                                      *
! Project away the space which is spanned by the constraints.

call FZero(T,nInter**2)

call dcopy_(nInter,[One],0,T,1+nInter)
#ifdef _DEBUGPRINT_
call RecPrt('T(orig)',' ',T,nInter,nInter)
#endif

! Form 1 - P

do iLambda=1,nLambda
  do iInter=1,nInter
    do jInter=1,nInter
      T(iInter,jInter) = T(iInter,jInter)-drdq(iInter,iLambda)*drdq(jInter,iLambda)
    end do
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('1-P',' ',T,nInter,nInter)
#endif

! Orthonormalize

call GS_(T,nInter,nInter,Thr)

! Set the trailing vectors to null vectors.
! We might have noise here.

if (nLambda /= 0) then
  iStart = nInter-nLambda+1
  call FZero(T(1,iStart),nInter*nLambda)
end if
#ifdef _DEBUGPRINT_
call RecPrt('1-P(GS)',' ',T,nInter,nInter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Restore dRdQ

if (.not. RD) then
  call dcopy_(nInter*nLambda,Temp,1,drdq,1)
end if
call mma_deallocate(Temp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Reorder, place the null vectors at the start.

j = nInter
do i=nInter,1,-1
  XX = DDot_(nInter,T(1,i),1,T(1,i),1)
  if ((XX > Zero) .and. (i /= j)) then
    call dcopy_(nInter,T(1,i),1,T(1,j),1)
    j = j-1
  else if (XX > Zero) then
    j = j-1
  end if
end do

! Put drdq at the start

call dcopy_(nInter*nLambda,drdq,1,T,1)
#ifdef _DEBUGPRINT_
call RecPrt('T(ReOrdered)',' ',T,nInter,nInter)
#endif

!call GS_(T,nInter,nInter,Thr)
if (Swap) call dswap_(nInter,T(1,1),1,T(1,3),1)

return

end subroutine GS
