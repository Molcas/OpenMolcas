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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: nLambda
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(inout) :: drdq(nInter,nLambda)
real(kind=wp), intent(out) :: T(nInter,nInter)
logical(kind=iwp), intent(in) :: Swap, RD
integer(kind=iwp) :: i, iInter, iLambda, iStart, j, jLambda
real(kind=wp) :: XX
real(kind=wp), allocatable :: Temp(:,:)
real(kind=wp), parameter :: Thr = 1.0e-12_wp ! Be careful here so that noise is not converted to a basis!
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_

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
  XX = sqrt(DDot_(nInter,drdq(:,i),1,drdq(:,i),1))
  !write(u6,*) 'i,XX=',i,XX
  if (XX > Thr) then
    jLambda = jLambda+1
    ! RD = remove degeneracies
    if (RD .and. (jLambda /= i)) drdq(:,jLambda) = drdq(:,i)
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('GS: dRdQ(orth)',' ',drdq,nInter,nLambda)
#endif
if ((.not. RD) .and. (jLambda /= nLambda)) then
  write(u6,*) ' Constraints are linear dependent!'
  call Abend()
end if
nLambda = jLambda
!                                                                      *
!***********************************************************************
!                                                                      *
! Project away the space which is spanned by the constraints.

call unitmat(T,nInter)
#ifdef _DEBUGPRINT_
call RecPrt('T(orig)',' ',T,nInter,nInter)
#endif

! Form 1 - P

do iLambda=1,nLambda
  do iInter=1,nInter
    T(:,iInter) = T(:,iInter)-drdq(:,iLambda)*drdq(iInter,iLambda)
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
  T(:,iStart:) = Zero
end if
#ifdef _DEBUGPRINT_
call RecPrt('1-P(GS)',' ',T,nInter,nInter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Restore dRdQ

if (.not. RD) drdq(:,:) = Temp
call mma_deallocate(Temp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Reorder, place the null vectors at the start.

j = nInter
do i=nInter,1,-1
  XX = DDot_(nInter,T(:,i),1,T(:,i),1)
  if (XX > Zero) then
    if (i /= j) T(:,j) = T(:,i)
    j = j-1
  end if
end do

! Put drdq at the start

T(:,1:nLambda) = drdq(:,1:nLambda)
#ifdef _DEBUGPRINT_
call RecPrt('T(ReOrdered)',' ',T,nInter,nInter)
#endif

!call GS_(T,nInter,nInter,Thr)
if (Swap) call dswap_(nInter,T(:,1),1,T(:,3),1)

return

end subroutine GS
