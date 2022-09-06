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

subroutine Fix_Exponents(nP,mP,nC,Exp,CoeffC,CoeffP)

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
real*8, allocatable :: exp(:), CoeffC(:,:,:), CoeffP(:,:,:)
real*8, allocatable :: Scr(:,:,:)

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Fix_Exponents: Exp',' ',Exp,1,nP)
call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(1,1,1),nP,nC)
call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(1,1,2),nP,nC)
call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(1,1,1),nP,nP)
call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(1,1,2),nP,nP)
#endif

mP = nP
call Fix_Exp()

#ifdef _DEBUGPRINT_
write(6,*) 'After Fix_Exp'
call RecPrt('Fix_Exponents: Exp',' ',Exp,1,nP)
call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(1,1,1),nP,nC)
call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(1,1,2),nP,nC)
call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(1,1,1),nP,nP)
call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(1,1,2),nP,nP)
#endif

! Reallocate arrays if the number of primitives is reduced.

if (mP /= nP) then
  call mma_allocate(Scr,mP,1,1,Label='Scr')
  Scr(1:mP,1,1) = exp(1:mP)
  call mma_deallocate(Exp)
  call mma_allocate(Exp,mP,Label='Exp')
  exp(:) = Scr(:,1,1)
  call mma_deallocate(Scr)

  call mma_allocate(Scr,mP,nC,2,Label='Scr')
  Scr(1:mP,1:nC,:) = CoeffC(1:mP,1:nC,:)
  call mma_deallocate(CoeffC)
  call mma_allocate(CoeffC,mP,nC,2,Label='CoeffC')
  CoeffC(:,:,:) = Scr(:,:,:)
  call mma_deallocate(Scr)

  call mma_allocate(Scr,mP,mP,2,Label='Scr')
  Scr(1:mP,1:mP,:) = CoeffP(1:mP,1:mP,:)
  call mma_deallocate(CoeffP)
  call mma_allocate(CoeffP,mP,mP,2,Label='CoeffP')
  CoeffP(:,:,:) = Scr(:,:,:)
  call mma_deallocate(Scr)
end if

#ifdef _DEBUGPRINT_
write(6,*) 'After Reallocation'
call RecPrt('Fix_Exponents: Exp',' ',Exp,1,mP)
call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(1,1,1),mP,nC)
call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(1,1,2),mP,nC)
call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(1,1,1),mP,mP)
call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(1,1,2),mP,mP)
#endif

return
!                                                                      *
!***********************************************************************
!                                                                      *
contains
!                                                                      *
!***********************************************************************
!                                                                      *
subroutine Fix_Exp()

  ! First, put the exponents with all coefficients less than the
  ! threshold, Thr_Skip, at the end.

  Thr_Skip = 1.0D-13
  do iP=nP,1,-1

    iSkip = 1
    do iC=1,nC
      if (abs(CoeffC(iP,iC,1)) >= Thr_Skip) iSkip = 0
    end do

    if (iSkip == 1) then
      if (iP < mP) then
        Temp = exp(iP)
        exp(iP) = exp(mP)
        exp(mP) = Temp
        do i=1,2
          Temp = CoeffP(iP,iP,i)
          CoeffP(iP,iP,i) = CoeffP(mP,mP,i)
          CoeffP(mP,mP,i) = Temp
          do iC=1,nC
            Temp = CoeffC(iP,iC,i)
            CoeffC(iP,iC,i) = CoeffC(mP,iC,i)
            CoeffC(mP,iC,i) = Temp
          end do
        end do
      end if
      mP = mP-1
    end if

  end do

  ! Second, order from largest to smallest exponent

  do iP=1,mP-1
    do jP=iP+1,mP
      if (exp(jP) > exp(ip)) then
        Temp = exp(iP)
        exp(iP) = exp(jP)
        exp(jP) = Temp
        do i=1,2
          Temp = CoeffP(iP,iP,i)
          CoeffP(iP,iP,i) = CoeffP(jP,jP,i)
          CoeffP(jP,jP,i) = Temp
          do iC=1,nC
            Temp = CoeffC(iP,iC,i)
            CoeffC(iP,iC,i) = CoeffC(jP,iC,i)
            CoeffC(jP,iC,i) = Temp
          end do
        end do
      end if
    end do
  end do

  return

end subroutine Fix_Exp
!                                                                      *
!***********************************************************************
!                                                                      *
end subroutine Fix_Exponents
