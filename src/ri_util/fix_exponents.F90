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

! This subroutine should be in a module, to avoid explicit interfaces
#ifdef _IN_MODULE_

subroutine Fix_Exponents(nP,mP,nC,Expn,CoeffC,CoeffP)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nP, nC
integer(kind=iwp), intent(out) :: mP
real(kind=wp), allocatable, intent(inout) :: Expn(:), CoeffC(:,:,:), CoeffP(:,:,:)
integer(kind=iwp) :: i, iC, iP, jP, iSkip
real(kind=wp) :: Temp, Thr_Skip
real(kind=wp), allocatable :: Scr1(:), Scr2(:,:,:)

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Fix_Exponents: Expn',' ',Expn,1,nP)
call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(:,:,1),nP,nC)
call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(:,:,2),nP,nC)
call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(:,:,1),nP,nP)
call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(:,:,2),nP,nP)
#endif

mP = nP

! First, put the exponents with all coefficients less than the
! threshold, Thr_Skip, at the end.

Thr_Skip = 1.0e-13_wp
do iP=nP,1,-1

  iSkip = 1
  do iC=1,nC
    if (abs(CoeffC(iP,iC,1)) >= Thr_Skip) iSkip = 0
  end do

  if (iSkip == 1) then
    if (iP < mP) then
      Temp = Expn(iP)
      Expn(iP) = Expn(mP)
      Expn(mP) = Temp
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
    if (Expn(jP) > Expn(ip)) then
      Temp = Expn(iP)
      Expn(iP) = Expn(jP)
      Expn(jP) = Temp
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

#ifdef _DEBUGPRINT_
write(u6,*) 'After Fix_Exp'
call RecPrt('Fix_Exponents: Expn',' ',Expn,1,nP)
call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(:,:,1),nP,nC)
call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(:,:,2),nP,nC)
call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(:,:,1),nP,nP)
call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(:,:,2),nP,nP)
#endif

! Reallocate arrays if the number of primitives is reduced.

if (mP /= nP) then
  call mma_allocate(Scr1,mP,Label='Expn')
  Scr1(:) = Expn(1:mP)
  call mma_deallocate(Expn)
  call move_alloc(Scr1,Expn)

  call mma_allocate(Scr2,mP,nC,2,Label='CoeffC')
  Scr2(:,:,:) = CoeffC(1:mP,1:nC,:)
  call mma_deallocate(CoeffC)
  call move_alloc(Scr2,CoeffC)

  call mma_allocate(Scr2,mP,mP,2,Label='CoeffP')
  Scr2(:,:,:) = CoeffP(1:mP,1:mP,:)
  call mma_deallocate(CoeffP)
  call move_alloc(Scr2,CoeffP)

# ifdef _DEBUGPRINT_
  write(u6,*) 'After Reallocation'
  call RecPrt('Fix_Exponents: Expn',' ',Expn,1,mP)
  call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(:,:,1),mP,nC)
  call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(:,:,2),mP,nC)
  call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(:,:,1),mP,mP)
  call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(:,:,2),mP,mP)
else
  write(u6,*) 'No Reallocation'
# endif
end if

return

end subroutine Fix_Exponents

#endif
