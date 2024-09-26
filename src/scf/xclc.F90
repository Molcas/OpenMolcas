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

!#define _DEBUGPRINT_
subroutine XClc()
! Compute the x parameters value as a function of Iter_ref

use LnkLst, only: GetNod, iVPtr, LLx, LstPtr, PutVec, SCF_V
use InfSCF, only: Iter, Iter_Ref, Iter_Start, mOV
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, inode, jpgrd
real(kind=wp), allocatable :: Scr(:)

call mma_allocate(Scr,mOV,Label='Scr')

jpgrd = LstPtr(Iter_Ref,LLx)   ! Pointer to X_old(i_ref)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'iter=',iter
write(u6,*) 'iter_Start=',iter_Start
write(u6,*) 'iter_ref=',iter_ref
call NrmClc(SCF_V(jpgrd)%A(:),mOV,'XClc','X(i_ref)(:)')
write(u6,*)
#endif

! Loop over all iterations starting at Iter_Start+1

do i=Iter_Start,Iter

  !X_new(i) = X_old(i)-X_old(i_ref)

  call GetNod(i,LLx,inode)
  if (inode == 0) then
    write(u6,*) 'inode == 0'
    call Abend()
  end if
  call iVPtr(Scr,mOV,inode)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'X(i) before  i=',i
  call NrmClc(Scr(:),mOV,'XClc','Scr(:)')
# endif
  Scr(:) = Scr(:)-SCF_V(jpgrd)%A(:)
# ifdef _DEBUGPRINT_
  write(u6,*) 'X(i) after  i=',i
  call NrmClc(Scr(:),mOV,'XClc','Scr(:)')
# endif

  call PutVec(Scr,mOV,i,'OVWR',LLx)
end do

call mma_deallocate(Scr)

return

end subroutine XClc
