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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine prints the body part to the postscript figure of        *
! occupation numbers.                                                  *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine FigPrt(lu,Title,nLqn,nOrb,Occ)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lu, nLqn, nOrb(0:nLqn)
character(len=*), intent(in) :: Title
real(kind=wp), intent(in) :: Occ(*)
integer(kind=iwp) :: i, k0, k1, lqn0, lqn1, mxNqn, n, nOff, iOff(0:25)
real(kind=wp) :: O, Oa, Ob
character, parameter :: tp(0:5) = ['s','p','d','f','g','h']

write(lu,'(a)') 'HF setfont'
write(lu,'(a)') '2.25 X -7.0 Y moveto'
k0 = 0
k1 = 0
do i=len(Title),1,-1
  if (Title(i:i) /= ' ') k0 = i
end do
do i=1,len(Title)
  if (Title(i:i) /= ' ') k1 = i
end do
write(lu,'(3a)') '(',Title(k0:k1),') CenterLine'
nOff = 0
mxNqn = 0
do i=0,nLqn
  if (nOrb(i) > 0) mxNqn = max(mxNqn,nOrb(i)+i)
  iOff(i) = nOff
  nOff = nOff+nOrb(i)
end do
do n=1,mxNqn
  write(lu,'(a,i2)') '%--- shell n=',n
  lqn0 = nLqn
  lqn1 = 0
  do i=0,min(n-1,nLqn)
    if (n-i <= nOrb(i)) then
      if (i > lqn1) lqn1 = i
      if (i < lqn0) lqn0 = i
      O = log10(Occ(iOff(i)+n-i))
      write(lu,'(i2,a,f7.4,a)') i,' X ',O,' Y Draw'
    end if
  end do
  if (n < 10) then
    write(lu,'(a,i1,a,a,i1,a,a)') ' (',n,tp(lqn0),'\261',n,tp(lqn1),') Label'
  else
    write(lu,'(a,i2,a,a,i2,a,a)') ' (',n,tp(lqn0),'\261',n,tp(lqn1),') Label'
  end if
  do i=1,min(n-1,nLqn)
    if ((n-i <= nOrb(i)).and.(n-i+1 <= nOrb(i-1))) then
      Oa = log10(Occ(iOff(i-1)+n-i+1))
      Ob = log10(Occ(iOff(i)+n-i))
      write(lu,'(i2,a,f7.4,a,i2,a,f7.4,a)') i-1,' X ',Oa,' Y ',i,' X ',Ob,' Y Connect'
    end if
  end do
end do

return

end subroutine FigPrt
