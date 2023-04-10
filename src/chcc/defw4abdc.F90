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

subroutine DefW4abdc(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (ab|dc)

implicit none
integer dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
real*8 Wx(1:abLen,1:cdLen)
! help variables
integer a, b, c, d, ab, dc

if (aSGrp == bSGrp) then
  ! case (a=b|d,c)
  dc = 0
  do c=1,dimc
    do d=1,dimd
      dc = dc+1
      do a=2,dima
        ab = a*(a-1)/2
        do b=1,a-1
          ab = ab+1
          W(a,b,c,d) = W(a,b,c,d)+Wx(ab,dc)
          W(b,a,c,d) = W(b,a,c,d)+Wx(ab,dc)
        end do
      end do
      do a=1,dima
        ab = a*(a+1)/2
        W(a,a,c,d) = W(a,a,c,d)+Wx(ab,dc)
      end do
    end do
  end do

else
  ! case (a,b|d,c)
  dc = 0
  do c=1,dimc
    do d=1,dimd
      dc = dc+1
      ab = 0
      do b=1,dimb
        do a=1,dima
          ab = ab+1
          W(a,b,c,d) = W(a,b,c,d)+Wx(ab,dc)
        end do
      end do
    end do
  end do

end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(cSGrp)
  call Unused_integer(dSGrp)
end if

end subroutine DefW4abdc
