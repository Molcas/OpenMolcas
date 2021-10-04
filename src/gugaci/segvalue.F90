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

subroutine stermha4(w,ww,ind1,jbr)

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    ! case a&l
    ! case a&r
    w = fq
  case (2)
    w = One
  case (3)
    w = sqrt(b/(b+One))
  case (4)
    w = -fq*sqrt((b+Two)/(b+One))
    !if (abs(w) > 1.0e-13_wp) then
end select
ww = w

return

end subroutine stermha4

subroutine stermhd1(w,ww,ind1,jbr)

use gugaci_global, only: v_onevsqtwo, v_sqtwo
use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default !(1)
    ! d1: case d&r&l
    w = -fq*v_onevsqtwo
    ww = -fq*sqrt((b-One)/(b+b+Two))
  case (2)
    ww = -sqrt(b/(b+One))
    !if (dldr == 2101) ww=(b+Two)/(b+One)
  case (3)
    w = -fq*v_onevsqtwo
    ww = fq*sqrt((b+Three)/(b+b+Two))
  case (4)
    w = fq*v_sqtwo
end select

return

end subroutine stermhd1

subroutine stermhd5(w,ww)

use gugaci_global, only: v_sqtwo
use Constants, only: Zero
use Definitions, only: wp

implicit none
real(kind=wp), intent(out) :: w, ww

w = -v_sqtwo
ww = Zero

return

end subroutine stermhd5

subroutine stermla1(w,ww,ind1,jbr)
! case a^l

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    w = -sqrt((b+1)/(b+2))
  case (2)
    w = -One
  case (3)
    w = fq*sqrt((b+1)/b)
  case (4)
    w = fq
end select
ww = w

return

end subroutine stermla1

subroutine stermla2(w,ww,ind1,jbr)
! case a^r

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    w = -fq
  case (2)
    w = -fq*sqrt((b+One)/(b+Two))
  case (3)
    w = One
  case (4)
    w = sqrt((b+One)/b)
end select
ww = w

return

end subroutine stermla2

subroutine stermld2(w,ww,ind1,jbr)
! d2: case d^r^l

use gugaci_global, only: v_onevsqtwo, v_sqtwo
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case (1)
    ww = -sqrt((b+One)/(b+Two))
  case default ! (2)
    w = -fq*v_onevsqtwo
    ww = fq*sqrt((b+Two)/(b*Two))
  case (3)
    w = -fq*v_onevsqtwo
    ww = -fq*sqrt(b/(b+b+Four))
  case (4)
    w = -fq*v_sqtwo
  case (5)
    ww = sqrt((b+One)/b)
end select

return

end subroutine stermld2

subroutine stermld6(w,ww)
! d6: case d^r^r

use gugaci_global, only: v_sqtwo
use Constants, only: Zero
use Definitions, only: wp

implicit none
real(kind=wp), intent(out) :: w, ww

w = -v_sqtwo
ww = Zero

return

end subroutine stermld6

subroutine segmidc1(w,ww,ind1,jbr)
! case c1

use Constants, only: Zero, One, Three, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    w = One
  case (2)
    w = One
  case (3)
    w = fq/sqrt((b*b+Four*b+Four))
  case (4)
    w = -sqrt((b+One)*(b+Three)/(b*b+Four*b+Four))
  case (5)
    w = -One
  case (6)
    w = One
  case (7)
    w = sqrt((b+One)*(b-One)/(b*b))
  case (8)
    w = fq/b
  case (9)
    w = -One
  case (10)
    w = -One
end select
ww = w

return

end subroutine segmidc1

subroutine segmidc2(w,ww,ind1,jbr)
! case c2

use Constants, only: Zero, One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case (1)
    ww = One
  case (2)
    ww = -One
  case (3)
    ww = -fq*sqrt(Two/((b+Two)*(b+Three)))
  case (4)
    ww = -sqrt((b+One)*(b+Four)/((b+Two)*(b+Three)))
  case (5)
    ww = One
  case default ! (6)
    w = One
    ww = One
  case (7)
    w = -One
    ww = -sqrt((b-One)*(b+Two)/(b*b+b))
  case (8)
    ww = -fq*sqrt(Two/(b*b+Three*b+Two))
  case (9)
    ww = fq*sqrt(Two/(b*(b+One)))
  case (10)
    w = -One
    ww = -sqrt(b*(b+Three)/(b*b+Three*b+Two))
  case (11)
    w = One
    ww = One
  case (12)
    ww = One
  case (13)
    ww = -sqrt((b-Two)*(b+One)/(b*b-b))
  case (14)
    ww = fq*sqrt(Two/(b*b-b))
  case (15)
    ww = -One
  case (16)
    ww = One
end select

return

end subroutine segmidc2

!subroutine segmidc22(w,ww,ind1,jbr)
!! case c22
!
!use Constants, only: Zero, One, Two, Three, Four
!use Definitions, only: wp, iwp
!
!implicit none
!real(kind=wp), intent(out) :: w, ww
!integer(kind=iwp), intent(in) :: ind1, jbr
!real(kind=wp) :: b, fq
!
!w = Zero
!ww = Zero
!! calculate w,ww
!b = real(jbr,kind=wp)
!if (mod(jbr,2) == 0) fq = One
!if (mod(jbr,2) /= 0) fq = -One
!select case (ind1)
!  case (1)
!    ww = One
!  case (2)
!    ww = -One
!  case (3)
!    ww = -fq*sqrt(Two/((b+Two)*(b+Three)))
!  case (4)
!    ww = -sqrt((b+One)*(b+Four)/((b+Two)*(b+Three)))
!  case (5)
!    ww = One
!  case default ! (6)
!    w = One
!    ww = One
!  case (7)
!    w = -One
!    ww = -sqrt((b-One)*(b+Two)/(b*b+b))
!  case (8)
!    ww = -fq*sqrt(Two/(b*b+Three*b+Two))
!  case (9)
!    w = -One
!    ww = -sqrt(b*(b+Three)/(b*b+Three*b+Two))
!  case (10)
!    w = One
!    ww = One
!  case (11)
!    ww = One
!  case (12)
!    ww = -sqrt((b-Two)*(b+One)/(b*b-b))
!  case (13)
!    ww = -One
!  case (14)
!    ww = One
!end select
!
!return
!
!end subroutine segmidc22

subroutine segmidb3(w,ww,ind1,jbr)
! submid b3(b&l)

use gugaci_global, only: v_onevsqtwo
use Constants, only: Zero, One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    w = -v_onevsqtwo
    ww = -sqrt(b/(b+b+Four))
  case (2)
    ww = -fq*sqrt((b+Three)/(b+Two))
  case (3)
    ww = fq
  case (4)
    w = sqrt((b+One)/(b+b+Four))
    ww = sqrt((b+Three)/(b+b+Four))
  case (5)
    ww = -sqrt((b-One)/b)
  case (6)
    w = fq*v_onevsqtwo
    ww = -fq*sqrt((b+Two)/(b+b))
  case (7)
    w = -fq*sqrt((b+One)/(b+b))
    ww = fq*sqrt((b-One)/(b+b))
  case (8)
    ww = One
end select

return

end subroutine segmidb3

subroutine segmidb4(w,ww,ind1,jbr)
! segmid b4(b&r)

use gugaci_global, only: v_onevsqtwo
use Constants, only: Zero, One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case (1)
    ww = -One
  case (2)
    w = fq*sqrt((b+One)/(b+b+Four))
    ww = -fq*sqrt((b+Three)/(b+b+Four))
  case (3)
    w = -fq*v_onevsqtwo
    ww = fq*sqrt(b/(b+b+Four))
  case (4)
    ww = sqrt((b+Three)/(b+Two))
  case default ! (5)
    w = -sqrt((b+One)/(b+b))
    ww = -sqrt((b-One)/(b+b))
  case (6)
    ww = -fq
  case (7)
    ww = fq*sqrt((b-One)/b)
  case (8)
    w = v_onevsqtwo
    ww = sqrt((b+Two)/(b+b))
end select

return

end subroutine segmidb4

subroutine segmidd10(w,ww,ind1,jbr)
! segmid d10(d^r&l)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: fq

w = Zero
ww = Zero
! calculate w,ww
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case (1)
    w = fq
  case default ! (2)
    w = fq
end select

return

end subroutine segmidd10

subroutine segmidb2(w,ww,ind1,jbr)
! segmid b2(b^r)

use gugaci_global, only: v_onevsqtwo
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case (1)
    ww = One
  case (2)
    ww = -sqrt((b+One)/(b+Two))
  case default ! (3)
    w = -sqrt(b/(b+b+Two))
    ww = sqrt((b+Two)/(b+b+Two))
  case (4)
    w = fq*sqrt((b+Two)/(b+b+Two))
    ww = fq*sqrt(b/(b+b+Two))
  case (5)
    w = fq*v_onevsqtwo
    ww = fq*sqrt((b+Two)/(b+b))
  case (6)
    w = v_onevsqtwo
    ww = -sqrt(b/(b+b+Four))
  case (7)
    ww = fq
  case (8)
    ww = fq*sqrt((b+One)/b)
end select

return

end subroutine segmidb2

subroutine segmidb1(w,ww,ind1,jbr)
! segmid b1(b^l)

use gugaci_global, only: v_onevsqtwo
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case (1)
    ww = -fq*sqrt((b+One)/(b+Two))
  case (2)
    ww = -fq
  case default ! (3)
    w = -v_onevsqtwo
    ww = sqrt((b+Two)/(b+b))
  case (4)
    w = -fq*v_onevsqtwo
    ww = -fq*sqrt(b/(b+b+Four))
  case (5)
    w = -fq*sqrt(b/(b+b+Two))
    ww = -fq*sqrt((b+Two)/(b+b+Two))
  case (6)
    w = sqrt((b+Two)/(b+b+Two))
    ww = -sqrt(b/(b+b+Two))
  case (7)
    ww = sqrt((b+One)/b)
  case (8)
    ww = -One
end select

return

end subroutine segmidb1

!subroutine stermh(isq,w,ww,ind1,jbr)
!
!use gugaci_global, only: v_onevsqtwo, v_sqtwo
!use Constants, only: Zero, One, Two, Three
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(out) :: isq
!real(kind=wp), intent(out) :: w, ww
!integer(kind=iwp), intent(in) :: ind1, jbr
!real(kind=wp) :: b, fq
!
!isq = 0
!w = Zero
!ww = Zero
!! calculate w,ww
!b = real(jbr,kind=wp)
!if (mod(jbr,2) == 0) fq = One
!if (mod(jbr,2) /= 0) fq = -One
!select case (ind1)
!  case default ! (1)
!    ! case a&l
!    ! case a&r
!    w = fq
!    isq = 1
!    ww = w
!  case (2)
!    w = One
!    isq = 1
!    ww = w
!  case (3)
!    ! case d&l&l
!    ! case d&r&r
!    w = -v_sqtwo
!    isq = 3
!  case (4)
!    ! case d&r&l
!    w = -fq*v_onevsqtwo
!    ww = -fq*sqrt((b-One)/(b+b+Two))
!    isq = 2
!  case (5)
!    ww = -sqrt(b/(b+One))
!    !if (dldr == 2101) ww=(b+Two)/(b+One)
!    isq = 2
!  case (6)
!    w = sqrt(b/(b+One))
!    isq = 1
!    ww = w
!  case (7)
!    w = -fq*sqrt((b+Two)/(b+One))
!    !if (abs(w) > 1.0e-13_wp) then
!    isq = 1
!    ww = w
!  case (8)
!    w = -fq*v_onevsqtwo
!    ww = fq*sqrt((b+Three)/(b+b+Two))
!    isq = 2
!  case (9)
!    w = fq*v_sqtwo
!    isq = 2
!end select
!
!return
!
!end subroutine stermh

subroutine stmh(isq,w,ww,mw,ind1,jbr)

use gugaci_global, only: v_onevsqtwo, v_sqtwo
use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: isq, mw
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

mw = 0
isq = 0
w = Zero
ww = Zero
! calculate w,w
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    w = -fq*v_onevsqtwo
    ww = -fq*sqrt((b-One)/(b+b+Two))
  case (2)
    w = -fq*v_onevsqtwo
    ww = fq*sqrt((b+Three)/(b+b+Two))
  case (3)
    w = fq*v_sqtwo
end select
if (abs(ww) > 1.0e-13_wp) mw = 2
if (abs(w) > 1.0e-13_wp) mw = mw+1
isq = 401

return

end subroutine stmh

subroutine smidc2(isq,w,ww,mw,ind1,jbr)

use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: isq, mw
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b

w = Zero
ww = Zero
isq = 0
mw = 0
! calculate w,ww
b = real(jbr,kind=wp)
select case (ind1)
  case default ! (1)
    !write(6,*) 'this is case c2'
    w = One
    ww = One
  case (2)
    w = -One
    ww = -sqrt((b-One)*(b+Two)/(b*b+b))
  case (3)
    w = -One
    ww = -sqrt(b*(b+Three)/(b*b+Three*b+Two))
  case (4)
    w = One
    ww = One
end select
isq = 302

return

end subroutine smidc2

subroutine stml(isq,w,ww,mw,ind1,jbr)

use gugaci_global, only: v_onevsqtwo, v_sqtwo
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: isq, mw
real(kind=wp), intent(out) :: w, ww
integer(kind=iwp), intent(in) :: ind1, jbr
real(kind=wp) :: b, fq

mw = 0
isq = 0
w = Zero
ww = Zero
! calculate w,ww
b = real(jbr,kind=wp)
if (mod(jbr,2) == 0) fq = One
if (mod(jbr,2) /= 0) fq = -One
select case (ind1)
  case default ! (1)
    ! case d^r^l
    w = -fq*v_onevsqtwo
    ww = fq*sqrt((b+Two)/(b*Two))
  case (2)
    w = -fq*v_onevsqtwo
    ww = -fq*sqrt(b/(b+b+Four))
  case (3)
    w = -fq*v_sqtwo
end select
if (abs(ww) > 1.0e-13_wp) mw = 2
if (abs(w) > 1.0e-13_wp) mw = mw+1
isq = 402

return

end subroutine stml

subroutine neoc(kcoe,nocc,tcoe)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: kcoe
integer(kind=iwp), intent(out) :: nocc
real(kind=wp), intent(out) :: tcoe

nocc = 1
tcoe = kcoe
if (kcoe == 0) nocc = 0
if (kcoe == 100) tcoe = Zero
if (kcoe == 200) then
  nocc = 2
  tcoe = -Half
end if

return

end subroutine neoc
