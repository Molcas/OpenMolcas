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

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, two = 2.0d0

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case default ! (1)
    ! case a&l
    ! case a&r
    w = fq
  case (2)
    w = done
  case (3)
    w = sqrt(b/(b+done))
  case (4)
    w = -fq*sqrt((b+two)/(b+done))
    !if (abs(w) > 1.e-13) then
end select
ww = w

return

end subroutine stermha4

subroutine stermhd1(w,ww,ind1,jbr)

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, three = 3.0d0, two = 2.0d0

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case default !(1)
    ! d1: case d&r&l
    w = -fq/sqrt(two)
    ww = -fq*sqrt((b-done)/(b+b+two))
  case (2)
    ww = -sqrt(b/(b+done))
    !if (dldr == 2101) ww=(b+two)/(b+done)
  case (3)
    w = -fq/sqrt(two)
    ww = fq*sqrt((b+three)/(b+b+two))
  case (4)
    w = fq*sqrt(two)
end select

return

end subroutine stermhd1

subroutine stermhd5(w,ww)

implicit none
real*8 :: w, ww
real*8, parameter :: two = 2.0d0

w = 0.0d0
ww = 0.0d0
! calculate w,ww
! d5: case d&r&r
w = -sqrt(two)

return

end subroutine stermhd5

subroutine stermla1(w,ww,ind1,jbr)
! case a^l

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8 :: done = 1.0d0

w = 0.0d0
ww = 0.0d0
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case default ! (1)
    w = -sqrt((b+1)/(b+2))
  case (2)
    w = -done
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

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, two = 2.0d0

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case default ! (1)
    w = -fq
  case (2)
    w = -fq*sqrt((b+done)/(b+two))
  case (3)
    w = done
  case (4)
    w = sqrt((b+done)/b)
end select
ww = w

return

end subroutine stermla2

subroutine stermld2(w,ww,ind1,jbr)
! d2: case d^r^l

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, two = 2.0d0

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case (1)
    ww = -sqrt((b+done)/(b+two))
  case default ! (2)
    w = -fq/sqrt(two)
    ww = fq*sqrt((b+two)/(b*two))
  case (3)
    w = -fq/sqrt(two)
    ww = -fq*sqrt(b/(b+b+4.0d0))
  case (4)
    w = -fq*sqrt(two)
  case (5)
    ww = sqrt((b+done)/b)
end select

return

end subroutine stermld2

subroutine stermld6(w,ww)
! d6: case d^r^r

implicit none
real*8 :: w, ww
real*8, parameter :: dzero = 0.0d0, two = 2.0d0

w = -sqrt(two)
ww = dzero

return

end subroutine stermld6

subroutine segmidc1(w,ww,ind1,jbr)
! case c1

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0

b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
select case (ind1)
  case default ! (1)
    w = done
  case (2)
    w = done
  case (3)
    w = fq/sqrt((b*b+4.0d0*b+4.0d0))
  case (4)
    w = -sqrt((b+done)*(b+three)/(b*b+4.0d0*b+4.0d0))
  case (5)
    w = -done
  case (6)
    w = done
  case (7)
    w = sqrt((b+done)*(b-done)/(b*b))
  case (8)
    w = fq/b
  case (9)
    w = -done
  case (10)
    w = -done
end select
ww = w

return

end subroutine segmidc1

subroutine segmidc2(w,ww,ind1,jbr)
! case c2

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0

b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
select case (ind1)
  case (1)
    ww = done
  case (2)
    ww = -done
  case (3)
    ww = -fq*sqrt(two/((b+two)*(b+three)))
  case (4)
    ww = -sqrt((b+done)*(b+4.0d0)/((b+two)*(b+three)))
  case (5)
    ww = done
  case default ! (6)
    w = done
    ww = done
  case (7)
    w = -done
    ww = -sqrt((b-done)*(b+two)/(b*b+b))
  case (8)
    ww = -fq*sqrt(two/(b*b+three*b+two))
  case (9)
    ww = fq*sqrt(two/(b*(b+done)))
  case (10)
    w = -done
    ww = -sqrt(b*(b+three)/(b*b+three*b+two))
  case (11)
    w = done
    ww = done
  case (12)
    ww = done
  case (13)
    ww = -sqrt((b-two)*(b+done)/(b*b-b))
  case (14)
    ww = fq*sqrt(two/(b*b-b))
  case (15)
    ww = -done
  case (16)
    ww = done
end select

return

end subroutine segmidc2

!subroutine segmidc22(w,ww,ind1,jbr)
!! case c22
!
!implicit none
!real*8 :: w, ww
!integer :: ind1, jbr
!real*8 :: b, fq
!real*8, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0
!
!b = dble(jbr)
!if (mod(jbr,2) == 0) fq = done
!if (mod(jbr,2) /= 0) fq = -done
!w = dzero
!ww = dzero
!select case (ind1)
!  case (1)
!    ww = done
!  case (2)
!    ww = -done
!  case (3)
!    ww = -fq*sqrt(two/((b+two)*(b+three)))
!  case (4)
!    ww = -sqrt((b+done)*(b+4.0d0)/((b+two)*(b+three)))
!  case (5)
!    ww = done
!  case default ! (6)
!    w = done
!    ww = done
!  case (7)
!    w = -done
!    ww = -sqrt((b-done)*(b+two)/(b*b+b))
!  case (8)
!    ww = -fq*sqrt(two/(b*b+three*b+two))
!  case (9)
!    w = -done
!    ww = -sqrt(b*(b+three)/(b*b+three*b+two))
!  case (10)
!    w = done
!    ww = done
!  case (11)
!    ww = done
!  case (12)
!    ww = -sqrt((b-two)*(b+done)/(b*b-b))
!  case (13)
!    ww = -done
!  case (14)
!    ww = done
!end select
!
!return
!
!end subroutine segmidc22

subroutine segmidb3(w,ww,ind1,jbr)
! submid b3(b&l)

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0, vtwo = 0.5d0

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
select case (ind1)
  case default ! (1)
    w = -sqrt(vtwo)
    ww = -sqrt(b/(b+b+4.0d0))
  case (2)
    ww = -fq*sqrt((b+three)/(b+two))
  case (3)
    ww = fq
  case (4)
    w = sqrt((b+done)/(b+b+4.0d0))
    ww = sqrt((b+three)/(b+b+4.0d0))
  case (5)
    ww = -sqrt((b-done)/b)
  case (6)
    w = fq*sqrt(vtwo)
    ww = -fq*sqrt((b+two)/(b+b))
  case (7)
    w = -fq*sqrt((b+done)/(b+b))
    ww = fq*sqrt((b-done)/(b+b))
  case (8)
    ww = done
end select

return

end subroutine segmidb3

subroutine segmidb4(w,ww,ind1,jbr)
! segmid b4(b&r)

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0, vtwo = 0.5d0

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
select case (ind1)
  case (1)
    ww = -done
  case (2)
    w = fq*sqrt((b+done)/(b+b+4.0d0))
    ww = -fq*sqrt((b+three)/(b+b+4.0d0))
  case (3)
    w = -fq*sqrt(vtwo)
    ww = fq*sqrt(b/(b+b+4.0d0))
  case (4)
    ww = sqrt((b+three)/(b+two))
  case default ! (5)
    w = -sqrt((b+done)/(b+b))
    ww = -sqrt((b-done)/(b+b))
  case (6)
    ww = -fq
  case (7)
    ww = fq*sqrt((b-done)/b)
  case (8)
    w = sqrt(vtwo)
    ww = sqrt((b+two)/(b+b))
end select

return

end subroutine segmidb4

subroutine segmidd10(w,ww,ind1,jbr)
! segmid d10(d^r&l)

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0

! calculate w,ww
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
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

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, two = 2.0d0, vtwo = 0.5d0

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
select case (ind1)
  case (1)
    ww = done
  case (2)
    ww = -sqrt((b+done)/(b+two))
  case default ! (3)
    w = -sqrt(b/(b+b+two))
    ww = sqrt((b+two)/(b+b+two))
  case (4)
    w = fq*sqrt((b+two)/(b+b+two))
    ww = fq*sqrt(b/(b+b+two))
  case (5)
    w = fq*sqrt(vtwo)
    ww = fq*sqrt((b+two)/(b+b))
  case (6)
    w = sqrt(vtwo)
    ww = -sqrt(b/(b+b+4.0d0))
  case (7)
    ww = fq
  case (8)
    ww = fq*sqrt((b+done)/b)
end select

return

end subroutine segmidb2

subroutine segmidb1(w,ww,ind1,jbr)
! segmid b1(b^l)

implicit none
real*8 :: w, ww
integer :: ind1, jbr
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, two = 2.0d0, vtwo = 0.5d0

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
select case (ind1)
  case (1)
    ww = -fq*sqrt((b+done)/(b+two))
  case (2)
    ww = -fq
  case default ! (3)
    w = -sqrt(vtwo)
    ww = sqrt((b+two)/(b+b))
  case (4)
    w = -fq*sqrt(vtwo)
    ww = -fq*sqrt(b/(b+b+4.0d0))
  case (5)
    w = -fq*sqrt(b/(b+b+two))
    ww = -fq*sqrt((b+two)/(b+b+two))
  case (6)
    w = sqrt((b+two)/(b+b+two))
    ww = -sqrt(b/(b+b+two))
  case (7)
    ww = sqrt((b+done)/b)
  case (8)
    ww = -done
end select

return

end subroutine segmidb1

!subroutine stermh(isq,w,ww,ind1,jbr)
!
!implicit none
!integer :: isq, ind1, jbr
!real*8 :: w, ww
!real*8 :: b, fq
!real*8, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0, vtwo = 0.5d0
!
!isq = 0
!w = dzero
!ww = dzero
!! calculate w,ww
!b = dble(jbr)
!if (mod(jbr,2) == 0) fq = done
!if (mod(jbr,2) /= 0) fq = -done
!select case (ind1)
!  case default ! (1)
!    ! case a&l
!    ! case a&r
!    w = fq
!    isq = 1
!    ww = w
!  case (2)
!    w = done
!    isq = 1
!    ww = w
!  case (3)
!    ! case d&l&l
!    ! case d&r&r
!    w = -sqrt(two)
!    isq = 3
!  case (4)
!    ! case d&r&l
!    w = -fq*sqrt(vtwo)
!    ww = -fq*sqrt((b-done)/(b+b+two))
!    isq = 2
!  case (5)
!    ww = -sqrt(b/(b+done))
!    !if (dldr == 2101) ww=(b+two)/(b+done)
!    isq = 2
!  case (6)
!    w = sqrt(b/(b+done))
!    isq = 1
!    ww = w
!  case (7)
!    w = -fq*sqrt((b+two)/(b+done))
!    !if (abs(w) > 1.e-13) then
!    isq = 1
!    ww = w
!  case (8)
!    w = -fq*sqrt(vtwo)
!    ww = fq*sqrt((b+three)/(b+b+two))
!    isq = 2
!  case (9)
!    w = fq*sqrt(two)
!    isq = 2
!end select
!
!return
!
!end subroutine stermh

subroutine stmh(isq,w,ww,mw,ind1,jbr)

implicit none
integer :: isq, mw, ind1, jbr
real*8 :: w, ww
real*8 :: b, fq
real*8 :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0, vtwo = 0.5d0

mw = 0
isq = 0
w = dzero
ww = dzero
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case default ! (1)
    w = -fq*sqrt(vtwo)
    ww = -fq*sqrt((b-done)/(b+b+two))
  case (2)
    w = -fq*sqrt(vtwo)
    ww = fq*sqrt((b+three)/(b+b+two))
  case (3)
    w = fq*sqrt(two)
end select
if (abs(ww) > 1.d-13) mw = 2
if (abs(w) > 1.d-13) mw = mw+1
isq = 401

return

end subroutine stmh

subroutine smidc2(isq,w,ww,mw,ind1,jbr)

implicit none
integer :: isq, mw, ind1, jbr
real*8 :: w, ww
real*8 :: b
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, three = 3.0d0, two = 2.0d0

! calculate w,ww
b = dble(jbr)
w = dzero
ww = dzero
isq = 0
mw = 0
select case (ind1)
  case default ! (1)
    !write(6,*) 'this is case c2'
    w = done
    ww = done
  case (2)
    w = -done
    ww = -sqrt((b-done)*(b+two)/(b*b+b))
  case (3)
    w = -done
    ww = -sqrt(b*(b+three)/(b*b+three*b+two))
  case (4)
    w = done
    ww = done
end select
isq = 302

return

end subroutine smidc2

subroutine stml(isq,w,ww,mw,ind1,jbr)

implicit none
integer :: isq, mw, ind1, jbr
real*8 :: w, ww
real*8 :: b, fq
real*8, parameter :: done = 1.0d0, dzero = 0.0d0, two = 2.0d0, vtwo = 0.5d0

mw = 0
isq = 0
w = dzero
ww = dzero
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
select case (ind1)
  case default ! (1)
    ! case d^r^l
    w = -fq*sqrt(vtwo)
    ww = fq*sqrt((b+two)/(b*two))
  case (2)
    w = -fq*sqrt(vtwo)
    ww = -fq*sqrt(b/(b+b+4.0d0))
  case (3)
    w = -fq*sqrt(two)
end select
if (abs(ww) > 1.d-13) mw = 2
if (abs(w) > 1.d-13) mw = mw+1
isq = 402

return

end subroutine stml

subroutine neoc(kcoe,nocc,tcoe)

implicit none
integer :: kcoe, nocc
real*8 :: tcoe

nocc = 1
tcoe = kcoe
if (kcoe == 0) nocc = 0
if (kcoe == 100) tcoe = 0.d0
if (kcoe == 200) then
  nocc = 2
  tcoe = -0.5d0
end if

return

end subroutine neoc
