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

implicit real*8(a-h,o-z)
data done,two/1.d0,2.d0/

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(34,35,40,44),ind1
! case a&l
! case a&r
34 continue
w = fq
goto 1
35 continue
w = done
goto 1
40 continue
w = sqrt(b/(b+done))
goto 1
44 continue
w = -fq*sqrt((b+two)/(b+done))
!if (abs(w) > 1.e-13) then
1 continue
ww = w

return

end subroutine stermha4

subroutine stermhd1(w,ww,ind1,jbr)

implicit real*8(a-h,o-z)
data done,two,three/1.d0,2.d0,3.d0/

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(38,39,43,48),ind1
! d1: case d&r&l
38 continue
w = -fq/sqrt(two)
ww = -fq*sqrt((b-done)/(b+b+two))
goto 2
43 continue
w = -fq/sqrt(two)
ww = fq*sqrt((b+three)/(b+b+two))
goto 2
48 continue
w = fq*sqrt(two)
goto 2
39 continue
ww = -sqrt(b/(b+done))
!if (dldr == 2101) ww=(b+two)/(b+done)
2 continue

return

end subroutine stermhd1

subroutine stermhd5(w,ww)

implicit real*8(a-h,o-z)
data two/2.d0/

w = 0.0d0
ww = 0.0d0
! calculate w,ww
! d5: case d&r&r
w = -sqrt(two)

return

end subroutine stermhd5

subroutine stermla1(w,ww,ind1,jbr)
! case a^l

implicit real*8(a-h,o-z)
data done/1.d0/

w = 0.0d0
ww = 0.0d0
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(19,24,50,60),ind1
19 continue
w = -sqrt((b+1)/(b+2))
goto 2
24 continue
w = -done
goto 2
50 continue
w = fq*sqrt((b+1)/b)
goto 2
60 continue
w = fq
2 continue
ww = w

return

end subroutine stermla1

subroutine stermla2(w,ww,ind1,jbr)
! case a^r

implicit real*8(a-h,o-z)
data done,two/1.d0,2.d0/

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(21,31,57,62),ind1
21 continue
w = -fq
goto 3
57 continue
w = done
goto 3
62 continue
w = sqrt((b+done)/b)
goto 3
31 continue
w = -fq*sqrt((b+done)/(b+two))
3 continue
ww = w

return

end subroutine stermla2

subroutine stermld2(w,ww,ind1,jbr)
! d2: case d^r^l

implicit real*8(a-h,o-z)
data done,two/1.d0,2.d0/

w = 0.0d0
ww = 0.0d0
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(7,38,43,48,74),ind1
38 continue
w = -fq/sqrt(two)
ww = fq*sqrt((b+two)/(b*two))
goto 4
43 continue
w = -fq/sqrt(two)
ww = -fq*sqrt(b/(b+b+4.0d0))
goto 4
48 continue
w = -fq*sqrt(two)
goto 4
7 continue
ww = -sqrt((b+done)/(b+two))
goto 4
74 continue
ww = sqrt((b+done)/b)
4 continue

return

end subroutine stermld2

subroutine stermld6(w,ww)
! d6: case d^r^r

implicit real*8(a-h,o-z)
data dzero,two/0.d0,2.d0/

w = -sqrt(two)
ww = dzero

return

end subroutine stermld6

subroutine segmidc1(w,ww,ind1,jbr)
! case c1

implicit real*8(a-h,o-z)
data dzero,done,three/0.d0,1.d0,3.d0/

b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(17,22,23,27,32,49,54,58,59,64),ind1
17 continue
w = done
goto 2
22 continue
w = done
goto 2
23 continue
w = fq/sqrt((b*b+4.0d0*b+4.0d0))
goto 2
27 continue
w = -sqrt((b+done)*(b+three)/(b*b+4.0d0*b+4.0d0))
goto 2
32 continue
w = -done
goto 2
49 continue
w = done
goto 2
54 continue
w = sqrt((b+done)*(b-done)/(b*b))
goto 2
58 continue
w = fq/b
goto 2
59 continue
w = -done
goto 2
64 continue
w = -done
2 continue
ww = w

return

end subroutine segmidc1

subroutine segmidc2(w,ww,ind1,jbr)
! case c2

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/

b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(1,6,7,11,16,33,38,39,42,43,48,65,70,74,75,80),ind1
33 continue
w = done
ww = done
goto 2
38 continue
w = -done
ww = -sqrt((b-done)*(b+two)/(b*b+b))
goto 2
43 continue
w = -done
ww = -sqrt(b*(b+three)/(b*b+three*b+two))
goto 2
48 continue
w = done
ww = done
goto 2
39 continue
ww = -fq*sqrt(two/(b*b+three*b+two))
goto 2
7 continue
ww = -fq*sqrt(two/((b+two)*(b+three)))
goto 2
1 continue
ww = done
goto 2
6 continue
ww = -done
goto 2
11 continue
ww = -sqrt((b+done)*(b+4.0d0)/((b+two)*(b+three)))
goto 2
16 continue
ww = done
goto 2
65 continue
ww = done
goto 2
70 continue
ww = -sqrt((b-two)*(b+done)/(b*b-b))
goto 2
75 continue
ww = -done
goto 2
80 continue
ww = done
goto 2
42 continue
ww = fq*sqrt(two/(b*(b+done)))
goto 2
74 continue
ww = fq*sqrt(two/(b*b-b))
2 continue

return

end subroutine segmidc2

subroutine segmidc22(w,ww,ind1,jbr)
! case c22

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/

b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(1,6,7,11,16,33,38,39,43,48,65,70,75,80),ind1
33 continue
w = done
ww = done
goto 2
38 continue
w = -done
ww = -sqrt((b-done)*(b+two)/(b*b+b))
goto 2
43 continue
w = -done
ww = -sqrt(b*(b+three)/(b*b+three*b+two))
goto 2
48 continue
w = done
ww = done
goto 2
39 continue
ww = -fq*sqrt(two/(b*b+three*b+two))
goto 2
7 continue
ww = -fq*sqrt(two/((b+two)*(b+three)))
goto 2
1 continue
ww = done
goto 2
6 continue
ww = -done
goto 2
11 continue
ww = -sqrt((b+done)*(b+4.0d0)/((b+two)*(b+three)))
goto 2
16 continue
ww = done
goto 2
65 continue
ww = done
goto 2
70 continue
ww = -sqrt((b-two)*(b+done)/(b*b-b))
goto 2
75 continue
ww = -done
goto 2
80 continue
ww = done
2 continue

return

end subroutine segmidc22

subroutine segmidb3(w,ww,ind1,jbr)
! submid b3(b&l)

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/vtwo/0.5d0/

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(21,25,30,31,53,57,62,63),ind1
21 continue
w = -sqrt(vtwo)
ww = -sqrt(b/(b+b+4.0d0))
goto 2
57 continue
w = fq*sqrt(vtwo)
ww = -fq*sqrt((b+two)/(b+b))
goto 2
62 continue
w = -fq*sqrt((b+done)/(b+b))
ww = fq*sqrt((b-done)/(b+b))
goto 2
31 continue
w = sqrt((b+done)/(b+b+4.0d0))
ww = sqrt((b+three)/(b+b+4.0d0))
goto 2
25 continue
ww = -fq*sqrt((b+three)/(b+two))
goto 2
30 continue
ww = fq
goto 2
53 continue
ww = -sqrt((b-done)/b)
goto 2
63 continue
ww = done
2 continue

return

end subroutine segmidb3

subroutine segmidb4(w,ww,ind1,jbr)
! segmid b4(b&r)

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/vtwo/0.5d0/

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(18,19,24,28,50,51,56,60),ind1
50 continue
w = -sqrt((b+done)/(b+b))
ww = -sqrt((b-done)/(b+b))
goto 6
19 continue
w = fq*sqrt((b+done)/(b+b+4.0d0))
ww = -fq*sqrt((b+three)/(b+b+4.0d0))
goto 6
24 continue
w = -fq*sqrt(vtwo)
ww = fq*sqrt(b/(b+b+4.0d0))
goto 6

60 continue
w = sqrt(vtwo)
ww = sqrt((b+two)/(b+b))
goto 6
18 continue
ww = -done
goto 6
28 continue
ww = sqrt((b+three)/(b+two))
goto 6
51 continue
ww = -fq
goto 6
56 continue
ww = fq*sqrt((b-done)/b)
6 continue

return

end subroutine segmidb4

subroutine segmidd10(w,ww,ind1,jbr)
! segmid d10(d^r&l)

implicit real*8(a-h,o-z)
data dzero,done/0.d0,1.d0/

! calculate w,ww
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(29,61),ind1
61 continue
w = fq
goto 100
29 continue
w = fq
100 continue

return

end subroutine segmidd10

subroutine segmidb2(w,ww,ind1,jbr)
! segmid b2(b^r)

implicit real*8(a-h,o-z)
data dzero,done,two/0.d0,1.d0,2.d0/vtwo/0.5d0/

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(5,15,37,41,46,47,73,78),ind1
37 continue
w = -sqrt(b/(b+b+two))
ww = sqrt((b+two)/(b+b+two))
goto 2
41 continue
w = fq*sqrt((b+two)/(b+b+two))
ww = fq*sqrt(b/(b+b+two))
goto 2
46 continue
w = fq*sqrt(vtwo)
ww = fq*sqrt((b+two)/(b+b))
goto 2
47 continue
w = sqrt(vtwo)
ww = -sqrt(b/(b+b+4.0d0))
goto 2
5 continue
ww = done
goto 2
73 continue
ww = fq
goto 2
78 continue
ww = fq*sqrt((b+done)/b)
goto 2
15 continue
ww = -sqrt((b+done)/(b+two))
2 continue

return

end subroutine segmidb2

subroutine segmidb1(w,ww,ind1,jbr)
! segmid b1(b^l)

implicit real*8(a-h,o-z)
data dzero,done,two/0.d0,1.d0,2.d0/vtwo/0.5d0/

! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
w = dzero
ww = dzero
goto(3,8,34,35,40,44,66,76),ind1
34 continue
w = -sqrt(vtwo)
ww = sqrt((b+two)/(b+b))
goto 4
35 continue
w = -fq*sqrt(vtwo)
ww = -fq*sqrt(b/(b+b+4.0d0))
goto 4
40 continue
w = -fq*sqrt(b/(b+b+two))
ww = -fq*sqrt((b+two)/(b+b+two))
goto 4
44 continue
w = sqrt((b+two)/(b+b+two))
ww = -sqrt(b/(b+b+two))
goto 4
66 continue
ww = sqrt((b+done)/b)
goto 4
3 continue
ww = -fq*sqrt((b+done)/(b+two))
goto 4
8 continue
ww = -fq
goto 4
76 continue
ww = -done
4 continue

return

end subroutine segmidb1

subroutine stermh(isq,w,ww,ind1,jbr)

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/vtwo/0.5d0/

isq = 0
w = dzero
ww = dzero
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(34,35,36,38,39,40,43,44,48),ind1
! case a&l
! case a&r
34 continue
w = fq
goto 1
35 continue
w = done
goto 1
40 continue
w = sqrt(b/(b+done))
goto 1
44 continue
w = -fq*sqrt((b+two)/(b+done))
!if (abs(w) > 1.e-13) then
1 continue
isq = 1
ww = w
goto 100
! case d&r&l
38 continue
w = -fq*sqrt(vtwo)
ww = -fq*sqrt((b-done)/(b+b+two))
goto 2
43 continue
w = -fq*sqrt(vtwo)
ww = fq*sqrt((b+three)/(b+b+two))
goto 2
48 continue
w = fq*sqrt(two)
goto 2
39 continue
ww = -sqrt(b/(b+done))
!if (dldr == 2101) ww=(b+two)/(b+done)
2 continue
isq = 2
goto 100
! case d&l&l
! case d&r&r
36 continue
w = -sqrt(two)
isq = 3

100 continue
return

end subroutine stermh

subroutine stmh(isq,w,ww,mw,ind1,jbr)

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.0d0,1.0d0,2.0d0,3.0d0/vtwo/0.5d0/

mw = 0
isq = 0
w = dzero
ww = dzero
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(38,43,48),ind1
38 continue
w = -fq*sqrt(vtwo)
ww = -fq*sqrt((b-done)/(b+b+two))
goto 2
43 continue
w = -fq*sqrt(vtwo)
ww = fq*sqrt((b+three)/(b+b+two))
goto 2
48 continue
w = fq*sqrt(two)
2 continue
if (abs(ww) > 1.d-13) mw = 2
if (abs(w) > 1.d-13) mw = mw+1
isq = 401
goto 100

100 continue
return

end subroutine stmh

subroutine smidc2(isq,w,ww,mw,ind1,jbr)

implicit real*8(a-h,o-z)
data dzero,done,two,three/0.0d0,1.0d0,2.0d0,3.0d0/

! calculate w,ww
b = dble(jbr)
w = dzero
ww = dzero
isq = 0
mw = 0
goto(33,38,43,48),ind1
!write(6,*) 'this is case c2'
33 continue
w = done
ww = done
goto 2
38 continue
w = -done
ww = -sqrt((b-done)*(b+two)/(b*b+b))
goto 2
43 continue
w = -done
ww = -sqrt(b*(b+three)/(b*b+three*b+two))
goto 2
48 continue
w = done
ww = done
goto 2
2 continue
isq = 302

return

end subroutine smidc2

subroutine stml(isq,w,ww,mw,ind1,jbr)

implicit real*8(a-h,o-z)
data dzero,done,two/0.0d0,1.0d0,2.0d0/vtwo/0.5d0/

mw = 0
isq = 0
w = dzero
ww = dzero
! calculate w,ww
b = dble(jbr)
if (mod(jbr,2) == 0) fq = done
if (mod(jbr,2) /= 0) fq = -done
goto(38,43,48),ind1
! case d^r^l
38 continue
w = -fq*sqrt(vtwo)
ww = fq*sqrt((b+two)/(b*two))
goto 4
43 continue
w = -fq*sqrt(vtwo)
ww = -fq*sqrt(b/(b+b+4.0d0))
goto 4
48 continue
w = -fq*sqrt(two)
4 continue
if (abs(ww) > 1.d-13) mw = 2
if (abs(w) > 1.d-13) mw = mw+1
isq = 402
goto 100

100 continue
return

end subroutine stml

subroutine neoc(kcoe,nocc,tcoe)

real*8 tcoe

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
