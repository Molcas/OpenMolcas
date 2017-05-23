************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
c
      subroutine stermha4(w,ww,ind1,jbr)
      implicit real*8 (a-h,o-z)
      data done,two/1.d0,2.d0/
      isq=0
      w=0.0d0
      ww=0.0d0
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(34,35,40,44),ind1
c     case a&l
c     case a&r
34    w=fq
      goto 1
35    w=done
      goto 1
40    w=sqrt(b/(b+done))
      goto 1
44    w=-fq*sqrt((b+two)/(b+done))
c      if(abs(w).gt.1.e-13) then
1     isq=1
      ww=w
      return
      end
c
      subroutine stermhd1(w,ww,ind1,jbr)
      implicit real*8 (a-h,o-z)
      data done,two,three/1.d0,2.d0,3.d0/
      isq=0
      w=0.0d0
      ww=0.0d0
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(38,39,43,48),ind1
c   d1: case d&r&l
38    w=-fq/sqrt(two)
      ww=-fq*sqrt((b-done)/(b+b+two))
      goto 2
43    w=-fq/sqrt(two)
      ww=fq*sqrt((b+three)/(b+b+two))
      goto 2
48    w=fq*sqrt(two)
      goto 2
39    ww=-sqrt(b/(b+done))
c      if(dldr.eq.2101) ww=(b+two)/(b+done)
2     isq=2
      return
      end
c
      subroutine stermhd5(w,ww)
      implicit real*8 (a-h,o-z)
      data two/2.d0/
      isq=0
      w=0.0d0
      ww=0.0d0
c     calculate w,ww
c   d5: case d&r&r
      w=-sqrt(two)
      isq=3
      return
      end
c
      subroutine stermla1(w,ww,ind1,jbr)
c     case a^l
      implicit real*8 (a-h,o-z)
      data done/1.d0/
      w=0.0d0
      ww=0.0d0
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(19,24,50,60),ind1
19    w=-sqrt((b+1)/(b+2))
      goto 2
24    w=-done
      goto 2
50    w=fq*sqrt((b+1)/b)
      goto 2
60    w=fq
2     ww=w
      return
      end
c
      subroutine stermla2(w,ww,ind1,jbr)
c     case a^r
      implicit real*8 (a-h,o-z)
      data done,two/1.d0,2.d0/
      w=0.0d0
      ww=0.0d0
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(21,31,57,62),ind1
21    w=-fq
      goto 3
57    w=done
      goto 3
62    w=sqrt((b+done)/b)
      goto 3
31    w=-fq*sqrt((b+done)/(b+two))
3     ww=w
      return
      end
c
      subroutine stermld2(w,ww,ind1,jbr)
c    d2: case d^r^l
      implicit real*8 (a-h,o-z)
      data done,two/1.d0,2.d0/
      w=0.0d0
      ww=0.0d0
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(7,38,43,48,74),ind1
38    w=-fq/sqrt(two)
      ww=fq*sqrt((b+two)/(b*two))
      goto 4
43    w=-fq/sqrt(two)
      ww=-fq*sqrt(b/(b+b+4.0d0))
      goto 4
48    w=-fq*sqrt(two)
      goto 4
7     ww=-sqrt((b+done)/(b+two))
      goto 4
74    ww=sqrt((b+done)/b)
4     return
      end
c
      subroutine stermld6(w,ww)
c   d6: case d^r^r
      implicit real*8  (a-h,o-z)
      data dzero,two/0.d0,2.d0/
      w=-sqrt(two)
      ww=dzero
      return
      end
c
      subroutine segmidc1(w,ww,ind1,jbr)
c  case c1
      implicit real*8 (a-h,o-z)
      data dzero,done,three/0.d0,1.d0,3.d0/
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(17,22,23,27,32,49,54,58,59,64),ind1
17    w=done
      goto 2
22    w=done
      goto 2
23    w=fq/sqrt((b*b+4.0d0*b+4.0d0))
      goto 2
27    w=-sqrt((b+done)*(b+three)/(b*b+4.0d0*b+4.0d0))
      goto 2
32    w=-done
      goto 2
49    w=done
      goto 2
54    w=sqrt((b+done)*(b-done)/(b*b))
      goto 2
58    w=fq/b
      goto 2
59    w=-done
      goto 2
64    w=-done
2     ww=w
      return
      end
c
      subroutine segmidc2(w,ww,ind1,jbr)
c  case c2
      implicit real*8 (a-h,o-z)
      data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(1,6,7,11,16,33,38,39,42,43,48,65,70,74,75,80),ind1
33    w=done
      ww=done
      goto 2
38    w=-done
      ww=-sqrt((b-done)*(b+two)/(b*b+b))
      goto 2
43    w=-done
      ww=-sqrt(b*(b+three)/(b*b+three*b+two))
      goto 2
48    w=done
      ww=done
      goto 2
39    ww=-fq*sqrt(two/(b*b+three*b+two))
      goto 2
7     ww=-fq*sqrt(two/((b+two)*(b+three)))
      goto 2
1     ww=done
      goto 2
6     ww=-done
      goto 2
11    ww=-sqrt((b+done)*(b+4.0d0)/((b+two)*(b+three)))
      goto 2
16    ww=done
      goto 2
65    ww=done
      goto 2
70    ww=-sqrt((b-two)*(b+done)/(b*b-b))
      goto 2
75    ww=-done
      goto 2
80    ww=done
      goto 2
42    ww=fq*sqrt(two/(b*(b+done)))
      goto 2
74    ww=fq*sqrt(two/(b*b-b))
2     return
      end
c
      subroutine segmidc22(w,ww,ind1,jbr)
c  case c22
      implicit real*8 (a-h,o-z)
      data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(1,6,7,11,16,33,38,39,43,48,65,70,75,80),ind1
33    w=done
      ww=done
      goto 2
38    w=-done
      ww=-sqrt((b-done)*(b+two)/(b*b+b))
      goto 2
43    w=-done
      ww=-sqrt(b*(b+three)/(b*b+three*b+two))
      goto 2
48    w=done
      ww=done
      goto 2
39    ww=-fq*sqrt(two/(b*b+three*b+two))
      goto 2
7     ww=-fq*sqrt(two/((b+two)*(b+three)))
      goto 2
1     ww=done
      goto 2
6     ww=-done
      goto 2
11    ww=-sqrt((b+done)*(b+4.0d0)/((b+two)*(b+three)))
      goto 2
16    ww=done
      goto 2
65    ww=done
      goto 2
70    ww=-sqrt((b-two)*(b+done)/(b*b-b))
      goto 2
75    ww=-done
      goto 2
80    ww=done
2     return
      end
c
      subroutine segmidb3(w,ww,ind1,jbr)
c  submid b3(b&l)
      implicit real*8     (a-h,o-z)
      data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/vtwo/0.5d0/
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(21,25,30,31,53,57,62,63),ind1
21    w=- sqrt(vtwo)
      ww=-sqrt(b/(b+b+4.0d0))
      goto 2
57    w=fq*sqrt(vtwo)
      ww=-fq*sqrt((b+two)/(b+b))
      goto 2
62    w=-fq*sqrt((b+done)/(b+b))
      ww=fq*sqrt((b-done)/(b+b))
      goto 2
31    w=sqrt((b+done)/(b+b+4.0d0))
      ww=sqrt((b+three)/(b+b+4.0d0))
      goto 2
25    ww=-fq*sqrt((b+three)/(b+two))
      goto 2
30    ww=fq
      goto 2
53    ww=-sqrt((b-done)/b)
      goto 2
63    ww=done
2     return
      end
c
      subroutine segmidb4(w,ww,ind1,jbr)
c   segmid b4(b&r)
      implicit real*8 (a-h,o-z)
      data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/vtwo/0.5d0/
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(18,19,24,28,50,51,56,60),ind1
50    w=-sqrt((b+done)/(b+b))
      ww=-sqrt((b-done)/(b+b))
      goto 6
19    w=fq*sqrt((b+done)/(b+b+4.0d0))
      ww=-fq*sqrt((b+three)/(b+b+4.0d0))
      goto 6
24    w=-fq*sqrt(vtwo)
      ww=fq*sqrt(b/(b+b+4.0d0))
      goto 6

60    w=sqrt(vtwo)
      ww=sqrt((b+two)/(b+b))
      goto 6
18    ww=-done
      goto 6
28    ww=sqrt((b+three)/(b+two))
      goto 6
51    ww=-fq
      goto 6
56    ww=fq*sqrt((b-done)/b)
6     return
      end
c
      subroutine segmidd10(w,ww,ind1,jbr)
c   segmid d10(d^r&l)
      implicit real*8 (a-h,o-z)
      data dzero,done/0.d0,1.d0/
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(29,61),ind1
61    w=fq
      goto 100
29    w=fq
100   return
      end

      subroutine segmidb2(w,ww,ind1,jbr)
c   segmid b2(b^r)
      implicit real*8 (a-h,o-z)
      data dzero,done,two/0.d0,1.d0,2.d0/vtwo/0.5d0/
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(5,15,37,41,46,47,73,78),ind1
37    w=-sqrt(b/(b+b+two))
      ww=sqrt((b+two)/(b+b+two))
      goto 2
41    w=fq*sqrt((b+two)/(b+b+two))
      ww=fq*sqrt(b/(b+b+two))
      goto 2
46    w=fq*sqrt(vtwo)
      ww=fq*sqrt((b+two)/(b+b))
      goto 2
47    w=sqrt(vtwo)
      ww=-sqrt(b/(b+b+4.0d0))
      goto 2
5     ww=done
      goto 2
73    ww=fq
      goto 2
78    ww=fq*sqrt((b+done)/b)
      goto 2
15    ww=-sqrt((b+done)/(b+two))
2     return
      end
c
      subroutine segmidb1(w,ww,ind1,jbr)
c   segmid b1(b^l)
      implicit real*8 (a-h,o-z)
      data dzero,done,two/0.d0,1.d0,2.d0/vtwo/0.5d0/
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      goto(3,8,34,35,40,44,66,76),ind1
34    w=-sqrt(vtwo)
      ww=sqrt((b+two)/(b+b))
      goto 4
35    w=-fq*sqrt(vtwo)
      ww=-fq*sqrt(b/(b+b+4.0d0))
      goto 4
40    w=-fq*sqrt(b/(b+b+two))
      ww=-fq*sqrt((b+two)/(b+b+two))
      goto 4
44    w=sqrt((b+two)/(b+b+two))
      ww=-sqrt(b/(b+b+two))
      goto 4
66    ww=sqrt((b+done)/b)
      goto 4
3     ww=-fq*sqrt((b+done)/(b+two))
      goto 4
8     ww=-fq
      goto 4
76    ww=-done
4     return
      end

      subroutine stermh(isq,w,ww,ind1,jbr)
      implicit real*8 (a-h,o-z)
      data dzero,done,two,three/0.d0,1.d0,2.d0,3.d0/vtwo/0.5d0/
      isq=0
      w=dzero
      ww=dzero
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(34,35,36,38,39,40,43,44,48),ind1
c     case a&l
c     case a&r
34    w=fq
      goto 1
35    w=done
      goto 1
40    w=sqrt(b/(b+done))
      goto 1
44    w=-fq*sqrt((b+two)/(b+done))
c      if(abs(w).gt.1.e-13) then
1     isq=1
      ww=w
      goto 100
c     case d&r&l
38    w=-fq*sqrt(vtwo)
      ww=-fq*sqrt((b-done)/(b+b+two))
      goto 2
43    w=-fq*sqrt(vtwo)
      ww=fq*sqrt((b+three)/(b+b+two))
      goto 2
48    w=fq*sqrt(two)
      goto 2
39    ww=-sqrt(b/(b+done))
c      if(dldr.eq.2101) ww=(b+two)/(b+done)
2     isq=2
      goto 100
c     case d&l&l
c     case d&r&r
36    w=-sqrt(two)
      isq=3
100   return
      end
c

      subroutine stmh(isq,w,ww,mw,ind1,jbr)
      implicit real*8    (a-h,o-z)
      data dzero,done,two,three/0.0d0,1.0d0,2.0d0,3.0d0/vtwo/0.5d0/
      mw=0
      isq=0
      w=dzero
      ww=dzero
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      goto(38,43,48),ind1
38    w=-fq*sqrt(vtwo)
      ww=-fq*sqrt((b-done)/(b+b+two))
      goto 2
43    w=-fq*sqrt(vtwo)
      ww=fq*sqrt((b+three)/(b+b+two))
      goto 2
48    w=fq*sqrt(two)
2     if(abs(ww).gt.1.d-13) mw=2
      if(abs(w).gt.1.d-13) mw=mw+1
      sq=401
      goto 100
100   return
      end
c
      subroutine smidc2(isq,w,ww,mw,ind1,jbr)
      implicit real*8 (a-h,o-z)
      data dzero,done,two,three/0.0d0,1.0d0,2.0d0,3.0d0/
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      w=dzero
      ww=dzero
      isq=0
      mw=0
      goto(33,38,43,48),ind1
c     write(6,*) 'this is case c2'
33    w=done
      ww=done
      goto 2
38    w=-done
      ww=-sqrt((b-done)*(b+two)/(b*b+b))
      goto 2
43    w=-done
      ww=-sqrt(b*(b+three)/(b*b+three*b+two))
      goto 2
48    w=done
      ww=done
      goto 2
2     sq=302
      return
      end

      subroutine stml(isq,w,ww,mw,ind1,jbr)
      implicit real*8 (a-h,o-z)
      data dzero,done,two/0.0d0,1.0d0,2.0d0/vtwo/0.5d0/
      mw=0
      isq=0
      w=dzero
      ww=dzero
c     calculate w,ww
      b=dble(jbr)
      if(mod(jbr,2).eq.0) fq=done
      if(mod(jbr,2).ne.0) fq=-done
      jf=1
      goto(38,43,48),ind1
c     case d^r^l
38    w=-fq*sqrt(vtwo)
      ww=fq*sqrt((b+two)/(b*two))
      goto 4
43    w=-fq*sqrt(vtwo)
      ww=-fq*sqrt(b/(b+b+4.0d0))
      goto 4
48    w=-fq*sqrt(two)
4     if(abs(ww).gt.1.d-13) mw=2
      if(abs(w).gt.1.d-13) mw=mw+1
      isq=402
      goto 100
100   return
      end

      subroutine neoc(kcoe,nocc,tcoe)
      real*8 tcoe
      nocc=1
      tcoe=kcoe
      if(kcoe.eq.0)   nocc=0
      if(kcoe.eq.100) tcoe=0.d0
      if(kcoe.eq.200) then
        nocc=2
        tcoe=-0.5d0
      endif
      return
      end
