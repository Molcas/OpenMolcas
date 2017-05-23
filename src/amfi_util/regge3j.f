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
      real*8 function  regge3j(
     *j1,     ! integer  2*j1
     *j2,     ! integer  2*j2
     *j3,     ! integer  2*j3
     *m1,     ! integer  2*m1
     *m2,     ! integer  2*m2
     *m3)     ! integer  2*m3
cbs   uses magic square of regge (see Lindner pp. 38-39)
cbs
cbs    ---                                            ---
cbs   |                                                  |
cbs   | -j1+j2+j3     j1-j2+j3         j1+j2-j3          |
cbs   |                                                  |
cbs   |                                                  |
cbs   |  j1-m1        j2-m2            j3-m3             |
cbs   |                                                  |
cbs   |                                                  |
cbs   |  j1+m1        j2+m2            j3+m3             |
cbs   |                                                  |
cbs    ---                                            ---
cbs
      implicit real*8(a-h,o-z)
      dimension MAT(3,3)
CBS   logical testup,testdown
#include "Regge.fh"
cbs  facul,   integer array (nprim,0:mxLinRE) prime-expansion of factorials
cbs  mxLinRE,    integer max. number for facul is given
cbs  nprim,   number of primes for expansion of factorials
cbs  prim,    integer array with the first nprim prime numbers
cbs  iwork)   integer array of size nprim
      regge3j=0d0
c     write(6,'(A24,6I3)') '3J to be calculated for ',
c    *j1,j2,j3,m1,m2,m3
cbs   quick check  if =/= 0 at all
      icheck=m1+m2+m3
      if (icheck.ne.0) then
c     write(6,*) 'sum over m =/= 0'
      return
      endif
cbs   check triangular relation (|j1-j2|<= j3 <= j1+j2 )
      imini=iabs(j1-j2)
      imaxi=j1+j2
      if (j3.lt.imini.or.j3.gt.imaxi) then
c     write(6,*) 'triangular relation not fulfilled'
      return
      endif
cbs   quick check  if =/= 0 at all  end
cbs
cbs   3J-symbol is not zero by simple rules
cbs
cbs   initialize MAT
      MAT(1,1) =-j1+j2+j3
      MAT(2,1) =j1-m1
      MAT(3,1) =j1+m1
      MAT(1,2) =j1-j2+j3
      MAT(2,2) =j2-m2
      MAT(3,2) =j2+m2
      MAT(1,3) =j1+j2-j3
      MAT(2,3) =j3-m3
      MAT(3,3) =j3+m3
      do I=1,3
      do J=1,3
cbs   check for even numbers (2*integer) and positive or zero
      if (mod(MAT(J,I),2).ne.0.or.MAT(J,I).lt.0)  then
c     write(6,*) 'J,I,MAT(J,I): ',J,I,MAT(J,I)
      return
      endif
      MAT(J,I)=MAT(J,I)/2
      if (Mat(j,i).gt.mxLinRE)
     *Call SysAbendMsg('regge3j','increase mxLinRE for regge3j',' ')
      enddo
      enddo
      Isigma=(j1+j2+j3)/2
cbs   check the magic sums
      do I=1,3
      IROW=0
      ICOL=0
      do J=1,3
      IROW=IROW+MAT(I,J)
      ICOL=ICOL+MAT(J,I)
      enddo
      if (IROW.ne.Isigma.or.ICOL.ne.Isigma) then
c     write(6,*) 'I,IROW,ICOL ',I,IROW,ICOL
      return
      endif
      enddo
cbs   if j1+j2+j3 is odd: check for equal rows or columns
      Isign=1
      if (iabs(mod(Isigma,2)).eq.1) then
      isign=-1
         do I=1,3
         do J=I+1,3
            if (MAT(1,I).eq.MAT(1,J).and.
     *         MAT(2,I).eq.MAT(2,J).and.
     *         MAT(3,I).eq.MAT(3,J)) return
            if (MAT(I,1).eq.MAT(J,1).and.
     *         MAT(I,2).eq.MAT(J,2).and.
     *         MAT(I,3).eq.MAT(J,3)) return
         enddo
         enddo
      endif
cbs   look for the lowest element indices: IFIRST,ISECOND
      imini=MAT(1,1)
      IFIRST=1
      ISECOND=1
      do I=1,3
      do J=1,3
      if (MAT(J,I).lt.imini) then
      IFIRST=J
      ISECOND=I
      imini=MAT(J,I)
      endif
      enddo
      enddo
c     write(6,*) 'Matrix before commuting vectors'
      do ibm=1,3
c     write(6,'(3I5)') (Mat(ibm,j),j=1,3)
      enddo
      if (IFIRST.ne.1) then  !interchange rows
c     write(6,*) 'IFIRST = ',ifirst
      do I=1,3
      IDUMMY=MAT(1,I)
      MAT(1,I)=MAT(IFIRST,I)
      MAT(IFIRST,I)=IDUMMY
      enddo
      endif
      if (ISECOND.ne.1) then  !interchange columns
c     write(6,*) 'ISECOND = ',isecond
      do I=1,3
      IDUMMY=MAT(I,1)
      MAT(I,1)=MAT(I,ISECOND)
      MAT(I,ISECOND)=IDUMMY
      enddo
      endif
cbs   lowest element is now on (1,1)
c     write(6,*) 'Matrix after commuting vectors'
c     do ibm=1,3
c     write(6,'(3I5)') (Mat(ibm,j),j=1,3)
c     enddo
cbs   begin to calculate Sum over s_n
cbs   first the simple cases
      if (Mat(1,1).eq.0) then
      isum=1
      elseif (Mat(1,1).eq.1) then
      isum=Mat(2,3)*Mat(3,2)-Mat(2,2)*Mat(3,3)
      elseif (Mat(1,1).eq.2) then
      isum=Mat(2,3)*(Mat(2,3)-1)*Mat(3,2)*(Mat(3,2)-1)-
     *2*Mat(2,3)*Mat(3,2)*Mat(2,2)*Mat(3,3)+
     *Mat(2,2)*(Mat(2,2)-1)*Mat(3,3)*(Mat(3,3)-1)
      else !  all the cases with Mat(1,1) >= 3
              Icoeff=1
              do Ibm=Mat(3,2)-Mat(1,1)+1,Mat(3,2)
                icoeff=icoeff*ibm
              enddo
              do Ibm=Mat(2,3)-Mat(1,1)+1,Mat(2,3)
                icoeff=icoeff*ibm
              enddo
              isum=icoeff
              do Icount=1,MAT(1,1)
                 icoeff=-icoeff*(Mat(1,1)+1-icount)*(Mat(2,2)+1-icount)*
     *           (Mat(3,3)+1-icount)
                 Idenom=icount*(Mat(2,3)-Mat(1,1)+icount)*
     *           (Mat(3,2)-Mat(1,1)+icount)
                 icoeff=icoeff/Idenom
                 isum=isum+icoeff
              enddo
      endif
cbs  additional sign from interchanging rows or columns
      if (ifirst.ne.1) isum=isum*isign
      if (isecond.ne.1) isum=isum*isign
c     write(6,*) 'isum = ',isum
cbs       Mat(2,3)+Mat(3,2)
cbs    (-)
      if (iabs(mod((Mat(2,3)+Mat(3,2)),2)).eq.1) isum=-isum
cbs   final factor
      LIMIT=ihigh(max(Mat(1,1),Mat(1,2),Mat(1,3),
     *Mat(2,1),Mat(2,2),Mat(2,3),Mat(3,1),Mat(3,2),
     *Mat(3,3),(Isigma+1)))
      do I=1,LIMIT
      iwork(I)=facul(I,Mat(1,2))+facul(I,Mat(2,1))+
     *facul(I,Mat(3,1))+facul(I,Mat(1,3))-
     *facul(I,Mat(1,1))-facul(I,Mat(2,2))-
     *facul(I,Mat(3,3))-facul(I,(Isigma+1))-
     *facul(I,Mat(2,3))-facul(I,Mat(3,2))
      enddo
c     write(6,*) 'Iwork: ',(iwork(i),i=1,LIMIT)
      factor=1d0
cbs   iup=1
CBS   idown=1
CBS   testup=.true.
CBS   testdown=.true.
CBS   do I=1,LIMIT
CBS   do J=1,iwork(I)
CBS   iup=iup*prim(i)
CBS   if (iup.lt.0) testup=.false. !check for Integer overflow
CBS   enddo
CBS   Enddo
CBS   up=DBLE(iup)
CBS   if(.not.testup) then ! if the integers did not run correctly
              up=1d0
              do I=1,LIMIT
              do J=1,iwork(I)
              up=up*DBLE(prim(i))
              enddo
              enddo
CBS   endif
CBS   do I=1,LIMIT
CBS   do J=1,-iwork(I)
CBS   idown=idown*prim(i)
CBS   if (idown.lt.0) testdown=.false.
CBS   enddo
CBS   enddo
CBS   down=DBLE(idown)
CBS   if(.not.testdown) then
              down=1d0
              do I=1,LIMIT
              do J=1,-iwork(I)
              down=down*DBLE(prim(i))
              enddo
              enddo
CBS   endif
c     if (.not.(testup.and.testdown)) then
c     write(6,*) 'j1,j2,j3,m1,m2,m3 ',j1,j2,j3,m1,m2,m3
c     write(6,*) 'iup,idown ',iup,idown,'up,down ',up,down
c     endif
      factor=factor*up/down
cbs   final result
      regge3j=sqrt(factor)*DBLE(isum)
      return
      end
