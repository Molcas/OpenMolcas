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
      Subroutine Pr2D(xyz2D,nT,nRys,la,lb,lc,ld,IfGrad,iPrint)
      Implicit Real*8 (a-h,o-z)
      Real*8 xyz2d(nT,nRys,0:la+1,0:lb+1,0:lc+1,0:ld+1,3)
      Logical IfGrad(3,4)
      Character Label*80, ch(3)*3
      Data ch/',x)',',y)',',z)'/
*
      Write (6,*)
      Write (6,*) ' Printing the 2d-integrals'
      Write (6,*)
*
      Label = ' '
      ja = 0
      If (IfGrad(1,1).or.IfGrad(2,1).or.IfGrad(3,1)) ja = 1
      jb = 0
      If (IfGrad(1,2).or.IfGrad(2,2).or.IfGrad(3,1)) jb = 1
      jc = 0
      If (IfGrad(1,3).or.IfGrad(2,3).or.IfGrad(3,3)) jc = 1
      jd = 0
      If (IfGrad(1,4).or.IfGrad(2,4).or.IfGrad(3,4)) jd = 1
      Do 10 ia = 0, la+ja
         kb = jb
         If (ia.gt.la) jb = 0
         Do 20 ib = 0, lb+jb
            kc = jc
            If (ia.gt.la.or.ib.gt.lb) jc = 0
            Do 30 ic = 0, lc+jc
               kd = jd
               If (ia.gt.la.or.ib.gt.lb.or.ic.gt.lc) kd = jd
               Do 40 id = 0, ld+jd
                  Do 50 iCar = 1, 3
                     If (ja.eq.1.and.ia.eq.la+ja.and.
     &                   .Not.IfGrad(iCar,1)) Go To 51
                     If (jb.eq.1.and.ib.eq.lb+jb.and.
     &                   .Not.IfGrad(iCar,2)) Go To 51
                     If (jc.eq.1.and.ic.eq.lc+jc.and.
     &                   .Not.IfGrad(iCar,3)) Go To 51
                     If (jd.eq.1.and.id.eq.ld+jd.and.
     &                   .Not.IfGrad(iCar,4)) Go To 51
                     Write (Label,'(A,4(I1,A))')
     &                ' xyz2D0(',ia,',',ib,',',ic,',',id,ch(iCar)
                     If (iPrint.ge.99) Then
                        Call RecPrt(Label,' ',
     &                  xyz2d(1,1,ia,ib,ic,id,iCar),nT,nRys)
                     Else
                        Write (6,'(A)') Label
                        Write (6,*) DDot_(nT*nRys,
     &                   xyz2d(1,1,ia,ib,ic,id,iCar),1,
     &                   xyz2d(1,1,ia,ib,ic,id,iCar),1)
                     End If
 51                  Continue
 50               Continue
 40            Continue
 30          Continue
 20       Continue
 10    Continue
*
       Return
       End
