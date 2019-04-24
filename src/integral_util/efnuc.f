************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1995, Roland Lindh                                     *
************************************************************************
      SubRoutine EFNuc(CoOP,Chrg,Coor,nAtm,ESIT,nOrdOp)
************************************************************************
*                                                                      *
* Object: to compute the electricstatic interaction tensor contribution*
*         from the nuclei. In the case that the test charge coincide   *
*         with a nucleau we will remove that center.                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, April '95.                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8 Chrg(nAtm), Coor(3,nAtm), ESIT((nOrdOp+1)*(nOrdOp+2)/2),
     &       CoOp(3)
*
*---- Statement function
*
      nElem(n)=(n+1)*(n+2)/2
*
      iRout = 185
      iPrint = nPrint(iRout)
      Call qEnter('EFNuc')
*
*     Compute the nuclear contribution to the electrostatic interation
*     tensor, ESIT.
*
      nComp=nElem(nOrdOp)
      call dcopy_(nComp,[Zero],0,ESIT,1)
*
      nTot=(nOrdOp+1)**6
      Call GetMem('ESIT','Allo','Inte',ipC,nTot)
      Call InitIA(iWork(ipC),nOrdOp)
*
      iPowR=2*nOrdOp+1
      Fact=One
      If (nOrdOp.ge.1) Fact=-One
      Do iAtom = 1, nAtm
         x = CoOp(1) - Coor(1,iAtom)
         y = CoOp(2) - Coor(2,iAtom)
         z = CoOp(3) - Coor(3,iAtom)
         r2 = x**2 + y**2 + z**2
         Thr=1.0D-12
         If (r2.gt.Thr) Then
            r  = Chrg(iAtom)/Sqrt(r2)**iPowR
            Do ix = nOrdOp, 0, -1
               Do iy = nOrdOp-ix, 0, -1
                  iz = nOrdOp - ix - iy
                  If (ix.eq.0) Then
                     EIx=One
                  Else
                     EIx=x**ix
                  End If
                  If (iy.eq.0) Then
                     EIy=One
                  Else
                     EIy=y**iy
                  End If
                  If (iz.eq.0) Then
                     EIz=One
                  Else
                     EIz=z**iz
                  End If
                  temp=Fact*EIx*EIy*EIz*r
*
                  Call ContEI(iWork(ipC),nOrdOp,ESIT,ix,iy,iz,temp)
*
               End Do
            End Do       ! End loop over cartesian combinations
         End If
      End Do             ! End loop over atoms
*
      Call GetMem('ESIT','Free','Inte',ipC,nTot)
*
      If (iPrint.ge.99) Call RecPrt(' The Electrostatic Interaction'
     &                 //' Tensor',' ',ESIT,nElem(nOrdOp),1)
*     Call GetMem(' Exit ESIT','CHECK','REAL',iDum,iDum)
      Call qExit('EFNuc')
      Return
      End
      Subroutine InitIA(I,mDeg)
      implicit integer (a-z)
      dimension I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg)
c
c Purpose: Express the interaction tensor, defined by the
c quantities T(a,b,c) as functions of the vector R=(x,y,z),
c where a,b, and c are nonnegative integers and
c T(a,b,c)=((d/dx)**a)((d/dy)**b)((d/dz)**c) 1/R, in terms
c of a polynomial:
c T(a,b,c)=
c  (sum over p,q,r of I(a,b,c,p,q,r) x**p y**q z**r)/(R**(2*n+1)),
c where n=a+b+c.
c The polynomial coefficients are integers, and are 0 unless
c p+q+r=n.
c Author: PAM
*
*----- Statement function
*
c      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
c initialize:
      do 10 a=0,mDeg
      do 10 b=0,mDeg
      do 10 c=0,mDeg
      do 10 p=0,mDeg
      do 10 q=0,mDeg
      do 10 r=0,mDeg
       I(a,b,c,p,q,r)=0
  10  continue
      I(0,0,0,0,0,0)=1
      If (mDeg.gt.0) Then
         I(1,0,0,1,0,0)=-1
         I(0,1,0,0,1,0)=-1
         I(0,0,1,0,0,1)=-1
      End If
      do 100 n=2,mDeg
      do 100 a=0,n
      do 100 b=0,n-a
      c=n-a-b
      do 100 p=0,n
      do 100 q=0,n-p
      r=n-p-q
      new=0
      if(a.gt.0) then
        if(p.gt.0) new=(p-(2*n))*I(a-1,b,c,p-1,q,r)
        if(q.gt.1) new=new+(p+1)*I(a-1,b,c,p+1,q-2,r)
        if(r.gt.1) new=new+(p+1)*I(a-1,b,c,p+1,q,r-2)
      else if(b.gt.0) then
        if(q.gt.0) new=(q-(2*n))*I(a,b-1,c,p,q-1,r)
        if(r.gt.1) new=new+(q+1)*I(a,b-1,c,p,q+1,r-2)
        if(p.gt.1) new=new+(q+1)*I(a,b-1,c,p-2,q+1,r)
      else
        if(r.gt.0) new=(r-(2*n))*I(a,b,c-1,p,q,r-1)
        if(p.gt.1) new=new+(r+1)*I(a,b,c-1,p-2,q,r+1)
        if(q.gt.1) new=new+(r+1)*I(a,b,c-1,p,q-2,r+1)
      end if
      I(a,b,c,p,q,r)=new
 100  continue
c
c write out only elements with a>=b>=c. The others are obtained
c by index permutation.
c This restriction has been removed! (Roland Lindh)
*     n=mDeg
*     do 200 a=n,0,-1
*     do 200 b=n-a,0,-1
*     c=n-a-b
*     write(*,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,
*    &     Ind(n,a,c)
*     do 150 p=n,0,-1
*     do 150 q=n-p,0,-1
*     r=n-p-q
*     coef=I(a,b,c,p,q,r)
*     if(coef.eq.0) goto 150
*     write(*,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)')
*    &  coef,p,q,r,Ind(n,p,r)
*150  continue
*200  continue
*
      Return
      End
      Subroutine ContEI(I,mDeg,ESIT,ix,iy,iz,temp)
      implicit integer (a-z)
      dimension I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg)
      Real*8 ESIT((mDeg+1)*(mDeg+2)/2), Temp
c
c Purpose: Express the interaction tensor, defined by the
c quantities T(a,b,c) as functions of the vector R=(x,y,z),
c where a,b, and c are nonnegative integers and
c T(a,b,c)=((d/dx)**a)((d/dy)**b)((d/dz)**c) 1/R, in terms
c of a polynomial:
c T(a,b,c)=
c  (sum over p,q,r of I(a,b,c,p,q,r) x**p y**q z**r)/(R**(2*n+1)),
c where n=a+b+c.
c The polynomial coefficients are integers, and are 0 unless
c p+q+r=n.
c
*
*----- Statement function
*
c      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
c
c write out only elements with a>=b>=c. The others are obtained
c by index permutation.
c This restriction has been removed! (Roland Lindh)
*     n=mDeg
*     do 200 a=n,0,-1
*     do 200 b=n-a,0,-1
*     c=n-a-b
*     write(*,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,
*    &     Ind(n,a,c)
*     do 150 p=n,0,-1
*     do 150 q=n-p,0,-1
*     r=n-p-q
*     coef=I(a,b,c,p,q,r)
*     if(coef.eq.0) goto 150
*     write(*,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)')
*    &  coef,p,q,r,Ind(n,p,r)
*150  continue
*200  continue
*
*     Write (*,*) ' Temp,ix,iy,iz=',temp,ix,iy,iz
      n=mDeg
      ip = 0
      Do a=n,0,-1
         Do b=n-a,0,-1
            c=n-a-b
            ip = ip + 1
*           Write (*,*) ip, I(a,b,c,ix,iy,iz)
            If (I(a,b,c,ix,iy,iz).ne.0)
     &      ESIT(ip)=ESIT(ip)+DBLE(I(a,b,c,ix,iy,iz))*temp
         End Do
      End Do
*
      Return
      End
