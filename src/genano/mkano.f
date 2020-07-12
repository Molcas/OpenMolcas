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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine MkAno
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
#include "symlab.fh"
      Dimension s(MxS*(MxS+1)/2),u(MxS*MxS),c(MxS*MxS),d(MxS*(MxS+1))
      Dimension q(MxS*(MxS+1)/2)
      Integer nContr(0:MxLqn)
      Dimension OccNo(MxBasX)
      call molcas_open(17,'ANO')
c      Open(Unit=17,File='ANO',Form='FORMATTED')
      iPtr=1
      Write(6,*)
      Call CollapseOutput(1,'   Contraction coefficients')
      Write(6,'(3X,A)')     '   ------------------------'
*
* Loop over angular quantum number
*
      nOccNo=0
      ChkSum=0.0
      Do 1000 iLqn=0,MxLqn
         n=nPrim(iLqn)
         nContr(iLqn)=0
         If(n.le.0) Go To 1000
         Write(6,*)
         Write(6,*) '*** Contraction coefficients for the ',
     &              type(iLqn*(iLqn+1)+1)(3:3),' shell ***'
         Write(6,*)
         nTri=n*(n+1)/2
         Do 100 i=1,n
         Do 101 j=1,n
           ind=j+(i-1)*n
           u(ind)=0.0d0
101      Continue
100      Continue
         Do 200 i=1,n
           ind=i+(i-1)*n
           u(ind)=1.0d0
200      Continue
         Do 300 i=1,nTri
            s(i)=Ssym(iPtr-1+i)
300      Continue
*        Call TriPrt('s matrix','(6F12.6)',s,n)
*        Write(*,*) 'u matrix'
*        Call sqprt(u,n)
         Call Jacob(s,u,n,n)
*        Call TriPrt('s matrix','(6F12.6)',s,n)
*        Write(*,*) 'u matrix'
*        Call sqprt(u,n)
*        Write(*,*)
         Do 400 i=1,n
            s(i)=sqrt(s(i*(i+1)/2))
400      Continue
         Do 500 i=1,n
         Do 501 j=1,n
            t=0.0d0
            p=0.0d0
            Do 510 k=1,n
               t=t+u(i+(k-1)*n)/s(k)*u(j+(k-1)*n)
               p=p+u(i+(k-1)*n)*s(k)*u(j+(k-1)*n)
510         Continue
            c(j+(i-1)*n)=t
            q(j+(i-1)*n)=p
501      Continue
500      Continue
         Do 600 i=1,n
         Do 601 j=1,i
            t=0.0d0
            Do 610 k=1,n
            Do 611 l=1,n
               ind=Min(k,l)+Max(k,l)*(Max(k,l)-1)/2
               p=q(i+(k-1)*n)*tDsym(ind+iPtr-1)*q(j+(l-1)*n)
*              Write(*,'(a,f12.6)') '   p:',p
               t=t+p
611         Continue
610         Continue
            d(j+i*(i-1)/2)=t
601      Continue
600      Continue
*        Call TriPrt('d matrix','(6F12.6)',d,n)
*        Write(*,*) 'c matrix'
*        Call sqprt(c,n)
         Call Jacob(d,c,n,n)
*        Call TriPrt('d matrix','(6F12.6)',d,n)
*        Write(*,*) 'c matrix'
*        Call sqprt(c,n)
*        Write(*,*)
         nBig=0
         Do 700 i=1,n
            d(i)=d(i*(i+1)/2)
            If(d(i).gt.thr) nBig=nBig+1
700      Continue
*        nBig=n
         nContr(iLqn)=nBig
         Call Sort_genano(d,c,n,n)
         Do 750 i=1,nContr(iLqn)
            OccNo(i+nOccNo)=d(i)
750      Continue
         nOccNo=nOccNo+nContr(iLqn)
*        Write(*,*) 'c matrix'
*        Call sqprt(c,n)
         Call NOphase(c,n)
*        Write(*,*) 'c matrix'
*        Call sqprt(c,n)
*        Write(*,*)
         Write(6,'(a,10(1x,f9.4))') 'occ    ',(d(i),i=1,nBig)
         Write(6,'(a,10(1x,f9.4))') 'lg(occ)',
     &      (Log10(d(i)+1.0d-12),i=1,nBig)
         Write(6,*)

         Do 800 j=1,n
            Write(6,'(7x,10(1x,f9.4))') (c(j+(i-1)*n),i=1,nBig)
800      Continue
         Do 850 i=1,nBig
            Do 851 j=1,n
               Do 852 k=1,n
                  ChkSum=ChkSum+d(i)*c(j+(i-1)*n)*c(k+(i-1)*n)
852            Continue
851         Continue
850      Continue
*        Write(17,*) 'Coefficients for the ',type(iLqn*(iLqn+1)+1),
*    &               ' shell.'
*        Write(17,*)
         Write(17,'(2i5)') n,nBig
         If(rowise) Then
            If(n*nBig.gt.0) Then
               Do 900 i=1,nBig
                  Write(17,'(10(1x,f15.10))') (c(j+(i-1)*n),j=1,n)
900            Continue
            End If
         Else
            If(n*nBig.gt.0) Then
               Do 950 j=1,n
                  Write(17,'(10(1x,f15.10))') (c(j+(i-1)*n),i=1,nBig)
950            Continue
            End If
         EndIf
*        Write(17,*)
*        Write(17,*)
         iPtr=iPtr+(2*iLqn+1)*nTri
1000  Continue
      Call CollapseOutput(0,'   Contraction coefficients')
      Write(6,*)
c      Write(6,'(a,f12.6)') 'Check sum is',ChkSum
      Call Add_Info('GENANO_CHKSUM',[ChkSum],1,5)
      Close(Unit=17)
      call molcas_open(18,'FIG')
c      Open(Unit=18,File='FIG',Form='FORMATTED')
      Call FigOpn(18)
      Call FigPrt(18,Title,MxLqn,nContr,OccNo)
      Call FigCls(18)
      Close(Unit=18)
      End
