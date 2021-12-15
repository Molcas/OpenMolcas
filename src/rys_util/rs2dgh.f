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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Rs2Dgh(xyz2D0,nT,nRys,la,lb,lc,ld,xyz2D1,xyz2D2,
     &                  IfHss,IndHss,IfGrad,IndGrd,IfG,
     &                  Coora,Alpha,Beta,Gamma,Delta,nZeta,
     &                  nEta,Scrtch,Scrtch2,Temp,Index1,Index2,
     &                  Index3,Index4,ng,nh,
     &                  ExpX,ExpY,mZeta,mEta,nIrrep,Tr)
************************************************************************
*                                                                      *
* Object:To compute the gradients and the Hessians of the 2D-integrals.*
*                                                                      *
*     Author: Anders Bernhardsson & Roland Lindh,                      *
*             Dept. of Theoretical Chemistry,                          *
*             University of Lund, SWEDEN                               *
*             Februar '95                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External ExpX, ExpY
#include "real.fh"
      Real*8 xyz2D0(nRys*nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3),
     &       xyz2D1(nRys*nT,0:la  ,0:lb  ,0:lc  ,0:ld  ,3,3),
     &       xyz2D2(nRys*nT,0:la  ,0:lb  ,0:lc  ,0:ld  ,3,6),
     &       Coora(3,4),
     &       Alpha(nZeta), Beta(nZeta), Gamma(nEta), Delta(nEta),
     &       Scrtch2(nRys*nT),Scrtch(nRys*nT), Temp(nT)
      Logical IfGrad(3,4),IfHss(4,3,4,3),IfG(4), EQ, Tr(4)
      Integer IndGrd(3,4,0:nIrrep-1), Ind1(3),
     &        Ind2(3),Ind3(3),Ind4(3),
     &        Index2(3,4,4),Index1(3,4),
     &        Index3(3,3),Index4(2,6,3),ng(3),nh(3),
     &        IndHss(4,3,4,3,0:nIrrep-1)
#ifdef NAGFOR
      Save Ind1, Ind2, Ind3, Ind4
#endif
*
      nx = 0
      ny = 0
      nz = 0
      mx=0
      my=0
      mz=0
      Call ICopy(12,[0],0,Index1,1)
      Call ICopy(48,[0],0,Index2,1)
*
*     Differentiate with respect to the first center
*
      If (IfG(1)) Then
        Call ExpX(Temp  ,mZeta,mEta,Alpha,Sqrt(Two))
        Call Exp_2(Scrtch,nRys,nT,Temp,Sqrt(Two))
        nVec = 0
        If (IfGrad(1,1)) Then
           nx = nx + 1
           nVec = nVec + 1
           Ind1(nVec) = nx
           Ind2(nVec) = 1
           Index1(1,1) =nx
           Index3(nx,1)=1
        End If
        If (IfGrad(2,1)) Then
           ny = ny + 1
           nVec = nVec + 1
           Ind1(nVec) = ny
           Ind2(nVec) = 2
           Index1(2,1) = ny
           Index3(ny,2)=1
        End If
        If (IfGrad(3,1)) Then
           nz = nz + 1
           nVec = nVec + 1
           Ind1(nVec) = nz
           Ind2(nVec) = 3
           Index1(3,1) = nz
           Index3(nz,3)=1
        End If
        Do i=nvec+1,3
         Ind1(i)=0
        End Do
*
        mvec=0
        If (IfHss(1,1,1,1)) Then
         mx=mx+1
         mvec=mvec+1
         Ind3(mVec) = mx
         Ind4(mVec) = 1
         Index2(1,1,1) = mx
         Index4(1,mx,1)=1
         Index4(2,mx,1)=1
        End If
        If (IfHss(1,2,1,2)) Then
         my=my+1
         mvec=mvec+1
         Ind3(mVec) = my
         Ind4(mVec) = 2
         Index2(2,1,1) = my
         Index4(1,my,2)=1
         Index4(2,my,2)=1
        End If
        If (IfHss(1,3,1,3)) Then
         mz=mz+1
         mvec=mvec+1
         Ind3(mVec) = mz
         Ind4(mVec) = 3
         Index2(3,1,1) = mz
         Index4(1,mz,3)=1
         Index4(2,mz,3)=1
        End If
        Do i=mvec+1,3
         Ind3(i)=0
        End Do
        nvecx=max(nvec,mvec)
        If (nVecx.ne.0) Then
*
*    Here we go with center 1
*
        Do n=1,nvec
         Do  id = 0, ld
          Do  ic = 0, lc
           Do  ib = 0, lb
            ra=-One
            Do  ia=0,la
             ra=ra+Two
             Do  iVec = 1, nRys*nT
              xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                    Scrtch(iVec) *
     &                   xyz2D0(iVec,ia+1,ib,ic,id,Ind2(n))
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do
        Do n=1,mvec
         Do  id = 0, ld
          Do  ic = 0, lc
           Do  ib = 0, lb
            ra=-One
            Do  ia=0,la
             ra=ra+Two
             Do  iVec = 1, nRys*nT
              xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n))=
     &                 Scrtch(iVec)**2*
     &               xyz2D0(iVec,ia+2,ib,ic,id,Ind4(n))-
     &                   ra*  Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do
        If (la.ge.1) Then
         Do n=1,nvec
          Do  id=0,ld
           Do  ic=0,lc
            Do  ib=0,lb
             ra=Zero
             Do  ia=1,la
              ra=ra+One
              Call DaXpY_inline(nRys*nT,-ra,
     &                   xyz2D0(1,ia-1,ib,ic,id,Ind2(n)),1,
     &                   xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
             End Do
            End Do
           End Do
          End Do
         End Do
        End If
        If (la.ge.2) Then
          Do n=1,mvec
              Do  id=0,ld
             Do  ic=0,lc
            Do  ib=0,lb
           ra=One
           Do  ia=2,la
            ra=ra+One
               Fact=ra*ra-ra
               Call Daxpy_inline(nRys*nT,Fact,
     &                    xyz2D0(1,ia-2,ib,ic,id,Ind4(n)),1,
     &                    xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)),1)
              End Do
             End Do
            End Do
           End Do
          End Do
        End If
        End If
      End If

*
*    Cross term center 1 and center 2
*
      If (IfG(2)) Then
         Call ExpX(Temp  ,mZeta,mEta,Beta,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
      End If
      If (IfG(2).and.Ifg(1)) Then
        nVec=0
        If (ifHss(2,1,1,1)) Then
           mx=mx+1
           nVec = nVec + 1
           Ind1(nvec)=mx
           Ind2(nVec)=1
           Index2(1,2,1)=mx
           Index4(1,mx,1)=2
           Index4(2,mx,1)=1
        End If
        If (ifHss(2,2,1,2)) Then
           my=my+1
           nVec = nVec + 1
           Ind1(nvec)=my
           Ind2(nVec)=2
           Index2(2,2,1)=my
           Index4(1,my,2)=2
           Index4(2,my,2)=1
        End If
        If (ifHss(2,3,1,3)) Then
           mz=mz+1
           nVec = nVec + 1
           Ind1(nvec)=mz
           Ind2(nVec)=3
           Index2(3,2,1)=mz
           Index4(1,mz,3)=2
           Index4(2,mz,3)=1
        End If
        Do i=nVec+1,3
          Ind1(i)=0
        End Do
        If (nvec.ne.0) Then
        Do n=1,nvec
         Do  id = 0, ld
          Do   ic = 0, lc
           Do  ib=0,lb
            Do  ia=0,la
             Do  iVec = 1, nRys*nT
              xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               Scrtch(iVec)  *  Scrtch2(iVec)*
     &               xyz2D0(iVec,ia+1,ib+1,ic,id,Ind2(n))
             End DO
            End DO
           End DO
          End DO
         End DO
        End DO
        If (la.ge.1) Then
         Do n=1,nvec
          Do  id = 0, ld
           Do  ic = 0, lc
            Do  ib=0,lb
             ra=Zero
             Do  ia=1,la
              ra=ra+One
              Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  ra  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia-1,ib+1,ic,id,Ind2(n))
              End Do
             End Do
            End Do
           End Do
          End Do
         End Do
        End If
        If (lb.ge.1) Then
         Do n=1,nvec
          Do   id = 0, ld
           Do  ic = 0, lc
            rb=Zero
            Do  ib=1,lb
             rb=rb+One
             Do  ia=0,la
              Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rb  *  Scrtch(iVec)*
     &                  xyz2D0(iVec,ia+1,ib-1,ic,id,Ind2(n))
              End Do
             End Do
            End Do
           End Do
          End Do
         End Do
        End If
        If ((la.ge.1).and.(lb.ge.1)) Then
         Do n=1,nvec
          Do  id = 0, ld
           Do  ic = 0, lc
            rb=Zero
            Do  ib=1,lb
             rb=rb+One
             ra=Zero
             Do  ia=1,la
              ra=ra+One
              Fact=ra*rb
              Call DaXpY_inline(nRys*nT,Fact,
     &                   xyz2D0(1,ia-1,ib-1,ic,id,Ind2(n)),1,
     &                   xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
             End Do
            End Do
           End Do
          End Do
         End Do
        End If
        End If
      End If
*
*     Differentiate with respect to the second center
*
*
      If (IfG(2) ) Then
         Call ExpX(Temp  ,mZeta,mEta,Beta,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
        nVec = 0
        If (IfGrad(1,2)) Then
         nx = nx + 1
         nVec = nVec + 1
         Ind1(nVec) = nx
         Ind2(nVec) = 1
         Index1(1,2) = nx
         Index3(nx,1)=2
        End If
        If (IfGrad(2,2)) Then
         ny = ny + 1
         nVec = nVec + 1
         Ind1(nVec) = ny
         Ind2(nVec) = 2
         Index3(ny,2)=2
         Index1(2,2) = ny
        End If
        If (IfGrad(3,2)) Then
         nz = nz + 1
         nVec = nVec + 1
         Ind1(nVec) = nz
         Ind2(nVec) = 3
         Index1(3,2) = nz
         Index3(nz,3)=2
        End If
        Do i=nvec+1,3
         Ind1(i)=0
        End Do
*
        mvec=0
        If (IfHss(2,1,2,1)) Then
          mx=mx+1
          mvec=mvec+1
          Ind3(mVec) = mx
          Ind4(mVec) = 1
          Index2(1,2,2) = mx
          Index4(1,mx,1)=2
          Index4(2,mx,1)=2
        End If
        If (IfHss(2,2,2,2)) Then
          my=my+1
          mvec=mvec+1
          Ind3(mVec) = my
          Ind4(mvec)=2
          Index2(2,2,2) = my
          Index4(1,my,2)=2
          Index4(2,my,2)=2
        End If
        If (IfHss(2,3,2,3)) Then
          mz=mz+1
          mvec=mvec+1
          Ind3(mVec) = mz
          Ind4(mVec) = 3
          Index2(3,2,2) = mz
          Index4(1,mz,3)=2
          Index4(2,mz,3)=2
        End If
        Do i=mvec+1,3
         Ind3(i)=0
        End Do
        nvecx=max(nvec,mvec)
        If (nVecx.ne.0) Then
*
        Do n=1,nVec
         Do  id = 0, ld
          Do  ic = 0, lc
           rb=-One
           Do  ib=0,lb
            rb=rb+Two
            Do  ia = 0, la
             Do  iVec = 1, nRys*nT
              xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                   Scrtch2(iVec) *
     &                 xyz2D0(iVec,ia,ib+1,ic,id,Ind2(n))
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do
        Do n=1,mVec
         Do  id = 0, ld
          Do  ic = 0, lc
           rb=-One
           Do  ib=0,lb
            rb=rb+Two
            Do  ia = 0, la
             Do  iVec = 1, nRys*nT
              xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n))=
     &                 Scrtch2(iVec)**2*
     &               xyz2D0(iVec,ia,ib+2,ic,id,Ind4(n))-
     &               rb * Scrtch2(iVec) *
     &               xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do

        If (lb.ge.1) Then
         Do n=1,nvec
             Do  id=0,ld
            Do  ic=0,lc
           rb=Zero
           Do  ib=1,lb
            rb=rb+One
          Do  ia=0,la
              Call DaXpy_inline(nRys*nT,-rb,
     &                   xyz2D0(1,ia,ib-1,ic,id,Ind2(n)),1,
     &                   xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
             End Do
            End Do
           End Do
          End Do
         End Do
        End If
        If (lb.ge.2) Then
         Do n=1,mvec
             Do  id=0,ld
            Do  ic=0,lc
           rb=One
           Do  ib=2,lb
            rb=rb+One
              Fact=rb*rb-rb
          Do  ia=0,la
              Call DaXpy_inline(nRys*nT,Fact,
     &                   xyz2D0(1,ia,ib-2,ic,id,Ind4(n)),1,
     &                   xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)),1)
             End Do
            End Do
           End Do
          End Do
         End Do
        End If
        End If
      End If
*   Cross Term center 2 and 3
*
      If (IfG(2).and.IfG(3)) Then
         Call ExpX(Temp  ,mZeta,mEta,Beta,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
         Call ExpY(Temp  ,mZeta,mEta,Gamma,Sqrt(Two))
         Call Exp_2(Scrtch,nRys,nT,Temp,Sqrt(Two))
         nVec=0
       If (ifHss(3,1,2,1)) Then
         mx=mx+1
         nVec = nVec + 1
         Ind1(nvec)=mx
         Ind2(nVec)=1
         Index2(1,3,2)=mx
         Index4(1,mx,1)=3
         Index4(2,mx,1)=2
       End If
       If (ifHss(3,2,2,2)) Then
         my=my+1
         nVec = nVec + 1
         Ind1(nvec)=my
         Ind2(nVec)=2
         Index2(2,3,2)=my
         Index4(1,my,2)=3
         Index4(2,my,2)=2
       End If
       If (ifHss(3,3,2,3)) Then
         mz=mz+1
         nVec = nVec + 1
         Ind1(nvec)=mz
         Ind2(nVec)=3
         Index2(3,3,2)=mz
         Index4(1,mz,3)=3
         Index4(2,mz,3)=2
       End If
       Do i=nVec+1,3
        Ind1(i)=0
       End Do
       If (nvec.ne.0) Then

       Do n=1,nvec
        Do  id = 0, ld
         Do  ic = 0, lc
          Do  ib = 0, lb
           Do  ia = 0, la
            Do  iVec = 1, nRys*nT
              xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               Scrtch(iVec)  *  Scrtch2(iVec)*
     &               xyz2D0(iVec,ia,ib+1,ic+1,id,Ind2(n))
            End Do
           End Do
          End Do
         End Do
        End Do
       End Do
       If (lb.ge.1) Then
        Do n=1,nvec
         Do  id = 0, ld
          Do ic = 0, lc
           rb=Zero
           Do ib=1,lb
            rb=rb+One
            Do  ia=0,la
             Do  iVec = 1, nRys*nT
              xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &        xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                rb  *  Scrtch(iVec)*
     &                xyz2D0(iVec,ia,ib-1,ic+1,id,Ind2(n))
             End DO
            End DO
           End DO
          End DO
         End DO
        End DO
       End If
       If (lc.ge.1) Then
        Do n=1,nvec
         Do   id = 0, ld
          rc=Zero
          Do  ic = 1, lc
           rc=rc+One
           Do  ib=0,lb
            Do  ia=0,la
             Do  iVec = 1, nRys*nT
               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                rc  * Scrtch2(iVec)*
     &                xyz2D0(iVec,ia,ib+1,ic-1,id,Ind2(n))
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do
       End If
       If ((lb.ge.1).and.(lc.ge.1)) Then
        Do n=1,nvec
         Do  id = 0, ld
          rc=Zero
          Do  ic = 1, lc
           rc=rc+One
           rb=Zero
           Do  ib=1,lb
            rb=rb+One
            Fact=rb*rc
            Do  ia=0,la
             Call DaxPy_inline(nRys*nT,Fact,
     &                  xyz2D0(1,ia,ib-1,ic-1,id,Ind2(n)),1,
     &                  xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
            End Do
           End Do
          End Do
         End Do
        End Do
       End If
       End If
      End If
*
*
*     Differentiate with respect to the third center
*
      If (IfG(3)) Then
         Call ExpY(Temp  ,mZeta,mEta,Gamma,Sqrt(Two))
         Call Exp_2(Scrtch,nRys,nT,Temp,Sqrt(Two))
       nvec=0
       If (IfGrad(1,3)) Then
         nx = nx + 1
         nVec = nVec + 1
         Ind1(nVec) = nx
         Ind2(nVec) = 1
         Index1(1,3) = nx
         Index3(nx,1)=3
       End If
       If (IfGrad(2,3)) Then
         ny = ny + 1
         nVec = nVec + 1
         Ind1(nVec) = ny
         Ind2(nVec) = 2
         Index1(2,3) = ny
         Index3(ny,2)=3
       End If
       If (IfGrad(3,3)) Then
         nz = nz + 1
         nVec = nVec + 1
         Ind1(nVec) = nz
         Ind2(nVec) = 3
         Index1(3,3) = nz
         Index3(nz,3)=3
       End If
       Do i=nvec+1,3
        Ind1(i)=0
       End Do
*
       mvec=0
       If (IfHss(3,1,3,1)) Then
        mx=mx+1
        mvec=mvec+1
        Ind3(mVec) = mx
        Ind4(mVec) = 1
        Index2(1,3,3) = mx
        Index4(1,mx,1)=3
        Index4(2,mx,1)=3
       End If
       If (IfHss(3,2,3,2)) Then
        my=my+1
        mvec=mvec+1
        Ind3(mVec) = my
        Ind4(mVec) = 2
        Index2(2,3,3) = my
        Index4(1,my,2)=3
        Index4(2,my,2)=3
       End If
       If (IfHss(3,3,3,3)) Then
        mz=mz+1
        mvec=mvec+1
        Ind3(mVec) = mz
        Ind4(mVec) = 3
        Index2(3,3,3) = mz
        Index4(1,mz,3)=3
        Index4(2,mz,3)=3
       End If
       Do i=mvec+1,3
        Ind3(i)=0
       End Do
       nvecx=max(nvec,mvec)
       If (nVecx.ne.0) Then
       Do n=1,nvec
        Do  id = 0, ld
         rc=-One
         Do  ic = 0, lc
          rc=rc+Two
          Do  ib=0,lb
           Do  ia = 0, la
            Do  iVec = 1, nRys*nT
             xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &            Scrtch(iVec) *
     &            xyz2D0(iVec,ia,ib,ic+1,id,Ind2(n))
            End Do
           End Do
          End Do
         End Do
        End Do
       End Do
       Do n=1,mvec
        Do  id = 0, ld
         rc=-One
         Do  ic = 0, lc
          rc=rc+Two
          Do  ib=0,lb
           Do  ia = 0, la
            Do  iVec = 1, nRys*nT
                      xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n))=
     &                    Scrtch(iVec)**2*
     &                    xyz2D0(iVec,ia,ib,ic+2,id,Ind4(n))-
     &                    Rc* Scrtch(iVec) *
     &                    xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
            End Do
           End Do
          End Do
         End Do
        End Do
       End Do

       If (lc.ge.1) Then
        Do n=1,nvec
            Do  id=0,ld
           rc=Zero
           Do ic=1,lc
            rc=rc+One
          Do  ib=0,lb
         Do  ia=0,la
             Call DaXpY_inline(nRys*nT,-rc,
     &                  xyz2D0(1,ia,ib,ic-1,id,Ind2(n)),1,
     &                  xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
            End Do
           End Do
          End Do
         End Do
        End Do
       End If
       If (lc.ge.2) Then
        Do n=1,mvec
            Do id=0,ld
           rc=One
           Do ic=2,lc
            rc=rc+One
            Fact=rc*rc-rc
          Do  ib=0,lb
         Do  ia=0,la
             Call DaXpY_inline(nt*nrys,Fact,
     &                  xyz2D0(1,ia,ib,ic-2,id,Ind4(n)),1,
     &                  xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)) ,1)
            End Do
           End Do
          End Do
         End Do
        End Do
       End If
       End If
      End If
*
*         Cross term 1 3
*
      If (IfG(1).and.IfG(3)) Then
         Call ExpX(Temp  ,mZeta,mEta,Alpha,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
         Call ExpY(Temp  ,mZeta,mEta,Gamma,Sqrt(Two))
         Call Exp_2(Scrtch,nRys,nT,Temp,Sqrt(Two))
         nVec = 0
         If (ifHss(3,1,1,1)) Then
           mx=mx+1
           nVec = nVec + 1
           Ind1(nvec)=mx
           Ind2(nVec)=1
           Index4(1,mx,1)=3
           Index4(2,mx,1)=1
           Index2(1,3,1)=mx
         End If
         If (ifHss(3,2,1,2)) Then
           my=my+1
           nVec = nVec + 1
           Ind1(nvec)=my
           Ind2(nVec)=2
           Index4(1,my,2)=3
           Index4(2,my,2)=1
           Index2(2,3,1)=my
         End If
         If (ifHss(3,3,1,3)) Then
           mz=mz+1
           nVec = nVec + 1
           Ind1(nvec)=mz
           Ind2(nVec)=3
           Index4(1,mz,3)=3
           Index4(2,mz,3)=1
           Index2(3,3,1)=mz
         End If
         Do i=nVec+1,3
           Ind1(i)=0
         End Do
         If (nVec.ne.0) Then
         Do n=1,nvec
          Do  id = 0, ld
           Do  ic = 0, lc
            Do  ib = 0, lb
             Do  ia = 0, la
              Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                   Scrtch(iVec)  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia+1,ib,ic+1,id,Ind2(n))
              End Do
             End Do
            End Do
           End Do
          End Do
         End Do
         If (la.ge.1) Then
          Do n=1,nvec
           Do  id = 0, ld
            Do  ic = 0, lc
             Do  ib=0,lb
              ra=Zero
              Do  ia=1,la
                ra=ra+One
                Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  ra  *  Scrtch(iVec)*
     &                  xyz2D0(iVec,ia-1,ib,ic+1,id,Ind2(n))
                End Do
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         If (lc.ge.1) Then
          Do n=1,nVec
           Do  id = 0, ld
            rc=Zero
            Do  ic = 1, lc
             rc=rc+One
             Do  ib=0,lb
              Do  ia=0,la
                Do  iVec = 1, nRys*nT
                   xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rc  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia+1,ib,ic-1,id,Ind2(n))
                End Do
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         If ((la.ge.1).and.(lc.ge.1)) Then
          Do n=1,nvec
           Do  id = 0, ld
            rc=Zero
            Do  ic = 1, lc
             rc=rc+One
             Do  ib=0,lb
              ra=Zero
              Do  ia=1,la
               ra=ra+One
               Fact=rc*ra
               Call DaXpy_inline(nt*nrys,Fact,
     &                    xyz2D0(1,ia-1,ib,ic-1,id,Ind2(n)),1,
     &               xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)) ,1)
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         End If
      End If
*                 1 4

      If (IfG(4)) Then
       Call ExpY(Temp  ,mZeta,mEta,Delta,Sqrt(Two))
       Call Exp_2(Scrtch,nRys,nT,Temp,Sqrt(Two))
       If (IfG(1)) Then
         Call ExpX(Temp  ,mZeta,mEta,Alpha,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
         nVec=0
         If (ifHss(4,1,1,1)) Then
          mx=mx+1
          nVec = nVec + 1
          Ind1(nvec)=mx
          Ind2(nVec)=1
          Index2(1,4,1)=mx
          Index4(1,mx,1)=4
          Index4(2,mx,1)=1
         End If
         If (ifHss(4,2,1,2)) Then
          my=my+1
          nVec = nVec + 1
          Ind1(nvec)=my
          Ind2(nVec)=2
          Index2(2,4,1)=my
          Index4(1,my,2)=4
          Index4(2,my,2)=1
         End If
         If (ifHss(4,3,1,3)) Then
          mz=mz+1
          nVec = nVec + 1
          Ind1(nvec)=mz
          Ind2(nVec)=3
          Index2(3,4,1)=mz
          Index4(1,mz,3)=4
          Index4(2,mz,3)=1
         End If
         Do i=nVec+1,3
          Ind1(i)=0
         End Do
*
         If (nVec.ne.0) Then
*
         Do n=1,nvec
          Do  id = 0, ld
           Do  ic = 0, lc
            Do  ib = 0, lb
             Do  ia = 0, la
              Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                   Scrtch(iVec)  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia+1,ib,ic,id+1,Ind2(n))
              End Do
             End Do
            End Do
           End Do
          End Do
         End Do

         If (la.ge.1) Then
          Do n=1,nvec
           Do  id = 0, ld
            Do  ic = 0, lc
             Do  ib=0,lb
              ra=Zero
              Do  ia=1,la
               ra=ra+One
               Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  ra  *  Scrtch(iVec)*
     &                  xyz2D0(iVec,ia-1,ib,ic,id+1,Ind2(n))
               End Do
              End Do
             End Do
            End Do
           End Do
          End Do
         End If


         If (ld.ge.1) Then
           Do n=1,nvec
           rd=Zero
            Do  id = 1, ld
             rd=rd+One
             Do  ic = 0, lc
              Do  ib=0,lb
               Do ia=0,la
                Do iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rd  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia+1,ib,ic,id-1,Ind2(n))
                End Do
               End Do
              End Do
             End Do
            End Do
           End Do
         End If
         If ((la.ge.1).and.(ld.ge.1)) Then
          Do n=1,nvec
           rd=Zero
           Do  id = 1, ld
            rd=rd+One
            Do  ic = 0, lc
             Do  ib=0,lb
              ra=Zero
              Do  ia=1,la
                ra=ra+One
                Fact=rd*ra
                      Call DaxPy_inline( nRys*nT,Fact,
     &                  xyz2D0(1,ia-1,ib,ic,id-1,Ind2(n)),1,
     &                  xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         End If
        End If
*
*        Cross terms between 2 4
*
        If (IfG(2).and.Ifg(4)) Then
         Call ExpX(Temp  ,mZeta,mEta,Beta,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
         nVec=0
         If (ifHss(4,1,2,1)) Then
            mx=mx+1
            nVec = nVec + 1
            Ind1(nvec)=mx
            Ind2(nVec)=1
            Index2(1,4,2)=mx
            Index4(1,mx,1)=4
            Index4(2,mx,1)=2
         End If
         If (ifHss(4,2,2,2)) Then
            my=my+1
            nVec = nVec + 1
            Ind1(nvec)=my
            Ind2(nVec)=2
            Index2(2,4,2)=my
            Index4(1,my,2)=4
            Index4(2,my,2)=2
         End If
         If (ifHss(4,3,2,3)) Then
            mz=mz+1
            nVec = nVec + 1
            Ind1(nvec)=mz
            Ind2(nVec)=3
            Index2(3,4,2)=mz
            Index4(1,mz,3)=4
            Index4(2,mz,3)=2
         End If
         Do i=nVec+1,3
           Ind1(i)=0
         End Do
*
         If (nvec.ne.0) Then
*

         Do n=1,nvec
          Do  id = 0, ld
           Do  ic = 0, lc
            Do  ib = 0, lb
             Do  ia = 0, la
                Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                   Scrtch(iVec)  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia,ib+1,ic,id+1,Ind2(n))
                End Do
             End Do
            End Do
           End Do
          End Do
         End Do

         If (lb.ge.1) Then
          Do n=1,nvec
           Do  id = 0, ld
            Do  ic = 0, lc
             rb=Zero
             Do  ib=1,lb
              rb=rb+One
              Do  ia=0,la
                Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rb  *  Scrtch(iVec)*
     &                  xyz2D0(iVec,ia,ib-1,ic,id+1,Ind2(n))
                End Do
              End Do
             End Do
            End Do
           End Do
          End Do
         End If


         If (ld.ge.1) Then
          Do n=1,nvec
           rd=Zero
           Do   id = 1, ld
            rd=rd+One
            Do  ic = 0, lc
             Do  ib=0,lb
              Do  ia=0,la
                Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rd  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia,ib+1,ic,id-1,Ind2(n))
                End do
              End do
             End do
            End do
           End do
          End do
         End If
         If ((lb.ge.1).and.(ld.ge.1)) Then
          Do n=1,nvec
           rd=Zero
           Do  id = 1, ld
            rd=rd+One
            Do  ic = 0, lc
             rb=Zero
             Do ib=1,lb
              rb=rb+One
              Do ia=0,la
                Fact=rb*rd
                     Call DaxPy_inline(nt*nrys,Fact,
     &                  xyz2D0(1,ia,ib-1,ic,id-1,Ind2(n)),1,
     &               xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         End If
        End If
*
*       Cross Term 3 4
*
        If (IfG(3).and.IfG(4)) Then
         Call ExpY(Temp  ,mZeta,mEta,Gamma,Sqrt(Two))
         Call Exp_2(Scrtch2,nRys,nT,Temp,Sqrt(Two))
         nVec=0
         If (ifHss(4,1,3,1)) Then
           mx=mx+1
           nVec = nVec + 1
           Ind1(nvec)=mx
           Ind2(nVec)=1
           Index2(1,4,3)=mx
           Index4(1,mx,1)=4
           Index4(2,mx,1)=3
         End If
         If (ifHss(4,2,3,2)) Then
           my=my+1
           nVec = nVec + 1
           Ind1(nvec)=my
           Ind2(nVec)=2
           Index2(2,4,3)=my
           Index4(1,my,2)=4
           Index4(2,my,2)=3

         End If
         If (ifHss(4,3,3,3)) Then
           mz=mz+1
           nVec = nVec + 1
           Ind1(nvec)=mz
           Ind2(nVec)=3
           Index2(3,4,3)=mz
           Index4(1,mz,3)=4
           Index4(2,mz,3)=3
         End If
         Do i=nVec+1,3
          Ind1(i)=0
         End Do
*
         If (nvec.ne.0) Then
*
         Do n=1,nvec
          Do  id = 0, ld
           Do  ic = 0, lc
            Do  ib = 0, lb
             Do  ia = 0, la
              Do iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                   Scrtch(iVec)  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia,ib,ic+1,id+1,Ind2(n))
              End Do
             End Do
            End Do
           End Do
          End Do
         End Do

         If (lc.ge.1) Then
          Do n=1,nvec
           Do  id = 0, ld
            rc=Zero
            Do  ic = 1, lc
             rc=rc+One
             Do  ib=0,lb
              Do  ia=0,la
                Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rc  *  Scrtch(iVec)*
     &                  xyz2D0(iVec,ia,ib,ic-1,id+1,Ind2(n))
                 End Do
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         If (ld.ge.1) Then
          Do n=1,nvec
           rd=Zero
           Do   id = 1, ld
            rd=rd+One
            Do  ic = 0, lc
             Do  ib=0,lb
              Do  ia=0,la
                Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &               xyz2D2(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) -
     &                  rd  *  Scrtch2(iVec)*
     &                  xyz2D0(iVec,ia,ib,ic+1,id-1,Ind2(n))
                End DO
              End DO
             End DO
            End DO
           End DO
          End DO
         End If
         If ((lc.ge.1).and.(ld.ge.1)) Then
          Do n=1,nvec
           rd=Zero
           Do  id = 1, ld
            rd=rd+One
            rc=Zero
            Do  ic = 1, lc
             rc=rc+One
             Do  ib=0,lb
              Do ia=0,la
                Fact=rc*rd
                Call DaxPy_inline(nt*nRys,fact,
     &                  xyz2D0(1,ia,ib,ic-1,id-1,Ind2(n)),1,
     &               xyz2D2(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
              End Do
             End Do
            End Do
           End Do
          End Do
         End If
         End If
        End If
*
*     Differentiate with respect to the fourth center
*
          nvec=0
          If (IfGrad(1,4)) Then
            nx = nx + 1
            nVec = nVec + 1
            Ind1(nVec) = nx
            Ind2(nVec) = 1
            Index1(1,4) = nx
            Index3(nx,1)=4
          End If
          If (IfGrad(2,4)) Then
             ny = ny + 1
             nVec = nVec + 1
             Ind1(nVec) = ny
             Ind2(nVec) = 2
             Index1(2,4) = ny
             Index3(ny,2)=4
          End If
          If (IfGrad(3,4)) Then
             nz = nz + 1
             nVec = nVec + 1
             Ind1(nVec) = nz
             Ind2(nVec) = 3
             Index1(3,4) = nz
             Index3(nz,3)=4
          End If
          Do i=nvec+1,3
           Ind1(i)=0
          End Do
*
          mvec=0
          If (IfHss(4,1,4,1)) Then
           mx=mx+1
           mvec=mvec+1
           Ind3(mVec) = mx
           Ind4(mVec) = 1
           Index2(1,4,4) = mx
           Index4(1,mx,1)=4
           Index4(2,mx,1)=4
          End If
          If (IfHss(4,2,4,2)) Then
            my=my+1
            mvec=mvec+1
            Ind3(mVec) = my
            Ind4(mVec) = 2
            Index2(2,4,4) = my
            Index4(1,my,2)=4
            Index4(2,my,2)=4
          End If
          If (IfHss(4,3,4,3)) Then
            mz=mz+1
            mvec=mvec+1
            Ind3(mVec) = mz
            Ind4(mVec) = 3
            Index2(3,4,4) = mz
            Index4(1,mz,3)=4
            Index4(2,mz,3)=4
          End If
          Do i=mvec+1,3
           Ind3(i)=0
          End Do
          nvecx=max(nvec,mvec)
*
          If (nVecx.ne.0) Then
*
          Do n=1,nvec
           rd=-One
           Do  id = 0, ld
            rd=rd+Two
            Do  ic = 0, lc
             Do ib=0,lb
              Do  ia = 0, la
               Do  iVec = 1, nRys*nT
                xyz2D1(iVec,ia,ib,ic,id,Ind2(n),Ind1(n)) =
     &                   Scrtch(iVec) *
     &                 xyz2D0(iVec,ia,ib,ic,id+1,Ind2(n))
               End do
              End do
             End do
            End do
           End do
          End do
          Do n=1,mvec
           rd=-One
           Do  id = 0, ld
            rd=rd+Two
            Do  ic = 0, lc
             Do ib=0,lb
              Do  ia = 0, la
               Do  iVec = 1, nRys*nT
                     xyz2D2(iVec,ia,ib,ic,id,Ind4(n),Ind3(n))=
     &                   Scrtch(iVec)**2*
     &                   xyz2D0(iVec,ia,ib,ic,id+2,Ind4(n))-
     &                   rd* Scrtch(iVec) *
     &                   xyz2D0(iVec,ia,ib,ic,id,Ind4(n))
               End do
              End do
             End do
            End do
           End do
          End do
          If (ld.ge.1) Then
            Do n=1,nvec
                rd=Zero
                Do  id=1,ld
                 rd=rd+One
               Do  ic=0,lc
              Do  ib=0,lb
             Do  ia=0,la
                 Call Daxpy_inline(nt*nrys,-rd,
     &                      xyz2D0(1,ia,ib,ic,id-1,Ind2(n)),1,
     &                      xyz2D1(1,ia,ib,ic,id,Ind2(n),Ind1(n)),1)
                End DO
               End DO
              End DO
             End DO
            End DO
           End If
           If (ld.ge.2) Then
            Do n=1,mvec
                rd=One
                Do  id=2,ld
                 rd=rd+One
                 Fact=rd*rd-rd
               Do  ic=0,lc
              Do ib=0,lb
             Do ia=0,la
                 Call Daxpy_inline(nt*nrys,Fact,
     &                      xyz2D0(1,ia,ib,ic,id-2,Ind4(n)),1,
     &                      xyz2D2(1,ia,ib,ic,id,Ind4(n),Ind3(n)),1)
                End Do
               End Do
              End Do
             End Do
            End Do
           End If

           End If
      End If
*
*-----Sum over common centers
*
      Do iCent = 1, 3
         If (IfG(iCent)) Then
         Do jCent =  iCent+1, 4
            If (EQ(Coora(1,iCent),Coora(1,jCent))) Then
              If  (IfG(jCent)) Then
               Do iCar = 1, 3
                     i1 = Index2(iCar,iCent,iCent)
                     i2 = Index2(iCar,jCent,jCent)
                     i3 = Index2(iCar,jCent,iCent)
                     j4=Index1(iCar,jCent)
                     j5=Index1(iCar,iCent)
                     If (IfHss(jCent,iCar,jCent,iCar).and.
     &                   IfHss(iCent,iCar,iCent,iCar)) Then
                       Call DaXpY_inline(
     &                          nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),
     &                          One,xyz2D2(1,0,0,0,0,iCar,i2),1,
     &                          xyz2D2(1,0,0,0,0,iCar,i1),1)
                     End If
                     If (IfHss(jCent,iCar,iCent,iCar).and.
     &                   IfHss(iCent,iCar,iCent,iCar)) Then
                        Call DaXpY_inline(
     &                          nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),
     &                          Two,xyz2D2(1,0,0,0,0,iCar,i3),1,
     &                          xyz2D2(1,0,0,0,0,iCar,i1),1)
                     End If
                     If ((j4.ne.0).and.(j5.ne.0).and.
     &               (ifgrad(iCar,iCent).and.ifgrad(iCar,jCent)))
     &               Call DaXpY_inline(
     &                          nRys*nT*(la+1)*(lb+1)*
     &                          (lc+1)*(ld+1),
     &                          One,xyz2D1(1,0,0,0,0,iCar,j4),1,
     &                          xyz2D1(1,0,0,0,0,iCar,j5),1)
                     Do kCent=1,4
                      If (IfG(kCent)) Then
                       If ((kcent.ne.iCent).and.(kcent.ne.jcent))
     &                 Then
                        If (ifHss(kCent,iCar,jCent,iCar).or.
     &                       ifHss(jCent,iCar,kCent,iCar)) Then
                          i4=Index2(iCar,Max(kCent,jCent),
     &                                    Min(jCent,kCent))
                          i5=Index2(iCar,Max(kCent,iCent),
     &                                    Min(iCent,kCent))
                          Call DaXpY_inline(
     &                          nRys*nT*(la+1)*(lb+1)*(lc+1)*(ld+1),
     &                          One,xyz2D2(1,0,0,0,0,iCar,i4),1,
     &                          xyz2D2(1,0,0,0,0,iCar,i5),1)
                        End If
                       End If
                      End If
               End Do ! kCent
               End Do ! iCar
*
               IfG(jCent)=.false.
               Tr(jCent)=.false.
               Do jCar=1,3
               IfGrad(jcar,jCent)=.false.
                Do iIrrep=0,nIrrep-1
                   IndGrd(jCar,jcent,iIrrep)=0
                End Do
                Do kCent=1,4
                 Do kCar=1,3
                  IfHss(jCent,jCar,kCent,kCar)=.false.
                  IfHss(kCent,kCar,jCent,jCar)=.false.
                  Do iIrrep=0,nIrrep-1
                   IndHss(jCent,jCar,kCent,kCar,iIrrep)=0
                   IndHss(kCent,kCar,jCent,jCar,iIrrep)=0
                  End Do
                 End Do
                End Do
               End Do
*
              End If
         End If ! end eq
         End Do    ! jCent
         End If
      End Do       ! iCent
*
      nh(1)=mx
      nh(2)=my
      nh(3)=mz
      ng(1)=nx
      ng(2)=ny
      ng(3)=nz

      Return
      End
      Subroutine Daxpy_inline(nt,r,A,inca,B,incB)
      Implicit Real*8 (A-H,O-Z)
      Real*8 A(*),B(*)
      Do i=1,nt
       B(i)=B(i)+r*A(i)
      end do
      return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(inca)
         Call Unused_integer(incB)
      End If
      end
