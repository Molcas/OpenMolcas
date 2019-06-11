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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Sphere(lMax)
************************************************************************
*                                                                      *
* Object: create the transformation matrices from cartesian gaussians  *
*         to spherical gaussians. By having these matricies being      *
*         defined dynamical we ensure that any extension of the        *
*         program to higher angular momentum is simply done by chang-  *
*         ing MxAng to the appropiate value.                           *
*         In addition, this will also allow us to have any order of    *
*         vectors in the matrix, i.e. we can have our own format or    *
*         any other odd order (MOLECULE).                              *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : GetMem                                                  *
*              Real_Sphere                                             *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      use Real_Spherical
      Implicit real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "status.fh"
*
      iAngMx=Max(iAngMx,lMax)
      If (iAngMx.gt.MxAng) Then
         Call WarningMessage(2,' Sphere: Increase MxAng!')
         Call Abend()
      End If
*
      If (Allocated(RSph)) Return
*
*     Make the labels
*
      Call Make_Labels(LblCbs,LblSbs,MxFnc,iAngMx)
*
*     Allocate memory for transformation matrices
*
      nSphr = 0
      Do iAng = 0, lMax
         nSphr = nSphr + (iAng*(iAng+1)/2 + iAng + 1)**2
      End Do
      Call mma_allocate(RSph,nSphr,label='RSph')
      Call mma_allocate(ipSph,[0,lMax],label='iSph')
      ipSph(0)=1
      Do 2 iAng = 0, lMax-1
         ipSph(iAng+1) = ipSph(iAng) + (iAng*(iAng+1)/2 + iAng + 1)**2
 2    Continue
*
      Call Real_Sphere(ipSph,lMax,RSph,nSphr)
*
*     Set up the symmetry properties of the spherical gaussians
*
      iii = 0
      jjj = 0
      Do 50 n = 0, lMax
         nElem = (n+1)*(n+2)/2
         ii = 0
         Do 55 m = n, 0, -2
            Do 60 l = -m, m
               iii = iii + 1
               Do 65 iElem = 1, nElem
                  If (RSph(iElem-1+ii+ipSph(n)).ne.Zero) Go To 66
65             Continue
66             iSphCr(iii) = iElem + jjj
               ii = ii + nElem
60          Continue
55       Continue
        jjj = jjj + nElem
50    Continue
*
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) ' Spherical Harmonic expansions '
      Write (6,*)
      iLbl=1
      Do n = 0, lMax
         nElem = (n+1)*(n+2)/2
         ii  = 0
         Write (6,*)
         Write (6,'(6X,31(2X,I1,I1,I1))') ((i,j,n-i-j,
     &         j=n-i,0,-1),i=n,0,-1)
         Write (6,*)
         Do m = n, 0, -2
            Do l = -m, m
               Write (6,'(1X,A4,1X,31F5.2)')
     &            LblSbs(iLbl),(RSph(i+ii+ipSph(n)),i=0,nElem-1)
               ii = ii + nElem
               iLbl = iLbl + 1
            End Do
            Write (6,*)
         End Do
         Write (6,*)
      End Do
#endif
*
      Return
      End
      Subroutine Real_Sphere(ipSph,lMax,RSph,nSphr)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 RSph(nSphr)
      Integer ipSph(0:lMax)
*
*     Call QEnter('Real_Sphere')
*
      i00 = ipSph(0)
      i10 = ipSph(0)
      Do i = 0, lMax
         i2  = ipSph(i)
         nElem = (i+1)*(i+2)/2
         i20= i2 + i*nElem
*        Write (*,*) i00,i10,i20,i
         Call Recurse(RSph(i00),RSph(i10),RSph(i20),i)
         Call Ladder(RSph(i2),i)
*
*        Now do the contaminant, by simply multiply with r**2
*
         j = i-2
         If (j.ge.0) Then
            iCont = i2 + (2*i+1)*nElem
            iOff=ipSph(j)
            mElem = (j+1)*(j+2)/2
            Do l = j, 0, -2
               Call Contaminant(RSph(iCont),i,RSph(iOff),j,l)
               iCont = iCont + (2*l+1)*nElem
               iOff  = iOff  + (2*l+1)*mElem
            End Do
         End If
*
         i00=i10
         i10=i20
      End Do
*
*.... Normalize
*
      Do i = 0, lMax
         Call NrmSph(RSph(ipSph(i)),i)
      End Do
*
      Return
      End
      Subroutine Recurse(P0,P1,P2,n2)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 P0((n2-1)*n2/2), P1(n2*(n2+1)/2),P2((n2+1)*(n2+2)/2)
*
      iad(ix,iy,iz)=(iz+iy)*(iz+iy+1)/2 +iz + 1
*
*     Call QEnter('Recurse')
*
      call dcopy_((n2+1)*(n2+2)/2,[Zero],0,P2,1)
*
*---- Use recurrence relation for Lagrange polynomials
*
*
*     (n+1) P_{n+1} = (2n+1) z P_n - n r^2 P_{n-1}
*
      If (n2.eq.0) then
         P2(1)=One
*
      Else
         Fact_1=DBLE(2*n2-1)/DBLE(n2)
         n1=n2-1
         Do ix = n1, 0, -1
           Do iy = n1-ix, 0, -1
             iz = n1-ix-iy
             P2(iad(ix,iy,iz+1)) = P2(iad(ix,iy,iz+1))
     &            + Fact_1*P1(iad(ix,iy,iz))
           End Do
         End Do
*
         Fact_2=DBLE(n2-1)/DBLE(n2)
         n0=n1-1
         Do ix = n0, 0, -1
            Do iy = n0-ix, 0, -1
               iz = n0-ix-iy
               P2(iad(ix+2,iy,iz)) = P2(iad(ix+2,iy,iz))
     &                             - Fact_2*P0(iad(ix,iy,iz))
               P2(iad(ix,iy+2,iz)) = P2(iad(ix,iy+2,iz))
     &                             - Fact_2*P0(iad(ix,iy,iz))
               P2(iad(ix,iy,iz+2)) = P2(iad(ix,iy,iz+2))
     &                             - Fact_2*P0(iad(ix,iy,iz))
            End Do
         End Do
      End if

*     m=(n2-1)*n2/2
*     If (m.ge.1) Call RecPrt('P0',' ',P0,m,1)
*     m=n2*(n2+1)/2
*     If (m.ge.1) Call RecPrt('P1',' ',P1,m,1)
*     m = (n2+1)*(n2+2)/2
*     If (m.ge.1) Call RecPrt('P2',' ',P2,m,1)
*
*     Call QExit('Recurse')
      Return
      End
      Subroutine Ladder(P0,n)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 P0((n+1)*(n+2)/2,-n:n)
*
      iad(ix,iy,iz)=(iz+iy)*(iz+iy+1)/2 +iz + 1
*
*     Call QEnter('Ladder')
*
*     Call RecPrt('Ladder (in)',' ',P0,(n+1)*(n+2)/2,2*n+1)
      Do m = 0, n-1
         m_p=m+1
         m_m=-(m+1)
         call dcopy_((n+1)*(n+2)/2,[Zero],0,P0(1,m_p),1)
         call dcopy_((n+1)*(n+2)/2,[Zero],0,P0(1,m_m),1)
         Fact=One/(Two*Sqrt(DBLE(n*(n+1)-m*(m-1))))
*
*        The spherical harmonic is a two component (real,imaginary)
*        function.
*
*....... Y(n,m) =(S(+,m),S(-,m)) and Y(n,-m)=(S(+,m),-S(-,m))
*
*        with S(-,0)=0
*
*        The ladder operator is subdivided in a similar way
*
*        L(+)=(Lr,Li),  L(-)=(Lr,-Li)
*
*        Hence
*
*        L(+) Y(n,m)= C x Y(n,m+1)
*
*        or
*
*        C(S(+,m+1),S(-,m+1))=(Lr S(+,m)-Li S(-,m),Li S(+,m)+Lr S(-,m))
*
         Do ix = n, 0, -1
            Do iy = n-ix, 0, -1
               iz = n-ix-iy
*
*............. Generating the real part
*
               If (iz.ge.1)
     &         P0(iad(ix+1,iy,iz-1),m_p)= P0(iad(ix+1,iy,iz-1),m_p)
     &                            + Fact*DBLE(iz)*P0(iad(ix,iy,iz),m)
               If (ix.ge.1)
     &         P0(iad(ix-1,iy,iz+1),m_p)= P0(iad(ix-1,iy,iz+1),m_p)
     &                            - Fact*DBLE(ix)*P0(iad(ix,iy,iz),m)
               If (m.ne.0) Then
                  If (iz.ge.1)
     &            P0(iad(ix,iy+1,iz-1),m_p)= P0(iad(ix,iy+1,iz-1),m_p)
     &                            - Fact*DBLE(iz)*P0(iad(ix,iy,iz),-m)
                  If (iy.ge.1)
     &            P0(iad(ix,iy-1,iz+1),m_p)= P0(iad(ix,iy-1,iz+1),m_p)
     &                            + Fact*DBLE(iy)*P0(iad(ix,iy,iz),-m)
               End If
*
*............. Generating the imaginary part
*
               If (iz.ge.1)
     &         P0(iad(ix,iy+1,iz-1),m_m)= P0(iad(ix,iy+1,iz-1),m_m)
     &                            + Fact*DBLE(iz)*P0(iad(ix,iy,iz),m)
               If (iy.ge.1)
     &         P0(iad(ix,iy-1,iz+1),m_m)= P0(iad(ix,iy-1,iz+1),m_m)
     &                            - Fact*DBLE(iy)*P0(iad(ix,iy,iz),m)
               If (m.ne.0) Then
                  If (iz.ge.1)
     &            P0(iad(ix+1,iy,iz-1),m_m)= P0(iad(ix+1,iy,iz-1),m_m)
     &                            + Fact*DBLE(iz)*P0(iad(ix,iy,iz),-m)
                  If (ix.ge.1)
     &            P0(iad(ix-1,iy,iz+1),m_m)= P0(iad(ix-1,iy,iz+1),m_m)
     &                            - Fact*DBLE(ix)*P0(iad(ix,iy,iz),-m)
               End If
*
            End Do
         End Do
      End Do
*
*     Call QExit('Ladder')
      Return
      End
      Subroutine Contaminant(P0,i,Px,j,l)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 P0((i+1)*(i+2)/2,-l:l), Px((j+1)*(j+2)/2,-l:l)
*
      iad(ix,iy,iz)=(iz+iy)*(iz+iy+1)/2 +iz + 1
*
*     Call QEnter('Contaminant')
*
      Do m = -l, l
         call dcopy_((i+1)*(i+2)/2,[Zero],0,P0(1,m),1)
         Do ix = j, 0, -1
            Do iy = j-ix, 0, -1
               iz = j-ix-iy
               P0(iad(ix+2,iy,iz),m)=P0(iad(ix+2,iy,iz),m)
     &                              +Px(iad(ix,iy,iz),m)
               P0(iad(ix,iy+2,iz),m)=P0(iad(ix,iy+2,iz),m)
     &                              +Px(iad(ix,iy,iz),m)
               P0(iad(ix,iy,iz+2),m)=P0(iad(ix,iy,iz+2),m)
     &                              +Px(iad(ix,iy,iz),m)
            End Do
         End Do
      End Do
*
*     Call QExit('Contaminant')
      Return
      End
      Subroutine NrmSph(P,n)
      Implicit Real*8 (a-h,o-z)
      Real*8 P((n+1)*(n+2)/2,(n+1)*(n+2)/2)
#include "real.fh"
*
      iad(ix,iy,iz)=(iy+iz)*(iy+iz+1)/2+iz+1
*
      Do m = 1, (n+1)*(n+2)/2
         rMax=Zero
         Do k = 1, (n+1)*(n+2)/2
            If (Abs(P(k,m)).gt.rMax) rMax=Abs(P(k,m))
         End Do
         Do k = 1, (n+1)*(n+2)/2
            If (Abs(P(k,m)).lt.1.0D-12*rMax) P(k,m)=Zero
         End Do
         tmp=Zero
         Do ijx = 2*n, 0, -2
            Do ijy = 2*n-ijx, 0, -2
               ijz=2*n-ijx-ijy
               DF=DblFac(ijx-1)*DblFac(ijy-1)*DblFac(ijz-1)
               temp=Zero
               Do ix = Min(n,ijx), Max(0,ijx-n), -1
                  jx=ijx-ix
                  Do iy = Min(n-ix,ijy), Max(0,ijy-n+jx), -1
                     jy=  ijy-iy
                     iz=n-ix-iy
                     jz=n-jx-jy
                     temp=temp+P(iad(ix,iy,iz),m)*P(iad(jx,jy,jz),m)
                  End Do
               End Do
               tmp=tmp+DF*temp
            End Do
         End Do
         Call DScal_((n+1)*(n+2)/2,One/Sqrt(tmp),P(1,m),1)
      End Do
      Return
      End
      Function DblFac(n)
************************************************************************
*                                                                      *
* Object: to compute the double factorial of n.                        *
*                                                                      *
* Called from: NrmSph                                                  *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 DblFac
*
      DblFac = One
      Do 20 i = n , 1, -2
         DblFac = DblFac * DBLE(i)
 20   Continue
      Return
      End
