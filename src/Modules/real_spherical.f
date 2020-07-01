************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
* Copyright (C) 1990,2020, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Module Real_Spherical
      Private
#include "stdalloc.fh"
      Public :: ipSph, RSph, Sphere, Sphere_Free,
     &          Condon_Shortley_phase_factor
      Integer, Dimension(:), Allocatable :: ipSph
      Integer :: lmax_internal=-1
      Real*8, Dimension(:), Allocatable :: RSph
      Logical :: Condon_Shortley_phase_factor=.False.
*
***********************************************************************
*
      Contains
*
***********************************************************************
*
      SubRoutine Sphere_Free()
      If (Allocated(RSph)) Call mma_deallocate(RSph)
      If (Allocated(ipSph)) Call mma_deallocate(ipSph)
      lmax_internal=-1
      End SubRoutine Sphere_Free
*
***********************************************************************
*
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
* Calling    : Real_Sphere                                             *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
*               Credits.                                               *
*               2020, R. Lindh; P. R. Taylor; L. Birnoschi; A. Dzubak; *
*                     M. Navarrete; C. Gonzalez-Espinoza; G. Raggi;    *
*                     N. F. Chilton at OpenMolcas2020                  *
************************************************************************
      Implicit real*8 (a-h,o-z)
*     find MxAng, limiting the highest ang mom, in itmax.fh
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "status.fh"
*     iAngMx is the largest ang mom in the current basis
      iAngMx=Max(iAngMx,lMax)
*     check if required ang mom is greater than hard-coded limit
      If (iAngMx.gt.MxAng) Then
         Call WarningMessage(2,' Sphere: Increase MxAng!')
         Call Abend()
      End If
*
      If (lmax.lt.0) Then
         Write (6,*) 'Sphere: lmax<0'
         Call Abend()
      End If
      If (lmax.gt.lmax_internal) Then
         Call Sphere_Free()
         lmax_internal=lMax
      Else
         Return
      End If
*
*     Make the labels
*     Gives info on basis function angular momenta
*     n, l, ml or assigns it as a diffuse/polarising function with '*'
*
      Call Make_Labels(LblCbs,LblSbs,MxFnc,iAngMx)
*
*     Allocate memory for transformation matrices
*     Here, ipSph are the pointers to memory locations in RSph, for the
*     transformation matrices of given ang mom
      nSphr = 0
      Do iAng = 0, lMax
         nSphr = nSphr + (iAng*(iAng+1)/2 + iAng + 1)**2
      End Do
      Call mma_allocate(RSph,nSphr,label='RSph')
      Call mma_allocate(ipSph,[0,lMax],label='ipSph')
      ipSph(0)=1
      Do 2 iAng = 0, lMax-1
         ipSph(iAng+1) = ipSph(iAng) + (iAng*(iAng+1)/2 + iAng + 1)**2
 2    Continue

*     Here the transformation matrices from cartesian to spherical are
*     made
      Call Real_Sphere(ipSph,lMax,RSph,nSphr)
*
*     Set up the symmetry properties of the spherical gaussians
*     We are not sure if this Condon and Shortley phase....
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
*#define _DEBUG_
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) ' Spherical Harmonic expansions '
      Write (6,*)
      iLbl=1
      Do n = 0, lMax
         nElem = (n+1)*(n+2)/2
         ii  = 0
         Write (6,*)
         Write (6,'(8X,31(2X,I1,I1,I1))') ((i,j,n-i-j,
     &         j=n-i,0,-1),i=n,0,-1)
         Write (6,*)
         Do m = n, 0, -2
            Do l = -m, m
               Write (6,'(1X,A6,1X,31F5.2)')
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
      Real*8 RSph(nSphr)
      Integer ipSph(0:lMax)
*
      i00 = ipSph(0)
      i10 = ipSph(0)
      Do i = 0, lMax
         i2  = ipSph(i)
         nElem = (i+1)*(i+2)/2
         i20= i2 + i*nElem
*        First generate the coefficients for Y(i,0) -- always real
         Call Recurse(RSph(i00),RSph(i10),RSph(i20),i)
*        Use ladder operators to generate Y(i,m)
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
***********************************************************************
*                                                                     *
*     The Legendre polynomial is identical to Y(l,0).                 *
*     Note that it is real and that there is no Condon-Shortley phase *
*     factor to consider.                                             *
*                                                                     *
***********************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 P0((n2-1)*n2/2), P1(n2*(n2+1)/2),P2((n2+1)*(n2+2)/2)
*     Define statement function:
      iad(ix,iy,iz)=(iz+iy)*(iz+iy+1)/2 +iz + 1
*
      P2(:)=Zero
*
*---- Use recurrence relation for Legendre polynomials
*
*     (n+1) P_{n+1} = (2n+1) z P_n - n r^2 P_{n-1}
*
      If (n2.eq.0) then
*
         P2(1)=One
*
      Else
*
*        P_{n+1} = (2n+1)/(n+1) z P_n
*
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
*        P_{n+1} = - n/(n+1) (x^2+y^2+z^2) P_{n-1}
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
*
      End if

      Return
      End
      Subroutine Ladder(P0,n)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 P0((n+1)*(n+2)/2,-n:n)
*     Define statement function:
      iad(ix,iy,iz)=(iz+iy)*(iz+iy+1)/2 +iz + 1
*
*     Generate Y(l,m) from Y(l,m-1), starting the process from Y(l,0)
*
      Do m = 0, n-1
         m_p=  m+1
         m_m=-(m+1)
         P0(:,m_p)=Zero
         P0(:,m_m)=Zero
         Fact=One/(Two*Sqrt(DBLE(n*(n+1)-m*(m-1))))
*
*        The spherical harmonic is a two component (real,imaginary)
*        function.
*
*....... Y(n, m) =(-1)**  m  x (S(+,m), S(-,m)) and
*        Y(n,-m) =(-1)**(-m) x (S(+,m),-S(-,m))
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
*
*        Up to this point we have been operating on the Legendre and
*        associated Legendre polynomials. Let us now put in the
*        Condon-Shortley phase factor
*
         If (Condon_Shortley_phase_factor .and.
     &       MOD(m+1,2).ne.0) Then
            P0(:,m_p)=-P0(:,m_p)
            P0(:,m_m)=-P0(:,m_m)
         End If
*
      End Do ! m
*
      Return
      End
      Subroutine Contaminant(P0,i,Px,j,l)
*     This subroutine generates the lower ang mom contaminants for the
*     given ang mom
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 P0((i+1)*(i+2)/2,-l:l), Px((j+1)*(j+2)/2,-l:l)
*     Declare statement function
      iad(ix,iy,iz)=(iz+iy)*(iz+iy+1)/2 +iz + 1
*     Call QEnter('Contaminant')
*
*     Px = (x^2+y^2+z^2) x P0
*
      Do m = -l, l
         P0(:,m)=Zero
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
*
***********************************************************************
*
      End Module Real_Spherical
