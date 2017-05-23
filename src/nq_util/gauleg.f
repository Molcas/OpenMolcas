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
      Subroutine gauleg(x1,x2,xw,n)
      Implicit None
#include "real.fh"
      Integer n
      Real*8 x1,x2,xw(2,n)
      Real*8 EPS
      Parameter(EPS=3.d-14)
      Integer i, j, m
      Real*8 p1, p2, p3, pp, xl, xm, z, z1
*
      m=(n+1)/2
      xm=Half*(x2+x1)
      xl=Half*(x2-x1)
c     Write (*,*) 'm=',m
      Do i = 1, m
         z=Cos(Pi*(DBLE(i)-0.25D0)/(DBLE(n)+Half))
 1       Continue
            p1=One
            p2=Zero
            Do j = 1, n
               p3=p2
               p2=p1
               p1=((Two*DBLE(j)-One)*z*p2-(DBLE(j)-One)*p3)/DBLE(j)
            End Do
            pp=DBLE(n)*(z*p1-p2)/(z*z-One)
            z1=z
            z=z1-p1/pp
            If (Abs(z-z1).gt.EPS) Go To 1
         xw(1,i    )=xm-xl*z
         xw(1,n+1-i)=xm+xl*z
         xw(2,i    )=Two*xl/((One-z*z)*pp*pp)
         xw(2,n+1-i)=xw(2,i)
         If (Abs(xw(1,i    )).lt.EPS) xw(1,i    )=Zero
         If (Abs(xw(1,n+1-i)).lt.EPS) xw(1,n+1-i)=Zero
         If (Abs(xw(2,i    )).lt.EPS) xw(1,i    )=Zero
         If (Abs(xw(2,n+1-i)).lt.EPS) xw(1,n+1-i)=Zero
      End Do
*
      Return
      End
      Subroutine gauleg_(x1,x2,xw,n)
      Implicit None
#include "real.fh"
      Integer n
      Real*8 x1,x2,xw(3,n)
      Real*8 EPS
      Parameter(EPS=3.d-14)
      Integer i, j, m
      Real*8 p1, p2, p3, pp, xl, xm, z, z1
*
      m=(n+1)/2
      xm=Half*(x2+x1)
      xl=Half*(x2-x1)
c     Write (*,*) 'm=',m
      Do i = 1, m
         z=Cos(Pi*(DBLE(i)-0.25D0)/(DBLE(n)+Half))
 1       Continue
            p1=One
            p2=Zero
            Do j = 1, n
               p3=p2
               p2=p1
               p1=((Two*DBLE(j)-One)*z*p2-(DBLE(j)-One)*p3)/DBLE(j)
            End Do
            pp=DBLE(n)*(z*p1-p2)/(z*z-One)
            z1=z
            z=z1-p1/pp
            If (Abs(z-z1).gt.EPS) Go To 1
         xw(1,i    )=xm-xl*z
         xw(1,n+1-i)=xm+xl*z
         xw(2,i    )=Two*xl/((One-z*z)*pp*pp)
         xw(2,n+1-i)=xw(2,i)
         If (Abs(xw(1,i    )).lt.EPS) xw(1,i    )=Zero
         If (Abs(xw(1,n+1-i)).lt.EPS) xw(1,n+1-i)=Zero
         If (Abs(xw(2,i    )).lt.EPS) xw(1,i    )=Zero
         If (Abs(xw(2,n+1-i)).lt.EPS) xw(1,n+1-i)=Zero
      End Do
*
      Return
      End
