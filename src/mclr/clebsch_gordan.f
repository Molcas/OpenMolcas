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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      Real*8 Function Clebsch_Gordan(j1,m1,j2,m2,j,m)
*
*     Calculates the Clebsch-Gordan coefficients
*
      implicit Real*8(a-z)
*
      integer z,zmax,zmin
      If (j1.lt.0.0d0.or.j2.lt.0.0d0.or.j.lt.0.0d0)
     &Then
         Write(6,*) 'Error J is lower than 0'
         Call Abend()
      End If
      r=abs(2.0d0*j1-dint(2.0d0*j1))+abs(2.0d0*j2-dint(2.0d0*j2))
     & +abs(2.0d0*j -dint(2.0d0*j ))+abs(2.0d0*m1-dint(2.0d0*m1))
     & +abs(2.0d0*m2-dint(2.0d0*m2))+abs(2.0d0*m -dint(2.0d0*m))
      If (r.gt.1.0d-6)
     &Then
        Write(6,*) 'CG provided with not half integer'
         Call Abend()
      End If
      If (m1+m2.eq.m) Then
       Fct1=(2.0D0*j+1.0D0)*Fact(j1+j2-j)*Fact(j1-j2+j)*Fact(-j1+j2+j)
       Fct2=Fact(j1+j2+j+1.0D0)
       Fct=sqrt(Fct1/fct2)
       Fct=Fct*sqrt(Fact(j1+m1)*Fact(j1-m1)*Fact(j2+m2)
     &             *Fact(j2-m2)*Fact(j+m)*FacT(j-m))
       zmax=nint(Min(j1+j2-j,j1-m1,j2+m2))
       zmin=-nint(Min(j-j2+m1,j-j1-m2))
       sum=0.0d0
       zmin=Max(0,zmin)
       Do z=zmin,zmax
        T=DBLE((-1)**z)
        N=Facti(z)*FacT(j1+j2-j-DBLE(z))*fact(j1-m1-DBLE(z))
     &   *fact(j2+m2-DBLE(z))*
     &    fact(j-j2+m1+DBLE(z))*fact(j-j1-m2+DBLE(z))
        sum=sum+T/N
       End do
       Clebsch_Gordan=sum*fct
      else
       Clebsch_Gordan=0.0d0
      End If
      return
      end
*
      Function Facti(R)
      integer n,i,j,R
      Real*8 Facti
      n=R
      i=1
      If (n.eq.0) Then
        Facti=1.0d0
        return
      end if
      Do j=1,n
      i=i*j
      end do
      Facti=DBLE(i)
      return
      end
      Function Fact(R)
      Real*8 R,Fact
      integer n,i,j
      n=nint(R)
      i=1
      If (n.eq.0) Then
        Fact=1.0d0
        return
      end if
      Do j=1,n
      i=i*j
      end do
      Fact=DBLE(i)
      return
      end
