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
      Subroutine Box_On_Sphere(x_Min_,x_Max_, y_Min_,
     &                         y_Max_, z_Min_,z_Max_,
     &                         xMin_, xMax_,  yMin_,
     &                         yMax_,  zMin_ ,zMax_ )
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 xyz(3,2), xyz0(3,2), Roots(3,3)
*
      delta=1.0D-15
*
      xyz(1,1)=x_min_
      xyz(1,2)=x_max_
      xyz(2,1)=y_min_
      xyz(2,2)=y_max_
      xyz(3,1)=z_min_
      xyz(3,2)=z_max_
c     Write (*,*)
c     Write (*,*) 'Box limitations'
c     Write (*,*) 'x:',xyz(1,1), xyz(1,2)
c     Write (*,*) 'y:',xyz(2,1), xyz(2,2)
c     Write (*,*) 'z:',xyz(3,1), xyz(3,2)
*
*     Set extremal values
*
      xyz0(1,1)= One
      xyz0(1,2)=-One
      xyz0(2,1)= One
      xyz0(2,2)=-One
      xyz0(3,1)= One
      xyz0(3,2)=-One
*
      Do ix = 1, 3
         iy=ix+1
         If (iy.gt.3) iy=1
         iz=iy+1
         If (iz.gt.3) iz=1
*
         xMax=xyz(ix,2)
         xMin=xyz(ix,1)
c        Write (*,*)
         Roots(1,iy)=xyz(iy,1)
         Roots(2,iy)=xyz(iy,2)
         If (xyz(iy,1)*xyz(iy,2).lt.Zero) Then
            ny_Roots=3
            Roots(3,iy)=Zero
         Else
            ny_Roots=2
         End If
         Roots(1,iz)=xyz(iz,1)
         Roots(2,iz)=xyz(iz,2)
         If (xyz(iz,1)*xyz(iz,2).lt.Zero) Then
            nz_Roots=3
            Roots(3,iz)=Zero
         Else
            nz_Roots=2
         End If
c        Call RecPrt('Roots','(3G25.12)',Roots,3,3)
*
         Do i = 1, ny_Roots
c           Write (*,*) 'i=',i,ny_Roots
c           Write (*,*)
            y = Roots(i,iy)
            Do j = 1, nz_Roots
c              Write (*,*) 'j=',j,nz_Roots
               z = Roots(j,iz)
*
               x=xMin
               r=sqrt(x**2+y**2+z**2)
c              Write (*,*) x/r
               If (r.eq.Zero) Then
                  x_r=Zero
               Else
                  x_r=x/r
               End If
               xyz0(ix,1)=Min(xyz0(ix,1),x_r)
               xyz0(ix,2)=Max(xyz0(ix,2),x_r)
*
               x=xMax
               r=sqrt(x**2+y**2+z**2)
c              Write (*,*) x/r
               If (r.eq.Zero) Then
                  x_r=Zero
               Else
                  x_r=x/r
               End If
               xyz0(ix,1)=Min(xyz0(ix,1),x_r)
               xyz0(ix,2)=Max(xyz0(ix,2),x_r)
*
            End Do
         End Do
      End Do
*
c     Write (*,*) 'xMin=',xyz0(1,1)
c     Write (*,*) 'xMax=',xyz0(1,2)
c     Write (*,*) 'yMin=',xyz0(2,1)
c     Write (*,*) 'yMax=',xyz0(2,2)
c     Write (*,*) 'zMin=',xyz0(3,1)
c     Write (*,*) 'zMax=',xyz0(3,2)
      xMin_=xyz0(1,1)-Delta
      xMax_=xyz0(1,2)+Delta
      yMin_=xyz0(2,1)-Delta
      yMax_=xyz0(2,2)+Delta
      zMin_=xyz0(3,1)-Delta
      zMax_=xyz0(3,2)+Delta
*
      Return
      End
