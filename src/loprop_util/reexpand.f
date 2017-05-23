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
      Subroutine ReExpand(rMP,nij,nElem,A,B,ij,lMax)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "itmax.fh"
#include "binom.fh"
      Real*8 rMP(nij,nElem), A(3), B(3)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      mElem(i)=(i+1)*(i+2)*(i+3)/6
*                                                                      *
************************************************************************
*                                                                      *
*
C     Call RecPrt('A',' ',A,1,3)
C     Call RecPrt('B',' ',B,1,3)
C     Call RecPrt('rMP',' ',rMP,nij,nElem)
      Do l = lMax, 0, -1
         iElem = mElem(l-1)
         Do ix = l, 0, -1
            ABx = A(1)-B(1)
            Do iy = l-ix, 0, -1
               ABy = A(2)-B(2)
               iz = l-ix-iy
               ABz = A(3)-B(3)
               iElem = iElem + 1
C              Write (*,*)
C              Write (*,*)
C              Write (*,*) 'ix,iy,iz=',ix,iy,iz
C              Write (*,*)
*
               temp = Zero
               Do jx = 0, ix
               Do jy = 0, iy
               Do jz = 0, iz
C              Write (*,*) 'jx,jy,jz=',jx,jy,jz
*
                  If (ix-jx.eq.0) then
                   ABx_=One
                  else
                   ABx_=ABx**(ix-jx)
                  endif
                  If (iy-jy.eq.0) then
                   ABy_=One
                  else
                   ABy_=ABy**(iy-jy)
                  endif
                  If (iz-jz.eq.0) then
                   ABz_=One
                  else
                   ABz_=ABz**(iz-jz)
                  endif
                  k=jx+jy+jz
                  jElem = mElem(k-1)+(jy+jz)*(jy+jz+1)/2+jz+1
C                 Write (6,*) 'jElem=',jElem
C                 Write (6,*)   binom(ix,jx),binom(iy,jy),binom(iz,jz),
C    &                   rMP(ij,jElem) ,
C    &                    ABx_ , ABy_ , ABz_
                  temp = temp + binom(ix,jx)*binom(iy,jy)*binom(iz,jz)
     &                  *rMP(ij,jElem)
     &                  * ABx_ * ABy_ * ABz_
*
               End Do
               End Do
               End Do
               rMP(ij,iElem)=temp
*
            End Do
         End Do
      End Do
C     Call RecPrt('rMP',' ',rMP,nij,nElem)
*
      Return
      End
