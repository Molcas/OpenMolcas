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
C
C----------------------------------------------------------------------|
C
      subroutine XDR_mkutls(n,TL,TS,Tr,Bk,A,B,R,UL,US,M1,M2,M3,M4)
C
C Evaluate transform matrices in non-orthogonal basis space
C
      implicit none
C Input
      integer n
      Real*8 TL(n,n),TS(n,n),Tr(n,n),Bk(n,n),A(n),B(n),R(n)
C Output
      Real*8 UL(n,n),US(n,n)
C Temp
      integer i,j
      Real*8 M1(n,n),M2(n,n),M3(n,n),M4(n,n)
C
      do i=1,n
        do j=1,n
          M1(j,i) = Tr(j,i) * A(i)
          M2(j,i) = Tr(j,i) * A(i) * R(i)
        end do
      end do
      call dmxma(n,'N','N',M1,TL,M3,1.d0)
      call dmxma(n,'N','N',M2,TS,M4,1.d0)
      do i=1,n
        do j=1,n
          M3(j,i) = M3(j,i) - M4(j,i)
        end do
      end do
      call dmxma(n,'N','N',M3,Bk,UL,1.d0)
C
      do i=1,n
        do j=1,n
          M1(j,i) = Tr(j,i) * B(i)
          M2(j,i) = Tr(j,i) * B(i) / R(i)
        end do
      end do
      call dmxma(n,'N','N',M1,TL,M3,1.d0)
      call dmxma(n,'N','N',M2,TS,M4,1.d0)
      do i=1,n
        do j=1,n
          M3(j,i) = M3(j,i) + M4(j,i)
        end do
      end do
      call dmxma(n,'N','N',M3,Bk,US,1.d0)
C
      return
      end
