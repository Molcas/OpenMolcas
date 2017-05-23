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
      subroutine XDR_fpFWprop(n,Tr,X,pXp,A,B,R,EL,ES,OL,OS,tmp)
C
C  Transform property operator (X,pXp) to fpFW picture
C
      implicit none
C Input
      integer n
      Real*8 Tr(n,n),X(n,n),pXp(n,n),A(n),B(n),R(n)
C Output
      Real*8 EL(n,n),ES(n,n),OL(n,n),OS(n,n)
C Temp
      integer i,j
      Real*8 tmp(n,n),av,aw
C
      call dmxma(n,'C','N',Tr,X,tmp,1.d0)
      call dmxma(n,'N','N',tmp,Tr,X,1.d0)
      call dmxma(n,'C','N',Tr,pXp,tmp,1.d0)
      call dmxma(n,'N','N',tmp,Tr,pXp,1.d0)
C
      do i=1,n
        do j=1,n
          av = X(j,i) * A(i) * A(j)
          aw = pXp(j,i) * B(i) * B(j)
          EL(j,i) = av + aw
          ES(j,i) = aw/R(i)/R(j) + av*R(i)*R(j)
          OL(j,i) = aw/R(i) - av*R(i)
          OS(j,i) = aw/R(j) - av*R(j)
        enddo
      enddo
C
      return
      end
