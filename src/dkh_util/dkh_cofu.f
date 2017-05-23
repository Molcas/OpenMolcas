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
      subroutine dkh_cofu(DKH_order,DKH_para,dkcof)
C
C Calculate coefficients of different parameterization of DKH unitary transformation
C
      implicit none
C Input
      integer DKH_order,DKH_para
C Output
      Real*8 dkcof(*)
C Temp
      integer i,k,n
      Real*8 c,d
C
      N=max(4,DKH_order)
      if(DKH_para.eq.2)then
C       ! Exponential
        dkcof(1)=1.0d0
        do i=2,N
          dkcof(i)=dkcof(i-1)/i
        end do
      else if(DKH_para.eq.3)then
C       ! Square root
        do i=1,N
          dkcof(i)=0.d0
        end do
        dkcof(1)=1.0d0
        dkcof(2)=0.5d0
        do i=4,N,2
          dkcof(i)=-dkcof(i-2)*(i-3)/i
        end do
      else if(DKH_para.eq.5)then
C       ! Cayley
        dkcof(1)=1.0d0
        do i=2,N
          dkcof(i)=dkcof(i-1)*0.5d0
        end do
      else if(DKH_para.eq.4)then
C       ! McWeeny
        dkcof(1)=1.0d0
        dkcof(2)=0.5d0
        dkcof(3)=0.5d0
        do i=4,N,2
          dkcof(i)=dkcof(i-2)*(i-1)/i
          if(i.lt.n)then
            dkcof(i+1)=dkcof(i)
          end if
        end do
      else if(DKH_para.eq.1)then
C       ! opt ( M. Reiher )
        dkcof(1)=1.0d0
        dkcof(2)=0.5d0
        dkcof(3)=(2.d0-sqrt(2.d0))/4.d0
        dkcof(4)=dkcof(3)-0.125d0
        do i=5,n,2
C         ! W(i+3) in terms of a(k),k<i
          d=0.d0
          do k=(i+3)/2,i-1
            c=dkcof(k)*dkcof(i+3-k)
            if(k.ne.(i+3)/2) c=c*2.d0
            if(mod(k,2).eq.1) c=-c
            d=d-c
          end do
          dkcof(i)=d*sqrt(2.d0)
          if(i.lt.n)then
            dkcof(i+1)=dkcof(i)
          end if
        end do
      end if
C
      return
      end
