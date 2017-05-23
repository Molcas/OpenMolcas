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
      subroutine dkh_geneu(n,m,xord,c,w,xl,xs,t1,t2,t3)
C
C Calculate the DKH unitary transformation truncated at [xord] order
C   U_{DKH}=U_{0}U_{1}U_{2}...U_{xord}
C   U_{k}=\sum_{i=0}^{xord/k}c_{i}W_{k}^{i}
C
      implicit none
C
C Input
C   n    dimension of matrix
C   m    =2*n
C   c    expansion coefficients of unitary transformation in terms of W
C   w    stored W matrices
C
      integer n,m,xord
      Real*8 c(*),w(n,n,2,xord)
C Output
C   xl   upper part
C   xs   lower part
      Real*8 xl(n,n),xs(n,n)
C Temp
      Real*8 t1(m,m),t2(m,m),t3(m,m)
      integer i,j,k,iord
C
      do iord=1,xord
C       ! initial unit matrix
        do i=1,m
          do j=1,m
            if(j.eq.i)then
              t2(j,i) = 1.d0
            else
              t2(j,i) = 0.d0
            end if
          end do
        end do
        do k=1,xord/iord
          if(mod(k,2).eq.1)then
            if(k.eq.1)then
              do i=1,n
                do j=1,n
C                 ! unitary transformation in right side, had {\dag} to original definition of U
                  xs(j,i) = -w(j,i,1,iord)
                end do
              end do
            else
              call dmxma(n,'N','N',xl,w(1,1,1,iord),xs,-1.d0)
            end if
            do i=1,n
              do j=1,n
                t2(j  ,i+n) = t2(j  ,i+n) + xs(j,i) * c(k)
                t2(j+n,i  ) = t2(j+n,i  ) - xs(i,j) * c(k)
              end do
            end do
          else
            call dmxma(n,'C','N',w(1,1,1,iord),xs,xl,1.d0)
            do i=1,n
              do j=1,n
                t2(j+n,i+n) = t2(j+n,i+n) + xl(j,i) * c(k)
              end do
            end do
            call dmxma(n,'N','C',xs,w(1,1,1,iord),xl,1.d0)
            do i=1,n
              do j=1,n
                t2(j  ,i  ) = t2(j  ,i  ) + xl(j,i) * c(k)
              end do
            end do
          end if
        end do
        if(iord.eq.1)then
          do i=1,m
            do j=1,m
              t1(j,i) = t2(j,i)
            end do
          end do
        else
C         ! multiply U_{iord} in right side
          call dmxma(m,'N','N',t1,t2,t3,1.d0)
          do i=1,m
            do j=1,m
              t1(j,i) = t3(j,i)
            end do
          end do
        end if
      end do
      do i=1,n
        do j=1,n
          xl(j,i) = t1(j  ,i)
          xs(j,i) = t1(j+n,i)
        end do
      end do
C
      return
      end
