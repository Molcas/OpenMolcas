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
      subroutine x2c_ts1e(n,s,t,v,w,ul,us,clight)
C
C Evaluate the X2C Hamiltonian matrix and store the transform matrices
C
      implicit none
#include "WrkSpc.fh"
C Input
C w   aka pVp
      integer n
      Real*8 s(n,n),t(n,n),v(n,n),w(n,n),clight
C Output
      Real*8 ul(n,n),us(n,n)
C Temp
      integer m,i,j,k
      Real*8 c_2c2,c_4c2,c_2c
      integer itF,itS,itX,itA,itB,itC,itSS
C
C Construct the full dimensional (2*n) Fock matrix and overlap matrix
C
      m = n+n
      c_2c  = clight*2.d0
      c_2c2 = clight*clight*2.d0
      c_4c2 = c_2c2*2.d0
C     ! rescale, in order to make X matrix close to unit matrix
      do i=1,n
        do j=1,n
          w(j,i) = w(j,i)/c_4c2
        end do
      end do
      call getmem('TmpF ','ALLOC','REAL', itF, m*m+4 )
      call getmem('TmpS ','ALLOC','REAL', itS, m*m+4 )
      do k=0,m*m-1
        Work(itS+k) = 0.d0
      end do
      do i=1,n
        do j=1,n
          Work(itS+(j-1)  +(i-1)*m)   = s(j,i)
          Work(itS+(j+n-1)+(i+n-1)*m) = t(j,i)/c_2c2
          Work(itF+(j-1)  +(i-1)*m)   = v(j,i)
          Work(itF+(j-1)  +(i+n-1)*m) = t(j,i)
          Work(itF+(j+n-1)+(i-1)*m)   = t(j,i)
          Work(itF+(j+n-1)+(i+n-1)*m) = w(j,i)-t(j,i)
        end do
      end do
C
C Call diagonalization routine to obtain the X matrix
C
      call getmem('TmpX ','ALLOC','REAL', itX, n*n+4 )
      call x2c_makx(m,n,Work(itF),Work(itS),Work(itX))
C
C Calculate transformed Hamiltonian matrix
C
      call getmem('TmpA ','ALLOC','REAL', itA, n*n+4 )
      call getmem('TmpB ','ALLOC','REAL', itB, n*n+4 )
      call getmem('TmpC ','ALLOC','REAL', itC, n*n+4 )
      call getmem('TmpSS','ALLOC','REAL', itSS, n*n+4 )
      call dmxma(n,'C','N',Work(itX),t,Work(itA),1.d0)
      call dmxma(n,'N','N',t,Work(itX),Work(itB),1.d0)
      call dmxma(n,'N','N',Work(itA),Work(itX),Work(itC),1.d0)
      k = 0
      do i=1,n
        do j=1,n
C          ! X-projected overlap matrix
          Work(itSS+k) = s(j,i) + Work(itC+k)/c_2c2
C          ! X-projected kinetic matrix
          t(j,i) = Work(itA+k) + Work(itB+k) - Work(itC+k)
          k = k + 1
        end do
      end do
      call XDR_dmatsqrt(s,n)
      call dmxma(n,'C','N',s,Work(itSS),Work(itA),1.d0)
      call dmxma(n,'N','N',Work(itA),s,Work(itB),1.d0)
      call XDR_dmatsqrt(Work(itB),n)
      call dmxma(n,'N','N',s,Work(itB),Work(itC),1.d0)
      call XDR_dmatinv(s,n)
C      ! renormalization matrix, also the upper part of transformation matrix
      call dmxma(n,'N','N',Work(itC),s,ul,1.d0)
C      ! lower part of the transformation matrix
      call dmxma(n,'N','N',work(itX),ul,us,1.d0)
C
C Apply transformation to kinetic and potential matrices
C
      call dmxma(n,'C','N',ul,t,work(itA),1.d0)
      call dmxma(n,'N','N',work(itA),ul,t,1.d0)
      call dmxma(n,'C','N',ul,v,work(itA),1.d0)
      call dmxma(n,'N','N',work(itA),ul,v,1.d0)
      call dmxma(n,'C','N',us,w,work(itA),1.d0)
      call dmxma(n,'N','N',work(itA),us,w,1.d0)
      do i=1,n
        do j=1,n
          v(j,i) = t(j,i) + v(j,i) + w(j,i)
C          ! since pVp was rescaled, the lower part of transformation matrix also need to be rescaled
          us(j,i) = us(j,i)/c_2c
        end do
      end do
C
C Free temp memories
C
      call getmem('TmpF ','FREE','REAL', itF, m*m+4 )
      call getmem('TmpS ','FREE','REAL', itS, m*m+4 )
      call getmem('TmpX ','FREE','REAL', itX, n*n+4 )
      call getmem('TmpA ','FREE','REAL', itA, n*n+4 )
      call getmem('TmpB ','FREE','REAL', itB, n*n+4 )
      call getmem('TmpC ','FREE','REAL', itC, n*n+4 )
      call getmem('TmpSS','FREE','REAL', itSS, n*n+4 )
      return
      end
