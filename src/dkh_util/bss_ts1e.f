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
      subroutine bss_ts1e(n,s,t,v,w,ul,us,clight)
C
C Evaluate the BSS Hamiltonian matrix and store the transform matrices
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
      integer m,i,j,k,nn,lwork,iW,info
      integer iTr,iBk,iEL,iES,iOL,iOS,iEp,iE0,iKC,itF,itmp
      integer iA,iB,iX,iM
C
C Transform to the free-particle Foldy-Wothuysen picture
C
      nn=n*n+4
      call getmem('Tr  ','ALLOC','REAL',iTr ,nn)
      call getmem('Back','ALLOC','REAL',iBk ,nn)
      call getmem('mEL ','ALLOC','REAL',iEL ,nn)
      call getmem('mES ','ALLOC','REAL',iES ,nn)
      call getmem('mOL ','ALLOC','REAL',iOL ,nn)
      call getmem('mOS ','ALLOC','REAL',iOS ,nn)
      call getmem('Ep  ','ALLOC','REAL',iEp ,n+4)
      call getmem('E0  ','ALLOC','REAL',iE0 ,n+4)
      call getmem('KC  ','ALLOC','REAL',iKC ,n*3+4)
      call XDR_fpFW(n,s,t,v,w,Work(iTr),Work(iBk),Work(iEL),Work(iES),
     &              Work(iOL),Work(iOS),Work(iEp),work(iE0),Work(iKC),
     &              Work(iKC+n),Work(iKC+2*n),clight)
C
C Diagonalize to get the X matrix
C
      do i=1,n
C       ! add (diagonal) kinetic matrix in fpFW picture
        Work(iEL+(i-1)*(n+1)) = Work(iEL+(i-1)*(n+1)) + Work(iE0+i-1)
        Work(iES+(i-1)*(n+1)) = Work(iES+(i-1)*(n+1)) - Work(iEp+i-1)
     &                          - clight*clight
      end do
      m = n+n
      lwork = 8*m
      call getmem('Fock ','ALLOC','REAL',itF ,m*m+4)
      call getmem('Eig  ','ALLOC','REAL',iW  ,m+4  )
      call getmem('Tmp  ','ALLOC','REAL',itmp,lwork  )
      k = 0
      do i=1,n
        do j=1,n
          Work(itF+(j-1)  +(i-1)*m  ) = Work(iEL+k)
          Work(itF+(j-1)  +(i+n-1)*m) = Work(iOL+k)
          Work(itF+(j-1+n)+(i-1)*m  ) = Work(iOS+k)
          Work(itF+(j-1+n)+(i+n-1)*m) = Work(iES+k)
          k = k + 1
        end do
      end do
      call dsyev_('V','L',m,Work(itF),m,Work(iW),Work(itmp),lwork,info)
      call getmem('tmpA','ALLOC','REAL',iA ,nn)
      call getmem('tmpB','ALLOC','REAL',iB ,nn)
      call getmem('matX','ALLOC','REAL',iX ,nn)
      k = 0
      do i=1,n
        do j=1,n
          Work(iA+k) = Work(itF+j-1+(i+n-1)*m)
          Work(iB+k) = Work(itF+j+n-1+(i+n-1)*m)
          k = k + 1
        end do
      end do
      call XDR_dmatinv(Work(iA),n)
      call dmxma(n,'N','N',Work(iB),Work(iA),Work(iX),1.d0)
C
C Apply decoupling transformation to the Fock matrix
C
      call dmxma(n,'C','N',Work(iX) ,Work(iOS),Work(iA) ,1.d0)
      call dmxma(n,'C','N',Work(iX) ,Work(iES),Work(iB) ,1.d0)
      call dmxma(n,'N','N',Work(iB) ,Work(iX) ,Work(iES),1.d0)
      call dmxma(n,'N','N',Work(iOL),Work(iX) ,Work(iB) ,1.d0)
      do i=0,n*n-1
C       ! X-projected Fock matrix
        Work(iEL+i)=Work(iEL+i)+Work(iA+i)+Work(iB+i)+Work(iES+i)
      end do
      call dmxma(n,'C','N',Work(iX),Work(iX),Work(iA),1.d0)
      do i=1,n
        Work(iA+(i-1)*(n+1)) = Work(iA+(i-1)*(n+1)) + 1.d0
      end do
C     ! renormalization matrix
      call XDR_dmatsqrt(Work(iA),n)
C     ! decoupled electron Fock matrix in moment space (eigenfunction of T matrix)
      call dmxma(n,'C','N',Work(iA),Work(iEL),Work(iB),1.d0)
      call dmxma(n,'N','N',Work(iB),Work(iA),v,1.d0)
C
C Back transform to non-orthogonal basis picture
C
      call dmxma(n,'C','N',Work(iBk),v,Work(iB),1.d0)
      call dmxma(n,'N','N',Work(iB),Work(iBk),v,1.d0)
C
C Calculate transform matrices in non-orthogonal basis space
C
      call dmxma(n,'N','N',Work(iX),Work(iA),Work(iB),1.d0)
C     ! iA/iB is the upper/lower part of transformation matrix in fpFW picture
      call getmem('TmpM','ALLOC','REAL',iM,nn*4+4)
      call XDR_mkutls(n,Work(iA),Work(iB),Work(iTr),Work(iBk),
     &     Work(iKC),Work(iKC+n),Work(iKC+2*n),ul,us,Work(iM),
     &     Work(iM+nn),Work(iM+nn*2),Work(iM+nn*3) )
      call getmem('TmpM','FREE','REAL',iM,nn*4+4)
C
C Free temp memories
C
      call getmem('Tr  ','FREE','REAL',iTr ,nn)
      call getmem('Back','FREE','REAL',iBk ,nn)
      call getmem('mEL ','FREE','REAL',iEL ,nn)
      call getmem('mES ','FREE','REAL',iES ,nn)
      call getmem('mOL ','FREE','REAL',iOL ,nn)
      call getmem('mOS ','FREE','REAL',iOS ,nn)
      call getmem('Ep  ','FREE','REAL',iEP ,n+4)
      call getmem('E0  ','FREE','REAL',iE0 ,n+4)
      call getmem('KC  ','FREE','REAL',iKC ,n*3+4)
      call getmem('Fock ','FREE','REAL',itF ,m*m+4)
      call getmem('Eig  ','FREE','REAL',iW  ,m+4  )
      call getmem('Tmp  ','FREE','REAL',itmp,lwork  )
      call getmem('tmpA','FREE','REAL',iA ,nn)
      call getmem('tmpB','FREE','REAL',iB ,nn)
      call getmem('matX','FREE','REAL',iX ,nn)
      return
      end
