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
      subroutine dkh_ts1e(n,s,t,v,w,ul,us,clight,dkord,xord,dkparam)
C
C Evaluate the arbitrary order DKH Hamiltonian ( and transform matrices ul & us )
C
      implicit none
#include "WrkSpc.fh"
C Input
      integer n,dkord,xord,dkparam
      Real*8 s(n,n),t(n,n),v(n,n),w(n,n),clight
C Ouput
C     ! v : store the transformed relativistic one-electron Hamiltonian
      Real*8 ul(n,n),us(n,n)
C Temp
      integer i,nn,word,vord,nz,m,n2
      integer iTr,iBk,iEL,iES,iOL,iOS,iEp,iE0,iKC,iCo,iSco
      integer iM,iZ,iW,iXL,iXS
C
C Transform Hamiltonian matrix to the free-particle Foldy-Wouthuysen picture
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
C Call DKH transformation routines
C
C     ! number of W matrices needed
      word = max( dkord/2, xord )
C     ! order of DKH needed, since high Xorder may need more transformation than DKHorder (for Hamiltonian)
      vord = max( dkord, word*2 )
      m = n*n
      nz = m*vord
      call getmem('Wsav','ALLOC','REAL', iW, m*xord*2+4)
      call getmem('Cof ','ALLOC','REAL',iCo, vord+8)
C     ! calculate expansion coefficient of general unitary transformation ( in terms of anti-Hermitian W )
      call dkh_cofu(vord,dkparam,Work(iCo))
C
      if(dkparam.eq.2) then
C       ! special routine for EXP parameterization ( with fewer matrix multiplication than general routine )
        call GetMem('NWork ','ALLOC','REAL',iM,m*5+4)
        call GetMem('NNWork','ALLOC','REAL',iZ,nz*3+4)
        do i=0,m*5
          Work(iM+i) = 0.d0
        end do
        do i=0,nz*3
          Work(iZ+i) = 0.d0
        end do
        call AODKHEXP(n,vord,xord,dkord,Work(iEp),Work(iE0),
     &                Work(iEL),Work(iES),Work(iOL),
     &                Work(iM),Work(iM+m),Work(iM+m*2),
     &                Work(iM+m*3),Work(iM+m*4),
     &                Work(iZ),Work(iZ+nz),Work(iZ+nz*2),
     &                Work(iW) )
        call GetMem('NWork ','FREE','REAL',iM,m*5+4)
        call GetMem('NNWork','FREE','REAL',iZ,nz*3+4)
      else
C       ! general parameterization routine
        call getmem('Cof2','ALLOC','REAL',iSCo,vord+8)
        call getmem('Mat ','ALLOC','REAL', iM, m*6+4)
        call getmem('Mat2','ALLOC','REAL', iZ, nz*10+4)
        call dkh_ham(n,dkord,xord,vord,Work(iEL),Work(iES),Work(iOL),
     &               Work(iOS),Work(iEp),Work(iE0),Work(iCo),Work(iSco),
     &               Work(iM),Work(iM+m),Work(iM+m*2),Work(iM+m*3),
     &               Work(iM+m*4),Work(iM+m*5),
     &               Work(iZ),Work(iZ+nz),Work(iZ+nz*2),Work(iZ+nz*3),
     &               Work(iZ+nz*4),Work(iZ+nz*5),Work(iZ+nz*6),
     &               Work(iZ+nz*7),Work(iZ+nz*8),Work(iZ+nz*9),
     &               Work(iW) )
        call getmem('Cof2','FREE','REAL',iSCo,vord+8)
        call getmem('Mat ','FREE','REAL',iM, m*6+4)
        call getmem('Mat2','FREE','REAL',iZ, nz*10+4)
      end if
C
C Calculate the transform matrices
C
      if(xord.gt.0)then
        call getmem('fpUL','ALLOC','REAL',iXL ,nn)
        call getmem('fpUS','ALLOC','REAL',iXS ,nn)
        n2 = n+n
        call getmem('TmpZ','ALLOC','REAL',iZ,n2*n2*3+4)
C       ! obtain transform matrices ( XL and XS )in fpFW picture
        call dkh_geneu(n,n2,xord,Work(iCo),Work(iW),Work(iXL),Work(iXS),
     &                 Work(iZ),Work(iZ+n2*n2),Work(iZ+n2*n2*2) )
        call getmem('TmpZ','FREE','REAL',iZ,n2*n2*3+4)
C
        call getmem('TmpM','ALLOC','REAL',iM,m*4+4)
C       ! convert to original basis picture
        call XDR_mkutls(n,Work(iXL),Work(iXS),Work(iTr),Work(iBk),
     &     Work(iKC),Work(iKC+n),Work(iKC+2*n),ul,us,Work(iM),
     &     Work(iM+m),Work(iM+m*2),Work(iM+m*3) )
        call getmem('TmpM','FREE','REAL',iM,m*4+4)
        call getmem('fpUL','FREE','REAL',iXL ,nn)
        call getmem('fpUS','FREE','REAL',iXS ,nn)
      end if
C
C Back transform Hamiltonian matrix to original non-orthogonal basis picture
C
      call dmxma(n,'C','N',Work(iBk),Work(iEL),Work(iES),1.d0)
      call dmxma(n,'N','N',Work(iES),Work(iBk),v,1.d0)
C
C Free temp memories
C
      call getmem('Cof ','FREE','REAL',iCo, vord+8)
      call getmem('Wsav','FREE','REAL',iW, m*xord*2+4)
      call getmem('Tr  ','FREE','REAL',iTr ,nn)
      call getmem('Back','FREE','REAL',iBk ,nn)
      call getmem('mEL ','FREE','REAL',iEL ,nn)
      call getmem('mES ','FREE','REAL',iES ,nn)
      call getmem('mOL ','FREE','REAL',iOL ,nn)
      call getmem('mOS ','FREE','REAL',iOS ,nn)
      call getmem('Ep  ','FREE','REAL',iEP ,n+4)
      call getmem('E0  ','FREE','REAL',iE0 ,n+4)
      call getmem('KC  ','FREE','REAL',iKC ,n*3+4)
C
      return
      end
