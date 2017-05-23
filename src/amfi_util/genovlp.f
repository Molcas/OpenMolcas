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
      subroutine genovlp(Lhigh,coulovlp)
      implicit real*8 (a-h,o-z)
cbs   generates overlap of normalized  primitives.
#include "para.fh"
#include "param.fh"
      dimension evecinv(MxprimL,MxprimL)
     *,coulovlp(MxprimL,MxprimL,-1:1,-1:1,0:Lmax,0:Lmax)
      do L=0,Lhigh
              do Jrun=1,nprimit(L)
              do Irun=1,nprimit(L)
        normovlp(Irun,Jrun,L)=coulovlp(irun,jrun,0,0,L,L)
              enddo
              enddo
cbs   invert the matrix, not very elegant, but sufficient
      ipnt=0
      do jrun=1,nprimit(L)
      do irun=1,jrun
      ipnt=ipnt+1
      scratchinv(ipnt)=normovlp(irun,jrun,L)
      enddo
      enddo
      do Jrun=1,nprimit(L)
      do Irun=1,MxprimL
      evecinv(Irun,Jrun)=0d0
      enddo
      enddo
      do Jrun=1,nprimit(L)
      evecinv(jrun,jrun)=1d0
      enddo
      call Jacob(scratchinv,evecinv,nprimit(L),MxprimL)
      do irun=1,nprimit(L)
      eval(irun)=sqrt(scratchinv((irun*irun+irun)/2))
      enddo
cbs   ensure normalization of the vectors.
      do IRUN=1,nprimit(L)
      fact=0d0
      do JRUN=1,nprimit(L)
      fact=fact+evecinv(JRUN,IRUN)*evecinv(JRUN,IRUN)
      enddo
      fact=1d0/sqrt(fact)
      do JRUN=1,nprimit(L)
      evecinv(JRUN,IRUN)=fact*evecinv(JRUN,IRUN)
      enddo
      enddo
cbs   now generate rootOVLP
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLP(irun,jrun,l)=0d0
      enddo
      enddo
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLP(irun,jrun,L)=rootOVLP(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo
      enddo
      enddo
cbs   now generate rootOVLPinv
      do irun=1,nprimit(L)
      eval(irun)=1d0/eval(irun)
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLPinv(irun,jrun,l)=0d0
      enddo
      enddo
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLPinv(irun,jrun,L)=rootOVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo
      enddo
      enddo
cbs   now generate OVLPinv
      do irun=1,nprimit(L)
      eval(irun)=eval(irun)*eval(irun)
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      OVLPinv(irun,jrun,l)=0d0
      enddo
      enddo
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      OVLPinv(irun,jrun,L)=OVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo
      enddo
      enddo
      enddo
      return
      end
