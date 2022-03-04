!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine genovlp(Lhigh,coulovlp)
!bs   generates overlap of normalized  primitives.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "para.fh"
integer(kind=iwp) :: Lhigh
real(kind=wp) :: coulovlp(MxprimL,MxprimL,-1:1,-1:1,0:Lmax,0:Lmax)
#include "param.fh"
integer(kind=iwp) :: ipnt, Irun, Jrun, krun, L
real(kind=wp) :: evecinv(MxprimL,MxprimL), fact !IFG

do L=0,Lhigh
  do Jrun=1,nprimit(L)
    do Irun=1,nprimit(L)
      normovlp(Irun,Jrun,L) = coulovlp(irun,jrun,0,0,L,L)
    end do
  end do
  !bs invert the matrix, not very elegant, but sufficient
  ipnt = 0
  do jrun=1,nprimit(L)
    do irun=1,jrun
      ipnt = ipnt+1
      scratchinv(ipnt) = normovlp(irun,jrun,L)
    end do
  end do
  do Jrun=1,nprimit(L)
    do Irun=1,MxprimL
      evecinv(Irun,Jrun) = Zero
    end do
  end do
  do Jrun=1,nprimit(L)
    evecinv(jrun,jrun) = One
  end do
  call Jacob(scratchinv,evecinv,nprimit(L),MxprimL)
  do irun=1,nprimit(L)
    eval(irun) = sqrt(scratchinv((irun*irun+irun)/2))
  end do
  !bs ensure normalization of the vectors.
  do IRUN=1,nprimit(L)
    fact = Zero
    do JRUN=1,nprimit(L)
      fact = fact+evecinv(JRUN,IRUN)*evecinv(JRUN,IRUN)
    end do
    fact = One/sqrt(fact)
    do JRUN=1,nprimit(L)
      evecinv(JRUN,IRUN) = fact*evecinv(JRUN,IRUN)
    end do
  end do
  !bs now generate rootOVLP
  do irun=1,nprimit(L)
    do jrun=1,nprimit(L)
      rootOVLP(irun,jrun,l) = Zero
    end do
  end do
  do jrun=1,nprimit(L)
    do irun=1,nprimit(L)
      do krun=1,nprimit(L)
        rootOVLP(irun,jrun,L) = rootOVLP(irun,jrun,L)+evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      end do
    end do
  end do
  !bs now generate rootOVLPinv
  do irun=1,nprimit(L)
    eval(irun) = One/eval(irun)
  end do
  do irun=1,nprimit(L)
    do jrun=1,nprimit(L)
      rootOVLPinv(irun,jrun,l) = Zero
    end do
  end do
  do jrun=1,nprimit(L)
    do irun=1,nprimit(L)
      do krun=1,nprimit(L)
        rootOVLPinv(irun,jrun,L) = rootOVLPinv(irun,jrun,L)+evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      end do
    end do
  end do
  !bs now generate OVLPinv
  do irun=1,nprimit(L)
    eval(irun) = eval(irun)*eval(irun)
  end do
  do irun=1,nprimit(L)
    do jrun=1,nprimit(L)
      OVLPinv(irun,jrun,l) = Zero
    end do
  end do
  do jrun=1,nprimit(L)
    do irun=1,nprimit(L)
      do krun=1,nprimit(L)
        OVLPinv(irun,jrun,L) = OVLPinv(irun,jrun,L)+evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      end do
    end do
  end do
end do

return

end subroutine genovlp
