!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2008, Jonas Bostrom                                    *
!***********************************************************************

subroutine Construct_WDensIII(Xaibj,LnPQRSprod,LiPQRSprod,iBatch,jBatch,nOccLeftI,nOccLeftJ)
!
! Jonas Bostrom, October 2008
!
! Purpose: Construct the piece of the energy-weighted density
!          usually labeled III.

use ChoMP2, only: iFirstS, LnBatOrb, LnPQprod, LiPQprod

implicit real*8(a-h,o-z)
integer LnPQRSprod
real*8 Xaibj(LnPQRSprod)
integer LiPQRSprod(8)
integer iBatch, jBatch
integer nOccLeftI(8), nOccLeftJ(8)
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
character*20 ThisNm
character*25 SecNam
parameter(SecNam='Construct_WDensIII',ThisNm='Construct_WDensIII')
! Statement functions
MulD2h(i,j) = ieor(i-1,j-1)+1
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j
iDensAllAll(i,j,k) = ip_Density(k)+j-1+(nOrb(k)+nDel(k))*(i-1)
iWDensOccOcc(i,j,k) = ip_WDensity(k)+j-1+(nOrb(k)+nDel(k))*(i-1)

iSymPQ = 1
iSymIJ = 1
do iSymJ=1,nSym
  LjOrb = min(LnBatOrb(iSymJ,jBatch),nOccLeftJ(iSymJ))
  do Lj=1,LjOrb
    iJ = iFirstS(iSymJ,jBatch)+Lj-1
    do iSymQ=1,nSym
      do iQ=1,nOrb(iSymQ)+nDel(iSymQ)
        iSymP = MulD2h(iSymPQ,iSymQ)
        iSymQJ = MulD2h(iSymQ,iSymJ)
        iSymPI = iSymQJ
        iSymP = MulD2h(iSymPQ,iSymQ)
        Ljq = LiPQProd(iSymQ,iSymJ,jBatch)+(nFro(iSymQ)+nOcc(iSymQ)+nVir(iSymQ)+nDel(iSymQ))*(Lj-1)+iQ
        LpOrb = LnBatOrb(iSymP,iBatch)
        do Lp=1,LpOrb
          iP = iFirstS(iSymP,iBatch)+Lp-1
          iSymI = MulD2h(iSymIJ,iSymJ)
          Lpq = LiPQProd(iSymQ,iSymP,iBatch)+(nFro(iSymQ)+nOcc(iSymQ)+nVir(iSymQ)+nDel(iSymQ))*(Lp-1)+iQ
          do iI=1,nFro(iSymI)+nOcc(iSymI)
            Lpi = LiPQProd(iSymI,iSymP,iBatch)+(nFro(iSymI)+nOcc(iSymI)+nVir(iSymI)+nDel(iSymI))*(Lp-1)+iI
            Lji = LiPQProd(iSymI,iSymJ,jBatch)+(nFro(iSymI)+nOcc(iSymI)+nVir(iSymI)+nDel(iSymI))*(Lj-1)+iI

            if (iBatch == jBatch) then
              ip_pqij = LiPQRSprod(iSymPQ)+iTri(Lji,Lpq)
              ip_piqj = LiPQRSprod(iSymPI)+iTri(Ljq,Lpi)
            else
              ip_pqij = LiPQRSprod(iSymPQ)+LnPQprod(iSymPQ,iBatch)*(Lji-1)+Lpq
              ip_piqj = LiPQRSprod(iSymPI)+LnPQprod(iSymPI,iBatch)*(Ljq-1)+Lpi
            end if

            X = 2.0d0*Xaibj(ip_pqij)-Xaibj(ip_piqj)
            Work(iWDensOccOcc(iI,iJ,iSymJ)) = Work(iWDensOccOcc(iI,iJ,iSymJ))-X*Work(iDensAllAll(iP,iQ,iSymQ))
            ! Debug Comments -------------------------------------------
            !write(6,*) 'IJPQ',iI,iJ,iP,iQ
            !write(6,*) 'Symm',iSymI,iSymJ,iSymP,iSymQ
            !write(6,*) 'pqij',Xaibj(ip_pqij)
            !write(6,*) 'piqj',Xaibj(ip_piqj)
            !write(6,*) 'Dens',Work(iDensAllAll(iP,iQ,iSymQ))
            !------------------------------------------------------------

          end do
        end do
      end do
    end do
  end do
end do

! Avoid unused argument warnings
if (.false.) call Unused_integer_array(nOccLeftI)

end subroutine Construct_WDensIII
