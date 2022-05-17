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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine Distg2(g2,Hess,nHess,IndGrd,IfHss,IndHss,iuvwx,kOp,nop,Tr,IfGr)
!***********************************************************************
!                                                                      *
! @parameter kop   operators for center generator                      *
!                                                                      *
! Object: trace the gradient of the ERI's with the second order        *
!         density matrix                                               *
!                                                                      *
!     Author: Anders Bernhardsson Dept. of Theoretical Chemistry,      *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl, iChBas

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 g2(78), Prmt(0:7), Hess(nHess)
logical IfHss(4,3,4,3), Tr(4), IfGr(4)
integer IndGrd(3,4,0:(nIrrep-1)), kOp(4), iuvwx(4), IndHss(4,3,4,3,0:(nIrrep-1)), nop(4)
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
! Statement Functions
xPrmt(i,j) = Prmt(iand(i,j))
ix(icn,icar,jcn,jcar) = (((icn-1)*3+icar)*((icn-1)*3+icar-1))/2+(jcn-1)*3+jcar

!                                                                      *
!***********************************************************************
!                                                                      *
!iRout = 239
!iPrint = nPrint(iRout)
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call recprt('Distg2: g2(raw) ',' ',g2,1,78)
call recprt('Distg2: Hess(raw) ',' ',Hess,1,nHess)
#endif

! Compute some of the contributions via the translational invariance

do iCn=1,4
  do iCar=1,3
    do jCn=1,iCn
      if (iCn == jCn) then
        iStop = iCar
      else
        iStop = 3
      end if
      do jCar=1,iStop
        if (Tr(iCn) .or. Tr(jCn)) then
          ij = ix(iCn,iCar,jCn,jCar)
          g2(ij) = zero
          !------------------------------------------------------------*
          !
          ! Both derivatives by translation!
          !
          !------------------------------------------------------------*
          if (tr(iCn) .and. tr(jCn)) then
            do kCn=1,4
              do lCn=1,kCn
                if (lCn == kCn) then
                  !iMax = Max(iCar,jCar)
                  iCa2 = min(iCar,jCar)
                  iCa1 = max(iCar,jCar)
                  if (IfHss(kCn,iCa1,lCn,iCa2)) then
                    k1 = Ix(kCn,iCa1,lCn,iCa2)
                    g2(ij) = g2(ij)+g2(k1)
                  end if
                else
                  if (IfHss(kCn,iCar,lCn,jCar)) then
                    k1 = Ix(kCn,iCar,lCn,jCar)
                    k2 = Ix(kCn,jCar,lCn,iCar)
                    g2(ij) = g2(ij)+g2(k1)+g2(k2)
                  end if
                end if
              end do
            end do
          else if (ifgr(iCn) .and. tr(jCn)) then
            !----------------------------------------------------------*
            !
            ! Centre jCn by translation
            !
            !----------------------------------------------------------*
            do kCn=1,4
              if (kCn > iCn) then
                iCn1 = kCn
                iCn2 = iCn
                iCa1 = jCar
                iCa2 = iCar
              else if (kCn < iCn) then
                iCn1 = iCn
                iCn2 = kCn
                iCa1 = iCar
                iCa2 = jCar
              else
                iCn1 = iCn
                iCn2 = kCn
                iCa1 = max(iCar,jCar)
                iCa2 = min(iCar,jCar)
              end if
              if (IfHss(iCn1,iCa1,iCn2,iCa2)) then
                kl = Ix(iCn1,iCa1,iCn2,iCa2)
                g2(ij) = g2(ij)-g2(kl)
              end if
            end do
          else if (IfGr(jCn) .and. tr(iCn)) then
            !----------------------------------------------------------*
            !
            ! Centre iCn by translation
            !
            !----------------------------------------------------------*
            do kCn=1,4
              if (kCn > jCn) then
                iCn1 = kCn
                iCn2 = jCn
                iCa1 = iCar
                iCa2 = jCar
              else if (kCn < jCn) then
                iCn1 = jCn
                iCn2 = kCn
                iCa1 = jCar
                iCa2 = iCar
              else
                iCn1 = jCn
                iCn2 = kCn
                iCa1 = max(iCar,jCar)
                iCa2 = min(iCar,jCar)
              end if
              if (IfHss(iCn1,iCa1,iCn2,iCa2)) then
                kl = Ix(iCn1,iCa1,iCn2,iCa2)
                g2(ij) = g2(ij)-g2(kl)
              end if
            end do
          end if
        end if
      end do
    end do
  end do
end do
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
call recprt('Distg2: g2 ',' ',g2,1,78)
#endif
! Distribute contribution to the hessian.

!----------------------------------------------------------------------*
do iIrrep=0,nIrrep-1
  do iCn=1,4
    do iCar=1,3
      do jCn=1,iCn
        if (iCn == jCn) then
          iStop = iCar
        else
          iStop = 3
        end if
        do jCar=1,istop
          if (IndHss(iCn,iCar,jCn,jCar,iIrrep) /= 0) then
            !----------------------------------------------------------*
            !
            ! Get indices
            !
            !----------------------------------------------------------*
            ij = Ix(iCn,iCar,jCn,jCar)
            iHess = abs(IndHss(iCn,iCar,jCn,jCar,iIrrep))
            !----------------------------------------------------------*
            !
            ! Sign due to integral direction
            !
            !----------------------------------------------------------*
            ps = dble(iChTbl(iIrrep,nOp(iCn))*iChTbl(iIrrep,nOp(jCn)))
            !----------------------------------------------------------*
            !
            ! If over & under triangular integrals are needed
            ! multiply by two instead!
            !
            !----------------------------------------------------------*
            if ((iCn /= jCn) .and. (iCar == jCar) .and. (abs(indgrd(iCar,iCn,iIrrep)) == abs(indgrd(jCar,jCn,iIrrep)))) then
              ps = ps*Two
            end if
            !----------------------------------------------------------*
            !
            ! Sign due to which symmetry group the translation is in
            !
            !----------------------------------------------------------*
            iCh1 = iChBas(iCar+1)
            iCh2 = iChBas(jCar+1)
            ps = ps*xPrmt(kOp(iCn),iCh1)*xPrmt(kOp(jCn),iCh2)
            !----------------------------------------------------------*
            !
            ! Multiply by number of stabilisers.
            !
            !----------------------------------------------------------*
            Fact = ps*dble(iuvwx(iCn))/dble(nIrrep*nirrep)*dble(iuvwx(jCn))
            !----------------------------------------------------------*
            !
            ! Add to hessian
            !
            !----------------------------------------------------------*
            Hess(iHess) = Hess(iHess)+Fact*g2(ij)
          end if
        end do
      end do
    end do
  end do
end do
#ifdef _DEBUGPRINT_
call recprt('Distg2: Hess ',' ',Hess,1,nHess)
#endif

return

end subroutine Distg2
