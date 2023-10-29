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
! @parameter kOp   operators for center generator                      *
!                                                                      *
! Object: trace the gradient of the ERI's with the second order        *
!         density matrix                                               *
!                                                                      *
!     Author: Anders Bernhardsson Dept. of Theoretical Chemistry,      *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl, iChBas
use Index_Functions, only: iTri
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nHess, IndGrd(3,4,0:(nIrrep-1)), IndHss(4,3,4,3,0:(nIrrep-1)), iuvwx(4), kOp(4), nOp(4)
real(kind=wp), intent(inout) :: g2(78), Hess(nHess)
logical(kind=iwp), intent(in) :: IfHss(4,3,4,3), Tr(4), IfGr(4)
integer(kind=iwp) :: iCa1, iCa2, iCar, iCn, iCn1, iCn2, iHess, iIrrep, ij, iStop, jCar, jCn, k1, k2, kCn, kl, lCn
real(kind=wp) :: Fact, ps
real(kind=wp), parameter :: Prmt(0:7) = [One,-One,-One,One,-One,One,One,-One]

!                                                                      *
!***********************************************************************
!                                                                      *
!iRout = 239
!iPrint = nPrint(iRout)
!#define _DEBUGPRINT_
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
          ij = iTri((iCn-1)*3+iCar,(jCn-1)*3+jCar)
          g2(ij) = Zero
          !------------------------------------------------------------*
          !
          ! Both derivatives by translation!
          !
          !------------------------------------------------------------*
          if (Tr(iCn) .and. Tr(jCn)) then
            do kCn=1,4
              do lCn=1,kCn
                if (lCn == kCn) then
                  !iMax = Max(iCar,jCar)
                  iCa2 = min(iCar,jCar)
                  iCa1 = max(iCar,jCar)
                  if (IfHss(kCn,iCa1,lCn,iCa2)) then
                    k1 = iTri((kCn-1)*3+iCa1,(lCn-1)*3+iCa2)
                    g2(ij) = g2(ij)+g2(k1)
                  end if
                else
                  if (IfHss(kCn,iCar,lCn,jCar)) then
                    k1 = iTri((kCn-1)*3+iCar,(lCn-1)*3+jCar)
                    k2 = iTri((kCn-1)*3+jCar,(lCn-1)*3+iCar)
                    g2(ij) = g2(ij)+g2(k1)+g2(k2)
                  end if
                end if
              end do
            end do
          else if (IfGr(iCn) .and. Tr(jCn)) then
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
                kl = iTri((iCn1-1)*3+iCa1,(iCn2-1)*3+iCa2)
                g2(ij) = g2(ij)-g2(kl)
              end if
            end do
          else if (IfGr(jCn) .and. Tr(iCn)) then
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
                kl = iTri((iCn1-1)*3+iCa1,(iCn2-1)*3+iCa2)
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
            ij = iTri((iCn-1)*3+iCar,(jCn-1)*3+jCar)
            iHess = abs(IndHss(iCn,iCar,jCn,jCar,iIrrep))
            !----------------------------------------------------------*
            !
            ! Sign due to integral direction
            !
            !----------------------------------------------------------*
            ps = real(iChTbl(iIrrep,nOp(iCn))*iChTbl(iIrrep,nOp(jCn)),kind=wp)
            !----------------------------------------------------------*
            !
            ! If over & under triangular integrals are needed
            ! multiply by two instead!
            !
            !----------------------------------------------------------*
            if ((iCn /= jCn) .and. (iCar == jCar) .and. (abs(IndGrd(iCar,iCn,iIrrep)) == abs(IndGrd(jCar,jCn,iIrrep)))) then
              ps = ps*Two
            end if
            !----------------------------------------------------------*
            !
            ! Sign due to which symmetry group the translation is in
            !
            !----------------------------------------------------------*
            ps = ps*Prmt(iand(kOp(iCn),iChBas(1+iCar)))*Prmt(iand(kOp(jCn),iChBas(1+jCar)))
            !----------------------------------------------------------*
            !
            ! Multiply by number of stabilisers.
            !
            !----------------------------------------------------------*
            Fact = ps*real(iuvwx(iCn),kind=wp)/real(nIrrep*nirrep,kind=wp)*real(iuvwx(jCn),kind=wp)
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
