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
! Copyright (C) 2013, Giovanni Li Manni                                *
!               2013, Dongxia Ma                                       *
!***********************************************************************

subroutine DmatDmat(Dmat,DDarray)
!***********************************************************************
! Purpose: To construct an array containing Dpq*Drs elements
!          (product of one-body density matrix elements) ordered
!          according to the ordering of the 2-electron integrals
!          elements g(pqrs).
!
! Author : Giovanni Li Manni and Dongxia Ma
! Date   : June 21st 2013 ... When Minnesotan summer seems to arrive!
!***********************************************************************

use Symmetry_Info, only: Mul
use Index_Functions, only: i_Tri => iTri, nTri_Elem
use rasscf_global, only: ISTORP
use general_data, only: NASH, NSYM
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: Dmat(*)
real(kind=wp), intent(_OUT_) :: DDarray(*)
integer(kind=iwp) :: indx1, indx2, indxDpq, indxDrs, iOffOrb(nSym), iOrbP, iOrbQ, iOrbR, iOrbS, iorp, iPsm, iQsm, iRsm, iSmPQ, &
                     iSsm, iSym, nRS
real(kind=wp) :: FACT

! Initialization of variables

iorp = 0
indx1 = 0
iOffOrb(1) = 0

do iSym=2,nSym
  iOffOrb(iSym) = iOffOrb(iSym-1)+NaSh(iSym-1)
end do

DDarray(1:istorp(nSym+1)) = Zero

do iPsm=1,nSym
  do iOrbP=1,NASh(iPsm)
    do iQsm=1,nSym
      if (NASh(iQsm) == 0) cycle  ! next iQsm
      iSmPQ = Mul(iPsm,iQsm)
      indx2 = 0
      do iRsm=1,nSym
        iSsm = Mul(iSmPQ,iRsm)
        if ((min(NASh(iRsm),NASh(iSsm)) /= 0) .and. (iSsm <= iRsm)) then
          nRS = NASh(iRsm)*NASh(iSsm)
          if (iSsm == iRsm) then
            nRS = nTri_Elem(NASh(iRsm))
          end if
          if ((iRsm /= iSsm) .or. (iQsm /= iPsm)) then
            iorp = iorp+nRS*NASh(iQsm)
          else
            do iOrbR=1,NASH(iRsm)
              do iOrbS=1,iOrbR
                if (iOrbR /= iOrbS) then
                  FACT = Two
                else
                  FACT = One
                end if
                do iOrbQ=1,NASH(iQsm)
                  iorp = iorp+1
                  indxDpq = indx1+i_Tri(iOrbP,iOrbQ)
                  indxDrs = indx2+i_Tri(iOrbR,iOrbS)
                  DDarray(iorp) = Dmat(indxDpq)*Dmat(indxDrs)*FACT
                end do ! Loop over iOrbQ
              end do ! Loop over iOrbS
            end do ! Loop over iOrbR
          end if
        end if
        indx2 = indx2+nTri_Elem(NASH(iRsm))
      end do ! loop over iRsm
    end do ! Loop over iQsm
  end do ! Loop over iOrbP
  indx1 = indx1+nTri_Elem(NASH(iPsm))
end do ! Loop over iPsm

end subroutine DmatDmat
