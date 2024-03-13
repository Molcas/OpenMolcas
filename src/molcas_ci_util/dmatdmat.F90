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
      SUBROUTINE DmatDmat(Dmat,DDarray)
! *******************************************************************
!     Purpose: To construct an array containing Dpq*Drs elements
!              (product of one-body density matrix elements) ordered
!              according to the ordering of the 2-electron integrals
!              elements g(pqrs).
!
!     Author : Giovanni Li Manni and Dongxia Ma
!     Date   : June 21st 2013 ... When Minnesotan summer seems to arrive!
! *******************************************************************

      IMPLICIT None
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
      Real*8 Dmat(*),DDarray(*)
      Integer iOffOrb(nSym)

      Integer :: iorp, iSym,iMaxPQ, iMaxRS, iMinPQ, iMinRS, indx1, indx2, indxdpq, indxdrs, iOrbP, &
                 iOrbQ, iOrbR, iOrbS, iPsm, iQsm, iRsm, iSmPQ, iSsm, nRS
      Real*8 FACT

! Initialization of variables

      iorp       = 0
      iOffOrb(1) = 0

      Do iSym = 2, nSym
        iOffOrb(iSym) = iOffOrb(iSym-1) + NaSh(iSym-1)
      End do

      call FZero(DDarray,istorp(nSym+1))

      indx1 = - NASH(1)*(NASH(1)+1)/2
      DO iPsm=1, nSym
        indx1 = indx1 + NASH(iPsm)*(NASH(iPsm)+1)/2
        Do iOrbP=1, NASh(iPsm)
          DO iQsm=1,nSym
            If(NASh(iQsm)==0) Cycle ! next iQsm
            iSmPQ = ieor(iPsm-1, iQsm -1)+1
            indx2=-NASH(1)*(NASH(1)+1)/2
            DO iRsm=1,nSym
              indx2=indx2+NASH(iRsm)*(NASH(iRsm)+1)/2
              iSsm = ieor(iSmPQ-1,iRsm-1)+1
              If(min(NASh(iRsm),NASh(iSsm))==0.OR.iSsm>iRsm) Cycle  !next iRsm
              nRS = NASh(iRsm)*NASh(iSsm)
              if(iSsm==iRsm) nRS=NASh(iRsm)*(NASh(iRsm)+1)/2
              If(iRsm.ne.iSsm.OR.iQsm.ne.iPsm) then
                iorp = iorp+nRS*NASh(iQsm)
                Cycle !next iRsm
              End If
              Do iOrbR=1,NASH(iRsm)
                DO iOrbS=1,iOrbR
                  FACT= 1.0d0
                  IF(iOrbR/=iOrbS) FACT= 2.0d0
                  Do iOrbQ=1,NASH(iQsm)
                    iorp = iorp +1
                    iMaxPQ=max(iOrbP,iOrbQ)
                    iminPQ=min(iOrbP,iOrbQ)
                    iMaxRS=max(iOrbR,iOrbS)
                    iminRS=min(iOrbR,iOrbS)
                    indxDpq = indx1+iMaxPQ*(iMaxPQ-1)/2 + iminPQ
                    indxDrs = indx2+iMaxRS*(iMaxRS-1)/2 + iminRS
                    DDarray(iorp)=Dmat(indxDpq)*Dmat(indxDrs)*FACT
                  End Do ! Loop over iOrbQ
                End Do ! Loop over iOrbS
              End Do ! Loop over iOrbR
            End Do ! loop over iRsm

          End Do ! Loop over iQsm
        End Do ! Loop over iOrbP
      End Do ! Loop over iPsm

      END
