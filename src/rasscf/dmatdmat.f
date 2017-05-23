************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2013, Giovanni Li Manni                                *
*               2013, Dongxia Ma                                       *
************************************************************************
      SUBROUTINE DmatDmat(Dmat,DDarray)
c *******************************************************************
c     Purpose: To construct an array containing Dpq*Drs elements
c              (product of one-body density matrix elements) ordered
c              according to the ordering of the 2-electron integrals
c              elements g(pqrs).
C
c     Author : Giovanni Li Manni and Dongxia Ma
c     Date   : June 21st 2013 ... When Minnesotan summer seems to arrive!
c *******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      PARAMETER (ROUTINE='DMATDMAT')
      DIMENSION Dmat(*),DDarray(*), iOffOrb(nSym)

      nTest = 000
      IF(nTest.ge.1000) THEN
        write(6,*)('*', i=1,70)
        WRITE(6,*)' Entering DmatDmat routine'
        write(6,*)('*', i=1,70)
        write(6,*) ' One body Dmat in INPUT:'
        CALL TRIPRT('Aver. 1-body Dmat in MO in DmatDmat',' ',Dmat,NAC)
      END IF
c Initialization of variables

      iorp       = 0
      iOffOrb(1) = 0
      indx1      = 0

      Do iSym = 2, nSym
        iOffOrb(iSym) = iOffOrb(iSym-1) + NaSh(iSym-1)
      End do

      IF(nTest.ge.1000) THEN
        write(6,*)('*', i=1,70)
        WRITE(6,*)' Offsets for Orbitals: Sym, Offset value'
        write(6,*)('*', i=1,70)
        do iSym = 1, nSym
          write(6,*) iSym, iOffOrb(iSym)
        end do
      END IF

      call FZero(DDarray,istorp(nSym+1))

      IF(nTest.ge.1000) THEN
        write(6,*)('*', i=1,70)
        WRITE(6,*)'array of Dpq*Drs elements size :',istorp(nSym+1)
      END IF

      If(nTest.ge.1000) then
        write(6,*)('*', i=1,70)
        write(6,*) '  iorp  iPsm  iQsm  iRsm  iSsm'//
     &             '  iOrbP iOrbQ iOrbR iOrbS'//
     &             '            Dpq           Drs'
      End If
      DO iPsm=1, nSym
        Do iOrbP=1, NASh(iPsm)
          DO iQsm=1,nSym
            If(NASh(iQsm).eq.0) goto 10  ! next iQsm
            iSmPQ = ieor(iPsm-1, iQsm -1)+1
            indx2      = 0
            DO iRsm=1,nSym
              iSsm = ieor(iSmPQ-1,iRsm-1)+1
              If(min(NASh(iRsm),NASh(iSsm)).eq.0.OR.
     &           iSsm.gt.iRsm)  goto 20  !next iRsm
              nRS = NASh(iRsm)*NASh(iSsm)
              if(iSsm.eq.iRsm) then
                nRS=NASh(iRsm)*(NASh(iRsm)+1)/2
              end if
              If(iRsm.ne.iSsm.OR.iQsm.ne.iPsm) then
                iorp = iorp+nRS*NASh(iQsm)
                goto 20 !next iRsm
              End If
              Do iOrbR=1,NASH(iRsm)
                DO iOrbS=1,iOrbR
                  FACT= 1.0d0
                  IF(iOrbR.ne.iOrbS) FACT= 2.0d0
                  Do iOrbQ=1,NASH(iQsm)
                    iorp = iorp +1
                    iMaxPQ=max(iOrbP,iOrbQ)
                    iminPQ=min(iOrbP,iOrbQ)
                    iMaxRS=max(iOrbR,iOrbS)
                    iminRS=min(iOrbR,iOrbS)
                    indxDpq = indx1+iMaxPQ*(iMaxPQ-1)/2 + iminPQ
                    indxDrs = indx2+iMaxRS*(iMaxRS-1)/2 + iminRS
                    DDarray(iorp)=Dmat(indxDpq)*Dmat(indxDrs)*FACT
                    If(nTest.ge.1000) then
                      write(6,'(9I6,5X,2F15.8)') iorp,
     &                           iPsm,iQsm,iRsm,iSsm,
     &                           iOrbP+iOffOrb(iPsm),
     &                           iOrbQ+iOffOrb(iQsm),
     &                           iOrbR+iOffOrb(iRsm),
     &                           iOrbS+iOffOrb(iSsm),
     &                           Dmat(indxDpq),
     &                           Dmat(indxDrs)
                    End If
                  End Do ! Loop over iOrbQ
                End Do ! Loop over iOrbS
              End Do ! Loop over iOrbR
20            CONTINUE
              indx2=indx2+NASH(iRsm)*(NASH(iRsm)+1)/2
            End Do ! loop over iRsm
10          CONTINUE
          End Do ! Loop over iQsm
        End Do ! Loop over iOrbP
        indx1 = indx1 + NASH(iPsm)*(NASH(iPsm)+1)/2
      End Do ! Loop over iPsm
      IF(nTest.ge.1000) THEN
        write(6,*)('*', i=1,70)
        WRITE(6,*)'array of Dpq*Drs elements in output from DmatDmat'
        call wrtmat(DDarray,1,istorp(nSym+1),1,istorp(nSym+1))
        write(6,*)('*', i=1,70)
        write(6,*) 'End of DmatDmat routine'
        write(6,*)('*', i=1,70)
      END IF

      RETURN
      END
