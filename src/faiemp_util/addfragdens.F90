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
! Copyright (C) Ben Swerts                                             *
!***********************************************************************
      SubRoutine AddFragDens(Array, nDens, nDens_Valence, nBas_Valence)
!***********************************************************************
!                                                                      *
! Input: Array(Size) filled with valence density at proper positions   *
!        for a density including fragments.                            *
! Output: Updated with fragment densities at their proper positions.   *
!                                                                      *
!***********************************************************************
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iOper
      Implicit None
#include "real.fh"
#include "WrkSpc.fh"
      Integer nDens,nDens_Valence
      Real*8  Array(nDens)
      Integer nBas_Valence(0:7)
      Logical EnergyWeight
      Integer iPrint,maxDens,iCnttp,ipFragDensAO,iDpos,ipFragDensSO
      Integer i,j,iCnt,iFpos,iFD,mdc,iIrrep,nBasC
      Real*8  rDummy(1)

      If(nIrrep.ne.1) Then
        write(6,*) 'AddFragDens: Symmetry not implemented yet'
        Call Abend()
      End If
!
      iPrint=0
!
! Each fragment needs it''s (symmetrized) density matrix added along the diagonal
! This density matrix first has to be constructed from the MO coefficients
! so allocate space for the largest possible density matrix
      maxDens = 0
      Do iCnttp = 1, nCnttp
        If (dbsc(iCnttp)%nFragType.gt.0) maxDens = Max(maxDens,         &
     &                        dbsc(iCnttp)%nFragDens                    &
     &                      *(dbsc(iCnttp)%nFragDens+1)/2)
      End Do
      Call GetMem('FragDSO','Allo','Real',ipFragDensSO,maxDens)
!     If(nIrrep.ne.1) Then
!       Call GetMem('FragDAO','Allo','Real',ipFragDensAO,maxDens)
!     Else
        ipFragDensAO = ipFragDensSO
!     End If

      iDpos = 1 ! position in the total density matrix
      Do iIrrep = 0, nIrrep - 1
        nBasC = nBas_Valence(iIrrep)
        iDpos = iDpos + nBasC*(nBasC+1)/2
        mdc = 0
        Do 1000 iCnttp = 1, nCnttp
          If(dbsc(iCnttp)%nFragType.le.0) Then
            mdc = mdc + dbsc(iCnttp)%nCntr
            Go To 1000
          End If

! construct the density matrix
          EnergyWeight = .false.
          Call MakeDens(dbsc(iCnttp)%nFragDens,                         &
     &                  dbsc(iCnttp)%nFragEner,                         &
     &                  dbsc(iCnttp)%FragCoef,                          &
     &                  rDummy,                                         &
     &                  EnergyWeight,Work(ipFragDensAO))
! create the symmetry adapted version if necessary
! (fragment densities are always calculated without symmetry)
!         If(nIrrep.ne.1) Call SymmDens(Work(ipFragDensAO),
!    &      Work(ipFragDensSO))
          If(iPrint.ge.99) Call TriPrt('Fragment density',' ',          &
     &      Work(ipFragDensSO),dbsc(iCnttp)%nFragDens)
          Do iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
! only add fragment densities that are active in this irrep
! => the following procedure still has to be verified thoroughly
!    but appears to be working
            If(iAnd(dc(mdc)%iChCnt,iIrrep).eq.iOper(iIrrep)) Then
! add it at the correct location in the large custom density matrix
              iFpos = 1
!              ! position in fragment density matrix
              Do i = 1, dbsc(iCnttp)%nFragDens
                iDpos = iDpos + nBasC
                Do j = 0, i-1
                  Array(iDpos + j) = Work(ipFragDensSO + iFpos + j - 1)
                End Do
                iDpos = iDpos + i
                iFpos = iFpos + i
              End Do
              nBasC = nBasC + dbsc(iCnttp)%nFragDens
            End If
          End Do
 1000   Continue
      End Do
      If(iPrint.ge.19) Then
        iFD = 1
        Do iIrrep = 0, nIrrep - 1
          Call TriPrt('Combined density',' ',Array(iFD),nBas(iIrrep))
          iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
        End Do
      End If
      Call GetMem('FragDSO','Free','Real',ipFragDensSO,maxDens)
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(nDens_Valence)
      End
