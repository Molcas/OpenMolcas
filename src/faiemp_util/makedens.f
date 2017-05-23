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
* Copyright (C) Ben Swerts                                             *
*               2016, Liviu Ungur                                      *
************************************************************************
      SubRoutine MakeDens(nBas,nOrb,Cff,OrbEn,EnergyWeight,Dens)
************************************************************************
*                                                                      *
*     purpose: Compute (energy weighted) density matrix in AO basis    *
*                                                                      *
*     input:                                                           *
*       nBas    : number of basis functions                            *
*       nOrb    : number of orbitals                                   *
*       Cff     : molecular orbital coefficients                       *
*       OrbEn   : molecular orbital energies                           *
*       EnergyWeight : if .true. do energy weighting                   *
*                                                                      *
*     output:                                                          *
*       Dens    : density matrix in triangular storage                 *
*                                                                      *
*     called from: Drv2El_FAIEMP                                       *
*                  FragPInt                                            *
*                                                                      *
*     written by: B. Swerts                                            *
*     modified by L. Ungur                                             *
*     simplified version of scf/done_scf.f                             *
*                                                                      *
************************************************************************
*
      Implicit  None
      Integer   nBas,nOrb
      Real*8    Cff(nBas*nOrb)
      Real*8    OrbEn(nOrb)
      Real*8    Dens(nBas*(nBas+1)/2)
      Logical   EnergyWeight
#include "real.fh"
      Integer i,j,Ind
      Integer ij,iRow,iCol
      Real*8  energy,Sum
      Logical DBG
*
*---- Statement function for triangular storage
      Ind(i,j) = i*(i - 1)/2 + j
      energy = One
*
      DBG=.false.
      if(DBG) write(6,*) 'MakeDens:  EnergyWeight',EnergyWeight
      if(DBG) call xFlush(6)
      if(DBG) write(6,*) 'MakeDens:  nBas=',nBas
      if(DBG) call xFlush(6)
      if(DBG) write(6,*) 'MakeDens:  nOrb=',nOrb
      if(DBG) call xFlush(6)
      if(DBG) write(6,*) 'MakeDens: OrbEn=',(OrbEn(i),i=1,nOrb)
      if(DBG) call xFlush(6)
      if(DBG) write(6,*) 'MakeDens:   Cff=',(Cff(i),i=1,nBas*nOrb)
      if(DBG) call xFlush(6)


      Do iRow = 1, nBas
        Sum = Zero
        ij  = -1
        Do i = 1, nOrb
          ij  = ij  + 1
          If(EnergyWeight) energy = OrbEn(i)
          Sum = Sum + energy * Cff(iRow + ij*nBas)
     &                       * Cff(iRow + ij*nBas)
        End Do
        Dens(Ind(iRow,iRow)) = Two*Sum
*
        Do iCol = 1, iRow - 1
          Sum = Zero
          ij  = -1
          Do i = 1, nOrb
             ij = ij  + 1
             If(EnergyWeight) energy = OrbEn(i)
             Sum = Sum + energy * Cff(iRow + ij*nBas)
     &                          * Cff(iCol + ij*nBas)
          End Do
          Dens(Ind(iRow,iCol)) = Two*Two*Sum
        End Do
      End Do
      if(DBG) call TriPrt('Dens in MakeDens',' ',Dens,nBas)
      if(DBG) call xFlush(6)
* Title,FmtIn,A,N
      Return
      End
