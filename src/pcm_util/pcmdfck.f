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
      Subroutine PCMDFck(nFck,PCMFck)
      use PCM_arrays
      Implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "Molcas.fh"
#include "rctfld.fh"
#include "periodic_table.fh"
      Dimension PCMFck(nFck,*)
      Character*2 Elements(MxAtom*8)
      Logical DoPot, DoFld
*
************************************************************************
*
*     Driver for the computation of PCM contributions to Fock matrix
*     derivatives
*
C
C     Retrieve atomic info
      Call Get_nAtoms_All(nAtoms)
      Call Allocate_Work(ipCoor,3*nAtoms)
      Call Get_Coord_All(Work(ipCoor),nAtoms)
      Call Get_Name_All(Elements)
      Call GetMem('ANr','Allo','Inte',ipANr,nAtoms)
      Do i = 1, nAtoms
         Do j = 0, Num_Elem
            If (PTab(j).eq.Elements(i)) iWork(ipANr+i-1)=j
         End Do
      End Do
      Call GetMem('Chrg','Allo','Real',ipChrg,nAtoms)
      Call Get_dArray('Nuclear charge',Work(ipChrg),nAtoms)
C
C     Allocate space for total charges
      Call GetMem('Qtot','Allo','Real',ip_Qtot,nTs)
C
C     Allocate space for the PCM matrix derivative
      Call GetMem('DerMat','Allo','Real',ip_DerMat,nTs*nTs)
C
C     Allocate space for the potential on tesserae
      Call GetMem('V','Allo','Real',ip_V,nTs)
C
C     Allocate space for the uncontracted potential on tesserae
      Call GetMem('VMN','Allo','Real',ip_VMN,nTs*nFck)
C
C     Allocate space for the potential derivatives
      nAt3 = nAtoms * 3
      Call GetMem('VDer','Allo','Real',ip_VDer,nAt3*nTs)
C
C     Allocate space for the uncontracted derivatives of the
C     potential on tesserae
      Call GetMem('VDerMN','Allo','Real',ip_VDerMN,nAt3*nTs*nFck)
C
C     Allocate space for electric field
      nComp=3
      Call GetMem('EF_n','Allo','Real',ip_EF_n,nComp*nTs)
      Call GetMem('EF_e','Allo','Real',ip_EF_e,nComp*nTs)
C
C     Allocate two scratch vectors
      Call GetMem('Temp1','Allo','Real',ip_Temp1,nTs)
      Call GetMem('Temp2','Allo','Real',ip_Temp2,nTs)

C
      Call FZero(PCMFck,nAt3*nFck)
C
C     Compute the potential and the electric field on tesserae
      DoPot = .True.
      DoFld = .True.
      Call V_EF_PCM(nAtoms,nTs,DoPot,DoFld,Work(ipCoor),Work(ip_Tess),
     &     Work(ip_V),Work(ip_EF_n),Work(ip_EF_e))
C
C     Compute the derivatives of the total potential on tesserae
      Call VDer_PCM(nAtoms,nTs,nS,Work(ipCoor),Work(ipChrg),
     &     Work(ip_EF_n),Work(ip_EF_e),Work(ip_Tess),iWork(ip_ISph),
     &     dTes,dPnt,Work(ip_DRad),Work(ip_DCntr),
     &     Work(ip_VDer))
C
C     Actually compute the PCM correction
      Call PCM_Der_Fock(nFck,nAtoms,nTs,nS,Eps,Work(ip_Sph),
     &     iWork(ip_ISph),iWork(ip_N),Work(ip_Tess),Work(ip_Q),
     &     Work(ip_Qtot),Work(ip_DM),Work(ip_DerMat),
     &     dTes,dPnt,Work(ip_DCntr),
     &     Work(ip_V),Work(ip_VMN),Work(ip_VDer),Work(ip_VDerMN),
     &     Work(ip_Temp1),Work(ip_Temp2),PCMFck)
C
C     Free the space
      Call GetMem('Qtot','Free','Real',ip_Qtot,nTs)
      Call GetMem('V','Free','Real',ip_V,nTs)
      Call GetMem('VMN','Free','Real',ip_VMN,nTs*nFck)
      Call GetMem('VDer','Free','Real',ip_VDer,nAt3*nTs)
      Call GetMem('VDerMN','Free','Real',ip_VDerMN,nAt3*nTs*nFck)
      Call GetMem('EF_n','Free','Real',ip_EF_n,nComp*nTs)
      Call GetMem('EF_e','Free','Real',ip_EF_e,nComp*nTs)
      Call GetMem('Temp1','Free','Real',ip_Temp1,nTs)
      Call GetMem('Temp2','Free','Real',ip_Temp2,nTs)
      Return
      End
*
************************************************************************
*
      Subroutine PCM_Der_Fock(nFck,nAt,nTs,nS,Eps,Sphere,ISphe,nOrd,
     &           Tessera,Q,Qtot,DM,DerDM,DerTes,DerPunt,DerCentr,
     &           V,VMN,VDer,VDerMN,Temp1,Temp2,PCMFck)
      Implicit real*8 (a-h,o-z)
#include "real.fh"
      Dimension Sphere(4,*),ISphe(*),nOrd(*),Temp1(*),Temp2(*)
      Dimension Tessera(4,*),Q(2,*),Qtot(*),DM(nTs,*),DerDM(nTs,*)
      DImension V(*),VMN(nTs,*)
      Dimension VDer(nTs,*),VDerMN(nTs,nFck,*),PCMFck(nFck,*)
      Dimension DerTes(nTs,nAt,3),DerPunt(nTs,nAt,3,3)
      Dimension DerCentr(nS,nAt,3,3)
C
      FPI = Four*PI
      Diag = - 1.0694d0 * Sqrt(FPI) / Two
      Sc_Cond = (Eps - One) / Eps
C
      Do iTs = 1, nTs
        Qtot(iTs) = Q(1,iTs) + Q(2,iTs)
      EndDo
C
C     Loop over the degrees of freedom
      Do 100 iAt = 1, nAt
        Do 100 iC = 1, 3
          Index = 3 * (iAt-1) + iC
C
C         Derivative of the PCM matrix for the conductor-like case
          Call DMat_CPCM(iAt,iC,Eps,nTs,nS,nAt,Diag,Tessera,
     &                   DerDM,DerTes,DerPunt,DerCentr,iSphe)
C
C         Solvation charges (weights) times the derivative of the PCM matrix
          Call PrMatVec(.True.,.True.,DerDM,-1.d0,nTs,nTs,
     &         QTot,Temp1)
C
C         The previous vector times the inverted PCM matrix
          Call PrMatVec(.True.,.True.,DM,1.d0,nTs,nTs,
     &         Temp1,Temp2)
C
C         Derivative of the potential times the inverted PCM matrix
          Call PrMatVec(.True.,.True.,DM,1.d0,nTs,nTs,
     &         VDer(1,Index),Temp1)
C
C         Loop over the Fock elements
          Do 200 iFck = 1, nFck
C
C------
C           First contribution: charges times the derivative of the
C           (uncontracted) electronic potential: Sum_i q_i V_mn^(i,x)
            Sum = Zero
            Do iTs = 1, nTs
              Sum = Sum + Qtot(iTs) * VDerMN(iTs,iFck,Index)
            EndDo
            PCMFck(iFck,Index) = PCMFck(iFck,Index) + Sum
C------
C           Second contribution: derivative of the potential times the
C           (symmetrized) PCM matrix times the uncontracted potential:
C           Sum_ij V_i^x Q_ij V_mn^j
            Sum = Zero
            Do iTs = 1, nTs
              Sum = Sum + Temp1(iTs) * VMN(iTs,iFck)
            EndDo
            PCMFck(iFck,Index) = PCMFck(iFck,Index) + Sum
C------
C           Third contribution: potential times the derivative of
C           the inverted PCM matrix times the uncontracted potential:
C           Sum_ij V_i Q_ij^x V_mn^j = Sum_ij V_i [-S^-1 S^x S^-1]_ij V_mn^j =
C           Sum_ij q_i [S^x S^-1]_ij V_mn^j
            Sum = Zero
            Do iTs = 1, nTs
              Do jTs = 1, nTs
                Sum = Sum + Temp2(iTs) * VMN(jTs,iFck)
              EndDo
            EndDo
            PCMFck(iFck,Index) = PCMFck(iFck,Index) + Sum
  200     Continue
  100 Continue
C
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Sphere)
        Call Unused_integer_array(nOrd)
        Call Unused_real_array(V)
      End If
      End
