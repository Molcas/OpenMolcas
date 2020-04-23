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
      Subroutine Read_Blocks(iTable,nBlocks,nBas,nIrrep,iOff,Buf,nBuf,
     &                       iSO2Shell,nSOs,Bin,nBin,nQuad,G_Toc,
     &                       iSO2cI,CutInt)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
#include "real.fh"
#include "aces_gamma.fh"
#include "mp2alaska.fh"
#include "pso.fh"
      Integer iTable(6,nBlocks), nBas(0:nIrrep-1),
     &        iOff(0:nIrrep-1), iSO2Shell(nSOs), iSO2cI(2,nSOs)
      Real*8 Buf(nBuf), Bin(2,nBin,nQuad), G_Toc(nQuad)
      Logical Triangular
*                                                                      *
************************************************************************
*                                                                      *
*---- Statement function for triangular indexation
*
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate table SO to contigous index
*
*     Write (*,*) 'nQuad=',nQuad
      Call SO2cI(iSO2cI,iSO2Shell,nSOs)
*                                                                      *
************************************************************************
*                                                                      *
*---- Initiate bins
*
      Do iQuad = 1, nQuad
         Bin(1,nBin,iQuad)=Zero  ! number of elements in bin
         Bin(2,nBin,iQuad)=-One  ! back chaining address
      End Do
      iDisk=0 ! initial disk address
      iWrite=1
*                                                                      *
************************************************************************
*                                                                      *
*---- Loop over symmetry blocks of gammas
*
*     Write (*,*) 'iSO2Shell',iSO2Shell
*     Write (*,*) 'iSO2cI(1)',(iSO2cI(1,i),i=1,nSOs)
*     Write (*,*) 'iSO2cI(2)',(iSO2cI(2,i),i=1,nSOs)
      iBlockAdr = 1
      Do iBlock = 1, nBlocks
         iType   =iTable(1,iBlock)
         iIrrep_A=iTable(2,iBlock)
         iIrrep_B=iTable(3,iBlock)
         iIrrep_C=iTable(4,iBlock)
         iIrrep_D=iTable(5,iBlock)
         IND     =iTable(6,iBlock)
*        Write (*,*)  'Irreps:',iIrrep_A, iIrrep_B, iIrrep_C, iIrrep_D
*
         nA=nBas(iIrrep_A)
         nB=nBas(iIrrep_B)
         nC=nBas(iIrrep_C)
         nD=nBas(iIrrep_D)
*        Write (*,*) 'nA,nB,nC,nD=',nA,nB,nC,nD
         If (iType.eq.1 .or. iType.eq.2) Then
            nAB = nA*(nA+1)/2
            nCD = nC*(nC+1)/2
            Triangular=.True.
         Else
            nAB=nA*nB
            nCD=nC*nD
            Triangular=.False.
         End If
         If (nAB*nCD.eq.0) Go To 999
*
         nAB_dist=Min(nAB,nBuf/nCD)
         iSO_A_r=1
         iSO_B_r=1
         Do iiAB = 1, nAB, nAB_dist
            iAB_s=iiAB
            iAB_e=Min(nAB,iiAB+nAB_dist-1)
*
* This call is to an ACES2 routine. This has to be replaced once this
* code is used by MOLCAS. RL 2007-10-18.
C           Call GetLst(Buf,iAB_s,nAB_dist,2,iType,IND)
************************************************************************
* Jonas B 2010. This code is used for reading nonseparable two-electron
*               density matrices when doing conventional mp2-gradients.
            If(Case_mp2) Then
               iSize = min(nAB_dist*nCD,(nAB-(iiAB-1))*nCD)
               iAdr = iBlockAdr + nCD*(iiAB-1)
               Call dDaFile(LuGam,2,Buf,iSize,iAdr)
            End If
*
************************************************************************
*           Call RecPrt('Aces_Gamma: from GetLst',' ',
*    &                  Buf,iAB_e-iAB_s+1,nCD)
*
            iBuf=0
            Do iAB = iAB_s, iAB_e
*
               iSO_A_a=iOff(iIrrep_A)+iSO_A_r
               iSO_B_a=iOff(iIrrep_B)+iSO_B_r
               iShell_A=iSO2Shell(iSO_A_a)
               iShell_B=iSO2Shell(iSO_B_a)
               iShell_AB=iTri(iShell_A,iShell_B)
*
               Index_A=iSO2cI(1,iSO_A_a)
               Index_B=iSO2cI(1,iSO_B_a)
               nDim_A=iSO2cI(2,iSO_A_a)
               nDim_B=iSO2cI(2,iSO_B_a)
               nDim_AB=nDim_A*nDim_B
               If (iShell_A.gt.iShell_B) Then
                  Index_AB=(Index_B-1)*nDim_A+Index_A
               Else If (iShell_A.eq.iShell_B) Then
                  Index_AB=iTri(Index_A,Index_B)
               Else
                  Index_AB=(Index_A-1)*nDim_B+Index_B
               End If
*
               iSO_C_r=1
               iSO_D_r=1
               Do iCD = 1, nCD
*
                  iBuf = iBuf + 1
*                 iBuf = (iCD-1)*nAB_dist+iAB
                  ABCD=Buf(iBuf)
*
*     Jonas B 2010. If calculating mp2-integrals the cutoff for small
*     gammas should not be used at the moment.
                  If(.not. Case_mp2) Then
                     If (Abs(ABCD).lt.CutInt) Go To 888
                  End If
                  iSO_C_a=iOff(iIrrep_C)+iSO_C_r
                  iSO_D_a=iOff(iIrrep_D)+iSO_D_r
                  iShell_C=iSO2Shell(iSO_C_a)
                  iShell_D=iSO2Shell(iSO_D_a)
                  iShell_CD=iTri(iShell_C,iShell_D)
                  iShell_ABCD=iTri(iShell_AB,iShell_CD)
*
*---------------- Compute canonical compound index
*
                  Index_C=iSO2cI(1,iSO_C_a)
                  Index_D=iSO2cI(1,iSO_D_a)
                  nDim_C=iSO2cI(2,iSO_C_a)
                  nDim_D=iSO2cI(2,iSO_D_a)
                  nDim_CD=nDim_C*nDim_D
                  If (iShell_C.gt.iShell_D) Then
                     Index_CD=(Index_D-1)*nDim_C+Index_C
                  Else If(iShell_C.eq.iShell_D) Then
                     Index_CD=iTri(Index_C,Index_D)
                  Else
                     Index_CD=(Index_C-1)*nDim_D+Index_D
                  End If
                  If (iShell_AB.gt.iShell_CD) Then
                     Index_ABCD=(Index_CD-1)*nDim_AB+Index_AB
                  Else If (iShell_AB.eq.iShell_CD) Then
                     Index_ABCD=iTri(Index_AB,Index_CD)
                  Else
                     Index_ABCD=(Index_AB-1)*nDim_CD+Index_CD
                  End If
*                                                                      *
************************************************************************
*                                                                      *
*---------------- Store information in the appropriate bin
*
*                 Write (*,*) 'ABCD, Index_ABCD=',ABCD, Index_ABCD
*                 Write (*,*) 'iSO:',iSO_A_a,iSO_B_a,iSO_C_a,iSO_D_a
*                 Write (*,*) 'iShell:',iShell_A,iShell_B,
*    &                                  iShell_C,iShell_D
*                 Write (*,*) 'cind:',Index_A,Index_B,Index_C,Index_D
*                 Write (*,*) 'ndims:',nDim_A,nDim_B,nDim_C,nDim_D
*                 Write (*,*) 'iShell_AB,iShell_CD=',iShell_AB,iShell_CD
*                 Write (*,*) 'nDim_AB,nDim_CD=',nDim_AB,nDim_CD
*                 Write (*,*) 'iShell_ABCD=',iShell_ABCD
                  kBin=Int(Bin(1,nBin,iShell_ABCD))+1
                  Bin(1,kBin,iShell_ABCD)=ABCD
                  Bin(2,kBin,iShell_ABCD)=DBLE(Index_ABCD)
                  Bin(1,nBin,iShell_ABCD)=DBLE(kBin)  ! Update counter
*                                                                      *
************************************************************************
*                                                                      *
*---------------- Write bin to disk if full.
*
                  If (kBin.eq.nBin-1) Then
                     BackChain=DBLE(iDisk)
                     Call dDaFile(LuGamma,iWrite,Bin(1,1,iShell_ABCD),
     &                            2*nBin,iDisk)
                     Bin(1,nBin,iShell_ABCD)=Zero
                     Bin(2,nBin,iShell_ABCD)=BackChain
                  End If
*                                                                      *
************************************************************************
*                                                                      *
*---------------- Increment indices for C and D
*
                  If (.Not.Triangular) Then
*                    iSO_C_r=iSO_C_r+1
*                    If (iSO_C_r.gt.nBas(iIrrep_C)) Then
*                       iSO_C_r=1
*                       iSO_D_r=iSO_D_r+1
*                    End If
                     iSO_D_r=iSO_D_r+1
                     If (iSO_D_r.gt.nBas(iIrrep_D)) Then
                        iSO_D_r=1
                        iSO_C_r=iSO_C_r+1
                     End If
                  Else
*Lower triangular
                     iSO_D_r=iSO_D_r+1
                     If (iSO_D_r.gt.iSO_C_r) Then
                        iSO_C_r=iSO_C_r+1
                        iSO_D_r=1
                     End If
*Upper triangular
*                    iSO_C_r=iSO_C_r+1
*                    If (iSO_C_r.gt.nBas(iIrrep_C)) Then
*                       iSO_D_r=iSO_D_r+1
*                       iSO_C_r=iSO_D_r
*                    End If
                  End If
*
 888              Continue
               End Do  ! iCD
*                                                                      *
************************************************************************
*                                                                      *
*------------- Increment indices for A and B
*
               If (.Not.Triangular) Then
*                 iSO_A_r=iSO_A_r+1
*                 If (iSO_A_r.gt.nBas(iIrrep_A)) Then
*                    iSO_A_r=1
*                    iSO_B_r=iSO_B_r+1
*                 End If
                  iSO_B_r=iSO_B_r+1
                  If (iSO_B_r.gt.nBas(iIrrep_B)) Then
                     iSO_B_r=1
                     iSO_A_r=iSO_A_r+1
                  End If
               Else
*Lower triangular
                  iSO_B_r=iSO_B_r+1
                  If (iSO_B_r.gt.iSO_A_r) Then
                     iSO_A_r=iSO_A_r+1
                     iSO_B_r=1
                  End If
*Upper triangular
*                 iSO_A_r=iSO_A_r+1
*                 If (iSO_A_r.gt.nBas(iIrrep_A)) Then
*                    iSO_B_r=iSO_B_r+1
*                    iSO_A_r=iSO_B_r
*                 End If
               End If
*
            End Do     ! iAB
*
         End Do        ! iiAB
 999     Continue
         iBlockAdr = iBlockAdr + nAB*nCD
*
      End Do           ! iBlock
*                                                                      *
************************************************************************
*                                                                      *
*---- Flush the bins to disk and save disk addresses
*
*     Write (*,*) 'Flush bins!'
      Do iQuad = 1, nQuad
         BackChain=DBLE(iDisk)
         Call dDaFile(LuGamma,iWrite,Bin(1,1,iQuad),2*nBin,iDisk)
         lxx=Int(Bin(1,nBin,iQuad))
*        Write (*,*) 'lxx=',lxx
*        Call RecPrt('Bins',' ',Bin(1,1,iQuad),2,lxx)
         G_Toc(iQuad)=BackChain
      End Do
*
*     Observe that backchain list is only stored in core!
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
