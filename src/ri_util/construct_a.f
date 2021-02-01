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
* Copyright (C) Jonas Bostrom                                          *
************************************************************************
      subroutine Construct_A(nFIorb,nIrrep)
************************************************************************
*     Author: Jonas Bostr{\"o}m
*
*     Purpose: Computation of the matrix A to contract with the 2-center
*              derivative in the Exchange part of the DF-gradients
*
*     Equation: A_KL = sum_ij = C^K_ij*C^L_ij
*
*     Input:
*            nFIOrb: Array of # of Inactive(Occupied) + frozen orbitals
*                    in each irrep.
*            nIrrep: Number of irreducible representations.
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer  nVec(8)
      Integer nFIorb(8)
      Character Fname*6
      Character*13 SECNAM
      Parameter (SECNAM = 'CONSTRUCT_A')

#include "cholesky.fh"
#include "WrkSpc.fh"
*
*************************
*     Define some indeces
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
*************************
      Do i = 1, nIrrep
         nVec(i) = NumCho(i)
      End Do
      Call GAIGOP(nVec,nIrrep,'+')
*

*
      Do iChSym = 1, nIrrep
         Nij = 0
         Do iSym = 1, nIrrep
            jSym = MulD2h(iChSym,iSym)
            If(iChSym.eq.1) Then
               Nij=Nij + nFIorb(iSym)*(nFIorb(iSym)+1)/2
            Else
               Nij=Nij + nFIorb(iSym)*nFIorb(jSym)
            End If
         End Do
         If(Nij.eq.0) Go To 1000
         Call GetMem('MemChk','Max','Real',iDum,MemMax)
         MemMax = 60
*     Setup batches in a way so that we maximize the number of full rows in A_KL we
*     can have in memory at one time and if zero row fits we maximize the number of
*     columns of one


         rNij = Real(Nij)
         rMemMax = Real(MemMax)
         rChoFull=Real(nVec(iChSym))
*--------------------------------------
*         Write(6,*) 'rNij', rNij
*         Write(6,*) 'rMemMax', rMemMax
*         Write(6,*) 'rChoFull', rChoFull
*---------------------------------------
         If(rChoFull*(rChoFull+1)/2+rNij*rChoFull.lt.rMemMax) Then
            nBatL = 1
            nBatK = 1
            nBatVecK = nVec(iChSym)
            nBatVecL = nVec(iChSym)
         Else
            rChoMaxK = (rMemMax - rNij*rChoFull)/rChoFull
            iChoMaxK = Int(rChoMaxK)
            nBatVecK = Min(iChoMaxK,nVec(iChSym))
*-------------------------------------------------
*            Write(6,*) 'rChoMaxK', rChoMaxK
*            Write(6,*) 'iChoMaxK', iChoMaxK
*            Write(6,*) 'nBatVecK', nBatVecK
*-------------------------------------------------
            If(nBatVecK.ne.0) Then
               nBatK = nVec(iChSym)/nBatVecK
               nFullK = nBatVecK*nBatK
               If(nFullK.lt.nVec(iChSym)) nBatK=nBatK+1
               nBatL = 1
               nBatVecL = nVec(iChSym)
            Else
               nBatK=0
            End If

*
            If(nBatK.le.0) Then
               nBatVecK = 1
               nBatK = nVec(iChSym)
               rChoMaxL = rMemMax/(1+2*Nij)
               iChoMaxL = Int(rChoMaxL)
               nBatVecL = min(iChoMaxL,nVec(iChSym))
*-------------------------------------------------
*               Write(6,*) 'rChoMaxL', rChoMaxL
*               Write(6,*) 'iChoMaxL', iChoMaxL
*               Write(6,*) 'nBatVecL', nBatVecL
*-------------------------------------------------
               nBatL = nVec(iChSym)/nBatVecL
               nFullL = nBatVecL*nBatL
               If(nFullL.lt.nVec(iChSym)) nBatL=nBatL+1
               If(nBatL.lt.0) Then
                  Call WarningMessage(2,'Error in Construct_A')
                  Write(6,*) 'Too little memory available'
                  Write(6,*) 'Consider increasing MOLCAS_MEM'
                  Call Abend
               End If
            Else
               nBatL = 1
            End If
         End If

*         Write(6,*) 'nBatK, nBatL', nBatK, nBatL
*         Write(6,*) 'nBatVecK,nBatVecL',nBatVecK,nBatVecL
         iSeed = 8
         LuCVec = IsFreeUnit(iSeed)
         Write(Fname,'(A4,I1)') 'CVEC',iChSym
         Call DaName_MF_WA(LuCVec,Fname)
*
         iSeed = 9
         LuAMat = IsFreeUnit(iSeed)
         Write(Fname,'(A4,I1)') 'AMAT',iChSym
         Call DaName_MF_WA(LuAMat,Fname)
*
         l_Cvec1 = nBatVecL*Nij
         l_Cvec2 = nBatVecK*Nij
         If((nBatL.eq.1).and.(nBatK.eq.1)) Then
            l_Amat = nBatVecK*(nBatVecK+1)/2
         Else
            l_Amat = nBatVecK*nBatVecL
         End If
*         Write(6,*) 'l_Amat', l_Amat
         Call GetMem('C_vec1','Allo','Real',ip_Cvec1,l_Cvec1)
         Call GetMem('C_vec2','Allo','Real',ip_Cvec2,l_Cvec2)
         Call GetMem('A_matrix','Allo','Real',ip_Amat,l_Amat)
         Call FZero(Work(ip_Amat),l_Amat)
*
         nVecK = nBatVecK
         Do iBatK = 1, nBatK
*            Write(6,*) 'iBatK', iBatK
            If(iBatK.eq.nBatK) nVecK = nVec(iChSym) - (nBatK-1)*nBatVecK

*            iAdr = (iBatK-1)*Nij*nBatVecK
*            Write(6,*) 'iAdr1', iAdr
*            Call dDaFile(LuCVec,2,Work(ip_Cvec1),nVecK*Nij,iAdr)
            nVecL = nBatVecL
            iBatL_Max = iBatK/nBatVecL
            iRest = mod(iBatK,nBatVecL)
            if(iRest.ne.0) iBatL_Max = iBatL_Max+1
            Do iBatL = 1, iBatL_Max
*               Write(6,*) 'iBatL', iBatL
               If(iBatL.eq.nBatL) nVecL =nVec(iChSym)-(nBatL-1)*nBatVecL
               iAdr = (iBatL-1)*Nij*nVecL
*               Write(6,*) 'iadr', iAdr
               Call dDaFile(LuCVec,2,Work(ip_Cvec1),nVecL*Nij,iAdr)
               If(nBatL.ne.1) Then
                  iAdr = (iBatK-1)*Nij*nVecK
                  Call dDaFile(LuCVec,2,Work(ip_Cvec2),nVecK*Nij,iAdr)
               End If

*------------------------------------------------
*               Write(6,*) 'C-Vector1'
*               Do i = 0, nVecL*Nij-1
*                  Write(6,*) Work(ip_Cvec1+i)
*               End Do
*               If(iBatL.ne.0) Then
*                  Write(6,*) 'C-Vector2'
*                  Do i = 0, nVecK*Nij-1
*                     Write(6,*) Work(ip_Cvec2+i)
*                  End Do
*               End If
*------------------------------------------------


               Do iK = 1, nVecK
*                  If(iBatK.eq.iBatL) nL = iK
                  iK_Real = iK + nVecK*(iBatK-1)
*                 The L-vectors already done in an earlier batch
                  iL_passed = (iBatL-1)*nBatVecL
                  iL_max = min(iK_Real - iL_passed,nBatVecL)
                  Do iL = 1, iL_max
                     A_KL = 0
*                     Write(6,*) 'iL,iK', iL,iK
                     Do ij = 1, Nij
                        If(nBatL.eq.1) Then
                           index_K = ij-1 + (iK_Real-1)*Nij
                        Else
                           index_K = ij-1 + (iK-1)*Nij
                        End If
                        index_L = ij-1 + (iL-1)*Nij
*                        Write(6,*) 'index_K', index_K
*                        Write(6,*) 'index_L', index_L
                        If(nBatL.eq.1) Then
                           A_KL = A_KL + Work(ip_Cvec1+index_L)*
     &                                   Work(ip_CVec1+index_K)
                        Else
                           A_KL = A_KL + Work(ip_CVec1+index_L)*
     &                                   Work(ip_CVec2+index_K)
                        End If
                     End Do
                     iSkip = (iBatK-1)*nBatVecK
                     index_A = iTri(iL,iK_Real) - 1 - iTri(iSkip,iSkip)
*-----------------------------------------------------------
*                     Write(6,*) 'iTri1', iTri(iL,iK_Real)
*                     Write(6,*) 'iTri2', iTri(iSkip,iSkip)
*                     Write(6,*) 'index_A', index_A
*-----------------------------------------------------------
                     Work(ip_Amat+index_A) = A_KL
                  End Do
               End Do
*----------------------------------------------
*               Write(6,*) 'A-matrixpiece'
*               Do i=0, l_Amat-1
*                  Write(6,*) Work(ip_Amat+i)
*               End Do
*----------------------------------------------

               iKmaxLast = nBatVecK*(iBatK-1)
               iKmaxThis = iKmaxLast + nVecK
               iAdr=iTri(iKmaxLast,iKmaxLast)+(iBatL-1)*nBatVecL
               l_Atemp = iTri(iKmaxThis,iKmaxThis)-
     &                   iTri(iKmaxLast,iKmaxLast)
               l_Atemp = l_Atemp - (iBatL-1)*nBatVecL
               l_Atemp = min(l_Atemp,l_Amat)
*               Write(6,*) 'iAdrA', iAdr
*               Write(6,*) 'l_Atemp', l_Atemp
               Call dDaFile(LuAMat,1,Work(ip_Amat),l_Atemp,iAdr)
            End Do
         End Do

*


         Call GetMem('C_vec1','Free','Real',ip_Cvec1,l_Cvec1)
         Call GetMem('C_vec2','Free','Real',ip_Cvec2,l_Cvec2)
         Call GetMem('A_matrix','Free','Real',ip_Amat,l_Amat)
#ifdef _DEBUG
         Call GetMem('A_mat','Allo','Real',ip_A,
     &               nVec(iChSym)*(nVec(iChSym)+1)/2)
         iAdr=0
         Call dDaFile(LuAMat,2,Work(ip_A),
     &                nVec(iChSym)*(nVec(iChSym)+1)/2,iAdr)
_
         Write(6,*) 'A-matrix in symm:',iChSym
         Do i = 0, nVec(iChSym)*(nVec(iChSym)+1)/2-1
            Write(6,*) Work(ip_A+i)
         End Do
         Call GetMem('A_mat','Free','Real',ip_A,
     &               nVec(iChSym)*(nVec(iChSym)+1)/2)
#endif
 1000    Continue
      End Do !iChSym



      Return
      End
