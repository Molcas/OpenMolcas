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
* Copyright (C) 1990,1991,1993,1998,2005, Roland Lindh                 *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Post_2Center_LDF(ipA_Diag,ipAB,MaxCntr,Lu_AB,ipLocal_A,
     &                            nLocal_A,ipSO2C,nSO_Aux)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals.                          *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              Timing                                                  *
*              Setup_Ints                                              *
*              Eval_Ints                                               *
*              Term_Ints                                               *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March 1990                                               *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August 1991                        *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. 1993         *
*             Modified driver. Jan. 1998                               *
*             Modified to 2-center ERIs for RI June 2005               *
*                                                                      *
************************************************************************
      use Basis_Info, only: nBas_Aux
      use Wrj12
      use Temporary_Parameters, only: force_out_of_core
      use RICD_Info, only: Thrshld_CD
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "nsd.fh"
      Character Name_Q*6
      Integer nQvec(0:7)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      Call QEnter('POST2LDF')
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the number of AB blocks
*
      nAB=MaxCntr*(MaxCntr+1)
      Call GetMem('A-blocks','Allo','Inte',ipAB,nAB)
      Call iZero(iWork(ipAB),nAB)
*
*     Find the size of the largest AB block!
*
      Max_AA=0
      Do iCenter = 1, MaxCntr
         nCenter=0
         Do i = 1, nSO_Aux
            jCenter=iWork(ipSO2C+i-1)
            If (jCenter.eq.iCenter) nCenter=nCenter+1
         End Do
         Max_AA = Max(Max_AA,nCenter)
      End Do
      nLocal_A=(2*Max_AA)**2
      Call GetMem('Local_A','Allo','Real',ipLocal_A,2*nLocal_A)
      ipLocal_AInv = ipLocal_A + nLocal_A
      Call GetMem('SO2lO','Allo','Inte',ipSO2lO,nSO_Aux)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Some setup for the generation of the Q-vectors
*
      nB=0
      nXZ=0
      nQm_Tot=0
      lScr=0
      Do iIrrep = 0, nIrrep-1
         iOffA(3,iIrrep)=nB
         mB=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) mB = mB - 1
         nB=nB+mB
         nQm_Tot=nQm_Tot + mB**2
         nXZ = Max(nXZ,mB)
      End Do
*
      nA_Diag=nB
      lScr=3*nXZ
      Call GetMem('iD_Diag','Allo','Inte',ip_iDiag,nA_Diag)
      Call GetMem('ip_Scr','Allo','Real',ip_Scr,lScr)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*     Process the atomic and diatomic blocks of the A matrix
*
      iAdr_AB=0
      Lu_AB=22
      Call DaName_MF_WA(Lu_AB,'ABLOCKS')
      Do iCenter = 1, MaxCntr
         niSO=0
         nlO=0
         Do iSO = 1, nSO_Aux
            kCenter = iWork(ipSO2C+iSO-1)
            If (kCenter.eq.iCenter) Then
               niSO=niSO+1
               nlO = nlO + 1
               iWork(ipSO2lO+iSO-1)=nlO
            End If
         End Do

         Do jCenter = 1, iCenter
            ijCenter = iCenter*(iCenter-1)/2 + jCenter
*
*           Process the atomic blocks
*
            iAddr=0
            iIrrep=0
            If (jCenter.eq.iCenter) Then
*
               Call FZero(Work(ipLocal_A),niSO**2)
*
               Do iSO = 1, nSO_Aux
                  Call dDaFile(Lu_A(iIrrep),2,Work(ip_Scr),
     &                         nSO_Aux,iAddr)
C                 Call RecPrt('Work(ip_Scr)','(6G23.15)',
C    &                         Work(ip_Scr),1,nSO_Aux)
                  kCenter = iWork(ipSO2C+iSO-1)
                  If (kCenter.eq.iCenter) Then
                      ilO=iWork(ipSO2lO+iSO-1)
                      Do jSO = 1, iSO
                         lCenter = iWork(ipSO2C+jSO-1)
                         If (lCenter.eq.iCenter) Then
                            AElement=Work(ip_Scr+jSO-1)
                            jlO=iWork(ipSO2lO+jSO-1)
                            ij = (ilO-1)*nlO + jlO
                            ji = (jlO-1)*nlO + ilO
                            AElement=Work(ip_Scr+jSO-1)
                            Work(ipLocal_A+ij-1)=AElement
                            Work(ipLocal_A+ji-1)=AElement
                         End If
                      End Do
                  End If
               End Do
C              Call RecPrt('Local_A','(6G23.15)',
C    &                     Work(ipLocal_A),nlO,nlO)
*
               Call CD_AInv(Work(ipLocal_A),nlO,Work(ipLocal_AInv),
     &                      Thrshld_CD)
*
C              Call RecPrt('Local_AInv','(6G23.15)',
C    &                     Work(ipLocal_AInv),nlO,nlO)
*
*              Store the address of the matrix, size and write the
*              matrix to *              disk.
*
               iWork(ipAB+(ijCenter-1)*2  )=iAdr_AB
               iWork(ipAB+(ijCenter-1)*2+1)=nlO
               Call dDaFile(Lu_AB,1,Work(ipLocal_AInv),nlO**2,iAdr_AB)
*
            Else
*
*              Process the diatomic blocks
*
               njSO=0
               mlO = 0
               Do jSO = 1, nSO_Aux
                  lCenter = iWork(ipSO2C+jSO-1)
                  If (lCenter.eq.jCenter) Then
                     njSO=njSO+1
                     mlO = mlO + 1
                     iWork(ipSO2lO+jSO-1)=mlO+nlO
                  End If
               End Do
*
               Call FZero(Work(ipLocal_A),(niSO+njSO)**2)
*
               Do iSO = 1, nSO_Aux
                  Call dDaFile(Lu_A(iIrrep),2,Work(ip_Scr),
     &                         nSO_Aux,iAddr)
                  kCenter = iWork(ipSO2C+iSO-1)
                  If (kCenter.eq.iCenter.or.
     &                kCenter.eq.jCenter) Then
                      ilO=iWork(ipSO2lO+iSO-1)
                      Do jSO = 1, iSO
                         lCenter = iWork(ipSO2C+jSO-1)
                         If (lCenter.eq.iCenter.or.
     &                       lCenter.eq.jCenter) Then
                            AElement=Work(ip_Scr+jSO-1)
                            jlO=iWork(ipSO2lO+jSO-1)
                            ij = (ilO-1)*(nlO+mlO) + jlO
                            ji = (jlO-1)*(nlO+mlO) + ilO
                            AElement=Work(ip_Scr+jSO-1)
                            Work(ipLocal_A+ij-1)=AElement
                            Work(ipLocal_A+ji-1)=AElement
                         End If
                      End Do
                  End If
               End Do
*
               Call CD_AInv(Work(ipLocal_A),nlO+mlO,Work(ipLocal_AInv),
     &                      Thrshld_CD)
*
*              Store the address of the matrix, size and write the
*               matrix to *              disk.
*
               iWork(ipAB+(ijCenter-1)*2  )=iAdr_AB
               iWork(ipAB+(ijCenter-1)*2+1)=nlO+mlO
               Call dDaFile(Lu_AB,1,Work(ipLocal_AInv),(nlO+mlO)**2,
     &                      iAdr_AB)
*
            End If
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Fill in the lower part of the A matrix as it is stored on disk.
*
      Call GetMem('MemMax','Max','Real',iDummy,MaxMem)
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB = nB - 1 ! subtract dummay af
         Call Square_A(Lu_A(iIrrep),nB,MaxMem,Force_Out_of_Core)
      End Do
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      ThrQ=1.0D-12
      ichk=0
      is=0
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB = nB - 1 ! subtract dummy aux. func.
         If (nB.gt.0) Then
*
            iSeed=55+iIrrep
            Lu_Q(iIrrep)=IsFreeUnit(iSeed)
            Write(Name_Q,'(A4,I2.2)') 'QVec',iIrrep
            Call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)
*
            Call get_pivot_idx(Work(ipA_Diag+is),nB,nQvec(iIrrep),
     &                         Lu_A(iIrrep),Lu_Q(iIrrep),
     &                         iWork(ip_iDiag+is),Work(ip_Scr),lScr,
     &                         ThrQ)
            nChV(iIrrep)=nQvec(iIrrep)
            ichk=ichk+Min(1,nB-nQvec(iIrrep))
            is=is+nB
         End If
         Call DaClos(Lu_A(iIrrep))
      End Do
*
      If (ichk.ne.0) Then
         Write(6,*)
         Write(6,*)'Detected lin. dependences in the auxiliary basis.'
         Write(6,'(A,8I6)')
     & ' # of AuxBas before l. d. removal: ',
     &   nBas_Aux(0)-1,(nBas_Aux(i),i=1,nIrrep-1)
         Write(6,'(A,8I6)')
     & ' # of AuxBas after  l. d. removal: ',(nQvec(i),i=0,nIrrep-1)
         Write(6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('ip_Scr','Free','Real',ip_Scr,lScr)
      Call GetMem('iD_Diag','Free','Inte',ip_iDiag,nA_Diag)
      Call GetMem('A_Diag','Free','Real',ipA_Diag,nA_Diag)
      Call GetMem('SO2lO','Free','Inte',ipSO2lO,nSO_Aux)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Post2LDF')
      Return
      End
