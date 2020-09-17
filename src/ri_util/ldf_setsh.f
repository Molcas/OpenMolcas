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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetSh(nShell_Val,nShell_Aux,Verbose,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Set data in localdf_bas.fh
C
      use Index_arrays, only: iSO2Sh
      use Basis_Info, only: nBas, nBas_Aux
      Implicit Real*8 (a-h,o-z) ! info.fh misses declarations
      Integer nShell_Val
      Integer nShell_Aux
      Logical Verbose
      Integer irc
#include "itmax.fh"
#include "info.fh"
#include "localdf_bas.fh"
#include "WrkSpc.fh"

      Integer nBT, nST
      Integer iB, ii, jj

      Integer i
      Integer iSOShl, iShlSO, nBasSh
      iSOShl(i)=iWork(ip_iSOShl-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      nBasSh(i)=iWork(ip_nBasSh-1+i)

      ! Init return code
      irc=0

      ! Set number of valence and auxiliary shells
      nShell_Valence=nShell_Val
      nShell_Auxiliary=nShell_Aux

      ! Get number of valence functions from info.fh
      nBas_Valence=nBas(0)

      ! Get number of auxiliary functions from info.fh
      ! Subtract the additional dummy function
      nBas_Auxiliary=nBas_Aux(0)-1

      ! Total number of basis functions (val+aux) + 1 dummy
      nBT=nBas_Valence+nBas_Auxiliary+1

      ! Total number of shells (val+aux) + 1 dummy
      nST=nShell_Valence+nShell_Auxiliary+1

      ! Get iSOShl from index_arrays
      l_iSOShl=nBT
      Call GetMem('LDF_iSOShl','Allo','Inte',ip_iSOShl,l_iSOShl)
      Call iCopy(l_iSOShl,iSO2Sh,1,iWork(ip_iSOShl),1)

      ! Get nBasSh
      l_nBasSh=nST
      Call GetMem('LDF_nBasSh','Allo','Inte',ip_nBasSh,l_nBasSh)
      Call iZero(iWork(ip_nBasSh),l_nBasSh)
      ii=ip_nBasSh-1
      Do iB=1,nBT
         jj=ii+iSOShl(iB)
         iWork(jj)=iWork(jj)+1
      End Do

      ! Get iShlSO
      l_iShlSO=l_iSOShl
      Call GetMem('LDF_iShlSO','Allo','Inte',ip_iShlSO,l_iShlSO)
      Call Cho_SetSh2(iWork(ip_iShlSO),iWork(ip_iSOShl),
     &                iWork(ip_nBasSh),nBT,nST)

      ! Print
      If (Verbose) Then
         Call Cho_Head('Info from LDF_SetSh','-',120,6)
         Write(6,'(/,A,I8)')
     &   'Number of valence shells:  ',nShell_Valence
         Write(6,'(A,I8)')
     &   'Number of auxiliary shells:',nShell_Auxiliary
         Write(6,'(A,I8)')
     &   'Number of valence BF:      ',nBas_Valence
         Write(6,'(A,I8)')
     &   'Number of auxiliary BF:    ',nBas_Auxiliary
         Write(6,'(/,A)')
     &   '      BF    Shell Index in Shell'
         Write(6,'(32A1)') ('-',ii=1,32)
         Do iB=1,nBT
            Write(6,'(I8,1X,I8,7X,I8)')
     &      iB,iSOShl(iB),iShlSO(iB)
         End Do
         Write(6,'(32A1)') ('-',ii=1,32)
         Write(6,'(/,A,/,A)')
     &   'Val Shell   Dimension',
     &   '---------------------'
         Do ii=1,nShell_Valence
            Write(6,'(1X,I8,4X,I8)') ii,nBasSh(ii)
         End Do
         Write(6,'(A)')
     &   '---------------------'
         Write(6,'(/,A,/,A)')
     &   'Aux Shell   Dimension',
     &   '---------------------'
         Do ii=nShell_Valence+1,nShell_Valence+nShell_Auxiliary
            Write(6,'(1X,I8,4X,I8)') ii,nBasSh(ii)
         End Do
         Write(6,'(A)')
     &   '---------------------'
         Call xFlush(6)
      End If

      End
