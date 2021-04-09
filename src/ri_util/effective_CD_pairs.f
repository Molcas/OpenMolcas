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
      Subroutine Effective_CD_Pairs(ij2,nij_Eff)
      use Basis_Info
      use Symmetry_Info, only: nIrrep
      use ChoArr, only: iSOShl
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
      Integer, Allocatable:: ij2(:,:)
#include "cholesky.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Integer, Allocatable :: SO_ab(:), ij3(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the max number of auxiliary shells in case of CD. Hence
*     we do not have any explicit auxiliary basis set!
*
      nSkal_Valence=0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%Aux) Then
            Do iAng = 0, dbsc(iCnttp)%nVal-1
               iShll = dbsc(iCnttp)%iVal + iAng
               If (.Not.Shells(iShll)%Aux) Then
                  nSkal_Valence = nSkal_Valence + dbsc(iCnttp)%nCntr
               End If
            End Do
         End If
      End Do
*
*     Max number of pairs
*
      nij=nSkal_Valence*(nSkal_Valence+1)/2
      Call mma_allocate(ij3,nij,Label='ij3')
      ij3(:)=0
C     Write (*,*) 'nij3=',nij
*                                                                      *
************************************************************************
*                                                                      *
      nAux_Tot=0
      nVal_Tot=0
      Do iIrrep = 0, nIrrep-1
         nAux_Tot=nAux_Tot+nBas_Aux(iIrrep)
         nVal_Tot=nVal_Tot+nBas(iIrrep)
      End Do

      Call mma_allocate(SO_ab,2*nAux_Tot,Label='SO_ab')
      SO_ab(:)=0
*
      iOff = 1
      jOff = 0
      nSym=nIrrep
      Do iSym = 1, nSym
         iIrrep=iSym-1
         ip_List_rs=ip_of_iWork(InfVec(1,1,iSym))
         Call CHO_X_GET_PARDIAG(iSym,ip_List_rs,SO_ab(iOff))
*
         Call Get_Auxiliary_Shells(SO_ab(iOff),nBas_Aux(iIrrep),
     &                             jOff,iSOShl,nVal_Tot,ij3,nij)
*
         jOff = jOff + nBas_Aux(iIrrep)
         iOff = iOff + 2*nBas_Aux(iIrrep)
      End Do
      Call mma_deallocate(SO_ab)
*                                                                      *
************************************************************************
*                                                                      *
      nij_Eff=0
      Do i = 1, nij
         nij_Eff = nij_Eff + ij3(i)
      End Do
C     Write (6,*) 'nij_Eff=',nij_Eff
      If (nij_Eff.gt.nij) Then
         Call WarningMessage(2,
     &               'Effective_CD_Pairs: nij_Eff.gt.nij')
         Call Abend()
      End If
*
      Call mma_allocate(ij2,2,nij_Eff,Label='ij2')
      ij = 0
      ij_Eff = 0
      Do i = 1, nSkal_Valence
         Do j = 1, i
            ij = ij + 1
C           Write (6,*) 'i,j,ij=',i,j,ij
            If (ij3(ij).eq.1) Then
               ij_Eff = ij_Eff + 1
C              Write (*,*) 'ij_Eff=',ij_Eff
               ij2(1,ij_Eff) = i
               ij2(2,ij_Eff) = j
            End If
         End Do
      End Do
      If (ij_Eff.ne.nij_Eff) Then
         Call WarningMessage(2,
     &               'Effective_CD_Pairs: ij_Eff.ne.nij_Eff')
         Call Abend()
      End If
      Call mma_deallocate(ij3)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine Get_Auxiliary_Shells(iSO,nSO,jOff,iSO2Shl,nSO2Shl,
     &                                iPair,nPair)
      Integer iSO(2,nSO), iSO2Shl(nSO2Shl), iPair(nPair)
*
C      Write (6,*) 'iSO'
C      Write (6,*) '==='
C      Do i = 1, nSO
C         Write (6,*) iSO(1,i), iSO(2,i)
C      End Do
*
C      Write (6,*) 'iSO2Shl'
C      Write (6,*) '======='
C      Do i = 1, nSO2Shl
C         Write (6,*) i, iSO2Shl(i)
C      End Do
       Do i = 1, nSO
          k=iSO(1,i) + jOff
          l=iSO(2,i) + jOff
          kSh=iSO2Shl(k)
          lSh=iSO2Shl(l)
C         Write (6,*) 'k,kSh=',k,kSh
C         Write (6,*) 'l,lSh=',k,lSh
          kl = Max(kSh,lSh)*(Max(kSh,lSh)-1)/2 + Min(kSh,lSh)
          iPair(kl)=1
      End Do
C     Write (6,*) 'iPairs'
C     Write (6,*) '======'
C     Do i = 1, nPair
C        Write (6,*) iPair(i)
C     End Do
*
      Return
      End
