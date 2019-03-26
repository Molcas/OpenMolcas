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
* Copyright (C) 2017, Roland Lindh                                     *
************************************************************************
      Subroutine Virt_Space(C_Occ,C_Virt,Ovrlp,nBas,nOcc,nVirt)
************************************************************************
*     The generation of starting orbitals suffer from the fact that    *
*     the virtual orbitals are not well defined. This routine is       *
*     supposed to generate a set well defined virtual orbitals from a  *
*     set of well defined occupid orbitals.                            *
************************************************************************
      Implicit None
#include "stdalloc.fh"
*define _DEBUG_
#ifdef _DEBUG_
      Integer i,j
      Integer ij, iTri
#endif
      Integer nOcc,nVirt,nBas
      Integer iBas, jBas, iOcc, kBas, iVirt, lBas
      Integer mVirt
      Logical Polished
      Real*8 C_Occ(nBas,nOcc), C_Virt(nBas,nVirt),
     &       Ovrlp(nBas*(nBas+1)/2)
      Real*8 tmp
      Real*8, Dimension(:,:), Allocatable:: P, Ovrlp_Sq, EVe, C_tmp
      Real*8, Dimension(:), Allocatable:: PNew, EVa
#ifdef _DEBUG_
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      Write (6,*) 'nBas,nOcc,nVirt=',nBas,nOcc,nVirt
      Call RecPrt('C_Occ',' ',C_Occ,nBas,nOcc)
      Call RecPrt('C_Virt',' ',C_Virt,nBas,nVirt)
      Call TriPrt('Ovrlp',' ',Ovrlp,nBas)
      Do iOcc = 1, nOcc
         tmp = 0.0D0
         Do iBas = 1, nBas
            Do jBas = 1, nBas
               ij = iTri(iBas,jBas)
               tmp = tmp + C_Occ(iBas,iOcc)*Ovrlp(ij)*
     &                     C_Occ(jBas,iOcc)
            End Do
         End Do
         Write (6,*) 'iOcc,tmp=',iOcc,tmp
      End Do
      Do iVirt = 1, nVirt
         tmp = 0.0D0
         Do iBas = 1, nBas
            Do jBas = 1, nBas
               ij = iTri(iBas,jBas)
               tmp = tmp + C_Virt(iBas,iVirt)*Ovrlp(ij)*
     &                     C_Virt(jBas,iVirt)
            End Do
         End Do
         Write (6,*) 'iVirt,tmp=',iOcc,tmp
      End Do
#endif
*
      If (nVirt.eq.0) Call Abend()
*
*     Compute S^{1/2}
*
      Call mma_allocate(Ovrlp_Sq,nBas,nBas,Label='Ovrlp_Sq')
      Call mma_allocate(EVa,nBas*(nBas+1)/2,Label='EVa')
      Call mma_allocate(EVe,nBas,nBas,Label='EVe')
      Call FZero(EVe,nBas**2)
      Call DCopy_(nBas,[1.0D0],0,EVe,nBas+1)
      Call DCopy_(nBas*(nBas+1)/2,Ovrlp,1,EVa,1)
      Call NIdiag(EVa,EVe,nBas,nBas,0)
*
      Do iBas = 2, nBas
         EVa(iBas) = EVa(iBas*(iBas+1)/2)
      End Do
*
      Do iBas = 1, nBas
         Do jBas = 1, nBas
            tmp = 0.0D0
            Do kBas = 1, nBas
               tmp = tmp +
     &               EVe(iBas,kBas)*Sqrt(EVa(kBas))*EVe(jBas,kBas)
            End Do
            Ovrlp_Sq(iBas,jBas)=tmp
         End Do
      End Do
*
*     Now transform the the CMOs of the occupied orbitals to an
*     orthonormal basis.
*
      Call mma_allocate(C_tmp,nBas,nOcc,Label='C_tmp')
      Call FZero(C_tmp,nBas*nOcc)
      Call DGEMM_('N','N',
     &            nBas,nOcc,nBas,
     &            1.0D0,Ovrlp_Sq,nBas,
     &                  C_Occ,nBas,
     &            0.0D0,C_tmp,nBas)
*
      Call mma_allocate(P,nBas,nBas,Label='P')
      Call mma_allocate(PNew,nBas,Label='PNew')
*
*     Form the P matrix = 1 - |C(occ)><C(occ)|
*
      Do iBas = 1, nBas
         Do jBas = 1, nBas
            tmp=0.0D0
            If (iBas.eq.jBas) tmp=1.0D0
            Do iOcc = 1, nOcc
               tmp = tmp - C_tmp(iBas,iOcc)*C_tmp(jBas,iOcc)
            End Do
            P(iBas,jBas) = tmp
         End Do
      End Do
#ifdef _DEBUG_
      Call RecPrt('P-mat',' ',P,nBas,nBas)
      Do iBas = 1, nBas
         Write (6,*) 'iBas,P(iBas,iBas)=',iBas,P(iBas,iBas)
      End Do
#endif
*
*     Now, use Cholesky Decomposition to reduce
*     the size of the P-matrix from a nBas x nBas size
*     down to something that is nBas x nVirt
*
      mVirt = 0
      Do kBas = 1, nBas
*
         tmp=0.0D0
         lBas = 0
         Do iBas = 1, nBas
            If (P(iBas,iBas).gt.tmp) Then
               tmp=P(iBas,iBas)
               lBas = iBas
            End If
         End Do
*
*        Pick up a vector and project it against the previous
*
         Call DCopy_(nBas,P(1,lBas),1,PNew,1)
*
*        Normalize PNew
*
         tmp=0.0D0
         Do iBas = 1, nBas
            tmp = tmp + PNew(iBas)**2
         End Do
*
*        Skip if this is a null vector.
*
         If (tmp.lt.1.0D-14) Cycle
*
         tmp=1.0D0/Sqrt(tmp)
         Call DScal_(nBas,tmp,PNew,1)
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) 'New Trial vector'
         Write (6,*) '================'
         Write (6,*) 'kBas,lBas,tmp=',kBas,lBas,tmp
         Call RecPrt('Normalized PNew',' ',PNew,nBas,1)
#endif
*
*        Polished = .False.
         Polished = .True.
 666     Continue
*
         Do iOcc = 1, nOcc
*
*           From the trial vector eliminate the occupied
*           space
*
            tmp=0.0D0
            Do iBas = 1, nBas
               tmp = tmp + PNew(iBas)*C_tmp(iBas,iOcc)
            End Do
#ifdef _DEBUG_
            Write (6,*) 'iOcc,tmp=',iOcc,tmp
#endif
*           Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)
*
            Call DaXpY_(nBas,-tmp,C_Occ(1,iOcc),1,PNew,1)
         End Do
         Do iVirt = 1, mVirt
*
*           From the trial vector eliminate parts which
*           already expressed.
*
            tmp=0.0D0
            Do iBas = 1, nBas
               tmp = tmp + PNew(iBas)*C_Virt(iBas,iVirt)
            End Do
#ifdef _DEBUG_
            Write (6,*) 'iVirt,tmp=',iVirt,tmp
#endif
*
*           Form PNew(2) = P(2) - <PNew(1)|Ovrlp|P(2)>PNew(1)
*
            Call DaXpY_(nBas,-tmp,C_Virt(1,iVirt),1,PNew,1)
         End Do
*
*        Test that it is not a null vector!
*
         tmp = 0.0D0
         Do iBas = 1, nBas
            tmp = tmp + PNew(iBas)**2
         End Do
#ifdef _DEBUG_
         Write (6,*) 'Norm after projection:',tmp
         Write (6,*)
#endif
         If (tmp.gt.1.0D-14) Then
            tmp=1.0D0/Sqrt(tmp)
            Call DScal_(nBas,1.0D0/Sqrt(tmp),PNew,1)
            If (.Not.Polished) Then
               Polished=.True.
               Go To 666
            End If
*
            If (tmp.gt.1.0D-14) Then
               If (mVirt+1.gt.nVirt) Then
                  Write (6,*) 'mVirt.gt.nVirt'
                  Write (6,*) 'mVirt=',mVirt
                  Write (6,*) 'nVirt=',nVirt
                  Call Abend()
               End If
               mVirt=mVirt+1
               Call DCopy_(nBas,PNew,1,C_Virt(1,mVirt),1)
*
*              Update the P-matrix.
*
               Do iBas = 1, nBas
                  Do jBas = 1, nBas
                     P(iBas,jBas)=P(iBas,jBas)-PNew(iBas)*PNew(jBas)
                  End Do
               End Do
            End If
         End If
*
         If (mVirt.eq.nVirt) Go To 667
*
      End Do
*
 667  Continue
      Call mma_deallocate(C_tmp)
*
      If (mVirt.ne.nVirt) Then
         Write (6,*) 'mVirt.ne.nVirt'
         Write (6,*) 'mVirt,nVirt=',mVirt,nVirt
         Call Abend()
      End If
*
      Call mma_deallocate(P)
      Call mma_deallocate(PNew)
*
*     Form S^(-1/2)
*
      Do iBas = 1, nBas
         Do jBas = 1, nBas
            tmp = 0.0D0
            Do kBas = 1, nBas
               tmp = tmp +
     &              (EVe(iBas,kBas)*EVe(jBas,kBas))/Sqrt(EVa(kBas))
            End Do
            Ovrlp_Sq(iBas,jBas)=tmp
         End Do
      End Do
      Call mma_deallocate(EVe)
      Call mma_deallocate(EVa)
*
      Call mma_allocate(C_tmp,nBas,nVirt,Label='C_tmp')
      Call DCopy_(nBas*nVirt,C_Virt,1,C_tmp,1)
      Call DGEMM_('N','N',
     &            nBas,nVirt,nBas,
     &            1.0D0,Ovrlp_Sq,nBas,
     &                  C_tmp,nBas,
     &            0.0D0,C_Virt,nBas)
      Call mma_deallocate(C_tmp)
      Call mma_deallocate(Ovrlp_Sq)
*
#ifdef _DEBUG_
      Call RecPrt('C_Virt(New)',' ',C_Virt,nBas,nVirt)
#endif
*
      Return
      End
