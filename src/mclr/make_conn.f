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
      Subroutine Make_Conn(F,Kappa,P,D)
      use Constants, only: Zero,Half,One
      use Arrays, only: F0SQMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only:  ipMat, n2Dens, nDens2
      use input_mclr, only: nSym,nBas
c
c kappa=\bar{kappa}
c P = \bar{d}
c D = \bar{D}
c
c
      Implicit None
      Real*8 P(*),D(*),F(*),Kappa(*)

      Real*8 dum(1)
      Real*8, Allocatable:: T1(:), T2(:), T3(:), T4(:)
      Real*8 Fact
      Integer iS, ijB, iB, ipTmp, ipTmp1, ipTmp2, jB
*
      Call mma_allocate(T1,n2dens,Label='MO')
      Call mma_allocate(T2,ndens2,Label='F1')
      Call mma_allocate(T3,ndens2,Label='F3')
      Call mma_allocate(T4,ndens2,Label='F2')
c
c OIT of the Fock matrix --> Ftilde the orbital part of the effective Fock
c F = Ftilde, the active part.
c T1 = Dtilde
c
      Call Rint_generic(kappa,T1,dum,T2,T3,T4,F,1,-One,0)
      Call DSCAL_(ndens2,Half,f,1)
c
c T2 = Fbar The ci part or the active Fock matrix
c
      Call FockGen(Zero,D,P,T2,T3,1)
      Fact=One
      Do iS=1,nsym
          If (nBas(is).ge.1) Then
             Call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),
     &                  Fact,kappa(ipMat(is,is)),nBas(is),
     &                       F0SQMO(ipMat(is,is)),nBas(is),
     &                  One,F(ipMat(is,is)),nBas(is))
             Call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),
     &                  -Fact,F0SQMO(ipMat(is,is)),nBas(is),
     &                        Kappa(ipmat(is,is)),nBas(is),
     &                  One,F(ipMat(is,is)),nBas(is))
         End If
      End Do
      Call DaxPy_(ndens2,One,F,1,T2,1)

      Call TCMO(T2,1,-2 )
      ijB=1
      Do iS=1,nsym
       Do iB=1,nBas(iS)
        iptmp=ipmat(iS,iS)
        iptmp1=iptmp-nbas(iS)+iB-1
        iptmp2=iptmp+nbas(iS)*(iB-1)-1
        Do jB=1,iB-1
         iptmp1=iptmp1+nbas(iS)
         iptmp2=iptmp2+1
         F(ijB)=(T2(iptmp1)+T2(iptmp2))
         ijB=ijB+1
        End Do
        F(ijB)=T2(iptmp2+1)
        ijB=ijB+1
       End Do
      End Do
*
      Call mma_deallocate(T4)
      Call mma_deallocate(T3)
      Call mma_deallocate(T2)
      Call mma_deallocate(T1)
      End Subroutine Make_Conn
