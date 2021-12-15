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
c
c kappa=\bar{kappa}
c P = \bar{d}
c D = \bar{D}
c
c
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 P(*),D(*),F(*),Kappa(*)
      Dimension dum(1)
*
      Call GetMem('MO','ALLO','REAL',ipT1,n2dens)
      Call GetMem('F1','ALLO','REAL',ipT2,ndens2)
      Call GetMem('F2','ALLO','REAL',ipT3,ndens2)
      Call GetMem('F3','ALLO','REAL',ipT4,ndens2)
c
c OIT of the Fock matrix --> Ftilde the orbital part of the effective Fock
c F = Ftilde, the active part.
c ipT1 = Dtilde
c
      Call Rint_generic(kappa,Work(ipT1),dum,Work(ipT2),
     &              Work(ipT3),
     &              Work(ipT4),F,1,-1.0d0,0)
      Call DSCAL_(ndens2,0.50d0,f,1)
c
c ipT2 = Fbar The ci part or the active Fock matrix
c
      Call FockGen(0.0d0,D,P,Work(ipT2),Work(ipT3),1)
      Fact=1.0d0
      Do iS=1,nsym
          ipF=ipF0SqMO+ipMat(is,is)-1
          If (nBas(is).ge.1) Then
             Call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),
     &                  Fact,kappa(ipMat(is,is)),nBas(is),
     &                  Work(ipF),nBas(is),1.0d0,
     &                  F(ipMat(is,is)),nBas(is))
             Call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),
     &                  -Fact,Work(ipF),nBas(is),
     &                  Kappa(ipmat(is,is)),nBas(is),1.0d0,
     &                  F(ipMat(is,is)),nBas(is))
         End If
      End Do
      Call DaxPy_(ndens2,1.0d0,F,1,Work(ipT2),1)

      Call TCMO(Work(ipT2),1,-2 )
      ijB=1
      Do iS=1,nsym
       Do iB=1,nBas(iS)
        iptmp=ipmat(iS,iS)+ipT2-1
        iptmp1=iptmp-nbas(iS)+iB-1
        iptmp2=iptmp+nbas(iS)*(iB-1)-1
        Do jB=1,iB-1
         iptmp1=iptmp1+nbas(iS)
         iptmp2=iptmp2+1
         F(ijB)=(Work(iptmp1)+Work(iptmp2))
         ijB=ijB+1
        End Do
        F(ijB)=Work(iptmp2+1)
        ijB=ijB+1
       End Do
      End Do
*
      Call GetMem('F3','Free','REAL',ipT4,ndens2)
      Call GetMem('F2','Free','REAL',ipT3,ndens2)
      Call GetMem('F1','Free','REAL',ipT2,ndens2)
      Call GetMem('MO','Free','REAL',ipT1,n2dens)
      Return
      End
