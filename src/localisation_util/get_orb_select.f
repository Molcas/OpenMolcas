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
      SUBROUTINE Get_Orb_Select(irc,CMO,XMO,Eorb,Smat,Saa,Name,NamAct,
     &                          nSym,nActa,mOrb,nBas,ortho,ThrSel,n_OK)

************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
      Real*8  CMO(*), XMO(*), Eorb(*), Smat(*), Saa(*), ThrSel
      Integer irc, nSym, nActa, mOrb(nSym), nBas(nSym), n_OK(nSym)
      Character*(LENIN8)  Name(*)
      Character*(LENIN) NamAct(nActa)
      Logical ortho
#include "WrkSpc.fh"
      Character*(LENIN) tmp
************************************************************************
      jD(i) = iWork(ip_iD-1+i)
******
      kD(i) = iWork(ip_iD+nBmx-1+i)
******
      lD(i) = iWork(ip_iD+nBmx+nOrbmx-1+i)
************************************************************************
*
      irc=0
*
      nBmx=0
      nOrbmx=0
      Do iSym=1,nSym
         nBmx=Max(nBmx,nBas(iSym))
         nOrbmx=Max(nOrbmx,mOrb(iSym))
      End Do
      Call GetMem('iD','Allo','Inte',ip_iD,2*nOrbmx+nBmx)
      Call GetMem('Smx','Allo','Real',iSQ,nBmx**2)
      Call GetMem('LCMO','Allo','Real',ip_C,(5*nBmx+nOrbmx+3)*nOrbmx)
      lScr=nBmx*nOrbmx
      ip_CC=ip_C+lScr
      ip_X=ip_CC+lScr
      iZ=ip_X+lScr
      ipScr=iZ+lScr
      ip_U=ipScr+lScr
      iQ=ip_U+nOrbmx**2
      ip_Fock=iQ+nOrbmx
      jp_Fock=ip_Fock+nOrbmx
      Call Fzero(Work(iQ),nOrbmx)
*
      iOff=0
      jOff=0
      kOff=0
      lOff=0
      Do iSym=1,nSym

         nBa=0
         Do ia=1,nBas(iSym)
            ja=ia+iOff
            tmp=Name(ja)(1:LENIN)
            Do j=1,nActa
               If(NamAct(j).eq.tmp) Then
                  iWork(ip_iD+nBa)=ia
                  nBa=nBa+1
               EndIf
            End Do
         End Do
         Do ia=1,nBa
            ifr=jOff+jD(ia)
            ito=ip_C+ia-1
            call dcopy_(mOrb(iSym),Xmo(ifr),nBas(iSym),Work(ito),nBa)
         End Do
         iS=kOff+1
         Do ia=1,nBa
            jb=jD(ia)
            jfr=iS+nBas(iSym)*(jb-1)
            jto=iSQ+nBas(iSym)*(ia-1)
            call dcopy_(nbas(iSym),Smat(jfr),1,Work(jto),1)
         End Do
         nBx=Max(1,nBas(iSym))
         nBax=Max(1,nBa)
         Call DGEMM_('T','N',nBa,mOrb(iSym),nBas(iSym),
     &                    1.0d0,Work(iSQ),nBx,
     &                          Xmo(jOff+1),nBx,
     &                    0.0d0,Work(iZ),nBax)
         Do i=0,mOrb(iSym)-1
            jQ=iQ+i
            jC=ip_C+nBa*i
            jZ=iZ+nBa*i
            Work(jQ)=ddot_(nBa,Work(jC),1,Work(jZ),1)**2
         End Do
         n_OK(iSym)=0
         n_KO=0
         Do i=1,mOrb(iSym)
            ThrS=ThrSel*Saa(lOff+i)
            jQ=iQ+i-1
            jfr=jOff+nBas(iSym)*(i-1)+1
            If (sqrt(Work(jQ)).ge.ThrS) Then
               jX=ip_X+nBas(iSym)*n_OK(iSym)
               call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jX),1)
               iWork(ip_iD+nBmx+n_OK(iSym))=i
               n_OK(iSym)=n_OK(iSym)+1
            Else
               jZ=iZ+nBas(iSym)*n_KO
               call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jZ),1)
               iWork(ip_iD+nBmx+nOrbmx+n_KO)=i
               n_KO=n_KO+1
            EndIf
         End Do
*
         mOx=Max(1,mOrb(iSym))

         If (.not.ortho) Then

            Call Ortho_orb(Work(ip_X),Smat(iS),nBas(iSym),
     &                     n_OK(iSym),2,.false.)
            Call Ortho_orb(Work(iZ),Smat(iS),nBas(iSym),
     &                     n_KO,2,.false.)
         EndIf
*
         Call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),
     &                           1.0d0,Cmo(jOff+1),nBx,
     &                                 Smat(iS),nBx,
     &                           0.0d0,Work(ip_CC),mOx)

         If (n_KO .gt. 0) Then
           Call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),
     &                             1.0d0,Work(ip_CC),mOx,
     &                                   Work(iZ),nBx,
     &                             0.0d0,Work(ip_U),mOx)

           Call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_KO,mOrb(iSym),
     &                     iWork(ip_iD+nBmx+nOrbmx),Work(ip_U),iSym)

           Call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,
     &                  1.0d0,Work(iZ),nBx,
     &                        Work(ip_U),n_KO,
     &                  0.0d0,Work(ipScr),nBx)
*
*  Reorder the final MOs such that those of the active site come first
           km=jOff+nBas(iSym)*n_OK(iSym)+1
           call dcopy_(nBas(iSym)*n_KO,Work(ipScr),1,Cmo(km),1)
           call dcopy_(nOrbmx,Work(ip_Fock),1,Work(jp_Fock),1)
         EndIf
*
         Call DGEMM_('N','N',mOrb(iSym),n_OK(iSym),nBas(iSym),
     &                           1.0d0,Work(ip_CC),mOx,
     &                                 Work(ip_X),nBx,
     &                           0.0d0,Work(ip_U),mOx)
*
         Call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_OK(iSym),
     &                     mOrb(iSym),
     &                     iWork(ip_iD+nBmx),Work(ip_U),iSym)
*
         nOx=Max(1,n_OK(iSym))
         Call DGEMM_('N','N',nBas(iSym),n_OK(iSym),n_OK(iSym),
     &                1.0d0,Work(ip_X),nBx,
     &                      Work(ip_U),nOx,
     &                0.0d0,Work(ipScr),nBx)
*
         Do i=1,n_OK(iSym)
            j=kD(i)
            k=lOff+i
            l=ip_Fock+j-1
            Eorb(k)=Work(l)
         End Do
         km=jOff+1
         call dcopy_(nBas(iSym)*n_OK(iSym),Work(ipScr),1,Cmo(km),1)
*
         Do i=1,n_KO
            j=lD(i)
            k=lOff+n_OK(iSym)+i
            l=jp_Fock+j-1
            Eorb(k)=Work(l)
         End Do
*
         iOff=iOff+nBas(iSym)
         jOff=jOff+nBas(iSym)*mOrb(iSym)
         kOff=kOff+nBas(iSym)**2
         lOff=lOff+mOrb(iSym)
      End Do

      Call GetMem('LCMO','Free','Real',ip_C,(5*nBmx+nOrbmx+3)*nOrbmx)
      CALL GetMem('Smx','Free','Real',iSQ,nBmx**2)
      Call GetMem('iD','Free','Inte',ip_iD,2*nOrbmx+nBmx)

      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Get_Can_Lorb(Ene,Fock,nO,nX,jOrb,Umat,iSym)

      Implicit Real*8 (a-h,o-z)
      Real*8  Ene(*), Fock(*), Umat(*)
      Integer nO, nX, jOrb(nO), iSym
#include "WrkSpc.fh"
*
*
      If (nO.lt.1)  Return
*
      Call GetMem('eta_ik','Allo','Real',ip_eta,2*nX**2+1)
      ip_Z=ip_eta+nX**2
      ip_ZZ=ip_Z+nX
      Call FZero(Work(ip_eta),nX**2)
      Do i=1,nX
         ii=ip_eta+nX*(i-1)+i-1
         Work(ii)=Ene(i)
      End Do
*
      nXx=Max(1,nX)
      nOx=Max(1,nO)
      Call DGEMM_('N','N',nX,nO,nX,1.0d0,Work(ip_eta),nXx,
     &                                  Umat(1),nXx,
     &                            0.0d0,Work(ip_Z),nXx)
      Call DGEMM_('T','N',nO,nO,nX,1.0d0,Umat(1),nXx,
     &                                  Work(ip_Z),nXx,
     &                            0.0d0,Work(ip_eta),nOx)

      Call Eigen_Molcas(nO,Work(ip_eta),Work(ip_Z),Work(ip_ZZ))

      call dcopy_(nO**2,Work(ip_eta),1,Umat,1)
      Do i=1,nO
         ii=ip_Z+i-1
         j=jOrb(i)
         Fock(j)=Work(ii)
         kk=ip_ZZ+i-1
      End Do
      Call GetMem('eta_ik','Free','Real',ip_eta,2*nX**2+1)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSym)
      End
************************************************************************
