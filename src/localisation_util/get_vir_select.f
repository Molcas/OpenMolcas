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
      SUBROUTINE Get_Vir_Select(irc,CMO,XMO,Eorb,Smat,Name,NamAct,
     &                   ind_V,nSym,nActa,mOrb,nBas,ortho,n_OK)

************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
      Real*8  CMO(*), XMO(*), Eorb(*), Smat(*)
      Integer irc, nSym, nActa, mOrb(nSym), nBas(nSym), n_OK(nSym)
      Integer ind_V(*)
      Character*(LENIN8)  Name(*)
      Character*(LENIN) NamAct(nActa)
      Logical ortho
#include "WrkSpc.fh"
      Character*(LENIN) tmp
************************************************************************
      kD(i) = iWork(ip_iD-1+i)
      lD(i) = iWork(ip_iD+nOrbmx-1+i)
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
      Call GetMem('iD','Allo','Inte',ip_iD,2*nOrbmx)
      Call GetMem('LCMO','Allo','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
      lScr=nBmx*nOrbmx
      ip_CC=ip_C+lScr
      ip_X=ip_CC+lScr
      iZ=ip_X+lScr
      ipScr=iZ+lScr
      ip_U=ipScr+lScr
      ip_Fock=ip_U+nOrbmx**2
      jp_Fock=ip_Fock+nOrbmx
*
      iOff=0
      jOff=0
      kOff=0
      lOff=0
      Do iSym=1,nSym

         iS=kOff+1
         n_OK(iSym)=0
         n_KO=0
         Do i=1,mOrb(iSym)
            ja=iOff+ind_V(i+iOff)
            tmp=Name(ja)(1:LENIN)
            jfr=jOff+nBas(iSym)*(i-1)+1
            kx=0
*
*            write(6,*)' We simulate Afreeze with all Vir'
*            kx=1
*
            Do j=1,nActa
               If(NamAct(j).eq.tmp) kx=kx+1
            End Do
            If (kx.gt.0) Then
                 jX=ip_X+nBas(iSym)*n_OK(iSym)
                 call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jX),1)
                 iWork(ip_iD+n_OK(iSym))=i
                 n_OK(iSym)=n_OK(iSym)+1
            Else
                 jZ=iZ+nBas(iSym)*n_KO
                 call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jZ),1)
                 iWork(ip_iD+nOrbmx+n_KO)=i
                 n_KO=n_KO+1
            EndIf
         End Do
*
         If (.not.ortho) Then

            Call Ortho_orb(Work(ip_X),Smat(iS),nBas(iSym),
     &                     n_OK(iSym),2,.false.)
            Call Ortho_orb(Work(iZ),Smat(iS),nBas(iSym),
     &                     n_KO,2,.false.)
         EndIf
*
         nBx=Max(1,nBas(iSym))
         mOx=Max(1,mOrb(iSym))
*
         Call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),
     &                           1.0d0,Cmo(jOff+1),nBx,
     &                                 Smat(iS),nBx,
     &                           0.0d0,Work(ip_CC),mOx)
*
         If (n_KO .gt.0 ) Then
           Call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),
     &                              1.0d0,Work(ip_CC),mOx,
     &                                    Work(iZ),nBx,
     &                              0.0d0,Work(ip_U),mOx)
*
           Call Get_Can_Lorb(Eorb(lOff+1),Work(ip_Fock),n_KO,mOrb(iSym),
     &                        iWork(ip_iD+nOrbmx),Work(ip_U),iSym)
*
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
     &                     iWork(ip_iD),Work(ip_U),iSym)
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
*
      Call GetMem('LCMO','Free','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
      Call GetMem('iD','Free','Inte',ip_iD,2*nOrbmx)
*
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
