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
* Copyright (C) 2008, Francesco Aquilante                              *
************************************************************************
      SUBROUTINE Loc_Nat_orb(irc,Cmo,Xmo,OccN,mOrb)
************************************************************************
*                                                                      *
*     Purpose: compute Localized Natural Orbitals (LNOs) to be used    *
*              for instance in Effective Bond Order (EBO) analysis     *
*                                                                      *
*              The density matrix in localized MO basis is             *
*              diagonalized separately in 2 subblocks defined by the   *
*              orbitals that extend within two distinct subregions     *
*              (atoms) of the molecule.                                *
*              The gross Mulliken population of each orbital on the    *
*              atoms of the two subregions determines the splitting    *
*              of the orbitals.                                        *
*              The atoms defining the "active subregion" are specified *
*              by the user with the keyword LOCN .                     *
*              The threshold used for the orbitals splitting criterium *
*              is also required within the keyword LOCN .              *
*                                                                      *
*     Author: F. Aquilante   (Geneva, Feb. 2008)                       *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8  Cmo(*), Xmo(*), OccN(*)
      Integer irc, mOrb(*)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "inflocal.fh"
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
*----------------------------------------------------------------------*
*     Read the overlap matrix                                          *
*----------------------------------------------------------------------*
      nnB=0
      nBmx=0
      nOrbmx=0
      Do iSym=1,nSym
         nnB=nnB+nBas(iSym)*(nBas(iSym)+1)/2
         nBmx=Max(nBmx,nBas(iSym))
         nOrbmx=Max(nOrbmx,mOrb(iSym))
      End Do
      Call GetMem('iD','Allo','Inte',ip_iD,2*nOrbmx+nBmx)
      Call GetMem('Smat','Allo','Real',ipS,nnB+nBmx**2)
      iSQ=ipS+nnB
      isymlbl=1
      Call RdOne(irc,6,'Mltpl  0',1,Work(ipS),isymlbl)
      If (irc.ne.0) Return
      Call GetMem('LCMO','Allo','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
      lScr=nBmx*nOrbmx
      ip_CC=ip_C+lScr
      ip_X=ip_CC+lScr
      iZ=ip_X+lScr
      ipScr=iZ+lScr
      ip_U=ipScr+lScr
      iQ=ip_U+nOrbmx**2
      ip_FOcc=iQ+nOrbmx
*
      iOff=0
      jOff=0
      kOff=0
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
         iS=ipS+kOff
         Do ia=1,nBa
            jb=jD(ia)
            jfr=iS+jb*(jb-1)/2
            jto=iSQ+nBas(iSym)*(ia-1)
            call dcopy_(jb,Work(jfr),1,Work(jto),1)
            jto=jto+jb
            Do ka=jb+1,nBas(iSym)
               iab=iS+ka*(ka-1)/2+jb-1
               Work(jto)=Work(iab)
               jto=jto+1
            End Do
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
         n_OK=0
         n_KO=0
         Do i=1,mOrb(iSym)
            jQ=iQ+i-1
            jfr=jOff+nBas(iSym)*(i-1)+1
            If (sqrt(Work(jQ)).ge.ThrSel) Then
               jX=ip_X+nBas(iSym)*n_OK
               call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jX),1)
               iWork(ip_iD+nBmx+n_OK)=i
               n_OK=n_OK+1
            Else
               jZ=iZ+nBas(iSym)*n_KO
               call dcopy_(nBas(iSym),Xmo(jfr),1,Work(jZ),1)
               iWork(ip_iD+nBmx+nOrbmx+n_KO)=i
               n_KO=n_KO+1
            EndIf
         End Do
         Call Square(Work(iS),Work(iSQ),1,nBas(iSym),nBas(iSym))
         mOx=Max(1,mOrb(iSym))
         Call DGEMM_('T','N',mOrb(iSym),nBas(iSym),nBas(iSym),
     &                           1.0d0,Cmo(jOff+1),nBx,
     &                                 Work(iSQ),nBx,
     &                           0.0d0,Work(ip_CC),mOx)

         Call DGEMM_('N','N',mOrb(iSym),n_OK,nBas(iSym),
     &                           1.0d0,Work(ip_CC),mOx,
     &                                 Work(ip_X),nBx,
     &                           0.0d0,Work(ip_U),mOx)
         jOcc=iOff+nFro(iSym)+1
         Call Get_Nat_Lorb(OccN(jOcc),Work(ip_FOcc),n_OK,mOrb(iSym),
     &                     iWork(ip_iD+nBmx),Work(ip_U),iSym)
         nOx=Max(1,n_OK)
         Call DGEMM_('N','N',nBas(iSym),n_OK,n_OK,
     &                1.0d0,Work(ip_X),nBx,
     &                      Work(ip_U),nOx,
     &                0.0d0,Work(ipScr),nBx)
         Do i=1,n_OK
            kl=ipScr+nBas(iSym)*(i-1)
            j=kD(i)
            km=jOff+nBas(iSym)*(j-1)+1
            call dcopy_(nBas(iSym),Work(kl),1,Xmo(km),1)
         End Do

         Call DGEMM_('N','N',mOrb(iSym),n_KO,nBas(iSym),
     &                           1.0d0,Work(ip_CC),mOx,
     &                                 Work(iZ),nBx,
     &                           0.0d0,Work(ip_U),mOx)
         Call Get_Nat_Lorb(OccN(jOcc),Work(ip_FOcc),n_KO,mOrb(iSym),
     &                     iWork(ip_iD+nBmx+nOrbmx),Work(ip_U),iSym)
         nOx=Max(1,n_KO)
         Call DGEMM_('N','N',nBas(iSym),n_KO,n_KO,
     &                1.0d0,Work(iZ),nBx,
     &                      Work(ip_U),nOx,
     &                0.0d0,Work(ipScr),nBx)
         Do i=1,n_KO
            kl=ipScr+nBas(iSym)*(i-1)
            j=lD(i)
            km=jOff+nBas(iSym)*(j-1)+1
            call dcopy_(nBas(iSym),Work(kl),1,Xmo(km),1)
            k=jOcc+j-1
            l=ip_FOcc+j-1
            OccN(k)=Work(l)
         End Do
         Do i=1,n_OK
            j=kD(i)
            k=jOcc+j-1
            l=ip_FOcc+j-1
            OccN(k)=Work(l)
         End Do

         iOff=iOff+nBas(iSym)
         jOff=jOff+nBas(iSym)*mOrb(iSym)
         kOff=kOff+nBas(iSym)*(nBas(iSym)+1)/2
      End Do

      Call GetMem('LCMO','Free','Real',ip_C,(5*nBmx+nOrbmx+2)*nOrbmx)
      CALL GetMem('Smat','Free','Real',ipS,nnB+nBmx**2)
      Call GetMem('iD','Free','Inte',ip_iD,2*nOrbmx+nBmx)

      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Get_Nat_Lorb(Occ,FOcc,nO,nX,jOrb,Umat,iSym)

      Implicit Real*8 (a-h,o-z)
      Real*8  Occ(*), FOcc(*), Umat(*)
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
         Work(ii)=Occ(i)
      End Do
      nXx=Max(1,nX)
      nOx=Max(1,nO)
      Call DGEMM_('N','N',nX,nO,nX,1.0d0,Work(ip_eta),nXx,
     &                                  Umat(1),nXx,
     &                            0.0d0,Work(ip_Z),nXx)
      Call DGEMM_('T','N',nO,nO,nX,1.0d0,Umat(1),nXx,
     &                                  Work(ip_Z),nXx,
     &                            0.0d0,Work(ip_eta),nOx)

      Call Eigen_Molcas(nO,Work(ip_eta),Work(ip_Z),Work(ip_ZZ))

      call dcopy_(nO**2,Work(ip_eta),1,Umat(1),1)
      Do i=1,nO
         ii=ip_Z+i-1
         j=jOrb(i)
         FOcc(j)=Work(ii)
         kk=ip_ZZ+i-1
      End Do
      Call GetMem('eta_ik','Free','Real',ip_eta,2*nX**2+1)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSym)
      End
************************************************************************
