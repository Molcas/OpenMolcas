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
      Subroutine Cho_ov_Loc(irc,Thrs,nSym,nBas,nFro,nIsh,
     &                               nAsh,nSsh,CMO,SMAT,iD_vir)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym), nAsh(nSym)
      Integer nIsh(nSym), nSsh(nSym), iD_vir(*)
      Real*8  Thrs, CMO(*), SMAT(*)
#include "WrkSpc.fh"

      irc=0
      l_Dens = 0
      Do iSym = 1,nSym
         l_Dens = max(l_Dens,nBas(iSym)**2)
      End Do
      Call GetMem('Density','Allo','Real',ip_Dens,2*l_Dens)
      ipD2=ip_Dens+l_Dens
      kOffC = 0
      jD = 1
      Do iSym = 1,nSym
         If (nIsh(iSym) .gt. 0) Then
            kOff1 = 1 + kOffC + nBas(iSym)*nFro(iSym)
            Call GetDens_Localisation(Work(ip_Dens),CMO(kOff1),
     &                                nBas(iSym),nIsh(iSym))
            Call FZero(CMO(kOff1),nBas(iSym)*nIsh(iSym))
            Call ChoLoc(irc,Work(ip_Dens),CMO(kOff1),Thrs,yNrm,
     &                      nBas(iSym),nIsh(iSym))
            If (irc .ne. 0) Then
               Call GetMem('Density','Free','Real',ip_Dens,2*l_Dens)
               irc  = 1
               Return
            End If
         End If
         Call izero(iD_vir(jD),nBas(iSym))
         If (nSsh(iSym) .gt. 0) Then
            kOff1= 1 + kOffC
            nOcc=nFro(iSym)+nIsh(iSym)+nAsh(iSym)
            Call GetDens_Localisation(Work(ip_Dens),CMO(kOff1),
     &                                nBas(iSym),nOcc)
            If ( nOcc+nSsh(iSym) .lt. nBas(iSym) ) Then  ! nDel > 0
               write(6,*) ' ******************************************'
               write(6,*) ' Cho_ov_Loc found Deleted orbitals in your '
               write(6,*) ' original MOs. She cannot properly handle  '
               write(6,*) ' this situation. The program may crash !! '
               write(6,*) ' ******************************************'
            EndIf
* compute -DS
            Call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),
     &                        -1.0d0,Work(ip_Dens),nBas(iSym),
     &                               SMAT(kOff1),nBas(iSym),
     &                         0.0d0,Work(ipD2),nBas(iSym))
* compute 1-DS = 1 + (-DS)
            Do i=0,nBas(iSym)-1
               ii=ipD2+nBas(iSym)*i+i
               Work(ii)=1.0d0+Work(ii)
            End Do
* compute (1-DS)*(1-DS)'
            Call GetDens_Localisation(Work(ip_Dens),Work(ipD2),
     &                                nBas(iSym),nBas(iSym))
            kOff2= kOff1+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
            Call FZero(CMO(kOff2),nBas(iSym)*nSsh(iSym))
            Call ChoLoc_xp(irc,Work(ip_Dens),CMO(kOff2),Thrs,yNrm,
     &                         nBas(iSym),nSsh(iSym),iD_vir(jD))
            If (irc .ne. 0) Then
               Call GetMem('Density','Free','Real',ip_Dens,2*l_Dens)
               irc  = 1
               Return
            End If
         End If
         kOffC = kOffC + nBas(iSym)**2
         jD = jD + nBas(iSym)
      End Do

      Call GetMem('Density','Free','Real',ip_Dens,2*l_Dens)
      Return
      End
