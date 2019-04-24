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
      Subroutine TCMO(A,isym,ictl)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Integer ip(8)
      Real*8 A(*)
      Call GetMem('Temp','ALLO','REAL',ipT,nDens2)
      Call ReLoad(A,isym,norb,nbas)
* irc used in call later must not be uninitialized
* since then MKL library gets upset...
      irc=0
      If (ictl.eq.-1) Then
       Do iS=1,nSym
        If (nbas(is).ne.0) Then
        Call GetMem('CMOINV','ALLO','REAL',ip(iS),
     &              nBas(is)**2+nbas(is))
        call dcopy_(nbas(is)**2,Work(ipCMO+ipcm(is)-1),1,Work(ip(is)),1)
        iip=ip_of_iWork_d(Work(ip(is)+nbas(is)**2))
        call dgetrf_(nBas(is),nBas(is),Work(ip(is)),
     &              nBas(is),iWork(iip),irc)
        If (irc.ne.0)
     &         Call SysAbendMsg('tcmo','DGETRF returns non zero', ' ')
        end if
       End Do
        Do iS=1,nSym
         js=ieor(is-1,isym-1)+1
        If (nbas(is)*nbas(js).ne.0) Then
         iip=ip_of_iWork_d(Work(ip(is)+nbas(is)**2))
         call dgetrs_('T',nbas(is),nbas(js),
     &                Work(ip(is)),nBas(is),
     &                iWork(iip),
     &                A(ipMat(is,js)),nBas(is),irc)
         if (irc.ne.0)
     &         Call SysAbendMsg('tcmo','DGETRS returns non zero', ' ')
         Call DGETMO(A(ipMat(is,js)),nBas(is),
     &                    nbas(is),nbas(js),
     &                    Work(ipT),nbas(js))
         iip=ip_of_iWork_d(Work(ip(js)+nbas(js)**2))
         call dgetrs_('T',nbas(js),nbas(is),
     &                Work(ip(js)),nBas(js),
     &                iWork(iip),
     &                Work(ipT),nBas(js),irc)
         if (irc.ne.0)
     &         Call SysAbendMsg('tcmo','DGETRS returns non zero', ' ')
         Call DGETMO(Work(ipT),nBas(js),
     &              nbas(js),nbas(is),
     &              A(ipMat(is,js)),nbas(is))
        End If
        End Do
       Do iS=1,nSym
        If (nbas(is).ne.0)
     &  Call GetMem('CMOINV','FREE','REAL',ip(iS),idum)
       End Do
      Else If (ictl.eq.1) Then
       Do iS=1,nSym
        js=ieor(is-1,isym-1)+1
        If (nBas(is)*nBas(js).ne.0) Then
           Call DGEMM_('T','N',
     &                 nOrb(iS),nBas(jS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(is),
     &                 A(ipmat(is,js)),nBas(iS),
     &                 0.0d0,Work(ipT),nOrb(iS))
           Call DGEMM_('N','N',
     &                 nOrb(is),nOrb(jS),nBas(jS),
     &                 1.0d0,Work(ipT),nOrb(iS),
     &                 Work(ipCMO+ipCM(jS)-1),nBas(jS),
     &                 0.0d0,A(ipMat(iS,jS)),nOrb(iS))
        End If
       End Do
      Else if (ictl.eq.-2) Then
       Do iS=1,nSym
        js=ieor(is-1,isym-1)+1
        If (nBas(is)*nBas(js).ne.0) Then
           Call DGEMM_('N','N',
     &                 nBas(iS),nOrb(jS),nOrb(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(is),
     &                 A(ipmat(is,js)),nOrb(iS),
     &                 0.0d0,Work(ipT),nBas(iS))
           Call DGEMM_('N','T',
     &                 nBas(is),nBas(jS),nOrb(jS),
     &                 1.0d0,Work(ipT),nBas(iS),
     &                 Work(ipCMO+ipCM(jS)-1),nBas(jS),
     &                 0.0d0,A(ipMat(iS,jS)),nBas(iS))
        End If
       End Do
      Else
        Write(6,*) 'Oink'
        Call SysHalt('tcmo')
      End If
      Call GetMem('Temp','FREE','REAL',ipT,nDens2)
      Return
      End
