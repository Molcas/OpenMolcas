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
      use Arrays, only: CMO
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer ip(8), iip(8)
      Real*8 A(*)
      Real*8, Allocatable:: Temp(:)
      Real*8, Allocatable:: CMOInv(:)
      Integer*8, Allocatable :: iCMOInv(:)

      Call mma_allocate(Temp,nDens2,Label='Temp')
      Call ReLoad(A,isym,norb,nbas)

* irc used in call later must not be uninitialized
* since then MKL library gets upset...
      irc=0

      If (ictl.eq.-1) Then

       nCMOInv=0
       niCMOInv=0
       ip(:)=0
       iip(:)=0
       Do iS=1,nSym
          If (nbas(is)==0) Cycle
          ip(iS)=1+nCMOInv
          iip(iS)=1+nICMOInv
          nCMOInv=nCMOInv + nBas(is)**2
          niCMOInv=niCMOInv + nBas(is)
       End Do
       Call mma_allocate(CMOInv,nCMOInv,Label='CMOINV')
       Call mma_allocate(iCMOInv,niCMOInv,Label='iCMOINV')

       Do iS=1,nSym
        If (nbas(is)==0) Cycle
        call dcopy_(nbas(is)**2,CMO(ipcm(is)),1,
     &                          CMOINV(ip(is)),1)
        call dgetrf_(nBas(is),nBas(is),CMOINV(ip(is)),
     &              nBas(is),iCMOINV(iip(is)),irc)
        If (irc.ne.0)
     &         Call SysAbendMsg('tcmo','DGETRF returns non zero', ' ')
       End Do

       Do iS=1,nSym
         js=ieor(is-1,isym-1)+1
         If (nbas(is)*nbas(js)==0) Cycle
         call dgetrs_('T',nbas(is),nbas(js),
     &                CMOINV(ip(is)),nBas(is),
     &                iCMOINV(iip(is)),
     &                A(ipMat(is,js)),nBas(is),irc)
         if (irc.ne.0)
     &         Call SysAbendMsg('tcmo','DGETRS returns non zero', ' ')
         Call DGETMO(A(ipMat(is,js)),nBas(is),
     &                    nbas(is),nbas(js),
     &                    Temp,nbas(js))
         call dgetrs_('T',nbas(js),nbas(is),
     &                CMOINV(ip(js)),nBas(js),
     &                iCMOINV(iip(js)),
     &                Temp,nBas(js),irc)
         if (irc.ne.0)
     &         Call SysAbendMsg('tcmo','DGETRS returns non zero', ' ')
         Call DGETMO(Temp,nBas(js),
     &              nbas(js),nbas(is),
     &              A(ipMat(is,js)),nbas(is))
       End Do

       Call mma_deallocate(CMOInv)
       Call mma_deallocate(iCMOInv)

      Else If (ictl.eq.1) Then

       Do iS=1,nSym
        js=ieor(is-1,isym-1)+1
        If (nBas(is)*nBas(js)==0) Cycle
           Call DGEMM_('T','N',
     &                 nOrb(iS),nBas(jS),nBas(iS),
     &                 1.0d0,CMO(ipCM(iS)),nBas(is),
     &                 A(ipmat(is,js)),nBas(iS),
     &                 0.0d0,Temp,nOrb(iS))
           Call DGEMM_('N','N',
     &                 nOrb(is),nOrb(jS),nBas(jS),
     &                 1.0d0,Temp,nOrb(iS),
     &                 CMO(ipCM(jS)),nBas(jS),
     &                 0.0d0,A(ipMat(iS,jS)),nOrb(iS))
       End Do

      Else if (ictl.eq.-2) Then

       Do iS=1,nSym
        js=ieor(is-1,isym-1)+1
        If (nBas(is)*nBas(js)==0) Cycle
           Call DGEMM_('N','N',
     &                 nBas(iS),nOrb(jS),nOrb(iS),
     &                 1.0d0,CMO(ipCM(iS)),nBas(is),
     &                 A(ipmat(is,js)),nOrb(iS),
     &                 0.0d0,Temp,nBas(iS))
           Call DGEMM_('N','T',
     &                 nBas(is),nBas(jS),nOrb(jS),
     &                 1.0d0,Temp,nBas(iS),
     &                 CMO(ipCM(jS)),nBas(jS),
     &                 0.0d0,A(ipMat(iS,jS)),nBas(iS))
       End Do

      Else

        Write(6,*) 'Oink'
        Call SysHalt('tcmo')

      End If

      Call mma_deallocate(Temp)

      Return
      End
