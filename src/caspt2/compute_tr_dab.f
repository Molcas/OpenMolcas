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
      Subroutine Compute_Tr_Dab(nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &                          CMO,OrbE,TrD)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym), nIsh(nSym)
      Integer nAsh(nSym), nSsh(nSym), nDel(nSym)
      Real*8  CMO(*), OrbE(*), TrD(nSym)
#include "WrkSpc.fh"
      Integer nAct(8), lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
*
*
      Call Izero(nAct,nSym)
      nVV=0
      nOrb=0
      Do iSym=1,nSym
         iE=1+nOrb+nFro(iSym)+nIsh(iSym)
         Do k=0,nAsh(iSym)-1
            If (OrbE(iE+k).lt.0.0d0) nAct(iSym)=nAct(iSym)+1
         End Do
         nVV=nVV+nSsh(iSym)**2
         nOrb=nOrb+nBas(iSym)
      End Do
*
      nBB=0
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnOrb(iSym)=nBas(iSym)
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)+nAct(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnDel(iSym)=nDel(iSym)
         nBB=nBB+nBas(iSym)**2
         nOA=nOA+lnOcc(iSym)
      End Do
*
      Call GetMem('EOV','Allo','Real',ipEorb,2*nOrb)
      kEOcc=ipEorb
      kEVir=kEOcc+nOrb
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff+nFro(iSym)
         ito=kEOcc+joff
         call dcopy_(lnOcc(iSym),OrbE(ifr),1,Work(ito),1)
         ifr=1+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),OrbE(ifr),1,Work(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+lnOcc(iSym)
         koff=koff+nSsh(iSym)
      End Do

      Call GetMem('Dmat','Allo','Real',ip_X,nVV+nOA)
      ip_Y=ip_X+nVV
      Call FZero(Work(ip_X),nVV+nOA)
*
      Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,
     &                           ip_Y,.true.)
      Call GetMem('CMON','Allo','Real',iCMO,nBB)
      Call FZero(Work(iCMO),nBB)
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=iCMO+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,Work(kto),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,Work(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
*
      Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,Work(iCMO),Work(kEOcc),Work(kEVir))
         If(irc.ne.0) then
           write(6,*) 'MP2 pseudodensity calculation failed !'
           Call Abend
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      Endif
      Call GetMem('CMON','Free','Real',iCMO,nBB)
*
      iV=ip_X
      Do iSym=1,nSym
        TrD(iSym)=ddot_(lnVir(iSym),Work(iV),1+lnVir(iSym),[1.0d0],0)
        iV=iV+lnVir(iSym)**2
      End Do
      Call GetMem('Dmat','Free','Real',ip_X,nVV+nOA)
      Call GetMem('EOV ','Free','Real',ipEorb,2*nOrb)
*
      Return
      End
