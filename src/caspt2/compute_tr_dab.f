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
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      integer(kind=iwp), intent(in):: nSym, nBas(nSym), nFro(nSym),
     &                                nIsh(nSym), nAsh(nSym),
     &                                nSsh(nSym), nDel(nSym)
      real(kind=wp), intent(in)::  CMO(*), OrbE(*)
      real(kind=wp), intent(out):: TrD(nSym)

      integer(kind=iwp) nAct(8), lnOrb(8), lnOcc(8), lnFro(8), lnDel(8),
     &                           lnVir(8)
      real(kind=wp), Allocatable:: EOrb(:), DMat(:), CMON(:)
      real(kind=wp) Dummy
      integer(kind=iwp) iE, ifr, ioff, ip_Y, irc, iSkip, iSym, ito, iV,
     &                  joff, k, kEOcc, kEVir, kfr, koff, kto, nBB, nOA,
     &                  nOrb, nVV
      real(kind=wp), external:: DDot_
*
      nAct(:)=0
      nVV=0
      nOrb=0
      Do iSym=1,nSym
         iE=1+nOrb+nFro(iSym)+nIsh(iSym)
         Do k=0,nAsh(iSym)-1
            If (OrbE(iE+k).lt.Zero) nAct(iSym)=nAct(iSym)+1
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
      Call mma_allocate(Eorb,2*nOrb,Label='EOrb')
      kEOcc=1
      kEVir=kEOcc+nOrb
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff+nFro(iSym)
         ito=kEOcc+joff
         call dcopy_(lnOcc(iSym),OrbE(ifr),1,EORb(ito),1)
         ifr=1+ioff+nFro(iSym)+nIsh(iSym)+nAsh(iSym)
         ito=kEVir+koff
         call dcopy_(nSsh(iSym),OrbE(ifr),1,EORb(ito),1)
         ioff=ioff+nBas(iSym)
         joff=joff+lnOcc(iSym)
         koff=koff+nSsh(iSym)
      End Do

      Call mma_allocate(Dmat,nVV+nOA,Label='DMat')
      ip_Y=1+nVV
      DMAT(:)=Zero
*
      Call LovCASPT2_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,.true.)
      Call mma_allocate(CMON,nBB,Label='CMON')
      CMON(:)=Zero
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=1+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr),1,CMON(kto),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym)+nAsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr),1,CMON(kto),1)
         iOff=iOff+nBas(iSym)**2
      End Do
*
      Call Check_Amp(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,CMON,EOrb(kEOcc),Eorb(kEVir),
     &                   DMAT(1:nVV),DMAT(ip_Y:))
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
      Call mma_deallocate(CMON)
*
      iV=1
      Do iSym=1,nSym
        TrD(iSym)=ddot_(lnVir(iSym),DMat(iV:),1+lnVir(iSym),[One],0)
        iV=iV+lnVir(iSym)**2
      End Do
      Call mma_deallocate(Dmat)
      Call mma_deallocate(Eorb)
*
      End Subroutine Compute_Tr_Dab
