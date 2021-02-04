!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This is a very preliminary routine for testing purposes. It relies   *
! on the basis set being of ANO type.                                  *
!                                                                      *
!                                                                      *
! Absolutely NOT to be used for production!!!!                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************
      Subroutine Fmod1n(StandAlone)
      use GuessOrb_Global, only: Label, MxBasis, MxSym, nBas, nSym, PrintMOs
      Implicit None
#include "stdalloc.fh"
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
      Logical StandAlone
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Logical Debug
      Logical Trace
      Integer iSym
      Integer iBas
      Integer jBas
      Integer kBas
      Integer iOff
      Integer ipCMO(MxSym)
      Integer ipFock(MxSym)
      Integer ipX
      Integer ipY
      Integer iSymlb
      Integer iRc
      Integer nBasMax
      Integer nBasTot
      Integer nTriTot
      Integer nSqrTot
!---
      Real*8  orbene(MxBasis)
      Real*8, Dimension(:), Allocatable :: CMO, Fock, EVec, Ovl, Nrm
      Real*8, Dimension(:), Allocatable :: SFk, Hlf, TFk, Aux1
!     Character*4 Name(MxAtom)
!     Character*4 Label(2,MxBasis)
!---
      Integer i,k
      Real*8  Sik
      Real*8  Sjk
      Real*8  eps
      real*8  sum
      Integer Lu
      Integer iDummy(7,8)
      Integer RC
      Character*80 Title
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
      If(StandAlone) Then
         Debug=.false.
         Trace=.false.
      Else
         Debug=.false.
         Trace=.false.
      End If
      If(Trace) Write(6,*) '>>> Entering fmod1n'
!----------------------------------------------------------------------*
! Setup various counters.                                              *
!----------------------------------------------------------------------*
      nBasMax=0
      nBasTot=0
      nTriTot=0
      nSqrTot=0
      Do iSym=1,nSym
         If(nBasMax.lt.nBas(iSym)) nBasMax=nBas(iSym)
         nTriTot=nTriTot+nBas(iSym)*(nBas(iSym)+1)/2
         nSqrTot=nSqrTot+nBas(iSym)*nBas(iSym)
         nBasTot=nBasTot+nBas(iSym)
      End Do
!----------------------------------------------------------------------*
! Make symmetric orthonormal orbital basis.                            *
!----------------------------------------------------------------------*
      Call mma_allocate(CMO,nSqrTot)
      ipCMO(1)=1
      Do iSym=1,nSym-1
         ipCMO(iSym+1)=ipCMO(iSym)+nBas(iSym)*nBas(iSym)
      End Do
      Call goLowdin(CMO)
!----------------------------------------------------------------------*
! Allocate Fock matrix.                                                *
!----------------------------------------------------------------------*
      Call mma_allocate(Fock,nTriTot)
      ipFock(1)=1
      Call mma_allocate(EVec,nBasTot)
      Do iSym=1,nSym-1
         ipFock(iSym+1)=ipFock(iSym)+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
      Call dCopy_(nTriTot,[0.0d0],0,Fock,1)
!----------------------------------------------------------------------*
! Create model Fock operator.                                          *
!----------------------------------------------------------------------*
      Call FockOper(RC,EVec)
      If(RC.ne.0) Then
         Call mma_deallocate(EVec)
         Call mma_deallocate(Fock)
         Call mma_deallocate(CMO)
         Return
      End If
!----------------------------------------------------------------------*
! Get overlap matrix and make normalization.                           *
!----------------------------------------------------------------------*
      Call mma_allocate(Ovl,nTriTot+4)
      Call mma_allocate(Nrm,nBasTot)
      iSymlb=1
      Call RdOne(irc,2,'Mltpl  0',1,Ovl,iSymlb)

      ipX=1
      ipY=1
      Do iSym=1,nSym
         Call goPickup(Ovl(ipX),Nrm(ipY),nBas(iSym))
         Do iBas=1,nBas(iSym)
            Nrm(ipY-1+iBas)=1.0d0/sqrt(Nrm(ipY-1+iBas))
         End Do
         ipX=ipX+nBas(iSym)*(nBas(iSym)+1)/2
         ipY=ipY+nBas(iSym)
      End Do
!----------------------------------------------------------------------*
! Build Fock matrix.                                                   *
!----------------------------------------------------------------------*
      ipX=1
      ipY=1
      iOff=0
      Do iSym=1,nSym
         Do iBas=1,nBas(iSym)
            Do jBas=1,iBas
               sum=0.0d0
               Do kBas=1,nBas(iSym)
                  eps=EVec(iOff+kBas)
                  i=Max(iBas,kBas)
                  k=Min(iBas,kBas)
                  Sik=Ovl(ipX-1+i*(i-1)/2+k)*                           &
     &               Nrm(ipY-1+iBas)*Nrm(ipY-1+kBas)
                  i=Max(jBas,kBas)
                  k=Min(jBas,kBas)
                  Sjk=Ovl(ipX-1+i*(i-1)/2+k)*                           &
     &               Nrm(ipY-1+jBas)*Nrm(ipY-1+kBas)
                  sum=sum+eps*Sik*Sjk
               End Do
               Fock(ipFock(iSym)-1+iBas*(iBas-1)/2+jBas)=sum
            End Do
         End Do
         If(Debug) Then
            Call TriPrt('Modified atomic Fock matrix','(12f12.6)',      &
     &                  Fock(ipFock(iSym)),nBas(iSym))
         End if
         ipX=ipX+nBas(iSym)*(nBas(iSym)+1)/2
         ipY=ipY+nBas(iSym)
         iOff=iOff+nBas(iSym)
      End Do
!----------------------------------------------------------------------*
! Scale Fock matrix.                                                   *
!----------------------------------------------------------------------*
      ipX=1
      Do iSym=1,nSym
         If(Debug) Then
            Write(6,*) '***'
            Write(6,*) '*** Symmetry',iSym
            Write(6,*) '***'
         End If
         Do iBas=1,nBas(iSym)
            Do jBas=1,iBas
               Fock(ipFock(iSym)-1+iBas*(iBas-1)/2+jBas)=               &
     &         Fock(ipFock(iSym)-1+iBas*(iBas-1)/2+jBas)*               &
     &         Sqrt(Ovl(ipX-1+iBas*(iBas+1)/2))*                        &
     &         Sqrt(Ovl(ipX-1+jBas*(jBas+1)/2))
            End Do
         End Do
         If(Debug) Then
            Call TriPrt('Scaled atomic Fock matrix','(12f12.6)',        &
     &                  Fock(ipFock(iSym)),nBas(iSym))
         End if
         ipX=ipX+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
!----------------------------------------------------------------------*
! Release overlap and norm arrays.                                     *
!----------------------------------------------------------------------*
      Call mma_deallocate(Nrm)
      Call mma_deallocate(Ovl)
!----------------------------------------------------------------------*
! Diagonalize Fock matrix.                                             *
!----------------------------------------------------------------------*
      Call mma_allocate(SFk,nBasMax*nBasMax)
      Call mma_allocate(Hlf,nBasMax*nBasMax)
      Call mma_allocate(TFk,nBasMax*(nBasMax+1)/2)
      iOff=0
      Do iSym=1,nSym
         If(Debug) Then
            Write(6,*) '***'
            Write(6,*) '*** Symmetry',iSym
            Write(6,*) '***'
         End If
         If(nBas(iSym).gt.0) Then
            Call Square(Fock(ipFock(iSym)), SFk,                        &
     &                  1,nBas(iSym), nBas(iSym))
            Call DGEMM_('N','N',                                        &
     &                  nBas(iSym),nBas(iSym),nBas(iSym),               &
     &                  1.0d0,SFk,nBas(iSym),                           &
     &                        CMO(ipCMO(iSym)),nBas(iSym),              &
     &                  0.0d0,Hlf,nBas(iSym))
            Call MxMt(CMO(ipCMO(iSym)),nBas(iSym),1,                    &
     &                Hlf,1,nBas(iSym),                                 &
     &                TFk,                                              &
     &                nBas(iSym),nBas(iSym))
            If(Debug) Then
               Call TriPrt('Transformed Fock matrix','(12f12.6)',       &
     &                     TFk,nBas(iSym))
            End If
         End If

!        Call Jacob(TFk,CMO(ipCMO(iSym)),nbas(iSym),nbas(iSym))
         Call NIdiag(TFk,CMO(ipCMO(iSym)),                              &
     &               nbas(iSym),nbas(iSym),0)
         If(Debug) Then
            Call TriPrt('Diagonalized atomic Fock matrix','(12f12.6)',  &
     &                  TFk,nBas(iSym))
         End If
         Call goPickup(TFk,orbene(iOff+1),nBas(iSym))
         Call goSort(orbene(iOff+1),CMO(ipCMO(iSym)),                   &
     &             nBas(iSym),nBas(iSym))
         iOff=iOff+nBas(iSym)
      End Do
      Call mma_deallocate(TFk)
      Call mma_deallocate(Hlf)
      Call mma_deallocate(SFk)
!----------------------------------------------------------------------*
! Present data.                                                        *
!----------------------------------------------------------------------*
      If(PrintMOs) then
         Call PriMO('Start orbitals',.false.,.true.,0.0d0,1.0d6,        &
     &              nSym,nBas,nBas,Label,orbene,orbene,CMO(ipCMO(1)),3)
      End If
      Call put_darray('Guessorb',CMO(ipCMO(1)),nSqrTot)
      Call put_darray('Guessorb energies',orbene,nBasTot)
      Call Put_iArray('nOrb',nBas,nSym)
      Call mma_allocate(Aux1,nBasTot)
      Do iBas=1,nBasTot
         Aux1(iBas)=0.0d0
      End Do
      Lu=20
      Title='Guess orbitals'
      Call WrVec('GSSORB',Lu,'COE',NSYM,NBAS,NBAS,CMO(ipCMO(1)),        &
     &           Aux1,orbene,iDummy,Title)
      Call mma_deallocate(Aux1)
!----------------------------------------------------------------------*
! Done, deallocate the rest.                                           *
!----------------------------------------------------------------------*
      Call mma_deallocate(Evec)
      Call mma_deallocate(Fock)
      Call mma_deallocate(CMO)
      If(trace) Write(6,*) '<<< Exiting fmod1n'
      Return
      End
