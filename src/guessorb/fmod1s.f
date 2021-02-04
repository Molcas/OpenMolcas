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
* Copyright (C) 2004, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* Absolutely NOT to be used for production!!!!                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: Oct 2004                                                    *
*                                                                      *
************************************************************************
      Subroutine Fmod1s(StandAlone)
      Implicit None
#include "Molcas.fh"
#include "stdalloc.fh"
#include "commgo.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Logical StandAlone
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Logical Debug
      Logical Trace
* Basis set indices
      Integer iSym
      Integer iBas
      Integer jBas
      Integer kBas
* Dimensions
      Integer nBasMax
      Integer nBasTot
      Integer nTriTot
      Integer nSqrTot
      Integer n2Full
* Pointers
      Integer ipTmp1
      Integer ipTmp2
      Integer ipTmp3
      Integer ipTmp4
* Matrix elements
      Real*8  Sii
      real*8  Sjj
      Real*8  Sik
      Real*8  Sjk
      Real*8  Skk
* Various variables
      Character*32 Line
      Integer indx
      Integer iSymlb
      Integer irc
      Real*8  Det
      Real*8  sum
      Real*8  eps
      Integer Lu
      Integer iDummy(7,8)
      Integer RC
      Character*80 Title
      Real*8, Dimension(:), Allocatable :: SmTr, DeTr, Esym, Edes, Smat,
     &        Ovl, Aux1, Sdes, Fdes, FSym, Fock, CMOs, Evec, Fmo, Aux2
*----------------------------------------------------------------------*
* Some setup                                                           *
*----------------------------------------------------------------------*
      If(StandAlone) Then
         Debug=.false.
         Trace=.false.
      Else
         Debug=.false.
         Trace=.false.
      End If
      If(Trace) Write(6,*) '>>> Entering fmod1s'
*----------------------------------------------------------------------*
* Should only be called if symmetry is used.                           *
*----------------------------------------------------------------------*
      If(nSym.eq.1) Then
         Call SysAbendMsg('fmod1s','internal error 001',' ')
      End If
*----------------------------------------------------------------------*
* Setup various counters.                                              *
*----------------------------------------------------------------------*
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
      n2Full=nBasTot*nBasTot
*----------------------------------------------------------------------*
* Get symmetry transformation matrix.                                  *
*----------------------------------------------------------------------*
      Call mma_allocate(SmTr,n2Full)
      Call mma_allocate(DeTr,n2Full)
      Call Get_dArray('SM',SmTr,n2Full)
      Call Minv(SmTr,DeTr,0,Det,nBasTot)
      If(Debug) Then
         Call RecPrt('Symmetrization transformation matrix',
     &               '(10f12.6)',
     &               SmTr,nBasTot,nBasTot)
         Call RecPrt('Desymmetrization transformation matrix',
     &               '(10f12.6)',
     &               DeTr,nBasTot,nBasTot)
      End If
*----------------------------------------------------------------------*
* Create model Fock operator.                                          *
*----------------------------------------------------------------------*
      Call mma_allocate(Esym,nBasTot)
      Call FockOper(RC,Esym)
      If(RC.ne.0) Then
         Call mma_deallocate(Esym)
         Call mma_deallocate(DeTr)
         Call mma_deallocate(SmTr)
         Return
      End If
      If(Debug) Then
         Call RecPrt('Diagonal of Fock','(10f12.6)',
     &               Esym,1,nBasTot)
      End If
*----------------------------------------------------------------------*
* Desymmetrize model Fock operator                                     *
*----------------------------------------------------------------------*
      Call mma_allocate(Edes,nBasTot)
      Do iBas=1,nBasTot
         Do kBas=1,nBasTot
            indx=kBas+nBasTot*(iBas-1)
            If(Abs(SmTr(indx)).gt.1.0d-3) Then
               Edes(kBas)=Esym(iBas)
            End if
         End Do
      End Do
      If(Debug) Then
         Call RecPrt('Desymmetrized diagonal of Fock','(10f12.6)',
     &               Edes,1,nBasTot)
      End If
*----------------------------------------------------------------------*
* Read and square overlap matrix.                                      *
*----------------------------------------------------------------------*
      Call mma_allocate(Smat,n2Full)
      Call mma_allocate(Ovl,nTriTot)
      iSymlb=1
      Call RdOne(irc,2,'Mltpl  0',1,Ovl,iSymlb)
      Call dCopy_(n2Full, [0.0d0],0, Smat,1)
      ipTmp1=1
      ipTmp2=1
      Do iSym=1,nSym
         Call Square(Ovl(ipTmp1),SMat(ipTmp2),1,nBasTot,nBas(iSym))
         ipTmp1=ipTmp1+nBas(iSym)*(nBas(iSym)+1)/2
         ipTmp2=ipTmp2+nBas(iSym)*nBasTot+nBas(iSym)
      End Do
      If(Debug) Then
         Call RecPrt('Full overlap matrix','(10f12.6)',
     &               Smat,nBasTot,nBasTot)
      End If
      Call mma_deallocate(Ovl)
*----------------------------------------------------------------------*
* Desymmetrize overlap matrix.                                         *
*----------------------------------------------------------------------*
      Call mma_allocate(Sdes,n2Full)
      Call mma_allocate(Aux1,n2Full)
      Call DGEMM_('N','N',
     &            nBasTot,nBasTot,nBasTot,
     &            1.0d0,Smat,nBasTot,
     &                  DeTr,nBasTot,
     &            0.0d0,Aux1,nBasTot)
      Call DGEMM_('T','N',
     &            nBasTot,nBasTot,nBasTot,
     &            1.0d0,DeTr,nBasTot,
     &                  Aux1,nBasTot,
     &            0.0d0,Sdes,nBasTot)
      If(Debug) Then
         Call RecPrt('Desymmetrized overlap matrix','(10f12.6)',
     &               Sdes,nBasTot,nBasTot)
      End If
      Call mma_deallocate(Aux1)
*----------------------------------------------------------------------*
* Build Fock matrix                                                    *
*----------------------------------------------------------------------*
      Call mma_allocate(Fdes,n2Full)
      Do iBas=1,nBasTot
         Sii=Sdes(iBas+nBasTot*(iBas-1))
         Do jBas=1,nBasTot
            Sjj=Sdes(jBas+nBasTot*(jBas-1))
            sum=0.0d0
            Do kBas=1,nBasTot
               eps=Edes(kBas)
               Skk=Sdes(kBas+nBasTot*(kBas-1))
               Sik=Sdes(iBas+nBasTot*(kBas-1))
               Sjk=Sdes(jBas+nBasTot*(kBas-1))
               sum=sum+eps*Sik*Sjk/Sqrt(Sii*Sjj*Skk*Skk)
            End Do
            Fdes(iBas+nBasTot*(jBas-1))=sum
         End Do
      End Do
      If(Debug) Then
         Call RecPrt('Desymmetrized Fock matrix','(10f12.6)',
     &               Fdes,nBasTot,nBasTot)
      End If
*----------------------------------------------------------------------*
* Symmetrize Fock matrix.                                              *
*----------------------------------------------------------------------*
      Call mma_allocate(Fsym,n2Full)
      Call mma_allocate(Aux1,n2Full)
      Call DGEMM_('N','N',
     &            nBasTot,nBasTot,nBasTot,
     &            1.0d0,Fdes,nBasTot,
     &                  SmTr,nBasTot,
     &            0.0d0,Aux1,nBasTot)
      Call DGEMM_('T','N',
     &            nBasTot,nBasTot,nBasTot,
     &            1.0d0,SmTr,nBasTot,
     &                  Aux1,nBasTot,
     &            0.0d0,Fsym,nBasTot)
      If(Debug) Then
         Call RecPrt('Symmetrized overlap matrix','(10f12.6)',
     &               Fsym,nBasTot,nBasTot)
      End If
      Call mma_deallocate(Aux1)
*----------------------------------------------------------------------*
* Compact the symmetrized Fock matrix.                                 *
*----------------------------------------------------------------------*
      Call mma_allocate(Fock,nTriTot)
      ipTmp1=1
      ipTmp2=1
      Do iSym=1,nSym
         indx=1
         Do iBas=1,nBas(iSym)
            Do jBas=1,iBas
               Fock(indx)=Fsym(iBas+nBasTot*(jBas-1))
               indx=indx+1
            End Do
         End Do
         If(Debug) Then
            Write(Line,'(a,i2)') 'Fock matrix, symmetry ',iSym
            Call TriPrt(Line,'(10f12.6)',Fock(ipTmp2),nBas(iSym))
         End If
         ipTmp1=ipTmp1+nBas(iSym)*nBasTot+nBas(iSym)
         ipTmp2=ipTmp2+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
*----------------------------------------------------------------------*
* Make ON basis                                                        *
*----------------------------------------------------------------------*
      Call mma_allocate(CMOs,nSqrTot)
      Call goLowdin(CMOs)
*----------------------------------------------------------------------*
* Transform Fock matrix to MO basis and diagonalize                    *
*----------------------------------------------------------------------*
      Call mma_allocate(Evec,nTriTot)
      Call mma_allocate(Fmo,nTriTot)
      Call mma_allocate(Aux1,nBasMax*nBasMax)
      Call mma_allocate(Aux2,nBasMax*nBasMax)

      ipTmp1=1
      ipTmp2=1
      ipTmp3=1
      ipTmp4=1
      Do iSym=1,nSym
         if(nBas(iSym).gt.0) Then
            Call Square(Fock(ipTmp1),Aux1,
     &                  1,nBas(iSym),nBas(iSym))
            Call DGEMM_('N','N',
     &                  nBas(iSym),nBas(iSym),nbas(iSym),
     &                  1.0d0,Aux1,nBas(iSym),
     &                        CMOs(ipTmp2),nBas(iSym),
     &                  0.0d0,Aux2,nBas(iSym))
            Call MxMt(CMOs(ipTmp2),nBas(iSym),1,
     &                Aux2,1,nBas(iSym),
     &                Fmo(ipTmp3),
     &                nBas(iSym),nBas(iSym))
            If(Debug) Then
               Write(Line,'(a,i2)')
     &            'MO Fock matrix, symmetry ',iSym
               Call TriPrt(Line,'(10f12.6)',Fmo(ipTmp3),nBas(iSym))
            End If
*           Call Jacob(Fmo(ipTmp3),CMOs(ipTmp2),nBas(iSym),nBas(iSym))
            Call NIdiag(Fmo(ipTmp3),CMOs(ipTmp2),
     &                  nBas(iSym),nBas(iSym),0)
            If(Debug) Then
               Write(Line,'(a,i2)')
     &            'Diagonal Fock matrix, symmetry ',iSym
               Call TriPrt(Line,'(10f12.6)',Fmo(ipTmp3),nBas(iSym))
            End If
            Call goPickup(Fmo(ipTmp3),Evec(ipTmp4),nBas(iSym))
            Call goSort(Evec(ipTmp4),CMOs(ipTmp2),nBas(iSym),nBas(iSym))
          End If
          ipTmp1=ipTmp1+nBas(iSym)*(nBas(iSym)+1)/2
          ipTmp2=ipTmp2+nBas(iSym)*nBas(iSym)
          ipTmp3=ipTmp3+nBas(iSym)*(nBas(iSym)+1)/2
          ipTmp4=ipTmp4+nBas(iSym)
      End Do
      Call mma_deallocate(Fmo)
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
*----------------------------------------------------------------------*
* Present data.                                                        *
*----------------------------------------------------------------------*
      If(PrintMOs) Then
         Call PriMO('Start orbitals',.false.,.true.,0.0d0,1.0d6,
     &              nSym,nBas,nBas,Label,Evec,Evec,CMOs,3)
      End If
      Call put_darray('Guessorb',CMOs,nSqrTot)
      Call put_darray('Guessorb energies',Evec,nBasTot)
      Call Put_iArray('nOrb',nBas,nSym)
      Call mma_allocate(Aux1,nBasTot)
      Do iBas=1,nBasTot
         Aux1(iBas)=0.0d0
      End Do
      Lu=20
      Title='Guess orbitals'
      Call WrVec('GSSORB',Lu,'COE',nSym,nBas,nBas,CMOs,
     &           Aux1,Evec,iDummy,Title)
      Call mma_deallocate(Aux1)
*----------------------------------------------------------------------*
* Done, deallocate the rest.                                           *
*----------------------------------------------------------------------*
      Call mma_deallocate(Evec)
      Call mma_deallocate(CMOs)
      Call mma_deallocate(Fock)
      Call mma_deallocate(Fsym)
      Call mma_deallocate(Fdes)
      Call mma_deallocate(Sdes)
      Call mma_deallocate(Smat)
      Call mma_deallocate(Edes)
      Call mma_deallocate(Esym)
      Call mma_deallocate(DeTr)
      Call mma_deallocate(SmTr)
      If(trace) Write(6,*) '<<< Exiting fmod1s'
      Return
      End
