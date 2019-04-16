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
* Copyright (C) 1995, Martin Schuetz                                   *
************************************************************************
************************************************************************
* This Module contains subroutines which are used to compute info on   *
* the size of the SO integral symmetry blocks for direct integral      *
* transformation                                                       *
*                                                                      *
* SubRoutine SOFSh1                                                    *
*  -> compute (1) # SO functions in irrep for all shells iShell        *
*             (2) position of 1st component of shell in irrep for      *
*                 all shells in all irreps                             *
*             (3) map vector between shells and psedoshells for each   *
*                 irrep (not any shell contributes to all irreps)      *
*             (4) map vector between SO indices and shells             *
* relevant data is declared and passed in common "inftra.fh"           *
*----------------------------------------------------------------------*
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
************************************************************************
      SubRoutine SOFSh1(nShBF,iShOff,iSh2Sh,iSO2Sh,icntr,nSkal,nSym,
     &                  nSOs,nSD,iSD,nShIrp,nShBFmx)
*
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
      Integer nShBF(0:nSym-1,nSkal), iShOff(0:nSym-1,nSkal),
     &        iSh2Sh(0:nSym-1,nSkal), iSO2Sh(nSOs), iSD(0:nSD,nSkal),
     &        nShIrp(0:nSym-1),nShOff(0:7),icntr(nSkal)
      Dimension iTmp(1)
*
*     Call QEnter('SOFSh1')
*
      iTmp=0 ! Use iTmp to get round compiler bug on some machines
      Call ICopy(nSkal*nSym,iTmp,0,nShBF,1)
*
      iTmp=9999999
      Call ICopy(nSkal*nSym,iTmp,0,iShOff,1)
*
      Do irp=0,nsym-1
        nShOff(irp)=1
      End Do
*
      Do iSkal = 1, nSkal
         iShell      = iSD(11,iSkal)
         iAO         = iSD( 7,iSkal)
         iCmp        = iSD( 2,iSkal)
         icntr(iSkal)= iSD(10,iSkal)
*
*        loop over components of shell...
*
         Do i=1, iCmp
*           loop over irreps...
            Do irp=0, nSym-1
               If (iAnd(IrrCmp(IndS(iShell)+i),2**irp).ne.0) Then
                  nShBF(irp,iSkal) = nShBF(irp,iSkal)+ iSD(3,iSkal)
*define _CHECK_
#ifdef _CHECK_
                  If (Basis_Mode.eq.Auxiliary_Mode) Then
                     iShOff(irp,iSkal)=Min(iShOff(irp,iSkal),
     &                                     iAOtSO(iAO+i,irp)-nBas(irp))
                  Else
                     iShOff(irp,iSkal)=Min(iShOff(irp,iSkal),
     &                                     iAOtSO(iAO+i,irp))
                  End If
#endif
               End If
            End Do
         End Do
         Do irp=0,nSym-1
#ifdef _CHECK_
            If(nShBF(irp,iskal).ne.0   .and.
     >         nShOff(irp).ne.iShOff(irp,iSkal)) Then
               Call WarningMessage(2,'PROGRAMMING ERROR IN SHELL_SIZES')
               Write (6,*) nShBF(irp,iskal)
               Write (6,*) nShOff(irp)
               Write (6,*) iShOff(irp,iSkal)
               write(6,*) 'PROGRAMMING ERROR IN SHELL_SIZES: ',
     >                    'SHELLS NOT CONTIGUOUS. IRP=',irp,
     >                    '  ISKAL=',iskal
               Call Abend()
            End if
#endif
            iShOff(irp,iSkal)=nShOff(irp)
            nShOff(irp)=nShOff(irp)+nShBF(irp,iSkal)
         End Do
*debug
c        Write(6,'(A)') 'nShBF'
c        Write(6,'(8I4)') iSkal,(nShBF(irp,iSkal),irp=0,nIrrep-1)
c        Write(6,'(A)') 'iShOff'
c        Write(6,'(8I4)') iSkal,(iShOff(irp,iSkal),irp=0,nIrrep-1)
*debug
      End Do
*
*     and now set up SO-Shell and Shell-Psudoshell index vectors...
*
      iTmp=0 ! Use iTmp to get round compiler bug on some machines
      Call ICopy(nSym,iTmp,0,nShIrp,1)
*
      iTmp=-9999999
      Call ICopy(nSOs,iTmp,0,iSO2Sh,1)
      Call ICopy(nSkal*nSym,iTmp,0,iSh2Sh,1)
*
*     Loop over irreps...
*
      iptr=0
      nShBFMx=0
      Do irp=0, nSym-1
         Do iSkal=1, nSkal
*
            nShBFi=nShBF(irp,iSkal)
            nShBFMx=Max(nShBFMx,nShBFi)
            iSOb=iShOff(irp,iSkal)
            Do iSO=iSOb, iSOb+nShBFi-1
               If (iSO.gt.nSOs) Then
                  Call WarningMessage(2,' Fucked again!')
                  Call Quit_OnUserError()
               End If
               iSO2Sh(iptr+iSO)=iSkal
            End Do
*
            If (nShBFi.gt.0) Then
              nShIrp(irp)=nShIrp(irp)+1
              iSh2Sh(irp,iSkal)=nShIrp(irp)
            End If
         End Do
         If (Basis_Mode.eq.Auxiliary_Mode) Then
            iptr=iptr+nBas_Aux(irp)
         Else
            iptr=iptr+nBas(irp)
         End If
      End Do
*debug
c     Write(6,'(A,I4)') 'max shell size:',nShBFMx
c     Write(6,'(A)') '# of shells contributing to each irrep:'
c     Write(6,'(8I4)') (nShIrp(irp),irp=0,nSym-1)
c     Write(6,'(A)') '# shell-psudoshell map vector:'
c     Do irp=0, nSym-1
c       Write(6,'(A4,2X,I4,2X,A4,2X,8I4)') 'irp=',irp,'map:',
c    &             (iSh2Sh(irp,iSkal),iSkal=1,nSkal)
c     End Do
c     Write(6,'(A)') 'SO-index to shell map vector:'
c     Write(6,'(50I4)') (iSO2Sh(iSO),iSO=1,nSOs)
*debug
*     Call QExit('SOFSh1')
      Return
      End
