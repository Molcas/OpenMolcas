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
* This routine creates a symmetric ON basis a la Lowdin                *
* Not true if orbitals are deleted.                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: Oct 2004                                                    *
*                                                                      *
************************************************************************
      Subroutine goLowdin(CMO)
      Implicit None
#include "Molcas.fh"
#include "stdalloc.fh"
#include "commgo.fh"
*----------------------------------------------------------------------*
* Dummy variables.                                                     *
*----------------------------------------------------------------------*
      Real*8 CMO(*)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Logical Debug
      Logical Trace
      Integer nBig
      Integer nTot
      Integer nTri
      Integer nTriTot
      Integer iSym
      Integer iBas
      Integer jBas
      Integer kBas
      Integer iOrb
      Integer ipOvl(8)
      Integer ipCMO
      Integer npSmat
      Integer irc
      Integer iSymlb
      Real*8 Temp, OrbPhase
      Real*8, Dimension(:), Allocatable :: Ovl, SMat, Vec, Eig
      Real*8, Dimension(:,:), Allocatable :: Tmp
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Debug=.false.
      Trace=.false.
      If(Trace) Write(6,*) '>>> Entering golowdin'
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      nBig=0
      nTot=0
      nTriTot=0
      Do iSym=1,nSym
         nTot=nTot+nBas(iSym)
         nTriTot=nTriTot+nBas(iSym)*(nBas(iSym)+1)/2
         If(nBig.lt.nBas(iSym)) nBig=nBas(iSym)
      End Do
*----------------------------------------------------------------------*
* Get overlap matrix                                                   *
*----------------------------------------------------------------------*
      npSmat=nBig*nBig
      Call mma_allocate(Ovl,nTriTot+4)
      ipOvl(1)=1
      Call mma_allocate(Smat,npSmat)
      iSymlb=1
      Call RdOne(irc,2,'Mltpl  0',1,Ovl,iSymlb)
      Do iSym=1,nSym-1
         ipOvl(iSym+1)=ipOvl(iSym)+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call mma_allocate(Vec,nBig**2)
      Call mma_allocate(Eig,nBig)

      ipCMO=1
      Do iSym=1,nSym
         nTri=nBas(iSym)*(nBas(iSym)+1)/2
         Call dCopy_(nTri,Ovl(ipOvl(iSym)),1,Smat,1)
*        If(Debug) Then
            Write(6,*)
            Write(6,*) '***'
            Write(6,*) '*** lowdin: symmetry',iSym
            Write(6,*) '***'
            Write(6,*)
            Call TriPrt('Overlap matrix','(12f18.12)',
     &                  Ovl(ipOvl(iSym)),nBas(iSym))
*        End If
         Call FZero(Vec,nBas(iSym)**2)
         Call DCopy_(nBas(iSym),1.0D0,0,Vec,nBas(iSym)+1)
         Call NIdiag_New(Ovl(ipOvl(iSym)),Vec,nBas(iSym),nbas(iSym),0)
*
         Do iBas = 1, nBas(iSym)
            temp=OrbPhase(Vec((iBas-1)*nBas(iSym)+1),nBas(iSym))
         End Do
*
*        If(Debug) Then
            Call RecPrt('Transformation','(12f18.12)',
     &                  Vec,nBas(iSym),nBas(iSym))
*        End If
         Call goPickUp(Ovl(ipOvl(iSym)),Eig,nBas(iSym))
*        If(Debug) Then
            Call RecPrt('Overlap eigenvalues before sort','(12f18.12)',
     &         Eig,1,nBas(iSym))
*        End If
         Do iBas=1,nBas(iSym)
            Eig(iBas)=-Eig(iBas)
         End Do
         Call goSort(Eig,Vec,nBas(iSym),nBas(iSym))
         Do iBas=1,nBas(iSym)
            Eig(iBas)=-Eig(iBas)
         End Do
*        If(Debug) Then
            Call RecPrt('Overlap eigenvalues after sort','(12f18.12)',
     &         Eig,1,nBas(iSym))
*        End If
         nDel(iSym)=0
         Do iBas=1,nBas(iSym)
            If(Eig(iBas).lt.SThr) nDel(iSym)=nDel(iSym)+1
         End Do
         Do iBas=1,nBas(iSym)
            Eig(iBas)=1.0d0/sqrt(Eig(iBas))
         End Do
         If(.false.) Then
            Do iBas=1,nBas(iSym)
               Do jBas=1,nBas(iSym)
                  Temp=0.0d0
                  Do kBas=1,nBas(iSym)
                     Temp=Temp + Eig(kBas)
     &                        *Vec(nBas(iSym)*(kBas-1)+iBas)
     &                        *Vec(nBas(iSym)*(kBas-1)+jBas)
                  End Do
                  CMO(ipCMO+nBas(iSym)*(jBas-1)+(iBas-1))=Temp
               End Do
            End Do
         Else
            Call dCopy_(nBas(iSym)*nBas(iSym),Vec,1,CMO(ipCMO),1)
            Do iOrb=1,nBas(iSym)
               Temp=Eig(iOrb)
               Do iBas=1,nBas(iSym)
                  CMO(ipCMO-1+iBas+nBas(iSym)*(iOrb-1))=
     &               Temp*CMO(ipCMO-1+iBas+nBas(iSym)*(iOrb-1))
               End Do
            End Do
         End If
*        If(Debug) Then
            Call RecPrt('Symmetric orbitals','(12f18.12)',
     &                  CMO(ipCMO),nBas(iSym),nBas(iSym))
*        End If
*        If(Debug) Then
            Call mma_allocate(Tmp,nBas(iSym),nBas(iSym))
            iSymlb=1
            Call RdOne(irc,2,'Mltpl  0',1,Ovl(ipOvl(1)),iSymlb)
            Do iBas=1,nBas(iSym)
               Do jBas=1,nBas(iSym)
                  Temp=0.0
                  Do kBas=1,nBas(iSym)
                     Temp=Temp+ CMO(ipCMO+nBas(iSym)*(kBas-1)+(iBas-1))
     &                       *CMO(ipCMO+nBas(iSym)*(jBas-1)+(kBas-1))
                  End Do
                  Tmp(iBas,jBas)=Temp
               End Do
            End Do
            Call RecPrt('Inverted overlap matrix','(12f18.12)',
     &                  Tmp,nBas(iSym),nBas(iSym))
            Call mma_deallocate(Tmp)
*        End If
         ipCMO=ipCMO+nBas(iSym)*nBas(iSym)
      End Do


      Call mma_deallocate(Eig)
      Call mma_deallocate(Vec)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call mma_deallocate(Smat)
      Call mma_deallocate(Ovl)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      If(Trace) Write(6,*) '<<< Exiting golowdin'
      Return
      End
