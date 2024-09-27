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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991,1993,1996, Markus P. Fuelscher                    *
************************************************************************
      Subroutine Fmat_m(CMO,PUVX,D,D1A,FI,FA)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Update the Fock matrix for the active orbitals and transform     *
*     it to MO basis as well as the matrix FI (Fock matrix) for        *
*     frozen and inactive orbitals).                                   *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : array of real*8                                        *
*               MO-coefficients                                        *
*     PUVX    : array of real*8                                        *
*               two-electron integrals (pu!vx)                         *
*     D       : array of real*8                                        *
*               averaged one-body density matrix                       *
*     D1A     : array of real*8                                        *
*               active one body density matrix in AO-basis             *
*     FI      : array of real*8                                        *
*               inactive Fock matrix                                   *
*     FA      : array of real*8                                        *
*               active Fock matrix                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     B.O. Roos and P.Aa. Malmqvist                                    *
*     University of Lund, Sweden, 1989                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - updated for MOLCAS version 2                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1991               *
*     - updated for MOLCAS version 3                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1993               *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      use printlevel, only: debug
      use mcpdft_output, only: lf, iPrLoc
      use stdalloc, only: mma_allocate, mma_deallocate

      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
      Character(LEN=16), Parameter:: ROUTINE='FMAT    '

      Real*8 CMO(*) , PUVX(*) , D(*) , D1A(*) , FI(*) , FA(*)

      Real*8, allocatable:: Tmp1(:), Tmp2(:)


C Local print level (if any)
      IPRLEV=IPRLOC(4)
      !iPrLev=DEBUG-1
      If ( iPrLev.ge.DEBUG ) then
        write(lf,*) repeat('*',65)
        write(lf,*) 'Entering FMAT routine called by MSCTL!'
        write(lf,*) repeat('*',65)
        write(lf,*) 'printing input matrices :'
        write(lf,*) repeat('*',65)
        Write(LF,*)
        Write(LF,*) ' CMOs in FMAT'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
           iBas = nBas(iSym)
           call wrtmat(CMO(ioff),iBas,iBas, iBas, iBas)
           iOff = iOff + iBas*iBas
         End Do

         Write(LF,*)
         Write(LF,*) ' PUVX in FMAT'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         call wrtmat(PUVX,1,nFint, 1, nFint)

         Write(LF,*)
         Write(LF,*) ' ---------------------'
        CALL TRIPRT('Averaged one-body density matrix D, in MO in FMAT',
     &              ' ',D,NAC)

         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in FMAT'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(D1A(iOff),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End DO

         Write(LF,*)
         Write(LF,*) ' FI in AO-basis in FMAT'
         Write(LF,*) ' --------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
           iOrb = nOrb(iSym)
           Call TriPrt(' ',' ',FI(iOff),iOrb)
           iOff = iOff + (iOrb*iOrb+iOrb)/2
         End Do

         Write(LF,*)
         Write(LF,*) ' FA in AO-basis in FMAT'
         Write(LF,*) ' --------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
           iOrb = nOrb(iSym)
           Call TriPrt(' ',' ',FA(iOff),iOrb)
           iOff = iOff + (iOrb*iOrb+iOrb)/2
         End Do
       End If

*************************************************************
* Here we should start the real work!
*************************************************************
C Local print level (if any)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*     create FA in AO basis
      Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
      Call Fold(nSym,nBas,D1A,Tmp1)

*     Inactive-active contribution to ECAS
      VIA=dDot_(nTot1,FI,1,Tmp1,1)
      ECAS=EMY+VIA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A,ES20.10)') ' Total core energy:            ',EMY
        Write(LF,'(A,ES20.10)') ' inactive-active interaction:  ',VIA
        Write(LF,'(A,ES20.10)') ' CAS energy (core+interaction):',ECAS
      End If
      Call mma_deallocate(Tmp1)

*     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in AO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(LF,*)
        Write(LF,*) ' FA in AO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

*     transform FI from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        If (iBas==0) Cycle
        iOrb = nOrb(iSym)
        If (iOrb==0) Cycle
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
        Call Square(FI(iOff1),Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
        Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                 1.0D0,Tmp2,iBas,
     &                       CMO(iOff2+(iFro*iBas)),iBas,
     &                 0.0D0,FI(iOff3),iOrb)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

*     transform FA from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        If (iBas==0) Cycle
        iOrb = nOrb(iSym)
        If (iOrb==0) Cycle
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
        Call Square(FA(iOff1),Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
        Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                 1.0D0,Tmp2,iBas,
     &                       CMO(iOff2+(iFro*iBas)),iBas,
     &                 0.0D0,FA(iOff3),iOrb)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

!***********************************************************************
!     update Fock matrix
      Call Upd_FA_m(PUVX,FA,D,ExFac)

!     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(LF,*)
        Write(LF,*) ' FA in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

      End Subroutine Fmat_m
