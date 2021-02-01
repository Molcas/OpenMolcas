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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine PrIte(QNR,CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)
************************************************************************
*                                                                      *
*     purpose: Print out informations in every iteration               *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(mBB,nD), Ovrlp(mBT), OccNo(mmB,nD)
      Logical QNR
      Real*8 :: Shift=0.0D0
      Logical :: Set_Shift=.False.
      Save Shift, Set_Shift
*
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
      character cEDiff, cDMOMax, cFMOMax,cDltNrm

      If(iterprlv.gt.0) Then
         Write(6,*)
         Write(6,'(a)') '*******************'
         Write(6,'(a,i3,a)') '** Iteration ',iter,' **'
         Write(6,'(a)') '*******************'
         Write(6,*)
         Write(6,'(a,f10.2)') 'Cpu time [sec]        ',CpuItr
         If(AccCon.eq.'None' .or. AccCon.eq.'NoneDa') Then
            Write(6,'(a)') 'No convergence acceleration'
         Else If(AccCon.eq.'EDIIS'.or.AccCon.eq.'ADIIS') Then
            Write(6,'(a)') 'Convergence is accelerated by damping'
         Else If(AccCon.eq.'QNRc1D') Then
            Write(6,'(2a)') 'Convergence is accelerated by QNR with ',
     &                     'c1-DIIS'
         Else If(AccCon.eq.'QNRc2D') Then
            Write(6,'(2a)') 'Convergence is accelerated by QNR with ',
     &                     'c2-DIIS'
         Else
            Write(6,'(2a)') 'Convergence accelerations is ',AccCon
         End If
         Write(6,*)
         Write(6,'(a,f16.8)') 'Total energy          ',EneV
         Write(6,'(a,f16.8)') 'One electron energy   ',E1V
         Write(6,'(a,f16.8)') 'Two electron energy   ',E2V
         If(Abs(Ediff).gt.Ethr .or. iter.le.1) Then
            Write(6,'(a,f16.8)') 'Energy difference     ',Ediff
         Else
            Write(6,'(a,f16.8,a)') 'Energy difference     ',Ediff,
     &                             ' is converged'
         End If
         If(QNR) Then
            If(DltNrm.gt.DltNth) Then
               Write(6,'(a,f16.8)') 'Delta norm            ',DltNrm
            Else
               Write(6,'(a,f16.8,a)') 'Delta norm            ',DltNrm,
     &                                ' is converged'
            End If
         Else
            If(Abs(DMOMax).gt.Dthr) Then
               Write(6,'(a,f16.8)') 'Max offdiagonal Dij   ',DMOmax
            Else
               Write(6,'(a,f16.8,a)') 'Max offdiagonal Dij   ',DMOmax,
     &                                ' is converged'
            End If
         End If
         If(Abs(FMOMax).gt.Fthr) Then
            Write(6,'(a,f16.8)') 'Max offdiagonal Fij   ',FMOmax
         Else
            Write(6,'(a,f16.8,a)') 'Max offdiagonal Fij   ',FMOmax,
     &                             ' is converged'
         End If
         Write(6,'(a,f16.8)') 'D-norm                ',Sqrt(Dnorm)
         Write(6,'(a,f16.8)') 'T-norm                ',Sqrt(Tnorm)
         If(iterprlv.ge.2) Then
            Call MulPop(CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)
         End If
      Else If(jPrint.ge.2) Then
         If (.Not.Set_Shift) Then
            If (ABS(EneV).gt.1.0D3) Then
               Shift=DBLE(INT(Abs(EneV)/1.0D3))*1.0D3
               Write(6,*)
               Write(6,'(1X,A,f10.0,A)') 'The total and one-electron'//
     &                                 ' energies are shifted'//
     &                                 ' by a value of ',Shift,' a.u.'
               Write(6,*)
            End If
            Set_Shift=.True.
         End If
         cEDiff=' '
         If (Abs(Ediff).gt.Ethr) cEDiff='*'
         cFMOMax=' '
         If (Abs(FMOMax).gt.Fthr) cFMOMax='*'
         If (QNR) Then
            cDltNrm=' '
            If (DltNrm.gt.DltNth) cDltNrm='*'
            Write(6,'(1X,i3,3f16.9,1x,3(e10.2,a1,1x),2e11.2,3x,A,f6.0)')
     &            Iter,EneV+Shift,E1V+Shift,
     &            E2V,EDiff,cEDiff,DltNrm,cDltNrm,
     &            FMOMax,cFMOMax,
     &            Sqrt(DNorm),Sqrt(TNorm),
     &            AccCon,CpuItr
         Else
            cDMOMax=' '
            If (Abs(DMOMax).gt.Dthr) cDMOMax='*'
            Write(6,'(1X,i3,3f16.9,1x,3(e10.2,a1,1x),2e11.2,3x,A,f6.0)')
     &            Iter,EneV+Shift,E1V+Shift,
     &            E2V,EDiff,cEDiff,DMOMax,cDMOMax,
     &            FMOMax,cFMOMax,
     &            Sqrt(DNorm),Sqrt(TNorm),
     &            AccCon,CpuItr

         End If
      End If

      End subroutine prite
