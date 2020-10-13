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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
C#define _DEBUGPRINT_
      Subroutine ChoLSOSMP2_Energy(irc,EMP2,EOcc,EVir,Sorted,DelOrig)
C
C     Thomas Bondo Pedersen, December 2012.
C
C     Compute Laplace-SOS-MP2 energy.
C
      Implicit None
      Integer irc
      Real*8  EMP2
      Real*8  EOcc(*)
      Real*8  EVir(*)
      Logical Sorted
      Logical DelOrig
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"

      Character*17 SecNam
      Parameter (SecNam='ChoLSOSMP2_Energy')

      Integer  CheckDenomRange
      Integer  TestMinimaxLaplace
      External TestMinimaxLaplace

      Logical Debug
#if defined (_DEBUGPRINT_)
      Parameter (Debug=.true.)
#else
      Parameter (Debug=.false.)
#endif

      Logical Verb, FermiShift

      Integer ip_w, l_w
      Integer ip_t, l_t

      Real*8 ELOMO, EHOMO
      Real*8 ELUMO, EHUMO
      Real*8 EFermi
      Real*8 xmin, xmax

      Integer iSym
      Integer i

C====================
C     Initializations
C====================

      ! for now, this code is restricted to usage by the author
      Call ThisIsRestrictedCode('Thomas Bondo Pedersen',
     &                          'Laplace-SOS-MP2',.false.)

      ! init return code
      irc=0

      ! init flag for Fermi shift done
      FermiShift=.false.

      ! check that Laplace is requested
      If (.not.Laplace) Then
         Call WarningMessage(1,SecNam//' was called - '
     &                       //'but this is not a Laplace calculation!')
         Call xFlush(6)
         irc=-1
         Return
      End If

      ! Debug: test minimax Laplace grid generation
      If (Debug) Then
         Verb=.false.
         irc=TestMinimaxLaplace(1.0d-7,Verb)
         If (irc.ne.0) Then
            Call WarningMessage(2,SecNam
     &         //': error detected in numerical Laplace transformation')
            Write(6,'(A,I6)') 'irc=',irc
            Call Abend()
         End If
      End If

C====================================================
C     Parameters for numerical Laplace transformation
C====================================================

      ! get max and min orbital energies
      ELOMO=0.0d0
      EHOMO=0.0d0
      ELUMO=0.0d0
      EHUMO=0.0d0
      i=0
      Do iSym=1,nSym
         If (nOcc(iSym).gt.0) Then
            If (i.eq.0) Then
               i=1
               ELOMO=EOcc(iOcc(iSym)+1)
               EHOMO=EOcc(iOcc(iSym)+nOcc(iSym))
            Else
               ELOMO=min(ELOMO,EOcc(iOcc(iSym)+1))
               EHOMO=max(EHOMO,EOcc(iOcc(iSym)+nOcc(iSym)))
            End If
         End If
      End Do
      If (i.eq.0) Then
         Call WarningMessage(2,
     &                        SecNam//': unable to determine LOMO,HOMO')
         Call Abend()
      End If
      i=0
      Do iSym=1,nSym
         If (nVir(iSym).gt.0) Then
            If (i.eq.0) Then
               i=1
               ELUMO=EVir(iVir(iSym)+1)
               EHUMO=EVir(iVir(iSym)+nVir(iSym))
            Else
               ELUMO=min(ELUMO,EVir(iVir(iSym)+1))
               EHUMO=max(EHUMO,EVir(iVir(iSym)+nVir(iSym)))
            End If
         End If
      End Do
      If (i.eq.0) Then
         Call WarningMessage(2,
     &                        SecNam//': unable to determine LUMO,HUMO')
         Call Abend()
      End If
C-tbp:
      write(6,*) 'ELOMO,EHOMO=',ELOMO,EHOMO
      write(6,*) 'ELUMO,EHUMO=',ELUMO,EHUMO

      ! compute "Fermi energy" as the midpoint between HOMO and LUMO.
      EFermi=0.5d0*(EHOMO+ELUMO)

      ! translate orbital energy origin to EFermi
      Do iSym=1,nSym
         Do i=1,nOcc(iSym)
            EOcc(iOcc(iSym)+i)=EOcc(iOcc(iSym)+i)-EFermi
         End Do
         Do i=1,nVir(iSym)
            EVir(iVir(iSym)+i)=EVir(iVir(iSym)+i)-EFermi
         End Do
      End Do
      ELOMO=ELOMO-EFermi
      EHOMO=EHOMO-EFermi
      ELUMO=ELUMO-EFermi
      EHUMO=EHUMO-EFermi
      FermiShift=.true.

      ! compute range of orbital energy denominator
      xmin=2.0d0*(ELUMO-EHOMO)
      xmax=2.0d0*(EHUMO-ELOMO)
      ! Debug: check range
      If (Debug) Then
         irc=CheckDenomRange(xmin,xmax,nSym,EOcc,Evir,iOcc,nOcc,
     &                                                iVir,nVir)
         If (irc.ne.0) Then
            Call WarningMessage(2,SecNam
     &         //': error detected in orbital energy denominator range')
            Write(6,'(A,I6)') 'irc=',irc
            Call Abend()
         End If
      End If

      ! get weights and grid points for numerical Laplace transform
      If (Laplace_nGridPoints.eq.0) Then
         l_w=Laplace_mGridPoints
      Else
         l_w=Laplace_nGridPoints
      End If
      l_t=l_w
      Call GetMem('Lap_w','Allo','Real',ip_w,l_w)
      Call GetMem('Lap_t','Allo','Real',ip_t,l_t)
      Call MinimaxLaplace(Verbose,Laplace_nGridPoints,xmin,xmax,
     &                    l_w,Work(ip_w),Work(ip_t),irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I6)') SecNam,': MinimaxLaplace returned',irc
         irc=1
         Go To 1 ! exit after cleanup actions
      End If

C======================================
C     Compute SOS-MP2 energy correction
C======================================

      If (Sorted) Then
         Call ChoLSOSMP2_Energy_Srt(Laplace_nGridPoints,
     &                              Work(ip_w),Work(ip_t),
     &                              EOcc,EVir,
     &                              DelOrig,EMP2,irc)
         If (irc .ne. 0) Then
            Write(6,'(A,A,I6)')
     &      SecNam,': ChoLSOSMP2_Energy_Srt returned',irc
            Go To 1 ! exit
         End If
      Else
         If (nBatch .eq. 1) Then
            Call ChoLSOSMP2_Energy_Fll(Laplace_nGridPoints,
     &                                 Work(ip_w),Work(ip_t),
     &                                 EOcc,EVir,
     &                                 DelOrig,EMP2,irc)
            If (irc .ne. 0) Then
               Write(6,'(A,A,I6)')
     &         SecNam,': ChoLSOSMP2_Energy_Fll returned',irc
               Go To 1 ! exit after cleanup
            End If
         Else
            Call WarningMessage(1,
     &                        SecNam//': unsorted case not implemented')
            irc=-2
            Go To 1 ! exit after cleanup
         End If
      End If

C============
C     Cleanup
C============
    1 Continue ! errors jump to this point

      ! translate orbital energy origin from Fermi back to original
      If (FermiShift) Then
         Do iSym=1,nSym
            Do i=1,nOcc(iSym)
               EOcc(iOcc(iSym)+i)=EOcc(iOcc(iSym)+i)+EFermi
            End Do
            Do i=1,nVir(iSym)
               EVir(iVir(iSym)+i)=EVir(iVir(iSym)+i)+EFermi
            End Do
         End Do
      End If

      ! deallocations
      Call GetMem('Lap_t','Free','Real',ip_t,l_t)
      Call GetMem('Lap_w','Free','Real',ip_w,l_w)

      End
      Integer Function CheckDenomRange(xmin,xmax,nSym,EOcc,Evir,
     &                                 iOcc,nOcc,iVir,nVir)
      Implicit None
      Real*8  xmin
      Real*8  xmax
      Integer nSym
      Real*8  EOcc(nSym)
      Real*8  EVir(nSym)
      Integer iOcc(nSym)
      Integer nOcc(nSym)
      Integer iVir(nSym)
      Integer nVir(nSym)

      Real*8 Tol
      Parameter (Tol=1.0d-12)

      Real*8 e, emin, emax

      Integer iSym, i
      Integer aSym, a
      Integer irc

      emin=9.9d15
      emax=-9.9d15
      Do iSym=1,nSym
         Do i=iOcc(iSym)+1,iOcc(iSym)+nOcc(iSym)
            Do aSym=1,nSym
               Do a=iVir(aSym)+1,iVir(aSym)+nVir(aSym)
                  e=EVir(a)-EOcc(i)
                  emin=min(emin,e)
                  emax=max(emax,e)
               End Do
            End Do
         End Do
      End Do
      emin=2.0d0*emin
      emax=2.0d0*emax

      irc=0
      If (abs(emin-xmin).gt.Tol) irc=irc+1
      If (abs(emax-xmax).gt.Tol) irc=irc+2

C-tbp:
      if (irc.ne.0) then
         write(6,'(A,1P,2D25.16)') 'xmin,xmax=',xmin,xmax
         write(6,'(A,1P,2D25.16)') 'emin,emax=',emin,emax
         write(6,'(A,1P,2D25.16)') 'diff=     ',xmin-emin,xmax-emax
      end if

      CheckDenomRange=irc

      End
