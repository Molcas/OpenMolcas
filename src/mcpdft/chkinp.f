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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine ChkInp_m()
************************************************************************
*                                                                      *
*     Check the input for obvious errors or violation of limits        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use mcpdft_output, only: lf
      use rasscf_global, only: lRoots, NAC, NIN, nRoots,
     &                         OutFmt2, iRoot
      use definitions,only:iwp
      implicit none

#include "rasdim.fh"
#include "general.fh"
#include "warnings.h"

      integer(kind=iwp) :: ierr, i, ia0, ib0, ic0, ierr1
      integer(kind=iwp) :: ierr2, iSym

      IERR=0
      If ( NTOT.gt.mxOrb ) Then
        Write(LF,*)
        Write(LF,*)          '************ ERROR ******************'
        Call WarningMessage(2,'Too many orbitals.')
        Write(LF,'(1X,A,I8)')'Too many orbitals NTOT=',NTOT
        Write(LF,'(1X,A,I8)')'Limit is MXORB=',MXORB
        Write(LF,*)          '*************************************'
        call Quit_OnUserError()
      End If
      If ( NAC.gt.mxAct ) then
        Write(LF,*)
        Write(LF,*)          '*************** ERROR ***************'
        Call WarningMessage(2,'Too many active orbitals.')
        Write(LF,'(1X,A,I8)')'Too many active orbitals NAC=',NAC
        Write(LF,'(1X,A,I8)')'Limit is MXACT=',MXACT
        Write(LF,*)          '*************************************'
        call Quit_OnUserError()
      Endif
      If ( NIN.gt.mxIna ) then
        Write(LF,*)
        Write(LF,*)          '*************** ERROR ***************'
        Call WarningMessage(2,'Too many inactive orbitals.')
        Write(LF,'(1X,A,I8)')'Too many inactive orbitals NIN=',NIN
        Write(LF,'(1X,A,I8)')'Limit is MXINA=',MXINA
        Write(LF,*)          '*************************************'
        call Quit_OnUserError()
      Endif
      If(NACTEL.gt.2*NAC) then
        Write(LF,*)
        Write(LF,*)'********************* ERROR **********************'
        Call WarningMessage(2,'Too many active electrons.')
        Write(LF,'(1X,A,I6)')'Too many active electrons NACTEL=',NACTEL
        Write(LF,'(1X,A,I6)')
     &             'Cannot be more than 2*Nr of active orbitals=',2*NAC
        Write(LF,*)'**************************************************'
        call Quit_OnUserError()
      Endif
      If(NHOLE1.gt.2*NRS1T) then
        Write(LF,*)
        Write(LF,*)'******************** WARNING *********************'
        Call WarningMessage(2,'Too many holes in Ras1.')
        Write(LF,'(1X,A,I6)')
     &             'You allow too many holes in Ras1 NHOLE1=',NHOLE1
        Write(LF,'(1X,A,I6)')
     &             'Cannot be more than 2*Nr of Ras1 orbitals=',2*NRS1T
        Write(LF,*)'**************************************************'
        call Quit_OnUserError()
      Endif
      If(NELEC3.gt.2*NRS3T) then
        Write(LF,*)
        Write(LF,*)'******************** WARNING *********************'
        Call WarningMessage(1,'Too many electrons in Ras3.')
        Write(LF,'(1X,A,I6)')
     &             'You allow too many electrons in Ras3 NELEC3=',NELEC3
        Write(LF,'(1X,A,I6)')
     &             'Cannot be more than 2*Nr of Ras3 orbitals=',2*NRS3T
        NELEC3=2*NRS3T
        Write(LF,'(1X,A,I6)')
     &             'NELEC3 has been reset to ',NELEC3
        Write(LF,*)'**************************************************'
      Endif
        If(NACTEL.gt.2*(NRS1T+NRS2T)+NELEC3) then
         Write(LF,*)
         Write(LF,*)'********************* ERROR **********************'
         Call WarningMessage(2,'Too many active electrons.')
         Write(LF,'(1X,A,I8)')
     &              'Too many active electrons NACTEL=',NACTEL
         Write(LF,'(1X,A,I8)')
     &              '(Incompatible with RAS restrictions).'
         Write(LF,*)'**************************************************'
          call Quit_OnUserError()
        End If
        If(NACTEL.lt.2*NRS1T-NHOLE1) then
         Write(LF,*)
         Write(LF,*)'********************* ERROR **********************'
         Call WarningMessage(2,'Too few active electrons.')
         Write(LF,'(1X,A,I8)')
     &              'Too few active electrons NACTEL=',NACTEL
         Write(LF,'(1X,A,I8)')
     &              '(Incompatible with RAS restrictions).'
         Write(LF,*)'**************************************************'
          call Quit_OnUserError()
        End If

      If (NSYM.ne.1 .and. NSYM.ne.2 .and.
     &    NSYM.ne.4 .and. NSYM.ne.8) Then
        Write(LF,*)
        Call WarningMessage(2,'Wrong nr of symmetries.')
        Write(LF,*)          '*************** ERROR ***************'
        Write(LF,'(1X,A,I8)')'Nr of symmetries NSYM=',NSYM
        Write(LF,*)          ' Only possible values are 1,2,4 or 8.'
        Write(LF,*)          '*************************************'
        call Quit_OnUserError()
      End If

      IERR1=0
      Do iSym=1,nSym
         If ( nBas(iSym).lt.0 ) IERR1=1
         If ( nFro(iSym).lt.0 ) IERR1=1
         If ( nDel(iSym).lt.0 ) IERR1=1
         If ( nOrb(iSym).lt.0 ) IERR1=1
         If ( nISh(iSym).lt.0 ) IERR1=1
         If ( nASh(iSym).lt.0 ) IERR1=1
         If ( nSSh(iSym).lt.0 ) IERR1=1
         If ( nRS1(iSym).lt.0 ) IERR1=1
         If ( nRS2(iSym).lt.0 ) IERR1=1
         If ( nRS3(iSym).lt.0 ) IERR1=1
      End Do
      IERR2=0
      Do iSym=1,nSym
         If ( nBas(iSym).gt.mxBas ) IERR2=1
         If ( nFro(iSym).gt.mxBas ) IERR2=1
         If ( nDel(iSym).gt.mxBas ) IERR2=1
         If ( nOrb(iSym).gt.mxBas ) IERR2=1
         If ( nISh(iSym).gt.mxBas ) IERR2=1
         If ( nASh(iSym).gt.mxBas ) IERR2=1
         If ( nSSh(iSym).gt.mxBas ) IERR2=1
         If ( nRS1(iSym).gt.mxBas ) IERR2=1
         If ( nRS2(iSym).gt.mxBas ) IERR2=1
         If ( nRS3(iSym).gt.mxBas ) IERR2=1
      End Do
      If (IERR1+IERR2.gt.0) Then
        Write(LF,*)
        Write(LF,*)'****************** ERROR *******************'
        Call WarningMessage(2,'Erroneous nr of orbitals.')
        Write(LF,*)'Inappropriate nr of orbitals. One or more of'
        Write(LF,*)'these orbital counts is wrong or too large.'
        If (IERR1.gt.0) Write(LF,*)' Negative values.'
        If (IERR2.gt.0) Write(LF,*)' Extremely large values.'
        Write(LF,'(1X,A,8I4)')
     &   '   All orbitals:',(NBAS(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '         Frozen:',(NFRO(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '       Inactive:',(NISH(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '         Active:',(NASH(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '          RAS-1:',(NRS1(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '          RAS-2:',(NRS2(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '          RAS-3:',(NRS3(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '      Secondary:',(NSSH(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   '        Deleted:',(NDEL(ISYM),ISYM=1,NSYM)
        Write(LF,'(1X,A,8I4)')
     &   'Basis functions:',(NBAS(ISYM),ISYM=1,NSYM)
        If (IERR1.gt.0) Then
         Write(LF,*)' Have you used a too small basis set?'
        End If
        Write(LF,*)'********************************************'
        Call Quit_OnUserError()
      End If

      IERR=0
      IB0=ISPIN-1
      IA0=(NACTEL-IB0)/2
      IC0=(NAC-IA0-IB0)
      If ( (2*IA0+IB0).ne.NACTEL ) IERR=1
      If ( IA0.lt.0 ) IERR=1
      If ( IB0.lt.0 ) IERR=1
      If ( IC0.lt.0 ) IERR=1
      If (IERR.eq.1) Then
        Write(LF,*)
        Write(LF,*)'******************** ERROR *********************'
        Call WarningMessage(2,'No such wave function.')
        Write(LF,*)'The following combined specifications are wrong.'
        Write(LF,'(1X,A,8I4)')'Nr of active electrons NACTEL=',NACTEL
        Write(LF,'(1X,A,8I4)')'Nr of active orbitals     NAC=',NAC
        Write(LF,'(1X,A,8I4)')'Spin degeneracy         ISPIN=',ISPIN
        Write(LF,*)'There can be no such wave function.'
        Write(LF,*)'************************************************'
        Call Quit_OnUserError()
      End If

      IERR=0
      If ( NROOTS.gt.mxRoot ) IERR=1
      If ( LROOTS.gt.mxRoot ) IERR=1
      If (IERR.eq.1) Then
        Write(LF,*)
        Write(LF,*) '***************** ERROR *****************'
        Call WarningMessage(2,'Max roots exceeded.')
        Write(LF,'(1X,A,I6)')'Input Error: Max roots exceeded.',
     &                        mxRoot
        Write(LF,'(1X,A,I4)')'Nr of CI roots        LROOTS=',LROOTS
        Write(LF,'(1X,A,I4)')'Nr of optimized roots NROOTS=',NROOTS
        Write(LF,*)'************************************************'
        Call Quit_OnUserError()
      End If
      Do i=1,NROOTS
         If ( IROOT(i).lt.0 .or. IROOT(i).gt.LROOTS ) IERR=1
      End Do
      If (IERR.eq.1) Then
        Write(LF,*)
        Write(LF,*) '***************** ERROR *****************'
        Call WarningMessage(2,'Root specifications are wrong.')
        Write(LF,*)'CHKINP Error: Root specifications are wrong.'
        Write(LF,'(1X,A,I4)')'Nr of CI roots        LROOTS=',LROOTS
        Write(LF,'(1X,A,I4)')'Nr of optimized roots NROOTS=',NROOTS
        Write(LF,*)'************************************************'
        call Quit_OnUserError()
      End If

      If (STSYM.GT.NSYM) then
        Write(LF,*)
        Write(LF,*) '***************** ERROR *****************'
        Call WarningMessage(2,'Wrong symmetry.')
        Write(LF,'(1X,A,I8)')'CHKINP Error: Wrong symmetry.'
        Write(LF,'(1X,A,I8)')'State symmetry   STSYM=',STSYM
        Write(LF,'(1X,A,I8)')'Point group order NSYM=',NSYM
        Write(LF,*)'************************************************'
        call Quit_OnUserError()
      End If

      IERR=0
      If ( NHOLE1.lt.0 .and. NRS1T.ne.0 ) IERR=1
      If ( NELEC3.lt.0 .and. NRS3T.ne.0 ) IERR=1
      If ( NACTEL.lt.0 .and. NRS2T.ne.0 ) IERR=1
      If (IERR.eq.1) Then
        Write(LF,*)
        Write(LF,*)          '***************** ERROR *****************'
        Call WarningMessage(2,'Error in RAS specification.')
        Write(LF,*)          'Error in RAS specification.'
        Write(LF,'(1X,A,I8)')'Max holes in Ras1,     NHOLE1=',NHOLE1
        Write(LF,'(1X,A,I8)')'Max electrons in Ras3, NELEC3=',NELEC3
        Write(LF,'(1X,A,I8)')'Nr of active electrons NACTEL=',NACTEL
        Write(LF,*)          '*****************************************'
        call Quit_OnUserError()
      End If


! Second question: How should they be printed?
! No user selection, so fall back on default choice.
      IF(NTOT.LT.256) THEN
         OutFmt2='FULL    '
      ELSE
         OutFmt2='COMPACT '
      END IF

      End
