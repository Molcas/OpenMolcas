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
* Copyright (C) 1994, Markus P. Fuelscher                              *
*               1994, Per Ake Malmqvist                                *
*               Pavel Neogrady                                         *
************************************************************************
       Subroutine PrInpPN
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     - echo the input parameters                                      *
*                                                                      *
*     calling parameters: none                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P.-AA. Malmqvist                              *
*     University of Lund, Sweden, 1994                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     Modified by P.N.                                                 *
*                                                                      *
************************************************************************
       Implicit Real*8 (A-H,O-Z)

#include "ccsort.fh"
#include "reorg.fh"

CLD    Character*8   Fmt1,Fmt2
CLD    Character*120  Line,BlLine,StLine
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
*       lPaper=132
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
*       lLine=Len(Line)
*       Do i=1,lLine
*       BlLine(i:i)=' '
*       StLine(i:i)='*'
*       End Do
*       lPaper=132
*       left=(lPaper-lLine)/2
*       Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
*       Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*
*     Print the project title                                          *
*----------------------------------------------------------------------*
*       If ( nTit.gt.0 ) then
*       Write(*,*)
*       nLine=nTit+5
*       Do i=1,nLine
*       Line=BlLine
*       If ( i.eq.1 .or. i.eq.nLine )
*     & Line=StLine
*       If ( i.eq.3 )
*     & Line='Project:'
*       If ( i.ge.4 .and. i.le.nLine-2 )
*     & Write(Line,'(18A4)')(Title(i-3,j),j=1,18)
*       Call Center(Line)
*       Write(*,Fmt1) '*'//Line//'*'
*       End Do
*       Write(*,*)
*      End If
* ---------------------------------------------------------------------*
*     Stop if NOOPeration key is used                                  *
*----------------------------------------------------------------------*
       If (noop.eq.1) then
         write(6,'(6X,A)') ' No operation is required'
         write(6,'(6X,A)') ' Happy Landing '
         Call Finish(0)
       end if
*----------------------------------------------------------------------*
*     Print iokey                                                      *
*----------------------------------------------------------------------*
       if (iokey.eq.1) then
         write(6,'(6X,A)') 'Standard Fortran IO handling used '
       end if
c
       if (iokey.eq.2) then
         write(6,'(6X,A)')
     &     'MOLCAS DA IO handling used '
       end if
*----------------------------------------------------------------------*
*     Print zrkey                                                      *
*----------------------------------------------------------------------*
       if (fullprint.eq.2) then
       if (zrkey.eq.1) then
         write(6,'(6X,A)') 'Separate V and Ind IO'
       end if
c
       if (zrkey.eq.0) then
         write(6,'(6X,A)') 'Simultanneous V and Ind IO'
       end if
       end if
*----------------------------------------------------------------------*
*     Print cckey and t3key                                            *
*----------------------------------------------------------------------*
       if (cckey.eq.1) then
         write(6,'(6X,A)') 'Integrals for CCSD will be produced'
       end if
c
       if (t3key.eq.1) then
         write(6,'(6X,A)')
     &     'Integrals for Noniterative T3 will be produced'
       end if
c
       if (clopkey.eq.1) then
         write(6,'(6X,A)') 'ROHF open shell reference function'
       else
         write(6,'(6X,A)') 'RHF closed shell reference function'
       end if
*----------------------------------------------------------------------*
*     Print allocation and printing parameters                         *
*----------------------------------------------------------------------*
c      if (maxspace.eq.0) then
c      write(6,'(6X,A)') ' Allocatable work space   : Unlimited'
c      else
c      write(6,'(6X,A,I10)') ' Allocatable work space   : ',maxspace
c      end if
c      if (fullprint.eq.0) then
c      write(6,'(6X,A)') ' Level of output printing : Minimal'
c      else if (fullprint.eq.1) then
c      write(6,'(6X,A)') ' Level of output printing : Medium '
c      else if (fullprint.eq.2) then
c      write(6,'(6X,A)') ' Level of output printing : Full'
c      end if
*----------------------------------------------------------------------*
*     Print actual frozen and deleted orbitals                         *
*----------------------------------------------------------------------*
       Write(6,*)
       Write(6,'(6X,A)')'Actual numbers of frozen and '//
     & 'deleted orbitals :'
       Write(6,'(6X,A)')'-----------------------------'//
     & '------------------'
       Write(6,*)
       Write(6,'(6X,A,T47,8I4)') 'Symmetry species',
     & (iSym,iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Frozen orbitals',
     & (nFror(iSym),iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Deleted orbitals',
     & (nDelr(iSym),iSym=1,nSym)
       Write(6,*)
c
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
       Write(6,*)
       Write(6,'(6X,A)')'Wave function specifications '//
     & 'from previous RASSCF:'
       Write(6,'(6X,A)')'-----------------------------'//
     & '---------------------'
       Write(6,*)
       Write(6,'(6X,A,T45,I6)')'Number of closed shell electrons',
     & 2*NISHT
       Write(6,'(6X,A,T45,I6)')'Number of electrons in active shells',
     & NACTEL
       Write(6,'(6X,A,T45,I6)')'Max number of holes in RAS1 space',
     & NHOLE1
       Write(6,'(6X,A,T45,I6)')'Max number of electrons in RAS3 '//
     & 'space',NELE3
       Write(6,'(6X,A,T45,I6)')'Number of inactive orbitals',
     & NISHT
       Write(6,'(6X,A,T45,I6)')'Number of active orbitals',
     & NASHT
       Write(6,'(6X,A,T45,I6)')'Number of secondary orbitals',
     & NSSHT
       Write(6,'(6X,A,T45,F6.1)')'Spin quantum number',
     & (dble(ISPIN-1))/2.
       Write(6,'(6X,A,T45,I6)')'State symmetry',
     & LSYM
       Write(6,'(6X,A,T45,I6)')'Number of configuration state fnc.',
     & NCONF
       Write(6,'(6X,A,T45,I6)')'Number of root(s) available',
     & NROOTS
       Write(6,'(6X,A,T45,5I6)')'CI root used',
     & LROOT
       If ( ISCF.eq.0 ) then
       Write(6,'(6X,A)')
     & 'This is a CASSCF reference function'
       Else If ( ISCF.eq.1 ) then
       Write(6,'(6X,A)')
     & 'This is a closed shell RHF reference function'
       Else
       Write(6,'(6X,A)')
     & 'This is a high spin open shell RHF reference function'
       End If
       Write(6,*)
       Write(6,*)
       Write(6,'(6X,A)')'Orbital specifications from '//
     & 'previous RASSCF:'
       Write(6,'(6X,A)')'----------------------------'//
     & '----------------'
       Write(6,*)
       Write(6,'(6X,A,T47,8I4)') 'Symmetry species',
     & (iSym,iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Frozen orbitals',
     & (nFro(iSym),iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Inactive orbitals',
     & (nIsh(iSym),iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Active orbitals',
     & (nAsh(iSym),iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Secondary orbitals',
     & (nSsh(iSym),iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Deleted orbitals',
     & (nDel(iSym),iSym=1,nSym)
       Write(6,'(6X,A,T47,8I4)') 'Number of basis functions',
     & (nBas(iSym),iSym=1,nSym)
       Write(6,*)
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
       Return
       End
