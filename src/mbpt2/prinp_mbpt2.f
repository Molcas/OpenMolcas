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
* Copyright (C) 1992, Markus P. Fuelscher                              *
************************************************************************
      Subroutine PrInp_MBPT2(Eocc,Eext)
************************************************************************
*                                                                      *
*     Print the program banner, date and time of execution             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
*
#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "cdtfaux.fh"
#include "print_mbpt2.fh"
*
*     declare local variables...
      Dimension Eocc(*),Eext(*)
      Character*8 Fmt1,Fmt2
      Character*120 Line,BlLine,StLine
      Character*102 XLine
      Character*3 lIrrep(8)
      Logical lFro,lDel
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      lPaper=132
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      lPaper=132
      left=(lPaper-lLine)/2
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*
*     Print the project title                                          *
*----------------------------------------------------------------------*
      If ( nTit.gt.0 ) then
         Write(6,*)
         nLine=nTit+5
         Do i=1,nLine
           Line=BlLine
           If ( i.eq.1 .or. i.eq.nLine )
     &     Line=StLine
           If ( i.eq.3 )
     &     Line='Project:'
           If ( i.ge.4 .and. i.le.nLine-2 ) then
              Line=Title(i-3)
           End If
           Call Center(Line)
           Write(6,Fmt1) '*'//Line//'*'
         End Do
         Write(6,*)
      End If
*----------------------------------------------------------------------*
*     Print the coordinates of the system                              *
*----------------------------------------------------------------------*
      If (iPL.ge.2) Call PrCoor
*----------------------------------------------------------------------*
*     Print contents of the runfile RUNFILE                            *
*----------------------------------------------------------------------*
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         Call RightAd(lIrrep(iSym))
      End Do
*
      If (iPL.ge.2) Then
         Write(6,*)
         Write(6,Fmt2//'A)') 'Contents of RUNFILE file:'
         Write(6,Fmt2//'A)') '-------------------------'
         Write(6,*)
         Write(6,Fmt2//'A,T47,8I4)') 'Symmetry species',
     &                               (iSym,iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8(1X,A))') '                ',
     &                               (lIrrep(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     &                               (nBas(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Frozen occupied orbitals',
     &                               (nFro(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Active occupied orbitals',
     &                               (nOcc(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Active external orbitals',
     &                               (nExt(iSym),iSym=1,nSym)
         Write(6,Fmt2//'A,T47,8I4)') 'Deleted external orbitals',
     &                            (nDel(iSym),iSym=1,nSym)
      End If
*----------------------------------------------------------------------*
*     ordering and orbital energies of the frozen occupied orbitals    *
*----------------------------------------------------------------------*
      lFro=.false.
      Do iSym=1,nSym
         If ( nFro1(iSym)+nFro2(iSym).ne.0 ) lFro=.true.
      End Do
      If (iPL.ge.2) Then
         If ( lFro ) then
            Write(6,*)
            Write(6,*)
            Write(6,Fmt2//'A,T47)')
     &         'Reference numbers of frozen occupied orbitals '//
     &         'according to the original input sequence'
            Do iSym=1,nSym
               If ( nFro1(iSym).ne.0 ) then
                  Write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))')
     &               'symmetry species',
     &               iSym,(iOrb,iOrb=1,nFro1(iSym))
               End If
               If ( nFro2(iSym).ne.0 ) then
                  Write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))')
     &               'symmetry species',
     &               iSym,(iFro(iSym,iOrb),iOrb=1,nFro2(iSym))
               End If
            End Do
         End If
         Write(6,*)
         Write(6,*)
         Write(6,Fmt2//'A,T47)')
     &      'Energies of the active occupied orbitals'
         ii=0
         Do iSym=1,nSym
            If ( nOcc(iSym).ne.0 ) then
               Write(6,*)
               Write(6,Fmt2//'A,I2,(T40,5F14.6))')
     &            'symmetry species',iSym,
     &            (Eocc(ii+k),k=1,nOcc(iSym))
               ii=ii+nOcc(iSym)
            End If
         End Do
      End If
*----------------------------------------------------------------------*
*     ordering and orbital energies of the deleted external orbitals   *
*----------------------------------------------------------------------*
      lDel=.false.
      Do iSym=1,nSym
         If ( nDel(iSym)+nDel1(iSym)+nDel2(iSym).ne.0 ) lDel=.true.
      End Do
      If (iPL.ge.2) Then
         If ( lDel ) then
            Write(6,*)
            Write(6,*)
            Write(6,Fmt2//'A,T47)')
     &         'Reference numbers of deleted external orbitals '//
     &         'according to the original input sequence'
            Do iSym=1,nSym
               If ( nDel1(iSym).ne.0 .or. nDsto(iSym).ne.0 ) then
                  Write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))')
     &               'symmetry species',iSym,
     &               (nBas(iSym)-nFro(iSym)-nOcc(iSym)-iOrb+1,
     &               iOrb=nDsto(iSym)+nDel1(iSym),1,-1)
               End If
               If ( nDel2(iSym).ne.0 ) then
                  Write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))')
     &               'symmetry species',
     &               iSym,(iDel(iSym,iOrb),iOrb=1,nDel2(iSym))
               End If
            End Do
         End If
         Write(6,*)
         Write(6,*)
         Write(6,Fmt2//'A,T47)')
     &      'Energies of the active external orbitals'
         ii=0
         Do iSym=1,nSym
            If ( nExt(iSym).ne.0 ) then
               Write(6,*)
               Write(6,Fmt2//'A,I2,(T40,5F14.6))')
     &            'symmetry species',iSym,
     &            (Eext(ii+k),k=1,nExt(iSym))
               ii=ii+nExt(iSym)
            End If
         End Do
      End If
*----------------------------------------------------------------------*
*     # of orbitals in 1st index to be skipped (restart)               *
*----------------------------------------------------------------------*
      If ( iRest.ne.0 .and. iPL.ge.2) then
       Write(6,*)
       Write(6,Fmt2//'A)')        'restart information...'
       Write(6,Fmt2//'A,T47,I8)') 'max # integral passes',nPass
       Write(6,Fmt2//'A,T47,8I4)')'symmetry species',
     &                            (iSym,iSym=1,nSym)
       Write(6,Fmt2//'A,T47,8I4)')'# of orbitals in 1st index skipped',
     &                            (MOSkip(iSym-1),iSym=1,nSym)
      End If
*----------------------------------------------------------------------*
*     print header for final results                                   *
*----------------------------------------------------------------------*
      If ( iTst.eq.0 .and. iPL.ge.2) then
         Write(6,*)
         Write(6,*)
         nLine=3
         Do i=1,nLine
           XLine=Trim(BlLine)
           If ( i.eq.1 .or. i.eq.nLine )
     &     XLine=Trim(StLine)
           If ( i.eq.2 )
     &     XLine='Results'
           Call Center(XLine)
           Write(6,Fmt1) '*'//XLine//'*'
         End Do
         Write(6,*)
      End If
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Return
      End
