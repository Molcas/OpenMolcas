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
* Copyright (C) 1991, Markus P. Fuelscher                              *
************************************************************************
      Subroutine PrInp(CMO)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Echo all input                                                   *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)

#include "motra_global.fh"
      Real*8 CMO(*)
      Character*120  Line,BlLine,StLine
      Character*8    Fmt
      Logical   PrOcc,PrEne
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
      left=(lPaper-lLine)/2
      Write(Fmt,'(A,I3.3,A)') '(',left,'X,A)'
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
           Write(6,Fmt) '*'//Line//'*'
         End Do
         Write(6,*)
      End If
*----------------------------------------------------------------------*
*     Print the integral file header                                   *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,'(6X,A)') 'Header of the integral files:'
      Write(Line,'(A)') Header( 1: 72)
      Write(6,'(6X,A)') Line(:mylen(Line))
      Write(Line,'(A)') Header(73:144)
      Write(6,'(6X,A)') Line(:mylen(Line))
      Write(6,*)
*----------------------------------------------------------------------*
*     Print the header of the source file of MO coefficients           *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,'(6X,A)') 'Header of MO coefficients source file:'
      Write(6,'(6X,A)') VecTit
      Write(6,*)
*----------------------------------------------------------------------*
*     Print coordinates of the system                                  *
*----------------------------------------------------------------------*
      Call PrCoor
*----------------------------------------------------------------------*
*     Print the orbital specifications                                 *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,'(6X,A)')'Orbital specifications:'
      Write(6,'(6X,A)')'-----------------------'
      Write(6,*)
      Write(6,'(6X,A,T35,8i4)')
     * 'Symmetry species:',(i,i=1,nSym)
      Write(6,'(6X,A,T35,8i4)')
     * 'Number of basis functions:',(nBas(i),i=1,nSym)
      Write(6,'(6X,A,T35,8i4)')
     * 'Frozen orbitals:',(nFro(i),i=1,nSym)
      Write(6,'(6X,A,T35,8i4)')
     * 'Deleted orbitals:',(nDel(i),i=1,nSym)
      Write(6,'(6X,A,T35,8i4)')
     * 'Number of orbitals used:',(nOrb(i),i=1,nSym)
      If ( iAutoCut.eq.1 ) Then
        Write(6,'(6X,A)')
     *   'Automatic orbital deletion is turned on'
        If ( iAutoCut.eq.1 ) Then
          Write(6,'(6X,A,T35,8F10.8)')
     *     'Cutting thresholds:',(CutThrs(i),i=1,nSym)
        End If
      End If
      If ( iRFpert.ne.0 ) then
         Write(6,*)
         Write(6,*)
         Write(6,'(6X,A)')'Reaction field specifications:'
         Write(6,'(6X,A)')'------------------------------'
         Write(6,*)
         Write(6,'(6X,A)')'The Reaction field is added as a '//
     &                    'perturbation and has been determined '//
     &                    'in a previos calculation'
         Write(6,*)
      End If

*----------------------------------------------------------------------*
*     Print MO coefficients                                            *
*----------------------------------------------------------------------*
      If ( iPrint.ge.2 .or. Debug.eq.1 ) Then
        PrEne=.false.
        PrOcc=.false.
        If ( iAutoCut.eq.1 ) PrOcc=.true.
        ThrOcc=0.0D0
        ThrEne=0.0D0
        Ene=0.0D0
        Line='Input orbitals after orthogonalization'
        Call PRIMO(Line,PrOcc,PrEne,ThrOcc,ThrEne,
     *             nSym,nBas,nBas,BsLbl,[Ene],Occ,Cmo,-1)
      End If
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
*
      Return
      End
