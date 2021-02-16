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
* Copyright (C) 1995, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine WriteDip(DipGrad,Modes,Title,nOsc)
C!
      Real*8 DipGrad(3,nOsc)
      Integer  Modes(nOsc)
      Character Title*(*)
#include "inout.fh"
C!
C!
      Write(6,*)
      Write(6,*)
      Write(6,'(a2,a)') ' ',Title
      Write(6,'(a2,a)') ' ',
     &       '=============================================='
      Write(6,'(a2,a)') ' ',
     &       ' mode          X           Y           Z      '
      Write(6,'(a2,a)') ' ',
     &       '----------------------------------------------'
      Do i = 1,nOsc
      Write(6,'(a3,i2,a1,a3,3f12.5)') ' ',Modes(i),'.',
     &               ' ',(DipGrad(j,i),j=1,3)
      End Do
      Write(6,'(a2,a)') ' ',
     &       '=============================================='
      Write(6,*)
C!
      End
C!
C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine IntCalcHeader
C!
#include "inout.fh"

      Write(6,*)
      Write(6,*)
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       '|             Intensity calculation               |'
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,*)
C!
      End
C!
C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine ExpPointHeader
C!
#include "inout.fh"
C!
      Write(6,*)
      Write(6,*)
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       '|            Expansion point geometry             |'
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,*)
C!
      End
C!
C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine ISCHeader
C!
#include "inout.fh"

      Write(6,*)
      Write(6,*)
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       '|      InterSystem Crossing rate calculation      |'
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,*)
C!
      End
C!
C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine WriteHeader(Title)
C!
C!  Purpose:
C!    Write header and title to logfile.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Character*80 Title
#include "inout.fh"
C!
      Write(6,*)
C!
C!---- Write title of project.
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,'(A,A)') '  Title : ',Title
      Write(6,'(A)') '  -------'
      Write(6,*)
C!
      End
C!
      Subroutine WrMold(FName,NumOfAt,AtomLbl,AtCoord,NumInt,
     &       HarmFreq,QMat)
C!
C! Open file named FName, and create a MOLDEN input file.
C!
c       Implicit None
#include "Constants_mula.fh"
      Character*(*) FName
      Integer NumOfAt,NumInt
      Character*(*) AtomLbl(NumOfAt)
      Real*8 AtCoord(3,NumOfAt),HarmFreq(NumInt),
     &       QMat(3,NumOfAt,NumInt)
      Real*8 DMax,D2,DispMx,Factor
      Integer i,iInt,iAtom
      Character*2 AtName
#include "inout.fh"
        call molcas_open(9,Fname)
c      Open (9,file=FName,status='UNKNOWN')
      Write (9,*) '[MOLDEN FORMAT]'

C! Harmonic frequencies:
      Write (9,*) '[N_FREQ]'
      Write (9,*) NumInt
      Write (9,*) '[FREQ]'
      Do iInt = 1, NumInt
      Write (9,'(1X,F10.3)') HarToRcm*HarmFreq(iInt)
      End Do
c VV: Not implemented???
      Write (9,*) '[INT]'
      Do iInt = 1, NumInt
      Write (9,'(1X,F10.3)') 0.0
      End Do
C! Atom coordinates:
      Write (9,*) '[NATOM]'
      Write (9,*) NumOfAt
      Write (9,*) '[FR-COORD]'
      Do iAtom = 1, NumOfAt
      AtName=AtomLbl(iAtom)(1:2)
      If(IChar('0').le.IChar(AtName(2:2)).and.
     &    IChar(AtName(2:2)).le.IChar('9')) Then
      AtName(2:2)=' '
      End If
      Write (9,'(1x,A2,3F16.8)') AtName,(AtCoord(i,iAtom),i=1,3)
      End Do

C! Cartesian displacement coordinates:
C! Scale such that max displacement is 0.2 A.U.
      DMax=0.2D0
      Write (9,*) '[FR-NORM-COORD]'
      Do iInt = 1, NumInt
      Write (9,*) 'vibration ', iInt
      DispMx=0.0D0
      Do iAtom = 1, NumOfAt
      D2=0.0D0
      Do i=1,3
      D2=D2+QMat(i,iAtom,iInt)**2
      End Do
      DispMx=Max(DispMx,sqrt(D2))
      End Do
      Factor=DMax/DispMx
      Do iAtom = 1, NumOfAt
      Write (9,'(1X,3F16.8)') (Factor*QMat(i,iAtom,iInt),i=1,3)
      End Do
      End Do
C!
      Close(9)

      End
