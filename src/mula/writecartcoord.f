************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine WriteCartCoord(AtomLbl,Coord,Mass,NumOfAt)
C!
C!  Purpose:
C!    Write cartesian coordinates to log file.
C!
      Character*4  AtomLbl(NumOfAt)
      Real*8 Coord( 3, NumOfAt )
      Real*8 Mass (NumOfAt)
#include "inout.fh"
C!
C!---- Format declarations.
C!
C!---- Initialize.
      nAt = NumOfAt
C!
C!---- Write labels, coordinates and masses to log file.
      Write(6,*)
      Write(6,*)
      Write(6,'(a1,a)') ' ',
     &       'Cartesian coordinates (in bohr) and masses (in u)'
      Write(6,*) ('====',i=1,17)
      Write(6,'(a2,a)')
     &         ' ',
     & 'Atom         x             y             z                Mass'
      Write(6,*) ('----',i=1,17)
      Do j = 1,NumOfAt
      Write(6,'(a2,a4,3f14.8,f20.8)')
     &        ' ',AtomLbl(j),(Coord(i,j),i=1,3),Mass(j)
      End Do
      Write(6,*) ('====',i=1,17)
      Write(6,*)
      Write(6,*)
C!
      End
C!
C!-----------------------------------------------------------------------!

C!-----------------------------------------------------------------------!
C!
      Subroutine WriteIntCoord(InterVec,AtomLbl,xvec,NumInt)
C!
C!  Purpose:
C!    Write internal coordinates to log file.
C!
#include "Constants_mula.fh"
c      Integer InterVec (NumInt)
      Integer InterVec (*)
      Character*4 AtomLbl(NumInt)
      Character*128 Line
      Real*8 xvec(NumInt)
      Real*8  const
#include "inout.fh"
C!
C!---- Internal coordinates at equilibrium.
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,'(a1,a)') ' ','Internal coordinates at equilibrium'
      Write(6,*) ('====',i=1,16)
      Write(6,'(a2,a)') ' ',
     & 'Distances :                            bohr          aangstrom'
      Write(6,'(a2,a)') ' ',
     & 'Angles    :                           radians         degrees'
      Write(6,*) ('----',i=1,16)
      k = 1
      iLine=32
      Do j = 1,NumInt
      IntType = InterVec(k)
      If(IntType.eq.1) Then
C!---- Bond Stretching.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      Write(Line,fmt='(A2,A,A,A,A,A)')
     &              ' ','Bond    ',
     &             AtomLbl(i1),'- ',AtomLbl(i2),'            '
      k = k+3
      Endif
      If(IntType.eq.2) Then
C!---- Valence Angle Bending.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      Write(Line,fmt='(A2,A,A,A,A,A,A,A)')
     &              ' ','Angle   ',
     &        AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),'      '
      k = k+4
      EndIf
      If(IntType.eq.3) Then
C!---- Linear Valence Angle.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      Write(Line,fmt='(A2,A,A,A,A,A,A,A)')
     &              ' ','LinAng  ',
     &         AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),'      '
      k = k+4
      EndIf
      If(IntType.eq.4) Then
C!---- Torsion.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      i4 = InterVec(k+4)
      Write(Line,fmt='(A2,A,A,A,A,A,A,A,A)')
     &              ' ','Torsion ',
     &               AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),
     &                         '- ',AtomLbl(i4)
      k = k+5
      EndIf
      If(IntType.eq.5) Then
C!---- Out of Plane Angle Bending.
      i1 = InterVec(k+1)
      i2 = InterVec(k+2)
      i3 = InterVec(k+3)
      i4 = InterVec(k+4)
      Write(Line,fmt='(A2,A,A,A,A,A,A,A,A)')
     &               ' ','OutOfPl ',
     &               AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),
     &                         '- ',AtomLbl(i4)
      k = k+5
      EndIf
      If ( intType.eq.1 ) Then
      const = Angstrom
      Else
      const = 180.0d0/rpi
      End If
      Write(6,'(A,A1,F15.8,F16.8)')
     &               Line(1:32),' ',xvec(j),xvec(j)*const
      End Do
      Write(6,*) ('====',i=1,16)
      Write(6,*)
      Write(6,*)
C!
      End
