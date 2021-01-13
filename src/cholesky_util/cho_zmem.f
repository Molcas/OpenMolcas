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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_ZMem(irc,l_Z,NVT,l_NVT,DoPrint,DoCheck)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: compute dimension of Z vector array and print it to
C              output (if requested through DoPrint). Check if there
C              is sufficient memory at this point (if DoCheck).
C
C     Return codes:
C
C     irc=-1 : Input error (debug only)
C     irc=0  : All ok
C     irc=1  : Negative length of Z vector array (could be integer
C              overflow)
C     irc=999: Insufficient memory for Z vectors (only if DoCheck)
C
      Implicit None
      Integer irc
      Integer l_Z
      Integer l_NVT
      Integer NVT(l_NVT)
      Logical DoPrint, DoCheck
#include "cholesky.fh"

#if !defined (_I8_) || defined (_DEBUGPRINT_)
      Character*8 SecNam
      Parameter (SecNam='Cho_ZMem')
#endif
      Character*2 Unt
      Integer iSym, ip_Mx, l_Mx
      Real*8 Byte, xl_Z
      Real*8 Word(8)

#if defined (_DEBUGPRINT_)
      If (l_NVT .lt. nSym) Then
         irc=-1
         l_Z=-999999
         Return
      End If
#endif

      irc=0

      xl_Z=0.0d0
      Do iSym=1,nSym
         Word(iSym)=DBLE(NVT(iSym))*(DBLE(NVT(iSym))+1.0d0)/2.0d0
         xl_Z=xl_Z+Word(iSym)
      End Do
      l_Z=INT(xl_Z)

      If (DoPrint) Then
         Call Cho_Head('Z Vector Storage Requirements','-',80,LuPri)
         Write(LuPri,*)
         Do iSym=1,nSym
            Call Cho_RWord2Byte(Word(iSym),Byte,Unt)
            Write(LuPri,'(A,I2,A,I8,A,F8.3,1X,A,A)')
     &      'Symmetry',iSym,':   ',INT(Word(iSym)),' words (',
     &      Byte,Unt,')'
         End Do
         Write(LuPri,'(A)')
     &   '------------------------------------------'
         Call Cho_RWord2Byte(xl_Z,Byte,Unt)
         Write(LuPri,'(A,I8,A,F8.3,1X,A,A)')
     &   'Total:        ',l_Z,' words (',Byte,Unt,')'
      End If

#if !defined (_I8_) || defined (_DEBUGPRINT_)
      If (l_Z .lt. 0) Then
         Write(Lupri,'(A,A)')
     &   SecNam,': dimension of Z vector array is negative!'
         Write(Lupri,'(A,I8)') 'l_Z=',l_Z
         If (xl_Z .gt. 0.0d0) Then
            Write(LuPri,'(A)') 'This seems to be an integer overflow!'
            Call Cho_RWord2Byte(xl_Z,Byte,Unt)
            Write(LuPri,'(A,1P,D15.6,A,D15.6,1X,A,A)')
     &      'In double precision, xl_Z=',xl_Z,
     &      ' words (',Byte,Unt,')'
         End If
         irc=1
         Return
      End If
#endif

      If (DoCheck) Then
         Call mma_maxDBLE(l_Mx)
         If (l_Z.gt.l_Mx) Then
            irc=999
            Return
         End If
      End If

      End
      SubRoutine Cho_RWord2Byte(Word,Byte,Unt)
      Implicit None
      Real*8  Word
      Real*8  Byte
      Character*2 Unt

      Byte = Word*8.0d0
      Unt  = 'b '
      If (ABS(Byte) .gt. 1.0d3) Then
         Byte = Byte/1.024d3
         Unt  = 'kb'
         If (ABS(Byte) .gt. 1.0d3) Then
            Byte = Byte/1.024d3
            Unt  = 'Mb'
            If (ABS(Byte) .gt. 1.0d3) Then
               Byte = Byte/1.024d3
               Unt  = 'Gb'
               If (ABS(Byte) .gt. 1.0d3) Then
                  Byte = Byte/1.024d3
                  Unt  = 'Tb'
               End If
            End If
         End If
      End If

      End
