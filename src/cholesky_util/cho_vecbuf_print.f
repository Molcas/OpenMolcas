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
      SubRoutine Cho_VecBuf_Print(Lupri,nSym)
C
C     Purpose: print allocation information of Cholesky vector buffer to
C              unit Lupri (if Lupri<1 nothing is printed).
C
      use ChoVecBuf
      Implicit None
      Integer Lupri, nSym

      Character*16 SecNam
      Parameter (SecNam = 'Cho_VecBuf_Print')

      Integer     iSym
      Real*8      xGb
      Character*2 Unt

      If (Lupri .lt. 1) Return
      If (nSym.lt.1 .or. nSym.gt.8) Then
         Call Cho_Quit('nSym error in '//SecNam,104)
      End If

      Call Cho_Head('Size of Cholesky vector buffer','-',80,Lupri)
      Write(Lupri,*)
      Do iSym = 1,nSym
         Call Cho_Word2Byte(l_ChVBuf_Sym(iSym),8,xGb,Unt)
         Write(Lupri,'(A,I2,A,I10,A,F8.2,A,A,A)')
     &   'Dimension, sym.',iSym,': ',l_ChVBuf_Sym(iSym),
     &   ' 8-byte words (',xGb,' ',Unt,')'
      End Do
      Call Cho_Word2Byte(l_ChVBuf,8,xGb,Unt)
      Write(Lupri,'(/,A,I10,A,F8.2,A,A,A)')
     & 'Total dimension  : ',l_ChVBuf,' 8-byte words (',xGb,' ',Unt,')'

      End
