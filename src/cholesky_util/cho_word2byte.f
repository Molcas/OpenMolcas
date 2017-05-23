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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine Cho_Word2Byte(iWord,n,Byte,Unt)
************************************************************
*
*   <DOC>
*     <Name>Cho\_Word2Byte</Name>
*     <Syntax>Call Cho\_Word2Byte(iWord,n,Byte,Unt)</Syntax>
*     <Arguments>
*       \Argument{iWord}{Number of n-byte words}{Integer}{in}
*       \Argument{n}{Number of byte per word}{Integer}{in}
*       \Argument{Byte}{iWord in bytes/kb/Mb/Gb/Tb}{Real*8}{out}
*       \Argument{Unt}{Unit of Byte ['b ','kb','Mb','Gb', or 'Tb'}
*                     {Character*2}{out}
*     </Arguments>
*     <Purpose>
*        Convert integer number of n-byte words to byte with a
*        reasonable prefix (-,k,M,G, or T).
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*         Convert integer number of n-byte words to byte, kilobyte,
*         Megabyte, Gigabyte, or Terabyte.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Integer     iWord, n
      Real*8      Byte
      Character*2 Unt

      Byte = DBLE(iWord)*DBLE(n)
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
