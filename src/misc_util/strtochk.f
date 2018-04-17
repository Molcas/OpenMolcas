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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
*  StrToChk
*
*> @brief
*>   Compute a checksum for a character string. Used for example to compute checksum for basis set
*> @author Per-Olof Widmark
*>
*> @details
*> This routine takes a character string and compute an
*> integer number. This number is not necessarily unique
*> but the chance that two different string will have the
*> same checksum is small. Whitespace is ignored so that
*> '``HelloDolly``' and '``Hello Dolly``'
*> will have the same checksum.
*> By default (\p iOpt = ``0``) the checksum is case insensitive,
*> but can be made case sensitive by setting \p iOpt = ``1``.
*>
*> @param[in]  String String to be checksummed
*> @param[out] Chk    Checksum
*> @param[in]  iOpt   Bitswitch, ``1`` for case sensitive checksum
************************************************************************
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine StrToChk(String,Chk,iOpt)
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) String
      Integer       Chk
      Integer       iOpt
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Logical sensitive
      Integer i,k,m,n
*----------------------------------------------------------------------*
* Process option switch                                                *
*----------------------------------------------------------------------*
      If(iAnd(iOpt,1).ne.0) Then
         sensitive=.true.
      Else
         sensitive=.false.
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      n=0
      k=1
      Do i=1,Len(String)
         k=mod(k+12,17)+1
         m=iChar(String(i:i))
         If(m.eq.32) GoTo 111
         If(m.eq.9)  GoTo 111
         If(.not.sensitive) Then
            If(m.ge.iChar('a') .and. m.le.iChar('z')) Then
               m=m-iChar('a')+iChar('A')
            End If
         End If
         n=n+k*m
  111    Continue
      End Do
      Chk=n
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
