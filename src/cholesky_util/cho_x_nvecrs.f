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
      SubRoutine Cho_X_nVecRS(iRed,iSym,iVec,nVec)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_nVecRS</Name>
*     <Syntax>Call Cho\_X\_nVecRS(iRed,iSym,iVec,nVec)</Syntax>
*     <Arguments>
*       \Argument{iRed}{Reduced set}{Integer}{in}
*       \Argument{iSym}{Symmetry block (1-8)}{Integer}{in}
*       \Argument{iVec}{First vector in red. set iRed, sym. iSym}
*                {Integer}{out}
*       \Argument{nVec}{Number of vectors in red. set iRed, sym. iSym}
*                {Integer}{out}
*     </Arguments>
*     <Purpose>Find first vector and number of vectors in
*              reduced set iRed, sym. block iSym
*     </Purpose>
*     <Dependencies>Call Cho\_X\_Init(rc)</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>This routine finds the first vector and number of
*                  vectors in reduced set iRed, sym. block iSym.
*                  Note that iVec=nVec=0 may be returned - this is
*                  perfectly acceptable: a given reduced set may be
*                  empty. However, if negative numbers are returned
*                  (iVec<0 and nVec<0), an error has ocurred. This
*                  should be tested by the caller!!
*     </Description>
*    </DOC>
*
************************************************************

      Implicit None
      Integer iRed, iSym, iVec, nVec
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'Cho_X_nVecRS')

      Logical Found

      Integer irc, LastRed, jVec, jRed

      Integer InfVec, i, j, k
      Integer N2
      Parameter (N2 = InfVec_N2)

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)

C     Check input.
C     ------------

      irc = 0
      If (iSym.lt.1 .or. iSym.gt.nSym) Then
         irc = -1
      End If
      If (NumCho(iSym).lt.0 .or. NumCho(iSym).gt.MaxVec) Then
         irc = -2
      End If
      LastRed = InfVec(NumCho(iSym),2,iSym)
      If (LastRed .lt. 1) Then
         irc = -3
      End If
      If (iRed .lt. 1) Then
         irc = -4
      End If
      If (irc .ne. 0) Then
         iVec = irc
         nVec = irc
         Return
      End If
      If (iRed .gt. LastRed) Then
         iVec = 0
         nVec = 0
         Return
      End If

C     Find first vector in reduced set iRed.
C     --------------------------------------

      Found = .false.
      jVec  = 0
      Do While (jVec.lt.NumCho(iSym) .and. .not.Found)
         jVec = jVec + 1
         jRed = InfVec(jVec,2,iSym)
         If (jRed .eq. iRed) Then
            iVec  = jVec
            Found = .true.
         Else If (jRed .gt. iRed) Then
            jVec = NumCho(iSym)  ! break loop
         End If
      End Do

C     No first vector <=> 0 vectors in reduced set iRed.
C     --------------------------------------------------

      If (.not. Found) Then
         iVec = 0
         nVec = 0
         Return
      End If

C     Count number of vectors in reduced set iRed.
C     --------------------------------------------

      nVec = 1
      jVec = iVec
      Do While (jVec .lt. NumCho(iSym))
         jVec = jVec + 1
         jRed = InfVec(jVec,2,iSym)
         If (jRed .eq. iRed) Then
            nVec = nVec + 1
         Else
            jVec = NumCho(iSym) ! break loop
         End If
      End Do

#if defined (_DEBUG_)
C     Debug: print result.
C     --------------------

      Write(6,*) SecNam,': there are ',nVec,' vectors in reduced set ',
     &           iRed,' (sym. block ',iSym,')'
      Write(6,*) SecNam,': first vector is: ',iVec
#endif

      End
