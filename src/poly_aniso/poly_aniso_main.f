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
      Subroutine poly_aniso(iReturn)

      Implicit None
#include "warnings.fh"
      Integer       :: nneq, exch, neqv, nmax, nLoc, nCenter,
     &                 nT, nH, nTempMagn, nDir, nDirZee,
     &                 nMult, nPair, MxRank1, MxRank2
      Integer       :: iReturn
      Logical       :: dbg

      iReturn=0
      dbg=.false.

      Write(6,'(A)') 'by:   Liviu Unugur       '//
     &               '(chmlu@nus.edu.sg)'
      Write(6,'(A)') 'and   Liviu F. Chibotaru '//
     &               '(Liviu.Chibotaru@kuleuven.be)'
      Write(6,'(A)') 'Last updated - 2 July 2018'

      ! initialize some important variables
      If(dbg) Write(6,*) 'Enter fetch_init_const'
      Call fetch_init_const( nneq, neqv, nmax, exch, nLoc,
     &                       nCenter, nT, nH, nTempMagn, nDir,
     &                       nDirZee, nMult, nPair, MxRank1, MxRank2,
     &                       iReturn )
      If(dbg) Write(6,*) 'Exit fetch_init_const'
      If(dbg) Write(6,*) '     nneq=',     nneq
      If(dbg) Write(6,*) '     neqv=',     neqv
      If(dbg) Write(6,*) '     nmax=',     nmax
      If(dbg) Write(6,*) '     exch=',     exch
      If(dbg) Write(6,*) '     nLoc=',     nLoc
      If(dbg) Write(6,*) '  nCenter=',  nCenter
      If(dbg) Write(6,*) '       nT=',       nT
      If(dbg) Write(6,*) '       nH=',       nH
      If(dbg) Write(6,*) 'nTempMagn=',nTempMagn
      If(dbg) Write(6,*) '     nDir=',     nDir
      If(dbg) Write(6,*) '  nDirZee=',  nDirZee
      If(dbg) Write(6,*) '    nMult=',    nMult
      If(dbg) Write(6,*) '    nPair=',    nPair
      If(dbg) Write(6,*) '  MxRank1=',  MxRank1
      If(dbg) Write(6,*) '  MxRank2=',  MxRank2
      If(dbg) Write(6,*) '  iReturn=',  iReturn
      If(dbg) Call xFlush(6)


      If (iReturn.ne.0) Then
         Write(6,*) 'ERROR: something went wrong during '//
     &              'initialization of main variables'
         Write(6,*) 'Have to quit now'
         Call Quit(_RC_GENERAL_ERROR_)
      End If

      ! main program
      If(dbg) Write(6,*) 'Enter poly_aniso_1'
      If(dbg) Call xFlush(6)
      Call POLY_ANISO_1( nneq, neqv, nmax, exch, nLoc,
     &                   nCenter, nT, nH, nTempMagn, nDir,
     &                   nDirZee, nMult, nPair, MxRank1, MxRank2,
     &                   iReturn )
      If(dbg) Write(6,*) 'Exit poly_aniso_1'

      If (iReturn.ne.0) Then
         Write(6,*) 'ERROR: something went wrong during '//
     &              'the execution of the POLY_ANISO_1'
         Call Quit(_RC_GENERAL_ERROR_)
      End If
      If(dbg) Call xFlush(6)

      Return
      End
