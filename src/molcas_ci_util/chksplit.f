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
* Copyright (C) 2010, Giovanni Li Manni                                *
************************************************************************

      Subroutine ChkSplit()
************************************************************************
*     SplitCAS Check for obvious errors or violation  of limits        *
*----------------------------------------------------------------------*
*     written by:                                                      *
*     Giovanni Li Manni (GLMJ)                                         *
*     University of Geneva, Switzerland, 2010                          *
*----------------------------------------------------------------------*
*     history: none                                                    *
************************************************************************

      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"

#include "output_ras.fh"
#include "general.fh"
#include "splitcas.fh"
#include "warnings.fh"
      Parameter ( ROUTINE='ChkSplit' )

      Call qEnter('ChkSplit')
*      IF(DoSplitCAS) THEN
        IERR=0
        If (lRootSplit.gt.nConf) IERR=1
        If (IERR.eq.1) then
          Write(LF,*)
          Write(LF,*)'******************** ERROR *********************'
          Write(LF,'(1X,A)')'Input Error: '
          write(LF,*)' Root you are looking for is not reachable within'
          write(LF,*)' the selected active space.'
          Write(LF,*)' Try to select a bigger active space!'
          write(LF,*)' Root selected by user = ', lRootSplit
          write(LF,*)' Root reachable = ', nConf
          Write(LF,*)'************************************************'
          Call Quit(_RC_INPUT_ERROR_)
        End If

*       If ( NumSplit ) then
        IERR=0
        If ( iDimBlockA.lt.lRootSplit) IERR=1
        If (IERR.eq.1) then
          Write(LF,*)
          Write(LF,*)'******************** ERROR **********************'
          Write(LF,*)'Input Error: AA-Block selected is too small!'
          Write(LF,'(1X,A,I5)')' Root to be optimized :',lRootSplit
          Write(LF,'(1X,A,I5)')' AA-Block dimension   :',iDimBlockA
        write(LF,*)'AA-Block must be always equal or greater than root.'
          write(LF,*) 'In a NUSP calculation increase iDimBlockA'
          write(LF,*) 'In a ENSP calculation increase the energy-gap'
          write(LF,*) 'In a PESP calculation increase the percentage'
          Write(LF,*)'*************************************************'
          Call Quit(_RC_INPUT_ERROR_)
*          iDimBlockA=lRootSplit
*          ircWar = 1
*          Write(LF,'(1X,A,I8)')
*     &    'iDimBlockA has been reset to =',iDimBlockA
        End If

        IERR=0
        If (iDimBlockA.gt.mxDimBlockA) IERR=1
        If (IERR.eq.1) then
          Write(LF,*)
          Write(LF,*) '***************** ERROR *****************'
          Write(LF,'(1X,A,I6)')'Input Error: Max dim. BlockA exceeded',
     &                        mxDimBlockA
          write(LF,*) 'iDimBlockA selected by user = ', iDimBlockA
          write(LF,*) 'If you are running a NUSP calculation, ',
     &    'please, decrease the value of iDimBlockA!'
          write(LF,*) 'If you are running a ENSP calculation, ',
     &    'please, decrease the energy-gap!'
          write(LF,*) 'If you are running a PESP calculation, ',
     &    'please, decrease the percentage!'
          Write(LF,*)'************************************************'
          Call Quit(_RC_INPUT_ERROR_)
        End If
*       End If
*       ^ End chk over NumSplit

      IERR=0
      If (ThrSplit.lt.min_ThrSplit) IERR = 1
      If (IERR.eq.1) then
        Write(LF,*)
        Write(LF,*) '***************** ERROR *****************'
        Write(LF,'(1X,A,I6)')'Input Error: ThrSplit too small'
        Write(LF,'(1X,A,I6)')'minimum value ThrSplit = ', min_ThrSplit
        write(LF,*) 'ThrSplit selected by user = ', ThrSplit
        Write(LF,*)'************************************************'
        Call Quit(_RC_INPUT_ERROR_)
      End If

*      END IF

      Call qExit('ChkSplit')
      Return
      End
