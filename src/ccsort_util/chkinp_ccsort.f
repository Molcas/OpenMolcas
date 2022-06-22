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
       Subroutine ChkInp_ccsort
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Check input for consistency                                      *
*                                                                      *
************************************************************************

       Implicit Real*8 (A-H,O-Z)

#include "ccsort.fh"
#include "motra.fh"
*
c
c      Just print warning...
       If ( IPT2.eq.0 ) then
         Write(6,*)
         Write(6,*) '       !!!!!WARNING!!!!!'
         Write(6,*)
         Write(6,*) '      *** input error ***'
         Write(6,*) '  The JOBIPH file does not include '//
     &                 'canonical orbitals'
         Write(6,*)
         Write(6,*) '       !!!!!WARNING!!!!!'
         Write(6,*)
c        Call Quit_OnUserError()
       End If
c
       If ( NCONF.ne.1 ) then
         Write(6,*)
         Write(6,*) '  *** input error ***'
         Write(6,*) '  The JOBIPH file does not include '//
     &                 'a RHF or ROHF wave function'
         Write(6,*)
         Call Quit_OnUserError()
       End If
*
       iErr = 0
       If ( nSym.ne.nSymX ) iErr = 1
       Do iSym = 1,nSym
         If ( nBas(iSym).ne.nBasX(iSym) ) iErr = 1
       End Do
       If ( iErr.ne.0 ) then
         Write(6,*)
         Write(6,*) '  *** input error ***'
         Write(6,*) '  The JOBIPH and the TRAONE files '//
     &                 'are inconsistent'
         Write(6,*)
         Call Quit_OnUserError()
       End If
*
       Return
       End
