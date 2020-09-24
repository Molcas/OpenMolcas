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
      SubRoutine Cho_Alaska_RdInp(LuSpool)
***********************************************************
*
*     Purpose: Read and process input for Cholesky section
*              in alaska.
*
***********************************************************

      Implicit Real*8 (A-H,O-Z)
#include "exterm.fh"
      Character*180 KWord, Key, Get_Ln
      External Get_Ln
      character*16 SECNAM
      parameter (SECNAM = 'CHO_ALASKA_INPUT')
      Real*8 dmpK
      Integer nScreen
      Logical timings
      COMMON  /CHOTIME /timings
*
*     Set defaults
*
      dmpK = 1.0d0
      dmpK_default = dmpK
      nScreen = 10
*
      iPrint=5

*     Process the input
 1000 continue
      Key=Get_Ln(LuSpool)
      Kword=Key
      Call UpCase(Kword)
      If(KWord(1:1).eq.'*')    Go To 1000
      If(KWord(1:4).eq.'')  Go To 1000
      If(KWord(1:4).eq.'DMPK') Go To 100
      If(KWord(1:4).eq.'SCRN') Go To 110
      If(KWord(1:4).eq.'TIMI') Go To 120
      If(KWord(1:4).eq.'ENDC') Go To 998
      If(KWord(1:4).eq.'END ') Go To 998
      If(KWord(1:4).eq.'ENDO') Go To 998
*                                                             *
**** DMPK *****************************************************
*                                                             *
 100  Continue
      Read(LuSpool,*, err=210, end=200) dmpK
      If(dmpK .lt. 0.0d0) Then
         Write(6,*) 'OBS! Specified DMPK value is negative.'
         Write(6,*) 'Restoring Default!'
         dmpK = dmpK_default
      End If
      Go To 1000
*                                                             *
**** SCRN *****************************************************
*                                                             *
 110  Continue
      Read(LuSpool,*, err=210, end=200) nScreen
      Go To 1000
*                                                             *
**** TIMI *****************************************************
*                                                             *
 120  Continue
      Timings=.True.
      Go To 1000
*                                                             *
*** ENDChoinput ***********************************************
*                                                             *
 998  Continue
      Return
*                                                             *
***************************************************************
*                                                             *
      Call ErrTra
 200  Write(6,*) SECNAM, 'Premature end of input file.'
      Call Quit_onUserError()
      Call ErrTra
 210  Write(6,*) SECNAM, 'Error while reading input file.'
      Call Quit_onUserError()
*
      End
