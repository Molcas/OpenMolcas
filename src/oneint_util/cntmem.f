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
      Subroutine CntMem(nHer,MemCnt,la,lb,lr)
      use Sizes_of_Seward, only: s
      Implicit Real*8 (A-H,O-Z)
*
      nHer=S%mCentr
      MemCnt = 3*(la+1)*nHer +
     &         3*(lb+1)*nHer +
     &         ((la+1)*(la+2)/2) *
     &         ((lb+1)*(lb+2)/2)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
