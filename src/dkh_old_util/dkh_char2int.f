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
* Copyright (C) 2004,2006, Alexander Wolf                              *
*               2006, Markus Reiher                                    *
************************************************************************
      integer function dkh_char2int (leng,string)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1)
c                  and dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 07.07.2006  (MR, ETH Zurich, name changed)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer leng,idum,j
      character*(*) string
c
      if (leng.lt.1) then
        write (stdout,1001)
1001    format (/2X,'ERROR in function char2int: length of string ',
     *          'is smaller than 1.',//2X,'STOP.',/)
        call Abend
      endif
c
      idum=0
      do 10 j=1,leng
        if (string(j:j).eq.'0') then
          idum=idum+0*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'1') then
          idum=idum+1*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'2') then
          idum=idum+2*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'3') then
          idum=idum+3*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'4') then
          idum=idum+4*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'5') then
          idum=idum+5*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'6') then
          idum=idum+6*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'7') then
          idum=idum+7*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'8') then
          idum=idum+8*(10**(leng-j))
          goto 222
        endif
        if (string(j:j).eq.'9') then
          idum=idum+9*(10**(leng-j))
          goto 222
        endif
        write (stdout,1004)
1004    format (/2X,'ERROR in function char2int: string ',
     *          'contains illegal character (no figure).',
     *           //2X,'STOP.',/)
        call Abend
 222    continue
 10   continue
c
      dkh_char2int = idum
c
      return
      end
