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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
      subroutine removeB1 (length,coeff,operator)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 12.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************

      implicit none
#include "dkhparameters.fh"

      integer length
      REAL*8 coeff
      character*(maxlength) operator

      integer k,l,m,dumleng,idum

      do 10 k=length,2,-1
        if (operator(k:k).eq.'B') then
          if (operator(k-1:k-1).eq.'B') then
            operator(k:k)=' '
            operator(k-1:k-1)=' '
          else if (operator(k-1:k-1).eq.'P') then
            coeff=-coeff
            operator(k-1:k-1)='B'
            operator(k:k)='P'
          else
            operator(k:k)=operator(k-1:k-1)
            operator(k-1:k-1)='B'
          endif
        endif
  10  continue

      dumleng=length
      idum=0
      do 20 l=1,dumleng
        idum=idum+1
        if (operator(idum:idum).eq.' ' .and. idum.le.length) then
          do 30 m=idum+1,length
            operator(m-1:m-1)=operator(m:m)
  30      continue
          if (operator(idum:idum).eq.' ') idum=idum-1
          operator(length:length)=' '
          length=length-1
        endif
  20  continue

      return
      end
