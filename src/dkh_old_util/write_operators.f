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
* Copyright (C) 2004, Alexander Wolf                                   *
*               2004,2007, Markus Reiher                               *
************************************************************************
      subroutine write_operators (nbas,isize,dkhorder,xorder,dkhscfflg,
     *                            scrno2,scr6,scr7,scr8,paramtype,
     *                            clight)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2)
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.06.2006 MR (ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer nbas,isize,dkhorder,xorder,scrno2
      REAL*8 scr6(isize,0:dkhorder),scr7(isize,0:xorder),
     *                 scr8(isize,scrno2),clight
      character*(3) paramtype
c
      integer i,j,k,idum
      intrinsic INT
#ifdef _MOLCAS_
      Integer IsFreeUnit
      outunit1=21
      outunit1=IsFreeUnit(outunit1)
      Call Molcas_Open(outunit1,'dkhout.21')
#else
      open (outunit1, file='dkhout.21', status='UNKNOWN',
     *        form='FORMATTED')
#endif
      rewind (outunit1)
      write (outunit1,1001) nbas,
     *                      dkhorder,xorder,
     *                      paramtype,dkhscfflg
1001  format (60('-'),/2X,
     *        'Basis:',1X,'No. of basis functions:',2X,I3,/2X,6('-'),
     *        /12X,
     *        'Relativity:',3X,'DKH order (H) :',2X,I2,/2X,11('-'),3X,
     *        'DKH order (X) :',2X,I2,/16X,'Paramtype',5X,':',2X,A3,
     *        /16X,'dkhscfflg',5X,':',2X,L1)
      write (outunit1,1008) clight
1008  format (16X,'Speed of light: ',D18.11,/,60('-'))
c
      iDum =INT((iSize+3)/4)
      write (outunit1,1033) iDum
1033  format ('***',/3X,'sinv(i,j) (j=1,nbas,  i=1,j)',3X,I5)
      do j=1,iSize,4
        write (outunit1,1035) (scr8(k,scrno2), k=j,Min(j+3,iSize))
1035    format (4(D25.18,1X))
      enddo
c
      do i=0,dkhorder
        write (outunit1,1243) i,idum
1243    format ('***',/2X,I2,3X,I5)
        do j=1,iSize,4
          write (outunit1,1245) (scr6(k,i), k=j,Min(j+3,iSize))
1245      format (4(D25.18,1X))
        enddo
      enddo
c
      close (outunit1)
c
      if (.not.dkhscfflg) then
#ifdef _MOLCAS_
        outunit2=22
        outunit2=IsFreeUnit(outunit2)
        Call Molcas_Open(outunit2,'dkhout.22')
#else
        open (outunit2, file='dkhout.22', status='UNKNOWN',
     *          form='FORMATTED')
#endif
        rewind (outunit2)
c
        write (outunit2,1001) nbas,
     *                        dkhorder,xorder,
     *                        paramtype,dkhscfflg
        write (outunit2,1008) clight
c
        iDum = INT((iSize+3)/4)
        write (outunit2,2033) idum
2033    format ('***',/3X,'sinv(i,j) (j=1,nbas,  i=1,j)',3X,I5)
*
        do j=1,iSize,4
          write (outunit2,2035) (scr8(k,scrno2), k=j,Min(iSize,j+3))
2035      format (4(D25.18,1X))
        enddo
        do i=0,xorder
          write (outunit2,2043) i,idum
2043      format ('***',/2X,I2,3X,I5)
          do j=1,iSize,4
            write (outunit2,2045) (scr7(k,i), k=j,Min(iSize,j+3))
2045        format (4(D25.18,1X))
          enddo
        enddo
c
        close (outunit2)
c
      endif
c
      return
      end
