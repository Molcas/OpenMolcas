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
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      FUNCTION ILEX_FOR_CONF_NEW(ICONF,NOCC_ORB,NORB,NEL,IARCW,
     &       IDOREO,IREO_new,nconf_op,ib_occls)
*
* A configuration ICONF of NOCC_ORB orbitals are given
* ICONF(I) = IORB implies  IORB is singly occupied
* ICONF(I) = -IORB  implies that IORB is doubly occupied
*
* Find lexical address
*
* IF IDOREO .ne. 0, IREO is used to reorder lexical number
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
*. Arcweights for single and doubly occupied arcs
      INTEGER IARCW(NORB,NEL,2)
*. Reorder array
      INTEGER IREO_new(*)
c     integer nconf_per_open(1,*)
*. Configuration
      INTEGER ICONF(NOCC_ORB)
*
      IEL = 0
      ILEX = 1

      DO IOCC = 1, NOCC_ORB
       IF(ICONF(IOCC).GT.0) THEN
         IEL = IEL + 1
         ILEX = ILEX + IARCW(ICONF(IOCC),IEL,1)
       ELSE IF(ICONF(IOCC).LT.0) THEN
         IEL = IEL + 2
         ILEX = ILEX + IARCW(-ICONF(IOCC),IEL,2)
       END IF
      END DO
*
      IF(IDOREO.NE.0) THEN
c
c     length=nconf_per_open(2*nocc_orb-nel+1,1)
      length=nconf_op
      ntest =  0
      if(ntest.ge.10) then
      write(6,*)'in ilex_for_conf_new, check if all is right'
      write(6,*)'iconf:',(iconf(i),i=1,nocc_orb)
      write(6,*)'ilex and offset:',ilex,ib_occls
      write(6,*)'nconf_op =', nconf_op
      write(6,*)'ireo_new:'
      call iwrtma(ireo_new,1,length,1,length)
      endif
c      call iwrtma(nconf_per_open,1,4,1,4)
c.. look for the position where the value is ilex
        if (ireo_new(1).eq.ilex+ib_occls-1) then
          n_fin=1
        elseif (ireo_new(length).eq.ilex+ib_occls-1) then
          n_fin=length
        else
          n_ini = 1
          n_end = length
666          continue
            n_ave = (n_ini + n_end)/2
            if (ireo_new(n_ave).eq.ilex+ib_occls-1) goto 100
            if (ireo_new(n_ave).lt.ilex+ib_occls-1) n_ini = n_ave
            if (ireo_new(n_ave).gt.ilex+ib_occls-1) n_end = n_ave
          goto 666
100     continue
        n_fin = n_ave
        end if
        ilex_for_conf_new = n_fin

      ELSE
       ILEX_FOR_CONF_NEW = ILEX
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,*) ' new Lexical number = ', ILEX
        IF(IDOREO.NE.0)
     &  WRITE(6,*) ' new Reordered number = ', ILEX_FOR_CONF_NEW
      END IF
*
      RETURN
      END
