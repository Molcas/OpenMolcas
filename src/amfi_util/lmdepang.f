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
* Copyright (C) 1996, Bernd Schimmelpfennig                            *
************************************************************************
      REAL*8 function LMdepang(
     *L,M,l1,l2,l3,l4,m1,m2,m3,m4,cheater)
cbs   l1-l4 and m1-m4 are already shifted !!
cbs   purpose: calculates the angular part of the
cbs   coulomb-type integrals. See documentation for details...
cbs   LMdepang= LM dependent angular factors
cbs   cheater included for a correcting signs, as there were some
cbs   signs (only signs!!!!) missing when compared to HERMIT
cbs                                        B.S.  08.10.96
      implicit REAL*8 (a-h,o-z)
#include "real.fh"
      LMdepang=0d0
cbs   some quick checks
      if (L.lt.abs(M)) return
      if (l1.lt.abs(m1)) return
      if (l2.lt.abs(m2)) return
      if (l3.lt.abs(m3)) return
      if (l4.lt.abs(m4)) return
cbs   prefactor
      fact1=4d0*pi/DBLE(L+L+1)
cbs   determining the sign
      isum=-l3-l1-l4-l2+2*(M+m3+m4)   !???? I am not sure
      if (mod(isum,4).eq.0) then
         isign=1
      elseif (iabs(mod(isum,4)).eq.2) then
         isign=-1
      else
         isign=0
         write(6,*) 'L,l1,l2,l3,l4,M,m1,m2,m3,m4'
         write(6,'(10I3)') L,l1,l2,l3,l4,M,m1,m2,m3,m4
         write(6,*) 'isum= ',isum,' mod = ',mod(isum,4)
         Call SysHalt( 'lmdepang' )
      endif
      fact2=couple3J(L,l3,l1,-M,m3,-m1)
      fact3=couple3J(L,l4,l2,M,m4,-m2)
C     write(6,*) 'fact2,fact3 ',fact2,fact3
      LMdepang=cheater*DBLE(isign)*fact1*fact2*fact3
      return
      end
