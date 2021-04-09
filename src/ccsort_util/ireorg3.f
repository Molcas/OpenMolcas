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
       subroutine ireorg3 (symp,typp,typpv1,paddv1,rc)
c
c*    this routine def. constants to be added to index from v2
c     to determine proper index in v1
c     N.B. typp and typpv1 must be compatible (this is testet
c     in this version)
c
c     symp  - irrep of p index (I)
c     typp  - typ of p index in v2 (I)
c     typpv1- typ of corresponding p index in v1 (I)
c     paddv1- constant to be added (O) pv1 = pv2+paddv1
c     rc    - return (error) code (O)
c
#include "ccsort.fh"
       integer symp,typp,typpv1,paddv1,rc
c
       rc=0
c
       if ((typp.eq.1).or.(typp.eq.2)) then
          if ((typpv1.eq.1).or.(typpv1.eq.2).or.(typpv1.eq.5)) then
             paddv1=0
          else
             rc=1
c            RC=1 : typp=1 or 2, incompatible typpv1 (Stup)
             return
          end if
       else if (typp.eq.3) then
       if (typpv1.eq.3) then
          paddv1=0
       else if (typpv1.eq.4) then
           paddv1=nvb(symp)-nva(symp)
       else if (typpv1.eq.5) then
           paddv1=noa(symp)
       else
           rc=2
c          RC=2 : typp=3, incompatible typpv1 (Stup)
           return
       end if
       else if (typp.eq.4) then
       if (typpv1.eq.4) then
       paddv1=0
       else if (typpv1.eq.5) then
       paddv1=nob(symp)
       else
       rc=3
c     RC=3 : typp=4, incompatible typpv1 (Stup)
       return
       end if
       else if (typp.eq.5) then
       if (typpv1.eq.5) then
       paddv1=0
       else
c     RC=4 : typp=5, incompatible typpv1 (Stup)
       return
       end if
       else
       rc=5
c     RC=5 : improper typp (Stup)
       return
       end if
c
       return
       end
