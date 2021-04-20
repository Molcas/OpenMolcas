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
       subroutine ireorg2 (symp,typp,pup,rc)
c
c*    this routine def. sumation limits for given symp and typp
c     i.e. number of indexes for this symmetry and typ
c
c     symp  - irrep of p index (I)
c     typp  - typ of p index in v2 (I)
c     pup   - sumation limit
c     rc    - return (error) code (O)
c
#include "ccsort.fh"
       integer symp,typp,pup,rc
c
       if (typp.eq.1) then
       pup=noa(symp)
       else if (typp.eq.2) then
       pup=nob(symp)
       else if (typp.eq.3) then
       pup=nva(symp)
       else if (typp.eq.4) then
       pup=nvb(symp)
       else if (typp.eq.5) then
       pup=norb(symp)
       else
       rc=1
c     RC=1 : bad typp (Stup)
       return
       end if
c
       return
       end
