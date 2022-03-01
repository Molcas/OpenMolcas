!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine mkampq (wrk,wrksize,                                  &
     & a,ammap)
!
!     this routine reconstruct #2 V2<_a,m|p,q> from corresponding TEMPDA2 file
!
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer a
       integer ammap(1:mbas,1:8,1:8)
!
!     help variables
!
       integer symm,symp
       integer iiv2,length,poss,irec0
!
!*    loops over symmetry combinations
       do 100 symm=1,nsym
       do 101 symp=1,nsym
!
!*    def initioal record possition in TEMPDA2
!     and corresponding possition and length in wrk (#2)
!
       irec0=ammap(a,symm,symp)
       iiv2=mapi2(symm,symp,1)
       poss=mapd2(iiv2,1)
       length=mapd2(iiv2,2)
!
       if (length.gt.0) then
       call daread (lunda2,irec0,wrk(poss),length,recl)
       end if
!
 101    continue
 100    continue
!
       return
       end
