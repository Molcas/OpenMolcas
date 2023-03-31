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
       subroutine contt12 (wrk,wrksize)
!
!     this routine do T12 contribution:
!     t1n(a,i) <- sum(e) [FI(a,e) . T1o(e,i)
!
!     N.B. use and destroy : M1
!     N.B. Parallel : in the case where idaaaa.ne.idbaab and
!          idaaaa.ne.idbaab this routine runs contributions to
!          T1n also on idaaaa and idbbbb nodes, since on these
!          nodes there is a specific
!          part of contributions F13 (see notes in sumoverb routine)
!           which are not presented on pilot nodes
!
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
!
!     help variables
!
       integer rc,ssc
!
!1    t1n(a,i)aa <- sum(e-a) [F1(a,e)aa . T1o(e,i)aa ]
!
!par
      if ((myRank.eq.idbaab).or.(myRank.eq.idaaaa)) then
!1.1  mult M1(a,i) <= F1(a,e)aa . T1o(e,i)aa
       call mult (wrk,wrksize,                                          &
     & 2,2,2,1,mapdf11,mapif11,1,mapdt11,mapit11,1,mapdm1,              &
     &            mapim1,ssc,possm10,rc)
!
!1.2  add t1n(a,i)aa <- M1(a,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt13,mapit13,1,rc)
      end if
!
!
!2    t1n(a,i)bb <- sum(e-a) [F1(a,e)bb . T1o(e,i)bb ]
!
!par
      if ((myRank.eq.idaabb).or.(myRank.eq.idbbbb)) then
!2.1  mult M1(a,i) <= F1(a,e)bb . T1o(e,i)bb
       call mult (wrk,wrksize,                                          &
     & 2,2,2,1,mapdf12,mapif12,1,mapdt12,mapit12,1,mapdm1,              &
     &            mapim1,ssc,possm10,rc)
!
!2.2  add t1n(a,i)bb <- M1(a,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdt14,mapit14,1,rc)
      end if
!
       return
       end
