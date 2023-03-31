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
       subroutine contf23 (wrk,wrksize)
!
!     this routine do:
!     f2(m,i) <- sum(e,n) [ <ie||mn> . T1o(e,n)]
!
       use Para_Info, only: MyRank
#include "ccsd2.fh"
#include "parallel.fh"
#include "wrk.fh"
!
!     help variables
!
       integer posst,rc,ssc
!
!1    f2(m,i)aa <- sum(e,n-aa) [ <ie||mn>aaaa . t1o(e,n)aa ]
!
!par
      if (myRank.eq.idbaab) then
!
!1.1  expand V1(i,e,m,n) <= <ie||mn>aaaa
       call expand (wrk,wrksize,                                        &
     & 4,3,mapdw11,mapiw11,1,possv10,mapdv1,mapiv1,rc)
!
!1.2  map V2(m,i,e,n) <= V1(i,e,m,n)
       call map (wrk,wrksize,                                           &
     & 4,2,3,1,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,           &
     &           rc)
!
!1.3  mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)aa
       call mult (wrk,wrksize,                                          &
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdm1,                &
     &            mapim1,ssc,possm10,rc)
!
!1.4  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
!
!
!
!2    f2(m,i)aa <- sum(e,n-bb) [ <ie||mn>abab . t1o(e,n)bb ]
!
!2.1  map V2(m,i,e,n) <= <ie||mn>abab
       call map (wrk,wrksize,                                           &
     & 4,2,3,1,4,mapdw13,mapiw13,1,mapdv2,mapiv2,possv20,               &
     &           posst,rc)
!
!2.2  mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)bb
       call mult (wrk,wrksize,                                          &
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdm1,                &
     &            mapim1,ssc,possm10,rc)
!
!2.3  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf21,mapif21,1,rc)
!
       end if
!
!
!3    f2(m,i)bb <- sum(e,n-bb) [ <ie||mn>bbbb . t1o(e,n)bb ]
!
!par
      if (myRank.eq.idaabb) then
!
!3.1  expand V1(i,e,m,n) <= <ie||mn>bbbb
       call expand (wrk,wrksize,                                        &
     & 4,3,mapdw12,mapiw12,1,possv10,mapdv1,mapiv1,rc)
!
!3.2  map V2(m,i,e,n) <= V1(i,e,m,n)
       call map (wrk,wrksize,                                           &
     & 4,2,3,1,4,mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,           &
     &           rc)
!
!3.3  mult M1(m,i) =  V2(m,i,e,n) . T2o(e,n)bb
       call mult (wrk,wrksize,                                          &
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt12,mapit12,1,mapdm1,                &
     &            mapim1,ssc,possm10,rc)
!
!3.4  add f2(m,i)bb <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
!
!
!
!4    f2(m,i)bb <- - sum(e,n-aa) [ <ie||nm>baab . t1o(e,n)aa ]
!
!4.1  map V2(m,i,e,n) <= <ie||nm>baab
       call map (wrk,wrksize,                                           &
     & 4,2,3,4,1,mapdw14,mapiw14,1,mapdv2,mapiv2,possv20,               &
     &           posst,rc)
!
!4.2  mult M1(m,i) = V2(m,i,e,n) . T2o(e,n)aa
       call mult (wrk,wrksize,                                          &
     & 4,2,2,2,mapdv2,mapiv2,1,mapdt11,mapit11,1,mapdm1,                &
     &            mapim1,ssc,possm10,rc)
!
!4.3  add f2(m,i)aa <- M1(m,i)
       call add (wrk,wrksize,                                           &
     & 2,2,0,0,0,0,1,1,-1.0d0,mapdm1,1,mapdf22,mapif22,1,rc)
!
       end if
!
       return
       end
