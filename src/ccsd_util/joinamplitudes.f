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
        subroutine joinamplitudes (wrk,wrksize)
c
c        this routine join contributions to amplitudes (Tn) from all
c        processors. Sum will be placed back to Tn on all machines
c
c        N.B. FREE ZONE in purpose of this routine is the space of
c       free working files V1-V4, .....
c        (free zone is not used in GA)
c        Since T2 amplitudes are no more then oovv, at most V1-V3 space
c        will be demaged. (acutally V1 and V2 only)
c
        use Para_Info, only: nProcs
        implicit none
c
#include "wrk.fh"
#include "ccsd1.fh"
#include "ccsd2.fh"
c
c        help variables
c
        integer ii,length
CLD        integer ii,length,rc,i
c
        if (nProcs.eq.1) return
c
c1        join t13,t14
c
c1.1    calc overall length of t13 and t14 (they must be one after the other)
        ii=mapdt14(0,5)
        length=mapdt14(ii,1)+mapdt14(ii,2)-posst130
c
c1.2        vanish required part in free zone
c        call mv0zero (length,length,wrk(possv10))
c
c1.3        allreduce t13 and t14 into free zone
c        call MPI_ALLREDUCE (wrk(posst130),wrk(possv10),length,
c    c  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
c
c1.4        put joined t13,t14 back to t13,t14 place from free zone
c        do i=0,length-1
c        wrk(posst130+i)=wrk(possv10+i)
c        end do
c
c1.om   allreduge t13,t14 together
        call gadgop (wrk(posst130),length,'+')
c
c
c2        join t21,t22,t23
c
c2.1    calc overall length of t21-t23 (they must be one after the other)
        ii=mapdt23(0,5)
        length=mapdt23(ii,1)+mapdt23(ii,2)-posst210
c
c2.2        vanish required part in free zone
c        call mv0zero (length,length,wrk(possv10))
c
c2.3        allreduce t13 and t14 into free zone
c        call MPI_ALLREDUCE (wrk(posst210),wrk(possv10),length,
c    c  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
c
c2.4        put joined t21-t23 back to t21-t23 place from free zone
c        do i=0,length-1
c        wrk(posst210+i)=wrk(possv10+i)
c        end do
c
c2.om   allreduge t21-t23 together
        call gadgop (wrk(posst210),length,'+')
c
        return
        end
