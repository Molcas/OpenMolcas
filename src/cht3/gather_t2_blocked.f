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
        subroutine gather_t2_blocked(
     & length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & t2,t2_tmp,tmp,switch)
c
c length1 = length of the 1st VO index (=< vblock)
c length2 = length of the 2nd VO index (=< vblock)
c
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer ngaf,ngal,ngbf,ngbl
        integer a,b,dima,dimb
cmp        integer length
        integer length
        integer lasta,lastb
        integer length1,length2
c
        real*8 t2(1:(length1*length2*no*no))
        real*8 tmp(1:(maxdim*maxdim*no*no))
        real*8 t2_tmp(1:(maxdim*maxdim*no*no))
        integer a_tmp,b_tmp
c
        logical sym
        logical switch
c
        sym=.false.
c
cmp        if (ngaf.eq.ngbf) sym=.true.
        if ((ngaf.eq.ngbf).and.(ngal.eq.ngbl)) sym=.true.
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=ngaf,ngal
        do b=ngbf,min0(a,ngbl)
c
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
c
         if (a.eq.b) then  ! a=b
          length=(dima*(dima+1)*no*no)/2
          call GetX_t3 (tmp,length,LunAux,T2Name(a,b),1,1)
         else ! a>b
          length=dima*dimb*no*no
          call GetX_t3 (tmp,length,LunAux,T2Name(a,b),1,1)
         end if
c
          lasta=0
        if (a.gt.ngaf) then
          do a_tmp=ngaf,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.ngbf) then
          do b_tmp=ngbf,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
cmp        write (6,'(A,2(i5,2x),2(i3,x))') 'lasta, lastb, a, b = ',
cmp     & lasta,lastb,a,b
c
        if (a.eq.b) then ! expand and map
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2_blocked(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b,sym,switch)
        else
        call grow_t2_blocked(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b,sym,switch)
        end if
c
        end do
        end do
c
        return
        end
