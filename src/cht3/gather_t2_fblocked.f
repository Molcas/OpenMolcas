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
        subroutine gather_t2_fblocked(
     & length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & t2,t2_tmp,tmp)
c
c length1 = length of the 1st VO index (nv)
c length2 = length of the 2nd VO index (=< vblock)
c
c This routine generates T2 amplitudes in this form :
c
c T2 = t2(a,b,j<i) - t2(b,a,j<i) finaly : <a,b,(i<j) >  ;  a in nv, b in vblock
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer ngaf,ngal,ngbf,ngbl
        integer a,b,dima,dimb
        integer length
        integer lasta,lastb
        integer length1,length2
c
        real*8 t2(*),tmp(*),t2_tmp(*)
        integer a_tmp,b_tmp
c
        logical switch
        integer aa,bb
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=1,NvGrp
        do b=ngbf,ngbl
c
          switch=.false.
        if (a.ge.b) then
          aa=a
          bb=b
        else
          aa=b
          bb=a
          switch=.true.
        end if
c
        dima=DimGrpaR(aa)
        dimb=DimGrpaR(bb)
c
         if (aa.eq.bb) then  ! aa=bb
          length=(dima*(dima+1)*no*no)/2
          call GetX_t3 (tmp,length,LunAux,T2Name(aa,bb),1,1)
         else ! aa>bb
          length=dima*dimb*no*no
          call GetX_t3 (tmp,length,LunAux,T2Name(aa,bb),1,1)
         end if
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
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
cmp        write (6,'(A,2(i5,2x),L,2(i3,x))') 'lasta, lastb, switch, a, b = ',
cmp     & lasta,lastb,switch,a,b
c
        if (aa.eq.bb) then ! expand tmp
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2_fblocked1(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
        else
         if (.not.switch) then
        call grow_t2_fblocked1(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
         else
        call grow_t2_fblocked2(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
        end if
        end if
c
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(ngaf)
        call Unused_integer(ngal)
      end if
        end
