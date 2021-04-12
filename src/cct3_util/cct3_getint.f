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
       subroutine cct3_getint (wrk,wrksize,
     & i,symi,possr0,mapdr,mapir,rc)
c
c     this routine read integrals R_i(a,bc) for given i in given symi
c
c     i      - number of orbital (I)
c     symi   - irrep of i (I)
c     possr0 - initial possition of R (I)
c     mapdr  - direct map of R (I)
c     mapir  - inverse map of R (I)
c     rc     - return (error) code (O)
c
       implicit none
#include "t31.fh"
#include "wrk.fh"
c
       integer i,symi,possr0,rc
c
       integer mapdr(0:512,1:6)
       integer mapir(1:8,1:8,1:8)
c
c     help variables
c
       integer iadd,lun,isym,num
c      integer rc1
       integer poss,length,im
c
c1    some tests
c
       if (i.gt.noa(symi)) then
c     RC=1 : i is higher than occipied in this irrep (Stup)
       rc=1
       return
       end if
c
       if (i.lt.1) then
c     RC=2 : i is less than 1 (Stup)
       rc=2
       return
       end if
c
c2    calc number for this orbital
c
       iadd=0
       if (symi.gt.1) then
       do 10 isym=1,symi-1
       iadd=iadd+noa(isym)
 10     continue
       end if
c
       num=iadd+i
c
c
c3    get R
c
       lun=1
       daddr(lun)=T3IntPoss(num)
c
       call daname (lun,t3nam)
c
       call idafile (lun,2,mapdr,513*6,daddr(lun))
       call idafile (lun,2,mapir,8*8*8,daddr(lun))
c
       poss=possr0
       length=0
       do im=1,mapdr(0,5)
         mapdr(im,1)=poss
         poss=poss+mapdr(im,2)
         length=length+mapdr(im,2)
c        write (*,99) ' MAP',(mapdr(im,k),k=1,6)
c99       format (a3,i8,2x,i8,4(2x,i2))
       end do
c
       if (length.gt.0) then
         call ddafile (lun,2,wrk(possr0),length,daddr(lun))
       end if
c      call cct3_getmediate (wrk,wrksize,
c    & lun,possr0,mapdr,mapir,rc1)

       call daclos (lun)
c
       return
       end
