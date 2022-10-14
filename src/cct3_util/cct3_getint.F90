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
       subroutine cct3_getint (wrk,wrksize,                             &
     & i,symi,possr0,mapdr,mapir,rc)
!
!     this routine read integrals R_i(a,bc) for given i in given symi
!
!     i      - number of orbital (I)
!     symi   - irrep of i (I)
!     possr0 - initial possition of R (I)
!     mapdr  - direct map of R (I)
!     mapir  - inverse map of R (I)
!     rc     - return (error) code (O)
!
       implicit none
#include "t31.fh"
#include "wrk.fh"
!
       integer i,symi,possr0,rc
!
       integer mapdr(0:512,1:6)
       integer mapir(1:8,1:8,1:8)
!
!     help variables
!
       integer iadd,lun,isym,num
!      integer rc1
       integer poss,length,im
!
!1    some tests
!
       if (i.gt.noa(symi)) then
!     RC=1 : i is higher than occipied in this irrep (Stup)
       rc=1
       return
       end if
!
       if (i.lt.1) then
!     RC=2 : i is less than 1 (Stup)
       rc=2
       return
       end if
!
!2    calc number for this orbital
!
       iadd=0
       if (symi.gt.1) then
       do 10 isym=1,symi-1
       iadd=iadd+noa(isym)
 10     continue
       end if
!
       num=iadd+i
!
!
!3    get R
!
       lun=1
       daddr(lun)=T3IntPos(num)
!
       call daname (lun,t3nam)
!
       call idafile (lun,2,mapdr,513*6,daddr(lun))
       call idafile (lun,2,mapir,8*8*8,daddr(lun))
!
       poss=possr0
       length=0
       do im=1,mapdr(0,5)
         mapdr(im,1)=poss
         poss=poss+mapdr(im,2)
         length=length+mapdr(im,2)
!        write (*,99) ' MAP',(mapdr(im,k),k=1,6)
!99       format (a3,i8,2x,i8,4(2x,i2))
       end do
!
       if (length.gt.0) then
         call ddafile (lun,2,wrk(possr0),length,daddr(lun))
       end if
!      call cct3_getmediate (wrk,wrksize,
!    & lun,possr0,mapdr,mapir,rc1)

       call daclos (lun)
!
       return
       end
