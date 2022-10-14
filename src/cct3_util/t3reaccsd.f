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
       subroutine t3reaccsd (wrk,wrksize,                               &
     & eccsd)
!
!     this routine read CCSD results, it T1 and T2 amplitudes
!     and CCSD energy from saverst file
!
!     eccsd - Converged CCSD energy (O)
!
#include "t31.fh"
#include "t32.fh"
#include "wrk.fh"
!
       real*8 eccsd,dum(1)
!
!     help parameters
!
       integer lunrst,rc1
!
!
!1    open file savename
       lunrst=1
!
       if (iokey.eq.1) then
!      Fortran IO
!       open (unit=lunrst,file=filerst,form='unformatted')
        call molcas_binaryopen_vanilla(lunrst,filerst)
!
       else
!      MOLCAS IO
       call daname (lunrst,filerst)
       daddr(lunrst)=0
       end if
!
!2    get T1aa
       call cct3_getmediate (wrk,wrksize,                               &
     & lunrst,posst110,mapdt11,mapit11,rc1)
!
!3    get T1bb
       call cct3_getmediate (wrk,wrksize,                               &
     & lunrst,posst120,mapdt12,mapit12,rc1)
!
!4    get T2aaaa
       call cct3_getmediate (wrk,wrksize,                               &
     & lunrst,posst210,mapdt21,mapit21,rc1)
!
!5    get T2bbbb
       call cct3_getmediate (wrk,wrksize,                               &
     & lunrst,posst220,mapdt22,mapit22,rc1)
!
!6    get T2abab
       call cct3_getmediate (wrk,wrksize,                               &
     & lunrst,posst230,mapdt23,mapit23,rc1)
!
!7    get energy,niter
       if (iokey.eq.1) then
!      Fortran IO
       read (lunrst,end=1) eccsd,rc1
       else
!
!      MOLCAS IO
       call ddafile (lunrst,2,dum,1,daddr(lunrst))
       eccsd=dum(1)
       end if
!
       goto 999
!
 1      write(6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'
       write(6,*) ' USE CCSD ENERGY FROM CCSD OUTPUT FILE'
       eccsd=0.0d0
!
 999   if (iokey.eq.1) then
!      Fortran IO
       close (lunrst)
!
       else
!      MOLCAS IO
       call daclos (lunrst)
       end if
!
       return
       end
