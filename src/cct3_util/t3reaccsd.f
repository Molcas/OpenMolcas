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
       subroutine t3reaccsd (wrk,wrksize,
     & eccsd)
c
c     this routine read CCSD results, it T1 and T2 amplitudes
c     and CCSD energy from saverst file
c
c     eccsd - Converged CCSD energy (O)
c
#include "t31.fh"
#include "t32.fh"
#include "wrk.fh"
c
       real*8 eccsd,dum(1)
c
c     help parameters
c
       integer lunrst,rc1
c
c
c1    open file savename
       lunrst=1
c
       if (iokey.eq.1) then
c      Fortran IO
c       open (unit=lunrst,file=filerst,form='unformatted')
        call molcas_binaryopen_vanilla(lunrst,filerst)
c
       else
c      MOLCAS IO
       call daname (lunrst,filerst)
       daddr(lunrst)=0
       end if
c
c2    get T1aa
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst110,mapdt11,mapit11,rc1)
c
c3    get T1bb
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst120,mapdt12,mapit12,rc1)
c
c4    get T2aaaa
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst210,mapdt21,mapit21,rc1)
c
c5    get T2bbbb
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst220,mapdt22,mapit22,rc1)
c
c6    get T2abab
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst230,mapdt23,mapit23,rc1)
c
c7    get energy,niter
       if (iokey.eq.1) then
c      Fortran IO
       read (lunrst,end=1) eccsd,rc1
       else
c
c      MOLCAS IO
       call ddafile (lunrst,2,dum,1,daddr(lunrst))
       eccsd=dum(1)
       end if
c
       goto 999
c
 1      write(6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'
       write(6,*) ' USE CCSD ENERGY FROM CCSD OUTPUT FILE'
       eccsd=0.0d0
c
 999   if (iokey.eq.1) then
c      Fortran IO
       close (lunrst)
c
       else
c      MOLCAS IO
       call daclos (lunrst)
       end if
c
       return
       end
