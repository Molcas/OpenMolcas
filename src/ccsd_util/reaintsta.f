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
       subroutine reaintsta (wrk,wrksize)
c
c     this routine read integral file INTSTA (reorg), which contains
c     following integrals: foka,fokb,
c     <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab
c     <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
c     <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
c
c     two electron integrals are readed to their fix files,
c     foka,fokb are readed to N,P help files
c
c     use and destroy : V1-3, N,P
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "filemgr.fh"
#include "wrk.fh"
c
c     help variables
c
       integer lunsta,rc,f_recl,f_iostat
       logical is_error
c
c*    open INTSTA file
       lunsta=1
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_open_ext2(lunsta,'INTSTA','sequential','unformatted',
     &                     f_iostat,.false.,f_recl,'unknown',is_error)
c       open (unit=lunsta,file='INTSTA',form='unformatted')
c
       else
c      MOLCAS IO
       call daname (lunsta,'INTSTA')
       daddr(lunsta)=0
       end if
c
c1    read foka to N
       call getmediate (wrk,wrksize,
     & lunsta,possn0,mapdn,mapin,rc)
c
c2    read fokb to P
       call getmediate (wrk,wrksize,
     & lunsta,possp0,mapdp,mapip,rc)
c
c
c3    read <kl||ij>aaaa to W01
       call getmediate (wrk,wrksize,
     & lunsta,possw010,mapdw01,mapiw01,rc)
c
c4    read <kl||ij>bbbb to W02
       call getmediate (wrk,wrksize,
     & lunsta,possw020,mapdw02,mapiw02,rc)
c
c5    read <kl||ij>abab to W03
       call getmediate (wrk,wrksize,
     & lunsta,possw030,mapdw03,mapiw03,rc)
c
c
c6    read <ie||mn>aaaa to W11
       call getmediate (wrk,wrksize,
     & lunsta,possw110,mapdw11,mapiw11,rc)
c
c7    read <ie||mn>bbbb to W12
       call getmediate (wrk,wrksize,
     & lunsta,possw120,mapdw12,mapiw12,rc)
c
c8    read <ie||mn>abab to W13
       call getmediate (wrk,wrksize,
     & lunsta,possw130,mapdw13,mapiw13,rc)
c
c9    read <ie||mn>baab to W14
       call getmediate (wrk,wrksize,
     & lunsta,possw140,mapdw14,mapiw14,rc)
c
c
c10   read <ab||ij>aaaa to V1
       call getmediate (wrk,wrksize,
     & lunsta,possv10,mapdv1,mapiv1,rc)
c
c11   read <ab||ij>bbbb to V2
       call getmediate (wrk,wrksize,
     & lunsta,possv20,mapdv2,mapiv2,rc)
c
c12   read <ab||ij>abab to V3
       call getmediate (wrk,wrksize,
     & lunsta,possv30,mapdv3,mapiv3,rc)
c
c*    close INTSTA file
c
       if (iokey.eq.1) then
c      Fortran IO
       close (lunsta)
c
       else
c      MOLCAS IO
       call daclos (lunsta)
       end if
c
       return
       end
