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
       subroutine mkintsta (wrk,wrksize,
     & foka,fokb)
c
c     this routine produces integral file INTSTA, which contains
c     following integrals: foka,fokb,
c     <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab
c     <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
c     <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
c
c     N.B. 1. work file #1 is used for <ij|pq> integrals, #2,3,4
c     must be free. possb0 must be defined
c     N.B. 2. this routine can be used only after definition of <ij|pq>
c     N.B. 3. this routine use followuing help routines:
c     expandfok
c     wrtmediate (from SYMM)
c
c
#include "wrk.fh"
#include "reorg.fh"
#include "files_ccsd.fh"
       real*8 foka(*)
       real*8 fokb(*)
c
c     help variables
c
       integer rc
c
c*    open INTSTA file
       if (iokey.eq.1) then
c      Fortarn IO
       call molcas_binaryopen_vanilla(lunsta,'INTSTA')
c       open (unit=lunsta,file='INTSTA',form='unformatted')
c
       else
c      MOLCAS IO
       call daname (lunsta,'INTSTA')
       daddr(lunsta)=0
       end if
c
c*    expand foka into work #2 and write to INTSTA
       call expandfok (wrk,wrksize,
     & foka)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    expand fokb into work #2 and write to INTSTA
       call expandfok (wrk,wrksize,
     & fokb)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c
c*    get #2 <kl||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,1,1,1,1,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <kl||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,2,2,2,2,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <kl| ij>abab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,1,2,1,2,1,0)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c
c*    get #2 <ka||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 3,1,3,1,1,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ka||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 3,2,4,2,2,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ka| ij>abab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,1,4,1,2,1,0)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ka||ij>baab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,2,3,1,2,0,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c
c*    get #2 <ab||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,3,3,1,1,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ab||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,4,4,2,2,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ab| ij>abab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,3,4,1,2,1,0)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
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
