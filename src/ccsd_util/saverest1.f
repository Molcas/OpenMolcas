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
       subroutine saverest1 (wrk,wrksize,
     & lunrst)
c
c     this routine save restart informations:
c     t13,t14,t21,t22,t23
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
c
       integer lunrst
c
c     help variables
c
       integer rc
c
c0    return if need
       if (keyrst.eq.0) return
c
c1    rewind tape
       call filemanager (2,lunrst,rc)
c
c2    write T1aa
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt13,mapit13,rc)
c
c3    write T1bb
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt14,mapit14,rc)
c
c4    write T2aaaa
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt21,mapit21,rc)
c
c5    write T2bbbb
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt22,mapit22,rc)
c
c6    write T2abab
       call wrtmediate (wrk,wrksize,
     & lunrst,mapdt23,mapit23,rc)
c
       return
       end
