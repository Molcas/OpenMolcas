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
      subroutine getAOs(lhigh)
      Implicit Real*8 (a-h,o-z)
cbs   get expansions of atomic orbitals in contracted functions
#include "para.fh"
#include "param.fh"
      character*12    occtext,occread
      character*18  textnorbmf,textnorbmf2
      logical EX
      occtext='OCCUPATION: '
      textnorbmf='Number of orbitals'
      Inquire(File='AO-expansion',exist=EX)
      if (.not.EX)  then
CBS      write(6,*) 'get occupations from DATA-block'
         call getAOs2(lhigh)
         return
      endif
      Lu_33=33
      Lu_33=IsFreeUnit(Lu_33)
      call molcas_open(Lu_33,'AO-expansion')
c      open(unit=Lu_33,file='AO-expansion',STATUS='UNKNOWN')
CBS   write(6,*) 'Orbitals for mean-field'
      do lrun=0,lhigh
CBS   write(6,'(A3,I3)') 'L= ',lrun
      read(Lu_33,'(A18,I3)') textnorbmf2,noccorb(lrun)
      if (textnorbmf.ne.textnorbmf2) call SysAbendMsg('getAOs',
     *'wrong keyword for number of orbitals in getAOs',' ')
CBS   write(6,*) 'number of orbitals ',noccorb(lrun)
      do iorbital=1,noccorb(lrun)
      read(Lu_33,'(A12,F5.3)')  occread,occup(iorbital,lrun)
CBS   write(6,'(A,F8.4)') occtext,occup(iorbital,lrun)
      if (occread.ne.occtext) call SysAbendMsg('getAOs',
     & 'error reading AOs',' ')
      read(Lu_33,*) (AOcoeffs(icont,iorbital,lrun),
     *icont=1,ncontrac(lrun))
CBS   write(6,'(8F10.4)') (AOcoeffs(icont,iorbital,lrun),
CBS  *icont=1,ncontrac(lrun))
CBS   write(6,*) ' '
      read(Lu_33,*)
      enddo
      enddo
      close(Lu_33)
      return
      end
