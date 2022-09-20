
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine cht3(ireturn)
! main driver for (T)

implicit none
#include "cht3_casy.fh"
#include "cht3_ccsd1.fh"
#include "files.fh"
#include "cht3_reord.fh"
#include "WrkSpc.fh"
integer ireturn
! DIRCC
!mp #include 'memvir_inc'
#include "ccsd_t3compat.fh"
integer wrksize
integer maxspace
integer isize
integer iOE, ioeh, ioep
!integer itmp,iW2,il1_1,ioff
integer i, nOrbE, nBas(8), nOrb(8)
character*24 Label
logical Found

!mp
! vynuluj hodiny
call CWTime(TCpu,TWall)
TCpu0 = TCpu
TWall0 = TWall
TCpu_l = TCpu
TCpu_l = TCpu
TWall_l = TWall
TWall_l = TWall
!mp
!.0 read input

call IniReord_t3(NvGrp,wrksize)

!.0.1 generate name convention for blocked integrals and T2 files

call DefParReord_t3(NvGrp,maxdim)
if (printkey >= 10) then
  write(6,*) 'Maxdim of virtual segment from CCSD = ',maxdim
end if

!.0.2 def commons for DIRCC

call defcommon(nfr,no,nv)

!.2.2 regenerate integrals from the Cholesky vectors (L1) L1(m,I ,A'')

!.2.2.1 (vo|vo) for testing purpose

!isize = nc*no*nv
!!write(6,*) 'size for l1_1 ',isize
!!write(6,*) 'size for l1_2 ',isize
!call GetMem('cht3_l1_1','Allo','Real',il1_1,isize)
!call GetMem('cht3_itmp','Allo','Real',itmp,isize)
!isize = nv*nv*no*no
!!write(6,*) 'size for W2 ',isize
!call GetMem('cht3_W2','Allo','Real',iW2,isize)
!call gen_vvoo(Work(iW2),Work(il1_1),Work(itmp))

!.2.2.2 get orbital energies

! Get Orbital energies

call Get_iArray('nBas',nBas,1)
call Get_iArray('nOrb',nOrb,1)

isize = nBas(1)

if (printkey >= 10) then
  write(6,*) 'Allocating memory for (tmp) OE files',isize
end if

call GetMem('cht3_oe','Allo','Real',iOE,isize)

Label = 'OrbE'
call qpg_dArray(Label,Found,nOrbE)
if (nOrbE /= nBas(1)) then
  write(6,*) 'Warning! in cht3 : (nOrbE /= nBas)!'
end if
if ((.not. Found) .or. (nOrbE == 0)) then
  call SysAbendMsg('get_orbe','Did not find:',Label)
end if
if (printkey >= 10) then
  write(6,*) 'isize = ',isize
  write(6,*) 'norbe = ',norbe
end if
call Get_dArray(Label,Work(iOE),nOrbE)

! write out the orbital energies

if (printkey >= 10) then
  write(6,*)
  write(6,*) 'Orbital energies for nfr+no+nv'
  write(6,*)
  do i=1,nfr+no+nv
    write(6,'(A,2x,i5,2x,f18.10)') 'Orbital Energy ',i,Work(iOE+i-1)
  end do
end if

!2.2.3 make OEH, OEP

isize = 2*no
call GetMem('cht3_oeh','Allo','Real',ioeh,isize)
isize = 2*nv
call GetMem('cht3_oeh','Allo','Real',ioep,isize)

call generate_juzekOE(Work(ioe+nfr),Work(ioeh),Work(ioep),no,nv)

!.2.3 Checkpoint. Calculate MP2 energy

!call calc_MP2(Work(iW2),Work(iOE+nfr),no,nv)
!call abend()

!.3 start (T) calculation

call GetMem('(T)','Max','Real',maxspace,maxspace)

write(6,*)
write(6,'(A,i13,A,f9.1,A,f5.1,A)') ' Memory available for (T) calc = ',maxspace-1,' in r*8 Words',(maxspace-1)*8.0d0/1024**2, &
                                   ' Mb',(maxspace-1)*8.0d0/1024**3,' Gb'

!mp call GetMem('t3_ampl_bti','Allo','Real',ioff,1)
!mp ioff = ioff+1
!mp !write(6,*) 'ioe   = ',ioe
!mp !write(6,*) 'ioeh  = ',ioeh
!mp !write(6,*) 'ioep  = ',ioep
!mp write(6,*) 'ioff volny  = ',ioff
!mp kvir1 = ioff+1
!mp !kvir2 = kvir1+(maxspace-1)-1

! toto sa da vyhodit a nahradit iba natvrdo definovanim kvir
!mp call adapt_mem(Work(1),Work(ioff),(maxspace-1),kvir1,kvir2)

!!call alloc_vm(WORK,maxspace,KVIR1,KVIR2)

!mp write(6,*) ' kvir1 = ',kvir1
!mp write(6,*) ' kvir2 = ',kvir2

!mp call T3AMPL_BTI(Work(ioff),Work(ioeh),Work(ioep))

call T3AMPL_BTI(Work(ioeh),Work(ioep))

!.2.4 Free the unnecessary memory

!mp isize = nfr+no+nv
!mp isize = nOrb(1)
isize = nBas(1)
call GetMem('cht3_oeh','Free','Real',iOE,isize)

isize = 2*no
call GetMem('cht3_oeh','Free','Real',ioeh,isize)
isize = 2*nv
call GetMem('cht3_oeh','Free','Real',ioep,isize)
!mp ioff = ioff-1
!mp call GetMem('t3_ampl_bti','Free','Real',ioff,1)

!!call GetMem('t3_ampl_bti','Free','Real',ioff,maxspace)

!Call EndGlb()

ireturn = 0

return

end subroutine cht3
