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

subroutine cht3(ireturn)
! main driver for (T)

use ChT3_global, only: DimGrpaR, L1Name, L2Name, maxdim, nfr, no, nv, NvGrp, printkey, T2Name, TCpu, TCpu_l, TCpu0, TWall, &
                       TWall_l, TWall0
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, maxspace, nBas(8), nOrb(8), nOrbE
character(len=24) :: Label
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: OE(:), oeh(:), oep(:)
real(kind=wp), parameter :: kb = 1024.0_wp

!mp
! vynuluj hodiny
call CWTime(TCpu,TWall)
TCpu0 = TCpu
TWall0 = TWall
TCpu_l = TCpu
TWall_l = TWall
!mp
!.0 read input

call IniReord_t3(NvGrp)

!.0.1 generate name convention for blocked integrals and T2 files

call DefParReord_t3(NvGrp,maxdim)
if (printkey >= 10) write(u6,*) 'Maxdim of virtual segment from CCSD = ',maxdim

!.0.2 def commons for DIRCC

call defcommon(no,nv)

!.2.2 regenerate integrals from the Cholesky vectors (L1) L1(m,I ,A'')

!.2.2.1 (vo|vo) for testing purpose

!isize = nc*no*nv
!!write(u6,*) 'size for l1_1 ',isize
!!write(u6,*) 'size for l1_2 ',isize
!call mma_allocate(l1_1,isize,label='cht3_l1_1')
!call mma_allocate(tmp,isize,label='cht3_itmp')
!isize = nv*nv*no*no
!!write(u6,*) 'size for W2 ',isize
!call mma_allocate(W2,isize,label='cht3_W2')
!call gen_vvoo(W2,l1_1,tmp)

!.2.2.2 get orbital energies

! Get Orbital energies

call Get_iArray('nBas',nBas,1)
call Get_iArray('nOrb',nOrb,1)

if (printkey >= 10) write(u6,*) 'Allocating memory for (tmp) OE files',nBas(1)

call mma_allocate(OE,nBas(1),label='cht3_oe')

Label = 'OrbE'
call qpg_dArray(Label,Found,nOrbE)
if (nOrbE /= nBas(1)) write(u6,*) 'Warning! in cht3 : (nOrbE /= nBas)!'
if ((.not. Found) .or. (nOrbE == 0)) call SysAbendMsg('get_orbe','Did not find:',Label)
if (printkey >= 10) then
  write(u6,*) 'nbas(1) = ',nBas(1)
  write(u6,*) 'norbe = ',norbe
end if
call Get_dArray(Label,OE,nOrbE)

! write out the orbital energies

if (printkey >= 10) then
  write(u6,*)
  write(u6,*) 'Orbital energies for nfr+no+nv'
  write(u6,*)
  do i=1,nfr+no+nv
    write(u6,'(A,2x,i5,2x,f18.10)') 'Orbital Energy ',i,OE(i)
  end do
end if

!2.2.3 make OEH, OEP

call mma_allocate(oeh,2*no,label='cht3_oeh')
call mma_allocate(oep,2*nv,label='cht3_oep')

call generate_juzekOE(OE(nfr+1),oeh,oep,no,nv)

!.2.3 Checkpoint. Calculate MP2 energy

!call calc_MP2(W2,OE(nfr+1),no,nv)
!call abend()

!.3 start (T) calculation

call mma_maxDBLE(maxspace)

write(u6,*)
write(u6,'(A,i13,A,f9.1,A,f5.1,A)') ' Memory available for (T) calc = ',maxspace-1,' in r*8 Words', &
                                    real((maxspace-1)*8,kind=wp)/kb**2,' Mb',real((maxspace-1)*8,kind=wp)/kb**3,' Gb'

!mp call GetMem('t3_ampl_bti','Allo','Real',ioff,1)
!mp ioff = ioff+1
!mp !write(u6,*) 'ioe   = ',ioe
!mp !write(u6,*) 'ioeh  = ',ioeh
!mp !write(u6,*) 'ioep  = ',ioep
!mp write(u6,*) 'ioff volny  = ',ioff
!mp kvir1 = ioff+1
!mp !kvir2 = kvir1+(maxspace-1)-1

! toto sa da vyhodit a nahradit iba natvrdo definovanim kvir
!mp call adapt_mem(Work(1),Work(ioff),(maxspace-1),kvir1,kvir2)

!!call alloc_vm(WORK,maxspace,KVIR1,KVIR2)

!mp write(u6,*) ' kvir1 = ',kvir1
!mp write(u6,*) ' kvir2 = ',kvir2

!mp call T3AMPL_BTI(Work(ioff),oeh,oep)

call T3AMPL_BTI(oeh,oep)

!.2.4 Free the unnecessary memory

!mp isize = nfr+no+nv
!mp isize = nOrb(1)
call mma_deallocate(OE)

call mma_deallocate(oeh)
call mma_deallocate(oep)
!mp ioff = ioff-1
!mp call GetMem('t3_ampl_bti','Free','Real',ioff,1)

call mma_deallocate(DimGrpaR)
call mma_deallocate(L1Name)
call mma_deallocate(L2Name)
call mma_deallocate(T2Name)

!Call EndGlb()

ireturn = 0

return

end subroutine cht3
