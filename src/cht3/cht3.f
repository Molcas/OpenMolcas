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
        subroutine cht3(ireturn)
c
c main driver for (T)
c
        implicit none
c
#include "cht3_casy.fh"
#include "cht3_ccsd1.fh"
#include "files.fh"
#include "cht3_reord.fh"
#include "WrkSpc.fh"
        integer ireturn
c DIRCC
cmp!        include 'memvir_inc'
#include "ccsd_t3compat.fh"
c
        integer wrksize
        integer maxspace
c
        integer isize
        integer iOE,ioeh,ioep
c       integer itmp,iW2,il1_1,ioff
c
        integer i,nOrbE,nBas(8),nOrb(8)
        character*24 Label
        logical Found
cmp
cmp
cmp
c       vynuluj hodiny
        Call CWTime(TCpu,TWall)
        TCpu0=TCpu
        TWall0=TWall
        TCpu_l=TCpu
        TCpu_l=TCpu
        TWall_l=TWall
        TWall_l=TWall
cmp
c.0        read input
c
        call IniReord_t3 (NvGrp,wrksize)
c
c.0.1        generate name convention for blocked integrals and T2 files
c
        call DefParReord_t3 (NvGrp,maxdim)
        if (printkey.ge.10) then
           write (6,*) 'Maxdim of virtual segment from CCSD = ',
     & maxdim
        end if
c
c.0.2        def commons for DIRCC
c
        call defcommon (nfr,no,nv)
c
c.2.2   regenerate integrals from the Cholesky vectors (L1) L1(m,I ,A'')
c
c.2.2.1 (vo|vo) for testing purpose
c
c        isize=nc*no*nv
c!        write (6,*) 'size for l1_1 ',isize
c!        write (6,*) 'size for l1_2 ',isize
c        Call GetMem('cht3_l1_1','Allo','Real',il1_1,isize)
c        Call GetMem('cht3_itmp','Allo','Real',itmp,isize)
c        isize=nv*nv*no*no
c!        write (6,*) 'size for W2 ',isize
c        Call GetMem('cht3_W2','Allo','Real',iW2,isize)
c        call gen_vvoo(Work(iW2),Work(il1_1),Work(itmp))
c
c.2.2.2 get orbital energies

c
c*      Get Oorital energies
c
        call Get_iArray('nBas',nBas,1)
        call Get_iArray('nOrb',nOrb,1)

        isize=nBas(1)

        if (printkey.ge.10) then
           write (6,*) 'Allocating memory for (tmp) OE files',
     &                  isize
        end if

        call GetMem('cht3_oe','Allo','Real',iOE,isize)
c
        Label='OrbE'
        Call qpg_dArray(Label,Found,nOrbE)
        if (nOrbE.ne.nBas(1)) then
          write (6,*) 'Warning! in cht3 : (nOrbE.ne.nBas)!'
        end if
        If(.not.Found .or. nOrbE.eq.0) Then
          Call SysAbendMsg('get_orbe','Did not find:',Label)
        End If
        if (printkey.ge.10) then
        write (6,*) 'isize = ',isize
        write (6,*) 'norbe = ',norbe
        end if
        call Get_dArray(Label,Work(iOE),nOrbE)
c
c        write out the orbital energies
c
        if (printkey.ge.10) then
           write (6,*)
           write (6,*) 'Orbital energies for nfr+no+nv'
           write (6,*)
           do i=1,nfr+no+nv
           write (6,'(A,2x,i5,2x,f18.10)') 'Orbital Energy ',
     &                                     i,Work(iOE+i-1)
           end do
        end if

c2.2.3        make OEH, OEP
c
        isize=2*no
         call GetMem('cht3_oeh','Allo','Real',ioeh,isize)
        isize=2*nv
         call GetMem('cht3_oeh','Allo','Real',ioep,isize)
c
        call generate_juzekOE (Work(ioe+nfr),
     & Work(ioeh),Work(ioep),no,nv)
c
c.2.3   Checkpoint. Calculate MP2 energy
c
c        call calc_MP2 (Work(iW2),Work(iOE+nfr),no,nv)
c        call abend()
c
c.3        start (T) calculation
c
        call GetMem('(T)','Max','Real',maxspace,maxspace)
c
        write (6,*)
        write (6,'(A,i13,A,f9.1,A,f5.1,A)')
     & ' Memory available for (T) calc = ',
     & (maxspace-1),' in r*8 Words',
     & ((maxspace-1)*8.0d0)/(1024*1024),' Mb',
     & ((maxspace-1)*8.0d0)/(1024*1024*1024),' Gb'
c
cmp        call GetMem('t3_ampl_bti','Allo','Real',ioff,1)
cmp        ioff=ioff+1
cmp!        write (6,*) 'ioe   = ',ioe
cmp!        write (6,*) 'ioeh  = ',ioeh
cmp!        write (6,*) 'ioep  = ',ioep
cmp        write (6,*) 'ioff volny  = ',ioff
cmp!        kvir1=ioff+1
cmp!        kvir2=kvir1+(maxspace-1)-1
c
c toto sa da vyhodit a nahradit iba natvrdo definovanim kvir
cmp        call adapt_mem(Work(1),Work(ioff),(maxspace-1),
cmp     & kvir1,kvir2)
c
c!        call alloc_vm(WORK, maxspace, KVIR1, KVIR2)
c
cmp        write (6,*) ' kvir1 = ',kvir1
cmp        write (6,*) ' kvir2 = ',kvir2
c
cmp        call T3AMPL_BTI(Work(ioff),Work(ioeh),Work(ioep))

        call T3AMPL_BTI(Work(ioeh),Work(ioep))
c
c.2.4        Free the unnecessary memory
c
cmp        isize=nfr+no+nv
cmp        isize=nOrb(1)
        isize=nBas(1)
        call GetMem('cht3_oeh','Free','Real',iOE,isize)
c
        isize=2*no
         call GetMem('cht3_oeh','Free','Real',ioeh,isize)
        isize=2*nv
         call GetMem('cht3_oeh','Free','Real',ioep,isize)
cmp        ioff=ioff-1
cmp        call GetMem('t3_ampl_bti','Free','Real',ioff,1)

c!        call GetMem('t3_ampl_bti','Free','Real',ioff,maxspace)
c
        !Call EndGlb

        ireturn = 0
        return
        end
c
c --------------
c
cmp!        subroutine adapt_mem(W1st,Wfree,maxsize,kvir1,kvir2)
cmp!c
cmp!        implicit none
cmp!        real*8 W1st(2)
cmp!        real*8 Wfree(2)
cmp!        integer locw
cmp!        integer maxsize
cmp!        integer kvir1,kvir2
cmp!        integer loc_1st,loc_free
cmp!#ifdef LOC
cmp!        integer loc
cmp!        external loc
cmp!#endif
cmp!        loc_1st=loc(W1st)
cmp!        loc_free=loc(Wfree)
cmp!        write(0,'(a,i21,2x,a)') 'FAKE ALLOC_VM,   W1st  =', loc_1st,
cmp!     &                          'byte address'
cmp!        write(0,'(a,i21,2x,a)') 'FAKE ALLOC_VM,   Wfree =', loc_free,
cmp!     &                          'byte address'
cmp!c
cmp!        kvir1= ( (loc_free - loc_1st)/8 + 1 )
cmp!        kvir2=kvir1+maxsize-1
cmp!c
cmp!      write(6,'(a,i21,2x,a)') 'FAKE ALLOC_VM, MAXSIZE =', maxsize,
cmp!     $     'real*8 words'
cmp!      write(6,'(a,i21,2x,a)') 'FAKE ALLOC_VM, KVIR1   =', KVIR1,
cmp!     $     'real*8 words'
cmp!      write(6,'(a,i21,2x,a)') 'FAKE ALLOC_VM, KVIR2   =', KVIR2,
cmp!     $     'real*8 words'
cmp!        return
cmp!        end
