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
      subroutine gugaci(ireturn)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "files_gugaci.fh"
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
      data istep_occ/0,1,1,2/
      data mul_tab/1,2,3,4,5,6,7,8,
     *             2,1,4,3,6,5,8,7,
     *             3,4,1,2,7,8,5,6,
     *             4,3,2,1,8,7,6,5,
     *             5,6,7,8,1,2,3,4,
     *             6,5,8,7,2,1,4,3,
     *             7,8,5,6,3,4,1,2,
     *             8,7,6,5,4,3,2,1/
!                           v   d  t  s dd tt
      data map_jplr/       25,23,17,10,24,18, !v
     :                     26,19,13, 6,22, 0, !d
     :                      0,14,11, 2, 0, 0, !t
     :                      0, 7, 3, 1, 9, 5, !s
     :                      0,21, 0, 8,20,15, !dd
     :                      0, 0, 0, 4,16,12/ !tt

      ireturn=100
      call version_info
      ! open scratch files and read from strach file
      call gugaciinit

      ! user input file
      logic_grad=.false.
      call mole_inf()

      call paras_calculate()
      call arrange_orbital()

      call allocate_casrst()
      call dbl_upwalk()         ! add by wyb 01.9.5
      call ext_downwalk()       ! add by wyb 01.9.5
      call active_drt()         ! add by wyb 01.9.5
      call value_of_pl_in_dbl()

      nc0=norb_all*(norb_all+1)/2
      nc1=nc0*(nc0+1)/2
      if(nc1.gt.2*max_vector) then
        write(6,*) "No enough space to store MO integrals! number of ",
     *             " orbital should be less than ",max_orb
#ifdef MOLPRO
#else
        call abend
#endif
#ifdef _XIANEST_
#endif
      endif

      call mem_intinnindex_alloc()
      lenvec=nc1
      allocate(vector1(nc1))
      vector1(1:nc1)=0.d0
      call int_sort()
      deallocate(vector1)

      mxvec=60000000*10
      if(mroot*nci_dim.le.mxvec) then
        nc1=mroot*nci_dim
        allocate(vector1(nc1))
!        allocate(vector2(nc1))
        lenvec=nc1
      else
        nc1=nci_dim
        allocate(vector1(nc1))
!        allocate(vector2(nc1))
        lenvec=nc1
      endif
      vector1(1:nc1)=0.d0

!      write(6,*)   '===================================================

      sc1=c_time()

      call allocate_subdrt(1,1)
      call allocate_subdrtl(1,1)

      call memcidiag_alloc()
      call diagonal_loop_wyb()
      call memcidiag_dealloc()
      sc2=c_time()
      !write(6,*) vector1(585)
      !stop 888

      write(6,*)
      write(6,*)   '==================================================='
      write(6,'(a30,i10,f14.2,a1)')
     :     '   end of h_diagonal, nci_dim=',nci_dim,sc2-sc1,'s'
      write(6,*)   '==================================================='
      write(6,*)


      call write_ml(lucidia,1,vector1,nci_dim,1)
      !do i=1,74902
      !  write(21,"(i8,f18.8)") i,vector1(i)
      !enddo
      !stop 888

      call allocate_vplp_memory()
      call allocate_int_memory()

      nc=nci_h0 !iw_sta(2,1)
      nc1=nc*(nc+1)/2

      allocate(vcm(mroot*nci_h0))
      if(nc1.le.lenvec) then
        allocate(vector2(lenvec))
      else
        deallocate(vector1)
        allocate(vector1(nc*(nc+1)/2))
        allocate(vector2(nc*(nc+1)/2))
        vector1=0.d0
        call read_ml(Lucidia,2,vector1,nci_dim,1)
      endif

      vector2(1:nc1)=0.d0

      call geth0()
      call xflush(6)

      if(nc1.gt.lenvec) then
        deallocate(vector1)
        deallocate(vector2)
        allocate(vector1(lenvec))
        allocate(vector2(lenvec))
      endif

      sc0=c_time()
      call guga_ploop(npl,maxplcon)

      call deallocate_subdrt()
      call deallocate_subdrtl()


      sc1=c_time()
      call xflush(6)

      write(6,'(a25,i10,f14.2,a1)')
     :          '  end of pl_serach, n_pl=',npl,sc1-sc0,'s'
      write(6,*)'=============================================='

      if(maxplcon.lt.max_node) maxplcon=max_node
      call allocate_subdrt(2,maxplcon)
      call allocate_subdrtl(2,maxplcon)

      call cidiagonalize(mxvec)

      sc2=c_time()

      call xflush(6)

      write(6,910)  sc2-sc0
      write(6,*)
      call deallocate_int_memory()
      !deallocate(vcm)
      deallocate(vector1)
      deallocate(vector2)

      if(logic_calpro) then
        logic_grad=.true.
        call memdengrad_alloc()

        nc0=norb_all*(norb_all+1)/2
        nc1=nc0*(nc0+1)/2
        allocate(vector1(nci_dim))
        allocate(vector2(nc1))
        vector1=0.d0
        vector2=0.d0

        call cidenmat()

! calculate properties
#ifdef MOLPRO
#else
        call cipro()
#endif
        deallocate(vector1)
        deallocate(vector2)
        call memdengrad_free()
      endif

      call deallocate_vplp_memory()
      call deallocate_subdrt()
      call deallocate_subdrtl()
      call deallocate_casrst()
      call mem_intinnindex_dealloc()
      call gugafinalize

      ireturn=0
910   format(/,1x,"end of ci energy calculation,
     *       takes ",f10.2,1x,"seconds"/)
      end
