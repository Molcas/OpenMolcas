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
c****************************************************
      subroutine arrange_orbital()
c****************************************************
c  arrange orbital for ci calculation, in meld and
c  molcas program, the orbital are arranged as the symmetry
c  block. we transfer them to ci order
c -----
c map_order_orbital    ab ---> ci
c
#include "drt_h.fh"
#include "intsort_h.fh"
      common /mcorb/ lsmorb(max_orb),noidx(8)
      dimension lsmorbcount(ng_sm),map_tmp(max_orb)
      logical logi_norb_inn(norb_all)

      logi_norb_inn(1:norb_all)=.false.
      iorb=norb_all
      do la=1,norb_all
         norb_number(la)=iorb
         iorb=iorb-1
      enddo

      nim=0
      lsmorbcount(1)=nim
      do im=2,ng_sm
        nim=nim+nlsm_all(im-1)
        lsmorbcount(im)=nim
      enddo

      if ( logic_assign_actorb ) then
        do lr=1,norb_inn
          lr_scf=map_orb_order(lr)
          logi_norb_inn(lr_scf)=.true.
        enddo
        goto 200
      else
        do lr=1,norb_inn
          lsmr=lsm_inn(lr)
          lsmorbcount(lsmr)=lsmorbcount(lsmr)+1
          lr_scf=lsmorbcount(lsmr)
          map_orb_order(lr)=lr_scf
          logi_norb_inn(lr_scf)=.true.
        enddo
      endif
200   lr_scf0=norb_all
      la=norb_inn+1
      do ms=ng_sm,1,-1
        lr_scf0=lr_scf0-nlsm_all(ms)
        lr_scf=lr_scf0
        do lra=1,nlsm_all(ms)
          lr_scf=lr_scf+1
          if(logi_norb_inn(lr_scf)) cycle
          map_orb_order(la)=lr_scf
          la=la+1
        enddo
      enddo

      isum2=0
      isum3=0
      do i=1,ng_sm
        jp2(i)=isum2
        isum2=isum2+i
        jp3(i)=isum3
        isum3=isum3+isum2
      enddo

      iccount=1
      do lrd = 1,norb_inn
        ipwt(lrd) = iccount
        iccount   = iccount+2
      enddo
      lsmorbcount =0
      do lrd=norb_dz,1,-1
        lsmid=lsm_inn(lrd)
        lsmorbcount(lsmid)=lsmorbcount(lsmid)+1
        ipws(lrd)=(lsmorbcount(lsmid)-1)*3+1
      enddo


      map_tmp(1:norb_all)=map_orb_order(1:norb_all)
      do 40 i=1,norb_all
        do j=1,norb_all
          if(map_tmp(j).eq.i) then
            map_orb_order(i)=j
            goto 40
          endif
        enddo
40    continue
      return
      end
