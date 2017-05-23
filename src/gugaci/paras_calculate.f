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
      subroutine paras_calculate()
#include "drt_h.fh"
#include "paraconstants_h.fh"
#include "intsort_h.fh"
c      data inlptb_new/
!         1  2  3  4  5  6  7  8  9  10  11  12 13
c a^r=1
c     *    -1, 0, 0, 0, 0, 0,-7, 0,-9,-10,  0,-12, 0,
c a^l=2
c     *     0, 0, 0, 0, 0,-6, 0,-8, 0,  0,-11,  0, 0,
c b_r=3
c     *    4, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c b_l=4
c     *    5, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c b^r=5
c     *    0, 7, 8,10,11, 0, 0, 0, 0,  0,  0,  0, 0,
c b^l=6
c     *    0, 0, 9, 0,12, 0, 0, 0, 0,  0,  0,  0, 9,
c c^'=7
c     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3,
c c^"=8
c     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,
c d_l^r=9
c     *    6, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c d^r^r=10
c     *    0,-2, 0,-4, 0, 0, 0, 0, 0,  0,  0,  0, 0,
c d^r^l=11
c     *    0, 0,-3, 0,-5, 0, 0, 0, 0,  0,  0,  0, 0,
c c^'" =12
c     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0/
!
!****************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl     =13 (*)
!****************************************************************

!     data inlptb_new
!    */ -1,  0,  4,  5,  0,  0,  1,  1,  6,  0,  0,  1,
!    *   0,  0,  0,  0,  7,  0,  2,  2,  0, -2,  0,  2,
!    *   0,  0,  0,  0,  8,  9,  3,  3,  0,  0, -3,  3,
!    *   0,  0,  0,  0, 10,  0,  4,  4,  0, -4,  0,  4,
!    *   0,  0,  0,  0, 11, 12,  5,  5,  0,  0, -5,  5,
!    *   0, -6,  0,  0,  0,  0,  6,  6,  0,  0,  0,  6,
!    *  -7,  0,  0,  0,  0,  0,  7,  7,  0,  0,  0,  7,
!    *   0, -8,  0,  0,  0,  0,  8,  8,  0,  0,  0,  8,
!    *  -9,  0,  0,  0,  0,  0,  9,  9,  0,  0,  0,  9,
!    * -10,  0,  0,  0,  0,  0, 10, 10,  0,  0,  0, 10,
!    *   0,-11,  0,  0,  0,  0, 11, 11,  0,  0,  0, 11,
!    * -12,  0,  0,  0,  0,  0, 12, 12,  0,  0,  0, 12,
!    *   0,  0,  0,  0,  0,  9,  3, 13,  0,  0,  0,  0/

!     data indord_plptype/0,0,0,1,1,3,3,1,1,2,2,2/   !severe_new_error_1
      dimension iwtsmsabtmp(maxgdm)
      ja_sys=int(n_electron*0.5d0-spin)-norb_dz
      jb_sys=int(spin+spin)
      jc_sys=norb_all-ja_sys-jb_sys
      jsm_sys=ns_sm
c      write(6,'(a9,3(1x,i3))') "ja,jb,jc ",ja_sys,jb_sys,jc_sys
!     if ( jb_sys.eq.jb_sys/2*2 ) then    ????? 11.2
!        w0_drl2_44=v_sqtwo
!        w0_drl2_1122=v_onevsqtwo
!     else
!        w0_drl2_44=-v_sqtwo
!        w0_drl2_1122=v_onevsqtwo
!     endif

      if ( jb_sys .eq.0 ) jroute_sys=1
      if ( jb_sys .eq.1 ) jroute_sys=2
      if ( jb_sys .gt.1 ) jroute_sys=3

      do i=1,ng_sm
      iwt_ext(i,1)=0
      iwt_ext(i,2)=nlsm_ext(i)
      enddo
      iwt_ext(1,1)=1
!
      do i=1,ng_sm
         ni=nlsm_ext(i)
      do j=i,ng_sm
         nj=nlsm_ext(j)
         if ( i.eq.j) then
            nij=ni*(ni-1)/2
         else
            nij=ni*nj
         endif
         ijsm=mul_tab(i,j)
         iwt_ext(ijsm,3)=iwt_ext(ijsm,3)+nij
         iwt_ext(ijsm,4)=iwt_ext(ijsm,4)+nij
      enddo
      enddo
      iwt_ext(1,4)=iwt_ext(1,4)+norb_ext

      do i=1,ng_sm
      iwt_dbl(i,1)=0
      iwt_dbl(i,2)=nlsm_dbl(mul_tab(i,jsm_sys))
      if ( jroute_sys .gt. 1 ) then
         iwt_dbl(i,3)=nlsm_dbl(mul_tab(i,jsm_sys))
      endif
      enddo
      iwt_dbl(ns_sm,1)=1
!
      do i=1,ng_sm
         ni=nlsm_dbl(i)
      do j=i,ng_sm
         nj=nlsm_dbl(j)
         if ( i.eq.j) then
            nij=ni*(ni-1)/2
         else
            nij=ni*nj
         endif
         ijsm=mul_tab(mul_tab(i,j),ns_sm)
         iwt_dbl(ijsm,4)=iwt_dbl(ijsm,4)+nij
      if ( jroute_sys .gt. 2 ) then
         iwt_dbl(ijsm,5)=iwt_dbl(ijsm,5)+nij
      endif
         iwt_dbl(ijsm,6)=iwt_dbl(ijsm,6)+nij
      if ( jroute_sys .gt. 1 ) then
         iwt_dbl(ijsm,6)=iwt_dbl(ijsm,6)+nij
      endif
      enddo
      enddo
      iwt_dbl(ns_sm,6)=iwt_dbl(ns_sm,6)+norb_dbl

!
      do ismb=1,ng_sm
         iwt_sm_sab(ismb)=0
      enddo
      do ismb=1,ng_sm
         ni=nlsm_ext(ismb)
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
         do isma=1,ismb
            nj=nlsm_ext(isma)
            if ( isma.eq.ismb) then
               nij=ni*(ni-1)/2
            else
               nij=ni*nj
            endif
            ijsm=mul_tab(isma,ismb)
!           iwt_sm_sab(ijsm)=iwt_sm_sab(ijsm)+nij
!           iwt_sm_sab fuction as a tmp array
            iwt_sm_ext(isma,ismb)=iwt_sm_sab(ijsm)

            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) then
               ibsta=iasta+1
            endif
            iwttmp=iwt_sm_sab(ijsm)
            do iborb=ibsta,ibend
         do iaorb=iasta,min(iaend,iborb-1)
            iwttmp=iwttmp+1
!     if ( iaorb.le.0.or.iaorb.gt.max_extorb ) then
!     write(6,*)
!     endif
!     if ( iborb.le.0.or.iborb.gt.max_extorb ) then
!     write(6,*)
!     endif
            iwt_orb_ext(iaorb,iborb)=iwttmp
         enddo
         enddo

         iwt_sm_sab(ijsm)=iwt_sm_sab(ijsm)+nij
      enddo
      enddo
      iwt_sm_s_ext=iwt_sm_sab(1)

      isumtmp=0
      do ismb=1,ng_sm
         icnttmp=iwt_sm_sab(ismb)
         iwt_sm_sab(ismb)=isumtmp
         isumtmp=isumtmp+icnttmp
      enddo
!     imap_revorbtoorb_ext()
      do ismb=1,ng_sm
         ibsta=ibsm_ext(ismb)
         ibend=iesm_ext(ismb)
      do isma=1,ismb
         iasta=ibsm_ext(isma)
         iaend=iesm_ext(isma)
         if ( isma .eq. ismb ) iaend=iaend-1

         ismab=mul_tab(isma,ismb)
         iwttmp=iwt_sm_ext(isma,ismb)+iwt_sm_sab(ismab)
         do iaorb=iasta,iaend
         do iborb=max(iaorb+1,ibsta),ibend
            iwttmp=iwttmp+1
            iwt_revorb_ext(iaorb,iborb)=iwttmp
            imap_revorbtoorb_ext(iwttmp)=iwt_orb_ext(iaorb,iborb)
         enddo
         enddo
      enddo
      enddo

      do ismb=1,ng_sm
         iwt_sm_sab(ismb)=0
      enddo
      do ismb=1,ng_sm
         ni=nlsm_dbl(ismb)
         ibsta=ibsm_dbl(ismb)
         ibend=iesm_dbl(ismb)
      do isma=1,ismb
         nj=nlsm_dbl(isma)
         if ( isma.eq.ismb) then
            nij=ni*(ni-1)/2
         else
            nij=ni*nj
         endif
         ijsm=mul_tab(isma,ismb)
         iwt_sm_dbl(isma,ismb)=iwt_sm_sab(ijsm)

         iasta=ibsm_dbl(isma)
         iaend=iesm_dbl(isma)
         if ( ismb .eq. isma ) then
            ibsta=iasta+1
         endif
         iwttmp=iwt_sm_sab(ijsm)
         do iborb=ibsta,ibend
         do iaorb=iasta,min(iaend,iborb-1)
            iwttmp=iwttmp+1
            iwt_orb_dbl(iaorb,iborb)=iwttmp
         enddo
         enddo

         iwt_sm_sab(ijsm)=iwt_sm_sab(ijsm)+nij
      enddo
      enddo

      isumtmp=0
      do ismb=1,ng_sm
         icnttmp=iwt_sm_sab(ismb)
         iwtsmsabtmp(ismb)=isumtmp
         isumtmp=isumtmp+icnttmp
      enddo
!     imap_revorbtoorb_dbl()
      do ismb=1,ng_sm
         ibsta=ibsm_dbl(ismb)
         ibend=iesm_dbl(ismb)
      do isma=1,ismb
         iasta=ibsm_dbl(isma)
         iaend=iesm_dbl(isma)
         if ( isma .eq. ismb ) iaend=iaend-1
         ismab=mul_tab(isma,ismb)
         iwttmp=iwt_sm_dbl(isma,ismb)+iwtsmsabtmp(ismab)
         do iaorb=iasta,iaend
         do iborb=max(iaorb+1,ibsta),ibend
            iwttmp=iwttmp+1
            iwt_revorb_dbl(iaorb,iborb)=iwttmp
            imap_revorbtoorb_dbl(iwttmp)=iwt_orb_dbl(iaorb,iborb)
         enddo
         enddo
      enddo
      enddo

      do ism=1,ng_sm
         jsm=mul_tab(ism,ns_sm)
         if ( ism .lt. jsm ) then
            itmp=iwt_sm_sab(jsm)
            iwt_sm_sab(jsm)=iwt_sm_sab(ism)
            iwt_sm_sab(ism)=itmp
         endif
      enddo

      iwt_sm_s_dbl=iwt_sm_sab(ns_sm)
      if ( jroute_sys .gt. 2 ) then
         iwt_sm_s_dbl=iwt_sm_s_dbl*2
      endif
      end
