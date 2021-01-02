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
      subroutine inn_ext_ss_drl_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) ' ss_test 2/2'
      ii=1
      ii0=ii
      if ( logic_g2g4a ) then
         ii=ii0
         mm0=iplplwei
         nn0=iplprwei+iwt_sm_s_ext
         mm=mm0
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1

               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii+1)
                  indexlp1=index_lpext1(ii+1)

!         if(indexlp.ne.0) then
           valuelp=value_lpext(ii+1)
           vector2(indexlp)=vector2(indexlp)
     :       +vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp
!         end if
!         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(ii+1)
           vector2(indexlp1)=vector2(indexlp1)
     :       +vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp1
!         end if
                  ii=ii+2
               enddo
            enddo
         enddo
      endif
      if ( logic_g2g4b ) then
         ii=ii0
         mm0=iplprwei
         nn0=iplplwei+iwt_sm_s_ext
         mm=mm0
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii+1)
                  indexlp1=index_lpext1(ii+1)
!         if(indexlp.ne.0) then
           valuelp=value_lpext(ii+1)
           vector2(indexlp)=vector2(indexlp)
     :       +vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp
!         end if
!         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(ii+1)
           vector2(indexlp1)=vector2(indexlp1)
     :       +vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp1
!         end if
                  ii=ii+2
               enddo
            enddo
         enddo
      endif
!     iaddii=(ii-ii0)/2
      ii0=ii-1

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(ii)
               indexlp1=index_lpext1(ii)
!         if(indexlp.ne.0) then
           valuelp=value_lpext(ii)
           vector2(indexlp)=vector2(indexlp)
     :       +vector1(mm)*vector1(nn)*valuelp
!         end if
!         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(ii)
           vector2(indexlp1)=vector2(indexlp1)
     :       +vector1(mm)*vector1(nn)*valuelp1
!         end if

               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
      else
            goto 935
      endif

         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
               indexlp1=index_lpext1(iij)
!         if(indexlp.ne.0) then
           valuelp=value_lpext(iij)
           vector2(indexlp)=vector2(indexlp)
     :       +vector1(mm)*vector1(nn)*valuelp
!         end if
!         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(iij)
           vector2(indexlp1)=vector2(indexlp1)
     :       +vector1(mm)*vector1(nn)*valuelp1
!         end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
935      continue
!     cycle

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
               indexlp1=index_lpext1(iij)
!         if(indexlp.ne.0) then
           valuelp=value_lpext(iij)
           vector2(indexlp)=vector2(indexlp)
     :       +vector1(mm)*vector1(nn)*valuelp
!         end if
!         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(iij)
           vector2(indexlp1)=vector2(indexlp1)
     :       +vector1(mm)*vector1(nn)*valuelp1
!         end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue

      enddo
      end


      subroutine inn_ext_st_drl_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) ' st_test 2/2'

      ii=1
      if ( logic_g1415 ) then
      mm=iplplwei
      nn=iplprwei
      do i=1,idownwei_g131415
         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
         indexlp1=index_lpext1(ii)
!         if(indexlp.ne.0) then
           valuelp=value_lpext(ii)
           vector2(indexlp)=vector2(indexlp)
     :       +vector1(mm)*vector1(nn)*valuelp
!         end if
!         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(ii)
           vector2(indexlp1)=vector2(indexlp1)
     :       +vector1(mm)*vector1(nn)*valuelp1
!         end if
         ii=ii+1
      enddo
      endif

      ii0=ii
      if ( logic_g2g4a ) then
         ii=ii0
         mm0=iplplwei
         nn0=iplprwei+iwt_sm_s_ext
            mm=mm0            !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1

               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1

                  indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!         end if
                  ii=ii+1

                  indexlp=index_lpext(ii)
!          if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!          end if
                  ii=ii+1
               enddo
            enddo
         enddo
      endif
      if ( logic_g2g4b ) then
         ii=ii0
         mm0=iplprwei
         nn0=iplplwei+iwt_sm_s_ext
            mm=mm0             !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1

                  indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!         end if
                  ii=ii+1

                  indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!         end if
                  ii=ii+1

               enddo
            enddo
         enddo
      endif
      ii0=ii-1
      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)

            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
               valuelp=value_lpext(ii)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!         end if
               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue
      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!          if(indexlp.ne.0) then
               valuelp=value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!          end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!           if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!           end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      endif

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!           if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!           end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue
      enddo

      end

      subroutine inn_ext_tt_drl_loop_unpack_g(iplplwei,iplprwei,n1415)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      logical logic_g14150,logic_g34b0,logic_g35b0,logic_g36b0
c     write(*,*) '  tt_test 2/2'
      logic_g14150=logic_g1415
      logic_g36b0=logic_g36b
      logic_g35b0=logic_g35b
      logic_g34b0=logic_g34b
      if(iplplwei.eq.iplprwei) then
         logic_g14150=.false.
         logic_g36b0=.false.
         logic_g35b0=.false.
         logic_g34b0=.false.
      endif
      ii=1
      if ( logic_g14150 ) then
      mm=iplplwei
      nn=iplprwei
!     ii=iista
      do i=1,idownwei_g131415
         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
         ii=ii+1
      enddo
      endif

      ii0=n1415         !severe_error

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b0 ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
               valuelp=value_lpext(ii)
               vector2(indexlp)=vector2(indexlp)
     :            +vector1(mm)*vector1(nn)*valuelp
!       end if
               indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
               valuelp1=value_lpext1(ii)
               vector2(indexlp1)=vector2(indexlp1)
     :            +vector1(mm)*vector1(nn)*valuelp1
!       end if
               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
      elseif (  icle.eq.2 .and. logic_g35b0 ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
      else
            goto 935
      endif

         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :            +vector1(mm)*vector1(nn)*valuelp
!       end if
               indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
               valuelp1=-value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :            +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
935      continue
!     cycle

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b0 ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
               valuelp=value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :            +vector1(mm)*vector1(nn)*valuelp
!       end if
               indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
               valuelp1=value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :            +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo

         enddo
934   continue

      enddo
      end

      subroutine inn_ext_ts_drl_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) '  ts_test 2/2'
      ii=1
c      logic_g1415=.false.
c      logic_g2g4a=.false.
c      logic_g2g4b=.false.
c      logic_g36a=.false.
c      logic_g36b=.false.
c      logic_g35a=.false.
c      logic_g35b=.false.
c      logic_g34a=.false.
c      logic_g34b=.false.
      if ( logic_g1415 ) then
      mm=iplplwei
      nn=iplprwei
      do i=1,idownwei_g131415
         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
!        if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!        end if
         indexlp1=index_lpext1(ii)
!        if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!        end if
         ii=ii+1
      enddo
      endif
      ii0=ii
      if ( logic_g2g4a ) then
         ii=ii0
         mm0=iplplwei
         nn0=iplprwei+iwt_sm_s_ext
            mm=mm0            !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!        if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!        end if
                  ii=ii+1
                  indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!         end if
                  ii=ii+1
               enddo
            enddo
         enddo
      endif
      if ( logic_g2g4b ) then
         ii=ii0
         mm0=iplprwei
         nn0=iplplwei+iwt_sm_s_ext
            mm=mm0             !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!         end if
                  ii=ii+1
                  indexlp=index_lpext(ii)
!         if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!         end if
                  ii=ii+1
               enddo
            enddo
         enddo
      endif

      ii0=ii-1         !severe_error

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(ii)
!          if(indexlp.ne.0) then
               valuelp=value_lpext(ii)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!          end if
               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!           if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!           end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!           if(indexlp.ne.0) then
               valuelp=value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!           end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      endif

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!           if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!           end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue
      enddo
      end


      subroutine inn_ext_ss_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) ' ss_test 1/2'
      ii=1
      if ( logic_g1415 ) then
      mm=iplplwei
      nn=iplprwei
      do i=1,idownwei_g131415

         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
         ii=ii+1

         indexlp=index_lpext(ii)
       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
        vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
       end if
         indexlp1=index_lpext1(ii)
       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
       end if
         ii=ii+1
      enddo
      endif

      ii0=ii
      if ( logic_g2g4a ) then
         ii=ii0
         mm0=iplplwei
         nn0=iplprwei+iwt_sm_s_ext
         mm=mm0
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!            if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!            end if
                  indexlp1=index_lpext1(ii)
!            if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nna)*valuelp1
!            end if
                  ii=ii+1
                  indexlp=index_lpext(ii)
!            if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!            end if
                  indexlp1=index_lpext1(ii)
!            if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nnb)*valuelp1
!            end if
                  ii=ii+1

               enddo
            enddo
         enddo
      endif
      if ( logic_g2g4b ) then
         ii=ii0
         mm0=iplprwei
         nn0=iplplwei+iwt_sm_s_ext
         mm=mm0
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!            if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!            end if
                  indexlp1=index_lpext1(ii)
!            if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nnb)*valuelp1
!            end if
                  ii=ii+1
                  indexlp=index_lpext(ii)
!            if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!            end if
                  indexlp1=index_lpext1(ii)
!            if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nna)*valuelp1
!            end if
                  ii=ii+1

               enddo
            enddo
         enddo
      endif
!     iaddii=(ii-ii0)/2
      ii0=ii-1

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(ii)
!            if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
                  indexlp1=index_lpext1(ii)
!            if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!            end if

               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
      else
            goto 935
      endif

         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!            if(indexlp.ne.0) then
               valuelp=value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
               indexlp1=index_lpext1(iij)
!            if(indexlp1.ne.0) then
               valuelp1=value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!            end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
935      continue

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!            if(indexlp.ne.0) then
               valuelp=value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
               indexlp1=index_lpext1(iij)
!            if(indexlp1.ne.0) then
               valuelp1=value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!            end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue
      ii0=ii0+nvalue_space_ss
      enddo
      end


      subroutine inn_ext_st_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) ' st_test 1/2'
      ii=1
      if ( logic_g1415 ) then
      mm=iplplwei
      nn=iplprwei
      do i=1,idownwei_g131415
         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
         ii=ii+1

         indexlp=index_lpext(ii)
       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
       end if
         indexlp1=index_lpext1(ii)
       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
       end if
         ii=ii+1
      enddo
      endif


      ii0=ii
      if ( logic_g2g4a ) then
         ii=ii0
         mm0=iplplwei
         nn0=iplprwei+iwt_sm_s_ext
            mm=mm0            !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nna)*valuelp1
!       end if
                  ii=ii+1

                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nnb)*valuelp1
!       end if
                  ii=ii+1
               enddo
            enddo
         enddo
      endif

      if ( logic_g2g4b ) then
         ii=ii0
         mm0=iplprwei
         nn0=iplplwei+iwt_sm_s_ext
            mm=mm0
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)

                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
                  valuelp1=value_lpext1(ii)
!       if(indexlp1.ne.0) then
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nnb)*valuelp1
!       end if
                  ii=ii+1

                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nna)*valuelp1
!       end if
                  ii=ii+1

               enddo
            enddo
         enddo
      endif
      ii0=ii-1

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

c      write(*,*)'st_g35a',iplplwei,iplprwei,vector2(137)
      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
               valuelp=value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
               indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
               valuelp1=value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
               indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
               valuelp1=-value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      endif

c      write(*,*)'st_g34',iplplwei,iplprwei,vector2(137)
      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
               indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
               valuelp=-value_lpext(iij)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
               indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
               valuelp1=-value_lpext1(iij)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue
      ii0=ii0+nvalue_space_ss
      enddo
      end


      subroutine inn_ext_ts_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) '  ts_test 1/2 '
      ii=1
      if ( logic_g1415 ) then
      mm=iplplwei
      nn=iplprwei
      do i=1,idownwei_g131415
         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
         ii=ii+1

         indexlp=index_lpext(ii)
       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
       end if
         indexlp1=index_lpext1(ii)
       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
       end if
         ii=ii+1
      enddo
      endif

      ii0=ii
      if ( logic_g2g4a ) then
         ii=ii0
         mm0=iplplwei
         nn0=iplprwei+iwt_sm_s_ext
            mm=mm0            !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nna)*valuelp1
!       end if
                  ii=ii+1

                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nnb)*valuelp1
!       end if
                  ii=ii+1

               enddo
            enddo
         enddo
      endif
      if ( logic_g2g4b ) then
         ii=ii0
         mm0=iplprwei
         nn0=iplplwei+iwt_sm_s_ext
            mm=mm0             !severe_error
         do ismb=1,ng_sm
            isma=mul_tab(ismb,ism_g2g4)
            if ( isma .gt. ismb ) cycle
            ibsta=ibsm_ext(ismb)
            ibend=iesm_ext(ismb)
            iasta=ibsm_ext(isma)
            iaend=iesm_ext(isma)
            if ( ismb .eq. isma ) ibsta=ibsta+1
            nna=nn0+ibsta-1
            do ib=ibsta,ibend
               nna=nna+1
               nnb=nn0+iasta-1
               do ia=iasta,min(iaend,ib-1)
                  mm=mm+1
                  nnb=nnb+1
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nnb)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nnb)*valuelp1
!       end if
                  ii=ii+1

                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nna)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nna)*valuelp1
!       end if
                  ii=ii+1

               enddo
            enddo
         enddo
      endif
      ii0=ii-1

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
                  valuelp=-value_lpext(iij)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
                  valuelp1=-value_lpext1(iij)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)
            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(iij)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(iij)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
      endif

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
                  valuelp=-value_lpext(iij)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
                  valuelp1=-value_lpext1(iij)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue
      ii0=ii0+nvalue_space_ss
      enddo
      end


      subroutine inn_ext_tt_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
c     write(*,*) '  tt_test 1/2'
      ii=1
      if ( logic_g1415 ) then
      mm=iplplwei
      nn=iplprwei

      do i=1,idownwei_g131415
         mm=mm+1
         nn=nn+1
         indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if

         ii=ii+1

         indexlp=index_lpext(ii)
       if(indexlp.ne.0) then
         valuelp=value_lpext(ii)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
       end if
         indexlp1=index_lpext1(ii)
       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(ii)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
       end if
         ii=ii+1
      enddo
      endif

      ii0=ii-1

      do icle=1,2
      if ( icle.eq.1 .and. logic_g36a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta36=lpsta36a
         lpend36=lpend36a
      elseif (  icle.eq.2 .and. logic_g36b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta36=lpsta36b
         lpend36=lpend36b
      else
            goto 936
      endif
         do iii=lpsta36,lpend36,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            ii=ii0+lpext_wei(iii+2)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(ii)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(ii)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(ii)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(ii)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if
               mm=mm+1
               nn=nn+1
            enddo
         enddo
936   continue

      if ( icle.eq.1 .and. logic_g35a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta35=lpsta35a
         lpend35=lpend35a
      elseif (  icle.eq.2 .and. logic_g35b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta35=lpsta35b
         lpend35=lpend35b
      else
            goto 935
      endif

         do iii=lpsta35,lpend35,4
            mm=mm0+lpext_wei(iii)
            nn=nn0+lpext_wei(iii+1)
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
                  valuelp=-value_lpext(iij)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
                  valuelp1=-value_lpext1(iij)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if

               iij=iij+1
               mm=mm+1
            enddo
         enddo
935      continue
!     cycle

      if ( icle.eq.1 .and. logic_g34a ) then
         mm0=iplplwei
         nn0=iplprwei
         lpsta34=lpsta34a
         lpend34=lpend34a
      elseif (  icle.eq.2 .and. logic_g34b ) then
         mm0=iplprwei
         nn0=iplplwei
         lpsta34=lpsta34b
         lpend34=lpend34b
      else
            goto 934
      endif
         do iii=lpsta34,lpend34,4
            ilwtmp=lpext_wei(iii)
            irwtmp=lpext_wei(iii+1)
            mm=mm0+ilwtmp
            nn=nn0+irwtmp
            iij=ii0+lpext_wei(iii+2)

            do i=1,lpext_wei(iii+3)
                  indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
                  valuelp=value_lpext(iij)
                  vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
                  indexlp1=index_lpext1(iij)
!       if(indexlp1.ne.0) then
                  valuelp1=value_lpext1(iij)
                  vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!       end if

               iij=iij+1
               mm=mm+1
            enddo
         enddo
934   continue
      ii0=ii0+nvalue_space_ss
      enddo
      end

      subroutine gsd_sequence_extspace_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      parameter (v_sqtwo=1.414213562373095d0 )
c      write(*,*) '  sd_test 1/2','  ds_test 1'

      ilpvalue=0
      if ( logic_g25a ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         do itmp=1,nint_g25
            ilpvalue=ilpvalue+1
            indexlp=index_lpext(ilpvalue)
            valuelp=value_lpext(ilpvalue)
            indexlp1=index_lpext1(ilpvalue)
            valuelp1=value_lpext1(ilpvalue)
            nn=nn0
            do i=1,nwei_g25
               mm=mm+1
               nn=nn+1
!            if(indexlp.ne.0) then
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
              if(indexlp1.ne.0) then
                vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
              end if
            enddo
         enddo
      elseif ( logic_g25b ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         ilpvalue=ilpvalue+1
         do itmp=2,nint_g25
            ilpvalue=ilpvalue+1
            indexlp=index_lpext(ilpvalue)
            valuelp=value_lpext(ilpvalue)
            indexlp1=index_lpext1(ilpvalue)
            valuelp1=value_lpext1(ilpvalue)
            nn=nn0
            do i=1,itmp-1
               mm=mm+1
               nn=nn+1
!            if(indexlp.ne.0) then
              vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
              if(indexlp1.ne.0) then
                vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
              end if
            enddo
         enddo

         mm=iplplwei+iweista_g28-1
         nn=iplprwei
         nn=nn+1
         do itmp=2,nwei_g28
            nn=nn+1
            ilpvalue=0
            do i=1,itmp-1
               ilpvalue=ilpvalue+1
               mm=mm+1
               indexlp=index_lpext(ilpvalue)
!            if(indexlp.ne.0) then
               valuelp=value_lpext(ilpvalue)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
              indexlp1=index_lpext1(ilpvalue)
              if(indexlp1.ne.0) then
                valuelp1=value_lpext1(ilpvalue)
                vector2(indexlp1)=vector2(indexlp1)
     :                    +vector1(mm)*vector1(nn)*valuelp1
              end if
           enddo
         enddo
      elseif ( logic_g28a ) then
         mm=iplplwei+iweista_g28-1
         nn0=iplprwei
         do nn=nn0+1,nn0+nwei_g28
            ilpvalue=0
            do i=1,nint_g28
               ilpvalue=ilpvalue+1
               mm=mm+1
               indexlp=index_lpext(ilpvalue)
!            if(indexlp.ne.0) then
               valuelp=value_lpext(ilpvalue)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
              indexlp1=index_lpext1(ilpvalue)
              if(indexlp1.ne.0) then
                valuelp1=value_lpext1(ilpvalue)
                vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
              end if
            enddo
          enddo
      endif

      if ( logic_g26 ) then
         ilpvalue=ivaluesta_g26
         mm=iplplwei+iweista_g26
         nn0=iplprwei
         do nn=nn0+1,nn0+nwei_g26
           ilpvalue=ilpvalue+1
           indexlp=index_lpext(ilpvalue)
!          if(indexlp.ne.0) then
           valuelp=value_lpext(ilpvalue)*v_sqtwo
           vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!          end if
           indexlp1=index_lpext1(ilpvalue)
           if(indexlp1.ne.0) then
             valuelp1=value_lpext1(ilpvalue)*v_sqtwo
             vector2(indexlp1)=vector2(indexlp1)
     :                  +vector1(mm)*vector1(nn)*valuelp1
           end if
           mm=mm+1
         enddo
      endif
c      print*, "out sd 1"

      end

      subroutine gtd_sequence_extspace_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"
      common /cont_tmp/icount_ext
      parameter (   v_sqtwo=1.414213562373095d0 )
c     write(*,*) ' td_test _1/2',' dt_test '
      ilpvalue=0
      if ( logic_g25a ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         do itmp=1,nint_g25
            ilpvalue=ilpvalue+1
            indexlp=index_lpext(ilpvalue)
            valuelp=value_lpext(ilpvalue)
            indexlp1=index_lpext1(ilpvalue)
            valuelp1=value_lpext1(ilpvalue)

            nn=nn0
            do i=1,nwei_g25
               mm=mm+1
               nn=nn+1
!            if(indexlp.ne.0) then
              vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
            if(indexlp1.ne.0) then
              vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
            end if
           enddo
         enddo
      elseif ( logic_g25b ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         ilpvalue=ilpvalue+1
         do itmp=2,nint_g25
            ilpvalue=ilpvalue+1
            indexlp=index_lpext(ilpvalue)
            valuelp=value_lpext(ilpvalue)
            indexlp1=index_lpext1(ilpvalue)
            valuelp1=value_lpext1(ilpvalue)

            nn=nn0
            do i=1,itmp-1
               mm=mm+1
               nn=nn+1
!            if(indexlp.ne.0) then
              vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!            end if
            if(indexlp1.ne.0) then
              vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
            end if
            enddo
         enddo
         mm=iplplwei+iweista_g28-1
         nn=iplprwei
         nn=nn+1
         do itmp=2,nwei_g28
            nn=nn+1
            ilpvalue=0
            do i=1,itmp-1
              ilpvalue=ilpvalue+1
              mm=mm+1
              indexlp=index_lpext(ilpvalue)
!          if(indexlp.ne.0) then
              valuelp=-value_lpext(ilpvalue)
              vector2(indexlp)=vector2(indexlp)
     :                 +vector1(mm)*vector1(nn)*valuelp
!          end if
             indexlp1=index_lpext1(ilpvalue)
             if(indexlp1.ne.0) then
               valuelp1=-value_lpext1(ilpvalue)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
             end if
            enddo
         enddo
      elseif ( logic_g28a ) then
         mm=iplplwei+iweista_g28-1
         nn0=iplprwei
         do nn=nn0+1,nn0+nwei_g28

            ilpvalue=0
            do i=1,nint_g28
               ilpvalue=ilpvalue+1
               mm=mm+1
               indexlp=index_lpext(ilpvalue)
!          if(indexlp.ne.0) then
               valuelp=-value_lpext(ilpvalue)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!          end if
               indexlp1=index_lpext1(ilpvalue)
          if(indexlp1.ne.0) then
               valuelp1=-value_lpext1(ilpvalue)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
          end if
            enddo
         enddo
      endif
      end


      subroutine gdv_sequence_extspace_g(ilw,irw)
#include "drt_h.fh"
#include "lpextmode_h.fh"

      mm=ilw
      nn=irw+1

      do iij=1,ilsegdownwei
         mm=mm+1
         indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(iij)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(iij)
         if(indexlp1.ne.0) then
           valuelp1=value_lpext1(iij)
           vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
         end if
      enddo
      end

      subroutine complete_sd_ar_ext_loop_g(ilwei,irwei,isdownwei)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      common /gext_sequence/icnt_base,icano_nnsta,icano_nnend,
     *   m_jc,m_jd, isegsta,isegupwei,isegdownwei
c      write(*,*) 'sd_test 2/2','  td_test_2/2 2'
      ilpvalue=0
      mm0=ilwei
      nn=irwei+icano_nnsta-1
      do nntmp=icano_nnsta,icano_nnend
       nn=nn+1
       mm=mm0
       do mmtmp=1,isdownwei
          ilpvalue=ilpvalue+1
          mm=mm+1
          indexlp=index_lpext(ilpvalue)
          valuelp=value_lpext(ilpvalue)
          vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
          indexlp1=index_lpext1(ilpvalue)
          if(indexlp1.ne.0) then

            valuelp1=value_lpext1(ilpvalue)
            vector2(indexlp1)=vector2(indexlp1)
     :                    +vector1(mm)*vector1(nn)*valuelp1
          endif
        enddo
      enddo
c      print*, "out sd 2"
      return
      end

      subroutine gdv_sequence_extspace1_g(ilw,irw,n)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "grad_h.fh"
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)

      mm=ilw
      nn=irw+1

      do ilpvalue=1,ilsegdownwei
         mm=mm+1
         indexlp2=index_lpext5(ilpvalue)
         valuelp2=value_lpext5(ilpvalue)
         dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2

         do nii=1,n
            indexlp3=index_lpext3(ilpvalue,nii)
!          if(indexlp3.ne.0) then
            valuelp3=value_lpext3(ilpvalue,nii)
            vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!          end if

            indexlp4=index_lpext4(ilpvalue,nii)
          if(indexlp4.ne.0) then
            valuelp4=value_lpext4(ilpvalue,nii)
            vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
          end if
         enddo
      enddo
      end


      subroutine gtd_sequence_extspace1_g(iplplwei,iplprwei,n)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "grad_h.fh"
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)
      common /cont_tmp/icount_ext
      parameter (   v_sqtwo=1.414213562373095d0 )
c     write(*,*) ' td_test _1/2',' dt_test '
      ilpvalue=0
      if ( logic_g25a ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         do itmp=1,nint_g25
            ilpvalue=ilpvalue+1
            indexlp2=index_lpext5(ilpvalue)
            valuelp2=value_lpext5(ilpvalue)
            nn=nn0
            do i=1,nwei_g25
               mm=mm+1
               nn=nn+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                 indexlp3=index_lpext3(ilpvalue,nii)
!          if(indexlp3.ne.0) then
                 valuelp3=value_lpext3(ilpvalue,nii)
                 vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!          end if
                 indexlp4=index_lpext4(ilpvalue,nii)
          if(indexlp4.ne.0) then
                 valuelp4=value_lpext4(ilpvalue,nii)
                 vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
          end if
              enddo
           enddo
         enddo
      elseif ( logic_g25b ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         ilpvalue=ilpvalue+1
         do itmp=2,nint_g25
            ilpvalue=ilpvalue+1
            indexlp2=index_lpext5(ilpvalue)
            valuelp2=value_lpext5(ilpvalue)

            nn=nn0
            do i=1,itmp-1
               mm=mm+1
               nn=nn+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                 indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                 valuelp3=value_lpext3(ilpvalue,nii)
                 vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                 indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
                 valuelp4=value_lpext4(ilpvalue,nii)
                 vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
              enddo
            enddo
         enddo
         mm=iplplwei+iweista_g28-1
         nn=iplprwei
         nn=nn+1
         do itmp=2,nwei_g28
            nn=nn+1
            ilpvalue=0
            do i=1,itmp-1
               ilpvalue=ilpvalue+1
               indexlp2=index_lpext5(ilpvalue)
               valuelp2=-value_lpext5(ilpvalue)
               mm=mm+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                 indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                 valuelp3=-value_lpext3(ilpvalue,nii)
                 vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                 indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
                 valuelp4=-value_lpext4(ilpvalue,nii)
                 vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
              enddo
            enddo
         enddo
      elseif ( logic_g28a ) then
         mm=iplplwei+iweista_g28-1
         nn0=iplprwei
         do nn=nn0+1,nn0+nwei_g28

            ilpvalue=0
            do i=1,nint_g28
               ilpvalue=ilpvalue+1
               indexlp2=index_lpext5(ilpvalue)
               valuelp2=-value_lpext5(ilpvalue)
               mm=mm+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                 indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                 valuelp3=-value_lpext3(ilpvalue,nii)
                 vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                 indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
                 valuelp4=-value_lpext4(ilpvalue,nii)
                 vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
               enddo
            enddo
         enddo
      endif
      end

      subroutine gsd_sequence_extspace1_g(iplplwei,iplprwei,n)
#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "grad_h.fh"
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)
      parameter (v_sqtwo=1.414213562373095d0 )
c      write(*,*) '  sd_test 1/2','  ds_test 0'

      ilpvalue=0
      if ( logic_g25a ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         do itmp=1,nint_g25
            ilpvalue=ilpvalue+1
            indexlp2=index_lpext5(ilpvalue)
            valuelp2=value_lpext5(ilpvalue)

            nn=nn0
            do i=1,nwei_g25
               mm=mm+1
               nn=nn+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                 indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                 valuelp3=value_lpext3(ilpvalue,nii)
                 vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                 indexlp4=index_lpext4(ilpvalue,nii)
                 if(indexlp4.ne.0) then
                   valuelp4=value_lpext4(ilpvalue,nii)
                   vector2(indexlp4)=vector2(indexlp4)
     :                     +vector1(mm)*vector1(nn)*valuelp4
                 end if
                enddo
            enddo
         enddo
      elseif ( logic_g25b ) then
         mm=iplplwei+iweista_g25-1
         nn0=iplprwei
         ilpvalue=ilpvalue+1
         do itmp=2,nint_g25
            ilpvalue=ilpvalue+1
            indexlp2=index_lpext5(ilpvalue)
            valuelp2=value_lpext5(ilpvalue)

            nn=nn0
            do i=1,itmp-1
               mm=mm+1
               nn=nn+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                 indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                 valuelp3=value_lpext3(ilpvalue,nii)
                 vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                 indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
                 valuelp4=value_lpext4(ilpvalue,nii)
                 vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
                enddo

            enddo
         enddo

         mm=iplplwei+iweista_g28-1
         nn=iplprwei
         nn=nn+1
         do itmp=2,nwei_g28
            nn=nn+1
            ilpvalue=0
            do i=1,itmp-1
               ilpvalue=ilpvalue+1
               indexlp2=index_lpext5(ilpvalue)
               valuelp2=value_lpext5(ilpvalue)
               mm=mm+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                  indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                  valuelp3=value_lpext3(ilpvalue,nii)
                  vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                  indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
                  valuelp4=value_lpext4(ilpvalue,nii)
                  vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
               enddo
            enddo
         enddo
      elseif ( logic_g28a ) then
         mm=iplplwei+iweista_g28-1
         nn0=iplprwei
         do nn=nn0+1,nn0+nwei_g28
            ilpvalue=0
            do i=1,nint_g28
               ilpvalue=ilpvalue+1
               indexlp2=index_lpext5(ilpvalue)
               valuelp2=value_lpext5(ilpvalue)
               mm=mm+1
               dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
               do nii=1,n
                  indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
                  valuelp3=value_lpext3(ilpvalue,nii)
                  vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
                  indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
                  valuelp4=value_lpext4(ilpvalue,nii)
                  vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
               enddo
            enddo
         enddo
      endif

      if ( logic_g26 ) then
         ilpvalue=ivaluesta_g26
         mm=iplplwei+iweista_g26
         nn0=iplprwei
         do nn=nn0+1,nn0+nwei_g26
            ilpvalue=ilpvalue+1
            indexlp2=index_lpext5(ilpvalue)
            valuelp2=value_lpext5(ilpvalue)*v_sqtwo

            dm1tmp(indexlp2)=dm1tmp(indexlp2)
     :                  +vector1(mm)*vector1(nn)*valuelp2
            do nii=1,n
               indexlp3=index_lpext3(ilpvalue,nii)
!             if(indexlp3.ne.0) then
               valuelp3=value_lpext3(ilpvalue,nii)*v_sqtwo
               vector2(indexlp3)=vector2(indexlp3)
     :                   +vector1(mm)*vector1(nn)*valuelp3
!             end if
               indexlp4=index_lpext4(ilpvalue,nii)
             if(indexlp4.ne.0) then
               valuelp4=value_lpext4(ilpvalue,nii)*v_sqtwo
               vector2(indexlp4)=vector2(indexlp4)
     :                   +vector1(mm)*vector1(nn)*valuelp4
             end if
            enddo
            mm=mm+1
         enddo
      endif
c      print*, "out ds 0"

      end


      subroutine inn_ext_dd_loop_unpack_g(iplplwei,iplprwei)
#include "drt_h.fh"
#include "lpextmode_h.fh"

c      write(nf2,*) 'logic_g49b',logic_g50,logic_g49a,logic_g49b

      ii=1
      if ( logic_g50 ) then
         if ( logic_g49b) then
            mm=iplplwei
            nn=iplprwei
            do i=1,ildownwei_segdd
               mm=mm+1
               nn=nn+1
              indexlp=index_lpext(ii)
           if(indexlp.ne.0) then
              valuelp=value_lpext(ii)
              vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
           end if
              indexlp1=index_lpext1(ii)
           if(indexlp1.ne.0) then
              valuelp1=value_lpext1(ii)
              vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
           end if
               ii=ii+1
            enddo
         endif

         ii=ii+int_dd_drl
         mm0=iplplwei
         nn=iplprwei+1
         do icle=1,2
         do j=2,ildownwei_segdd
            nn=nn+1
            mm=mm0
            do i=1,j-1
               mm=mm+1
               indexlp=index_lpext(ii)
           if(indexlp.ne.0) then
               valuelp=value_lpext(ii)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
           end if
               indexlp1=index_lpext1(ii)
           if(indexlp1.ne.0) then
               valuelp1=value_lpext1(ii)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
           end if
               ii=ii+1
            enddo
         enddo
         if ( .not. logic_g49b ) exit
         mm0=iplprwei
         nn=iplplwei+1
         enddo

      else               !drl=.true.


          ii=ii+int_dd_drl
         if ( logic_g49a ) then
            mm0=iplplwei
            nn=iplprwei
             ildownwei=ildownwei_segdd
             irdownwei=irdownwei_segdd
         else
            mm0=iplprwei
            nn=iplplwei
             ildownwei=irdownwei_segdd
             irdownwei=ildownwei_segdd
         endif
         do j=1,irdownwei
            nn=nn+1
            mm=mm0
            do i=1,ildownwei
              mm=mm+1
               indexlp=index_lpext(ii)
!           if(indexlp.ne.0) then
               valuelp=value_lpext(ii)
               vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!           end if
               indexlp1=index_lpext1(ii)
!           if(indexlp1.ne.0) then
               valuelp1=value_lpext1(ii)
               vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
!           end if
               ii=ii+1
            enddo
         enddo
      endif
      end

      subroutine inn_ext_sv_loop_unpack_g(ilw,irw)
#include "drt_h.fh"
#include "lpextmode_h.fh"

      mm=ilw
      nn=irw+1
      do iij=1,ilsegdownwei
         mm=mm+1
         indexlp=index_lpext(iij)
!       if(indexlp.ne.0) then
         valuelp=value_lpext(iij)
         vector2(indexlp)=vector2(indexlp)
     :                   +vector1(mm)*vector1(nn)*valuelp
!       end if
         indexlp1=index_lpext1(iij)
       if(indexlp1.ne.0) then
         valuelp1=value_lpext1(iij)
         vector2(indexlp1)=vector2(indexlp1)
     :                   +vector1(mm)*vector1(nn)*valuelp1
       end if
      enddo
      end


!  density matrix formation : vector2=ci*<i|epq,rs|j>*cj and dm1=ci*<i|e
!  for norb_act<>0         mg1,mg2,mg3,mg4,mg5:
!  idb=1  in dbl_space      ity_up=0-5               0 ,jpad,iwdl,iwdr,
!  idb=2  in act_space      ity_up=0-5,itdown=0,3      jph, jpe,iwal,iwa
!  idb=3  betwin dbl and act   ity_up=0-5,itdown=0,3      jpe,iwdl,iwdr,

!  this subroutine prodab_1 does the dm1 part, which corresponds to voin
      subroutine prodab_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)
#include "drt_h.fh"
      if(log_prod.eq.1) call prodab_h_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,
     *                                  mg6,mg7)
      if(log_prod.eq.2) call prodab_h0_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,
     *                                  mg6,mg7)
      return
      end

      subroutine prodab_h_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
#include "grad_h.fh"
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)

      goto(100,200,300),idb
! in dbl_space
100   jpad=mg2
      iwdl=mg3
      iwdr=mg4
      do ipae=1,25
        if(nu_ae(ipae).eq.0) cycle
        iwdown=iw_downwei(jpad,ipae)
        if(iwdown.eq.0) cycle
        lwnu=iseg_downwei(ipae)
        do iwa=0,iwdown-1
          iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
          mm=iwadl
          nn=iwadr
          do m=1,lwnu
            mm=mm+1
            nn=nn+1
c101        vector2(mm)=vector2(mm)+vector1(nn)*wl
c          vector2(nn)=vector2(nn)+vector1(mm)*wl
c          vector2(mm)=vector2(mm)+vector1(nn)*wl
c          vector2(nn)=vector2(nn)+vector1(mm)*wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
          enddo
        enddo
      enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do iwe=1,isegdownwei
                mm=iwadl+iwe
                nn=iwadr+iwe
c201        vector2(mm)=vector2(mm)+vector1(nn)*wl
c           vector2(nn)=vector2(nn)+vector1(mm)*wl
c              vector2(mm)=vector2(mm)+vector1(nn)*wl
c              vector2(nn)=vector2(nn)+vector1(mm)*wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
c             write(nf2,'(4i4,4f18.10)')mg6,mg7, mm,nn,wl,
c     :       vector1(mm),vector1(nn),dm1tmp(mg67)
              enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      isegdownwei=iseg_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          mm=iwalk_ad(jpadl,ipael,iwal,iwdl)
          nn=iwalk_ad(jpad,ipae,iwar,iwdr)
          do iwe=1,isegdownwei
            mm=mm+1                  ! iwl=iwalk_ad
            nn=nn+1                  ! iwl=iwalk_ad
c301        vector2(mm)=vector2(mm)+vector1(nn)*wl
c           vector2(nn)=vector2(nn)+vector1(mm)*wl
c           vector2(mm)=vector2(mm)+vector1(nn)*wl
c            vector2(nn)=vector2(nn)+vector1(mm)*wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
c            write(nf2,'(a9,3i8,3f18.10)') '1_dbl_act',mg6,mm,nn,wl,
c     :             vector1(mm),vector1(nn)
            enddo
        enddo
      enddo
      goto 1000
1000  return
      end

      subroutine prodab_h0_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
#include "grad_h.fh"
      common /iaib/ ican_a(max_orb),ican_b(mtmp+max_orb)

      goto(100,200,300),idb
! in dbl_space
100   jpad=mg2
      iwdl=mg3
      iwdr=mg4
      do ipae=1,25
        if(nu_ae(ipae).eq.0) cycle
        iwdown=iw_downwei(jpad,ipae)
        if(iwdown.eq.0) cycle
        lwnu=iseg_downwei(ipae)
        do iwa=0,iwdown-1
          iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
          mm=iwadl
          nn=iwadr
          do m=1,lwnu
            mm=mm+1
            nn=nn+1
c             if(mm.gt.nn) mntmp=mm*(mm-1)/2+nn
c             if(nn.gt.mm) mntmp=nn*(nn-1)/2+mm
c             vector2(mntmp)=vector2(mntmp)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
!              if(mntmp.eq.2) then
!                write(*,*)'  102',vector2(mntmp),wl
!           endif
          enddo
        enddo
      enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do iwe=1,isegdownwei
                mm=iwadl+iwe
                nn=iwadr+iwe
c             if(mm.gt.nn) mntmp=mm*(mm-1)/2+nn
c             if(nn.gt.mm) mntmp=nn*(nn-1)/2+mm
c             vector2(mntmp)=vector2(mntmp)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
!              if(mntmp.eq.2) then
!                write(*,*)'  202',vector2(mntmp),wl
!           endif
              enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      isegdownwei=iseg_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          mm=iwalk_ad(jpadl,ipael,iwal,iwdl)
          nn=iwalk_ad(jpad,ipae,iwar,iwdr)
          do iwe=1,isegdownwei
            mm=mm+1                  ! iwl=iwalk_ad
            nn=nn+1                  ! iwl=iwalk_ad
c             if(mm.gt.nn) mntmp=mm*(mm-1)/2+nn
c             if(nn.gt.mm) mntmp=nn*(nn-1)/2+mm
c             vector2(mntmp)=vector2(mntmp)+wl
            mg67=ican_a(mg7)+mg6
            dm1tmp(mg67)=dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
!              if(mntmp.eq.2) then
!                write(*,*)'  302',vector2(mntmp),wl
!           endif
          enddo
        enddo
      enddo
      goto 1000
1000  return
      end


!this subroutine prodab_2 does the dm2 part, which corresponds to vint_c
      subroutine prodab_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)
#include "drt_h.fh"
      if(log_prod.eq.1) call prodab_h_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,
     *                                  mg6)
      if(log_prod.eq.2) call prodab_h0_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,
     *                                  mg6)
      return
      end

      subroutine prodab_h_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
#include "grad_h.fh"

      goto(100,200,300),idb
! in dbl_space
100   jpad=mg2
      iwdl=mg3
      iwdr=mg4
      do ipae=1,25
        if(nu_ae(ipae).eq.0) cycle
        iwdown=iw_downwei(jpad,ipae)
        if(iwdown.eq.0) cycle
        lwnu=iseg_downwei(ipae)
        do iwa=0,iwdown-1
          iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
          mm=iwadl
          nn=iwadr
          do m=1,lwnu
            mm=mm+1
            nn=nn+1
c101        vector2(mm)=vector2(mm)+vector1(nn)*wl
c          vector2(nn)=vector2(nn)+vector1(mm)*wl
c          vector2(mm)=vector2(mm)+vector1(nn)*wl
c          vector2(nn)=vector2(nn)+vector1(mm)*wl
            vector2(mg6)=vector2(mg6)+vector1(nn)*wl*vector1(mm)
c       if(mg6.eq.29)
c     :     write(nf2,'(i8,2i4,4f18.10)')mg6,mm,nn,vector2(mg6),
c     :                          vector1(mm),wl,vector1(nn)
           enddo
        enddo
      enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do iwe=1,isegdownwei
                mm=iwadl+iwe
                nn=iwadr+iwe
c201        vector2(mm)=vector2(mm)+vector1(nn)*wl
c           vector2(nn)=vector2(nn)+vector1(mm)*wl
c              vector2(mm)=vector2(mm)+vector1(nn)*wl
c              vector2(nn)=vector2(nn)+vector1(mm)*wl
            vector2(mg6)=vector2(mg6)+vector1(nn)*wl*vector1(mm)
c       if(mg6.eq.15)  write(nf2,'(2i4,4f18.10)')
c     :          mm,nn,wl,vector1(mm),vector1(nn),vector2(mg6)
c            write(nf2,'(a3,2i4,i8,2f18.10)') 'act',mm,nn,mg6,
c     :       vector1(mm),vector1(nn)
            enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      isegdownwei=iseg_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          mm=iwalk_ad(jpadl,ipael,iwal,iwdl)
          nn=iwalk_ad(jpad,ipae,iwar,iwdr)
          do iwe=1,isegdownwei
            mm=mm+1                  ! iwl=iwalk_ad
            nn=nn+1                  ! iwl=iwalk_ad
c301        vector2(mm)=vector2(mm)+vector1(nn)*wl
c           vector2(nn)=vector2(nn)+vector1(mm)*wl
c           vector2(mm)=vector2(mm)+vector1(nn)*wl
c            vector2(nn)=vector2(nn)+vector1(mm)*wl
            vector2(mg6)=vector2(mg6)+vector1(nn)*wl*vector1(mm)
          enddo
        enddo
      enddo
      goto 1000
1000  return
      end

      subroutine prodab_h0_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      dimension lopu(4,loputmp)
#include "grad_h.fh"

      goto(100,200,300),idb
! in dbl_space
100   jpad=mg2
      iwdl=mg3
      iwdr=mg4
      mntmp=0
      do ipae=1,25
        if(nu_ae(ipae).eq.0) cycle
        iwdown=iw_downwei(jpad,ipae)
        if(iwdown.eq.0) cycle
        lwnu=iseg_downwei(ipae)
        do iwa=0,iwdown-1
          iwadl=iwalk_ad(jpad,ipae,iwa,iwdl)
          iwadr=iwalk_ad(jpad,ipae,iwa,iwdr)
          mm=iwadl
          nn=iwadr
          do m=1,lwnu
            mm=mm+1
            nn=nn+1
              if(mm.gt.nn) then
                mntmp=mm*(mm-1)/2+nn
              else
                mntmp=nn*(nn-1)/2+mm
              endif
c             vector2(mntmp)=vector2(mntmp)+wl
            vector2(mg6)=vector2(mg6)+vector1(mntmp)*wl*vector1(mntmp)

!              if(mntmp.eq.2) then
!                write(*,*)'  102',vector2(mntmp),wl
!           endif
          enddo
        enddo
      enddo
      goto 1000
! in act_space
200   if(jpad.ne.jpadl) return
      jph=mg1
      jpl=mg2
c     iwal=mg3
c     iwar=mg4
      iwupwei=jpad_upwei(jpad)
      isegdownwei=iseg_downwei(ipae)
      jpy=jphy(jph)
      in=ihy(jpy)

      call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
      do lp=1,mp
        iwl=lopu(1,lp)-1
        iwr=lopu(2,lp)-1
        jpe=lopu(3,lp)
        lwnu=iy(1,jpe)
        do jwu=jpy+1,jpy+in
          iwal=iwl+ihy(jwu)
          iwar=iwr+ihy(jwu)
          do jwd=1,lwnu
            iwal=iwal+1
            iwar=iwar+1
            do iwd=0,iwupwei-1
              iwadl=iwalk_ad(jpadl,ipael,iwal,iwd)
              iwadr=iwalk_ad(jpad,ipae,iwar,iwd)
              do iwe=1,isegdownwei
                mm=iwadl+iwe
                nn=iwadr+iwe
              if(mm.gt.nn) then
                mntmp=mm*(mm-1)/2+nn
              else
                mntmp=nn*(nn-1)/2+mm
              endif
c             vector2(mntmp)=vector2(mntmp)+wl
            vector2(mg6)=vector2(mg6)+vector1(mntmp)*wl*vector1(mntmp)

!              if(mntmp.eq.2) then
!                write(*,*)'  202',vector2(mntmp),wl
!           endif
              enddo
            enddo
          enddo
        enddo
      enddo
      goto 1000
! betwin act and dbl
300   jpl =mg1
      iwdl=mg2
      iwdr=mg3
      isegdownwei=iseg_downwei(ipae)

      call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
      do lp=1,mp
        iwal=lopu(1,lp)-1
        iwar=lopu(2,lp)-1
        jpe= lopu(3,lp)
        jwnu=iy(1,jpe)
        do ii=1,jwnu
          iwal=iwal+1
          iwar=iwar+1
          mm=iwalk_ad(jpadl,ipael,iwal,iwdl)
          nn=iwalk_ad(jpad,ipae,iwar,iwdr)
          do iwe=1,isegdownwei
            mm=mm+1                  ! iwl=iwalk_ad
            nn=nn+1                  ! iwl=iwalk_ad
              if(mm.gt.nn) then
                mntmp=mm*(mm-1)/2+nn
              else
                mntmp=nn*(nn-1)/2+mm
              endif
c             vector2(mntmp)=vector2(mntmp)+wl
            vector2(mg6)=vector2(mg6)+vector1(mntmp)*wl*vector1(mntmp)

!              if(mntmp.eq.2) then
!                write(*,*)'  302',vector2(mntmp),wl
!           endif
          enddo
        enddo
      enddo
      goto 1000
1000  return
      end
