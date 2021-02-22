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
!     juv,just(nost,nost),jud(nost)
!     |  \  1         |
!     | d,dd,s(i=i)   |
!     |    \ s,t,tt(i<j)|
!     |     \       1 2 |     deal with inner of dbl_space
!     |ss(i>j)\       |
!     |  2 1  \       |
      subroutine gugadrt_dbl_downwalk()
#include "gendrt.fh"
!     integer lsml(10,10)       !to del
      if(norb_dbl.ne.0) goto 200
!----------- norb_dbl=0 ------------------------------------------------
      do im=1,ng_sm
        nnd=iseg_sta(1+im)
        nnt=iseg_sta(9+im)
        nns=iseg_sta(17+im)
        do lri=norb_dz,norb_frz+1,-1
          ismi=lsm_inn(lri)
         if(ismi.ne.im) cycle
          jud(lri)=nnd
          nnd=nnd+iseg_downwei(1+im)
        enddo
        do lrj=norb_dz,norb_frz+1,-1
         ismj=lsm_inn(lrj)
          do lri=lrj,1,-1
           ismi=lsm_inn(lri)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lri,lrj)=nns
            nns=nns+iseg_downwei(17+im)
            if(lri.eq.lrj) cycle
           just(lrj,lri)=nnt
            nnt=nnt+iseg_downwei(9+im)
           enddo
        enddo
      enddo
!----------- norb_dbl=0 ------------------------------------------------
!----------- norb_dbl<>0 -----------------------------------------------
 200  continue
      do im=1,ng_sm
        nnd=0
        nns=0
        do lri=norb_frz+1,norb_dz
          ismi=mul_tab(lsm_inn(lri),ns_sm)
         if(ismi.ne.im) cycle
          jud(lri)=nnd
          nnd=nnd+1
        enddo
        do lri=norb_frz+1,norb_dz-1
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          do lrj=lri+1,norb_dz      !tmp
           ismj=lsm_inn(lrj)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lri,lrj)=nns
            nns=nns+1
         enddo
        enddo
        if(im.eq.ns_sm) then
          do lr0=norb_frz+1,norb_dz
            just(lr0,lr0)=nns
            nns=nns+1
          enddo
        endif
        do lri=norb_frz+1,norb_dz-1
          ismi=mul_tab(lsm_inn(lri),ns_sm)
          do lrj=lri+1,norb_dz      !tmp
           ismj=lsm_inn(lrj)
           ismij=mul_tab(ismi,ismj)
            if(ismij.ne.im) cycle
            just(lrj,lri)=nns
            nns=nns+1
         enddo
        enddo
      enddo

      return
      end
