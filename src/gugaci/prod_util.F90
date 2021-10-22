!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2007, Bingbing Suo                                     *
!***********************************************************************

subroutine jl_ne_jr(mp,jl,jr,jwl,jwr,lopu)

use gugaci_global, only: iy, iyl, jj_sub, jjl_sub, loputmp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: mp, lopu(4,loputmp)
integer(kind=iwp), intent(in) :: jl, jr, jwl, jwr
integer(kind=iwp) :: i, idlr, jlp, jrp, jwlp, jwrp, lpi, lpj, ml, mr, nlp
integer(kind=iwp), allocatable :: lopi(:,:), lopj(:,:)

!-----------------------------------------------------------------------
! on entry:
!     jl  - left DRT node
!     jr  - right DRT node
!     jwl - left weight
!     jwr - right weight
! output
!     mp  - number of partial loop tails
!     lopu(1,i) - left weight of parital loop i
!     lopu(2,i) - right weight of partial loop i
!     lopu(3,i) - loop tail node
!-----------------------------------------------------------------------

! lopu(1,*)=jwl lopu(2,*)=jwr,lopu(3,*)=jl,lopu(4,*)=jr
!write(u6,*) 'in subroutine jl_ne_jr',jl,jr
call mma_allocate(lopi,4,loputmp,label='lopi')
call mma_allocate(lopj,4,loputmp,label='lopj')
mp = 0
lpi = 1
lopi(1,1) = jwl
lopi(2,1) = jwr
lopi(3,1) = jl
lopi(4,1) = jr
do
  lpj = 0
  do nlp=1,lpi
    if (lopi(3,nlp) == lopi(4,nlp)) then
      mp = mp+1
      lopu(1:4,mp) = lopi(1:4,nlp)
      cycle
    end if
    jlp = lopi(3,nlp)
    jrp = lopi(4,nlp)
    do idlr=1,4
      ml = jjl_sub(idlr,jlp)
      mr = jj_sub(idlr,jrp)
      !write(u6,*) 'ml,mr',ml,mr
      if ((ml == 0) .or. (mr == 0)) cycle
      jwlp = lopi(1,nlp)
      jwrp = lopi(2,nlp)
      if (idlr /= 1) jwlp = jwlp+iyl(idlr,jlp)
      if (idlr /= 1) jwrp = jwrp+iy(idlr,jrp)
      lpj = lpj+1
      lopj(1,lpj) = jwlp
      lopj(2,lpj) = jwrp
      lopj(3,lpj) = ml
      lopj(4,lpj) = mr
    end do
  end do
  if (lpj == 0) exit
  lpi = lpj
  do i=1,lpj
    lopi(1:4,i) = lopj(1:4,i)
  end do
end do
call mma_deallocate(lopi)
call mma_deallocate(lopj)

return

end subroutine jl_ne_jr

!  idb=1  in dbl_space          ity_up=0-5                 0,  jpad,iwdl,iwdr,
!  idb=2  in act_space          ity_up=0-5,itdown=0,3      jph,jpe, iwal,iwa
!  idb=3  between dbl and act   ity_up=0-5,itdown=0,3      jpe,iwdl,iwdr,

subroutine prodab_h0(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
!***********************************************************************
! calculate h matrix elements which were contributed
! only by whole inner space loops and store them into
! vector2

use gugaci_global, only: ihy, ipae, ipael, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jpadl, jphy, loputmp, vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, isegdownwei, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, jpe, &
                     jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mm, mntmp, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

!write(u6,*) 'prodab_02 '

select case (idb)
  case default ! (1)
    ! in dbl_space
    !jpad = mg2
    iwdl = mg3
    iwdr = mg4
    ipae = 1
    jpad = 1
    mntmp = 0
    iwdown = iw_downwei(jpad,ipae)
    lwnu = iseg_downwei(ipae)
    do iwa=0,iwdown-1
      iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
      iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
      mm = iwadl
      nn = iwadr
      do m=1,lwnu
        mm = mm+1
        nn = nn+1
        if (mm > nn) mntmp = mm*(mm-1)/2+nn
        if (nn > mm) mntmp = nn*(nn-1)/2+mm
        vector2(mntmp) = vector2(mntmp)+wl
        !if (mntmp == 2) write(u6,*) '  102',vector2(mntmp),wl
      end do
    end do
    !end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do iwe=1,isegdownwei
              mm = iwadl+iwe
              nn = iwadr+iwe
              if (mm > nn) then
                mntmp = mm*(mm-1)/2+nn
              else
                mntmp = nn*(nn-1)/2+mm
              end if
              vector2(mntmp) = vector2(mntmp)+wl
              if (mntmp == 7) then
                write(u6,*) '  202',vector2(mntmp),wl
              end if
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    isegdownwei = iseg_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        mm = iwalk_ad(jpadl,ipael,iwal,iwdl)
        nn = iwalk_ad(jpad,ipae,iwar,iwdr)
        do iwe=1,isegdownwei
          mm = mm+1             ! iwl=iwalk_ad
          nn = nn+1             ! iwl=iwalk_ad
          if (mm > nn) then
            mntmp = mm*(mm-1)/2+nn
          else
            mntmp = nn*(nn-1)/2+mm
          end if
          vector2(mntmp) = vector2(mntmp)+wl
          !if (mntmp == 2) write(u6,*) '  302',vector2(mntmp),wl
        end do
      end do
    end do
    call mma_deallocate(lopu)
end select

return

end subroutine prodab_h0

subroutine prodab_h(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
!***********************************************************************
! whole inner space loops - h*c
! 26 feb 2007 - revised by suo bing for multi-root calculation

use gugaci_global, only: ihy, indx, ipae, ipael, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jpadl, jphy, loputmp, mcroot, &
                         nu_ae, vector1, vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, ipae_, ipaeend, irot, irtidx, isegdownwei, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, &
                     iwe, iwl, iwr, iwupwei, jpe, jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mm, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

!write(u6,*) 'prodab_02 ',mcroot,indx(1),iw_downwei(jpad,ipae)
select case (idb)
  case default ! (1)
    ! in dbl_space
    jpad = mg2
    iwdl = mg3
    iwdr = mg4
    ipaeend = 25
    do ipae_=1,ipaeend
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdown = iw_downwei(jpad,ipae)
      if (iwdown == 0) cycle
      lwnu = iseg_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
        do irot=1,mcroot
          irtidx = indx(irot)
          mm = iwadl+irtidx
          nn = iwadr+irtidx
          do m=1,lwnu
            mm = mm+1
            nn = nn+1
            !if ((mm > nci_dim) .or. (nn > nci_dim)) write(u6,*) jpad,ipae,iw_downwei(jpad,ipae),iseg_downwei(ipae), &
            !                                                    jpad_upwei(jpad),iw_sta(jpad,ipae),iwadl,iwadr
            vector2(mm) = vector2(mm)+vector1(nn)*wl
            vector2(nn) = vector2(nn)+vector1(mm)*wl
          end do
        end do
      end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do irot=1,mcroot
              irtidx = indx(irot)
              mm = iwadl+irtidx
              nn = iwadr+irtidx
              do iwe=1,isegdownwei
                mm = mm+1
                nn = nn+1
                vector2(mm) = vector2(mm)+vector1(nn)*wl
                vector2(nn) = vector2(nn)+vector1(mm)*wl
              end do
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    isegdownwei = iseg_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)

    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        iwadl = iwalk_ad(jpadl,ipael,iwal,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwar,iwdr)
        do irot=1,mcroot
          irtidx = indx(irot)
          mm = iwadl+irtidx
          nn = iwadr+irtidx
          do iwe=1,isegdownwei
            mm = mm+1             ! iwl=iwalk_ad
            nn = nn+1             ! iwl=iwalk_ad
            vector2(mm) = vector2(mm)+vector1(nn)*wl
            vector2(nn) = vector2(nn)+vector1(mm)*wl
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)
end select

return

end subroutine prodab_h

subroutine prodab(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)

use gugaci_global, only: log_prod
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr
real(kind=wp), intent(in) :: wl

select case (log_prod)
  case (1)
    call prodab_h(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
  case (2)  ! save H0 CI matrix, use for mrci
    call prodab_h0(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
  case (3)
    call prodab_h0_t(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
  case (4)  ! H0C, use for mrpt2
    call prodab_h0_d(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)
  case default ! HC, use for mrci
end select

return

end subroutine prodab

subroutine prodab_h0_d(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)

use gugaci_global, only: ihy, indx, ipae, ipael, iw_downwei, iy, jpad, jpad_upwei, jpadl, jpae_downwei, jphy, loputmp, mcroot, &
                         vector1, vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, irot, irtidx, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, &
                     jpaedownwei, jpe, jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mm, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

! log_prod=2:directly no_formh0
!write(u6,*) 'prodab_h0 '

select case (idb)
  case default ! (1)
    ! in dbl_space
    iwdl = mg3
    iwdr = mg4
    iwdown = iw_downwei(jpad,ipae)
    lwnu = jpae_downwei(ipae)
    do iwa=0,iwdown-1
      iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
      iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
      do irot=1,mcroot
        irtidx = indx(irot)
        mm = iwadl+irtidx
        nn = iwadr+irtidx
        do m=1,lwnu
          mm = mm+1
          nn = nn+1
          !if ((mm > nci_dim) .or. (nn > nci_dim)) write(u6,*) jpad,ipae,iw_downwei(jpad,ipae),iseg_downwei(ipae), &
          !                                                    jpad_upwei(jpad),iw_sta(jpad,ipae),iwadl,iwadr
          vector2(mm) = vector2(mm)+vector1(nn)*wl
          vector2(nn) = vector2(nn)+vector1(mm)*wl
        end do
      end do
      !do m=1,lwnu
      !  mm = mm+1
      !  nn = nn+1
      !  vector2(mm) = vector2(mm)+vector1(nn)*wl
      !  vector2(nn) = vector2(nn)+vector1(mm)*wl
      !end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    jpaedownwei = jpae_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do irot=1,mcroot
              irtidx = indx(irot)
              mm = iwadl+irtidx
              nn = iwadr+irtidx
              do iwe=1,jpaedownwei
                mm = mm+1
                nn = nn+1
                vector2(mm) = vector2(mm)+vector1(nn)*wl
                vector2(nn) = vector2(nn)+vector1(mm)*wl
              end do
            end do
            !do iwe=1,jpaedownwei
            !  mm = iwadl+iwe
            !  nn = iwadr+iwe
            !  vector2(mm) = vector2(mm)+vector1(nn)*wl
            !  vector2(nn) = vector2(nn)+vector1(mm)*wl
            !end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    jpaedownwei = jpae_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        iwadl = iwalk_ad(jpadl,ipael,iwal,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwar,iwdr)
        do irot=1,mcroot
          irtidx = indx(irot)
          mm = iwadl+irtidx
          nn = iwadr+irtidx
          do iwe=1,jpaedownwei
            mm = mm+1             ! iwl=iwalk_ad
            nn = nn+1             ! iwl=iwalk_ad
            vector2(mm) = vector2(mm)+vector1(nn)*wl
            vector2(nn) = vector2(nn)+vector1(mm)*wl
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)
end select

return

end subroutine prodab_h0_d

subroutine prodab_h0_t(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr)

use gugaci_global, only: ihy, ipae, ipael, iw_downwei, iy, jpad, jpad_upwei, jpadl, jpae_downwei, jphy, loputmp, vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, jpaedownwei, jpe, &
                     jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mm, mnh0, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

! log_prod=1:traditional formh0
!write(u6,*) 'prodab_h0 '

select case (idb)
  case default ! (1)
    ! in dbl_space
    iwdl = mg3
    iwdr = mg4
    iwdown = iw_downwei(jpad,ipae)
    lwnu = jpae_downwei(ipae)
    do iwa=0,iwdown-1
      mm = iwalk_ad(jpad,ipae,iwa,iwdl)
      nn = iwalk_ad(jpad,ipae,iwa,iwdr)
      do m=1,lwnu
        mm = mm+1
        nn = nn+1
        mnh0 = mm*(mm-1)/2+nn
        if (nn > mm) mnh0 = nn*(nn-1)/2+mm
        vector2(mnh0) = vector2(mnh0)+wl
        !if (((mm == 1) .and. (nn == 9)) .or. ((mm == 9) .and. (nn == 1))) write(nf2,'(a10,2i5,2f18.8)') ' in dbl ',mm,nn,wl, &
        !                                                                                                vector2(mnh
      end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    jpaedownwei = jpae_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do iwe=1,jpaedownwei
              mm = iwadl+iwe
              nn = iwadr+iwe
              mnh0 = mm*(mm-1)/2+nn
              if (nn > mm) mnh0 = nn*(nn-1)/2+mm
              vector2(mnh0) = vector2(mnh0)+wl
              !if (((mm == 1) .and. (nn == 9)) .or. ((mm == 9) .and. (nn == 1))) write(nf2,'(a10,2i5,2f18.8)') ' in act ',mm,nn,wl, &
              !                                                                                                vector2(mnh
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    jpaedownwei = jpae_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        mm = iwalk_ad(jpadl,ipael,iwal,iwdl)
        nn = iwalk_ad(jpad,ipae,iwar,iwdr)
        do iwe=1,jpaedownwei
          mm = mm+1             ! iwl=iwalk_ad
          nn = nn+1             ! iwl=iwalk_ad
          mnh0 = mm*(mm-1)/2+nn
          if (nn > mm) mnh0 = nn*(nn-1)/2+mm
          vector2(mnh0) = vector2(mnh0)+wl
          !if (((mm == 1) .and. (nn == 9)) .or. ((mm == 9) .and. (nn == 1))) write(nf2,'(a10,2i5,2f18.8)') ' act-dbl ',mm,nn,wl, &
          !                                                                                                vector2(mnh0)
        end do
      end do
    end do
    call mma_deallocate(lopu)
end select

return

end subroutine prodab_h0_t
