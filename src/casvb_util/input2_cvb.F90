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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine input2_cvb(iorbrel,mxdimrel,ifxorb,iorts,irots,izeta,orbs,irdorbs)

use casvb_global, only: i2s_fr, inputmode, mnion_fr, mxion_fr, nalf_fr, nbet_fr, nconf_fr, nconfion_fr, ndetvb_fr, ndetvb2_fr, &
                        nel_fr, nfrag, nMs_fr, nS_fr, nvbr_fr
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: mxdimrel, iorbrel(mxdimrel), ifxorb(mxorb_cvb), iorts(2,*), irots(2,*), izeta(*), irdorbs(mxorb_cvb)
real(kind=wp) :: orbs(mxaobf,mxorb_cvb)
#include "files_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i2s_min, ibase, iconf, iconf_add, idelstr, ifrag, ifrom, ifsc, ifxstr, ioffs, iorb, ip_cvb, ip_from, &
                     ip_iconfs, ip_symelm, ip_to, iS, isyme, ito, jconf, jorb, kbasiscvb_inp, need, nelcheck, nmov, noe1
real(kind=wp) :: swap
integer(kind=iwp), external :: ihlf_cvb, mheapiz_cvb, mheaprz_cvb, mstacki_cvb, nvb_cvb

ibase = mstacki_cvb(0)
ip_iconfs = mheapiz_cvb(0)
ip_cvb = mheaprz_cvb(0)
ip_symelm = mheaprz_cvb(0)
ifxstr = mheapiz_cvb(0)
idelstr = mheapiz_cvb(0)
noe = 2*mxorb_cvb

call hini_cvb()

call defs_cvb()
call casinfodef_cvb()

! Counters
nconf = 0
nvbinp = 0
nsyme = 0
norbrel = 0
ndimrel = 0
nijrel = 0
nfxorb = 0
nfxvb = 0
nzrvb = 0
nort = 0
ndrot = 0
lfxvb = 0
lzrvb = 0
call maxdims0_cvb()
call izero(ifxorb,mxorb_cvb)
call izero(izeta,mxsyme)
call fraginit_cvb()

call input3_cvb(iorbrel,mxdimrel,ifxorb,ifxstr,idelstr,iorts,irots,izeta,ip_iconfs,orbs,irdorbs,ip_cvb,ip_symelm,kbasiscvb_inp)

if (inputmode == 2) then
  ! Input parsing complete for this step ...
  ! ... Work out NORB, NEL, S ...
  noe1 = noe
  call casinfoset_cvb()
  ! ... Do ICONFS before others to get NVB and related info ...
  do iconf=1,nconf
    call imove_cvb(iwork((iconf-1)*noe1+ip_iconfs),iwork((iconf-1)*noe+ip_iconfs),noe)
  end do
  call mrealloci_cvb(ip_iconfs,noe*nconf)

  if (nfrag <= 1) then
    nMs_fr(1) = 1
    nalf_fr(1,1) = nalf
    nbet_fr(1,1) = nbet
  else
    do ifrag=1,nfrag
      nMs_fr(ifrag) = 1
      nalf_fr(1,ifrag) = (nel_fr(ifrag)+i2s_fr(1,ifrag))/2
      nbet_fr(1,ifrag) = nel_fr(ifrag)-nalf_fr(1,ifrag)
    end do
  end if
  if (nfrag == 0) then
    nfrag = 1
    nel_fr(1) = nel
    nconf_fr(1) = nconf
    nS_fr(1) = 1
    i2s_fr(1,1) = nalf-nbet
  end if
  do ifrag=1,nfrag
    if (nS_fr(ifrag) == 0) then
      nS_fr(ifrag) = 1
      i2s_fr(1,ifrag) = nalf-nbet
    end if
  end do
  iconf_add = 0
  do ifrag=1,nfrag
    if (nel_fr(ifrag) == 0) then
      nel_fr(ifrag) = nel
      nalf_fr(1,ifrag) = nalf
      nbet_fr(1,ifrag) = nbet
    end if
    if (nS_fr(ifrag) == 0) then
      nS_fr(ifrag) = 1
      i2s_fr(1,ifrag) = nalf-nbet
    end if
    if (nconf_fr(ifrag) == 0) then
      nconf = nconf+1
      nconf_fr(ifrag) = 1
      call mrealloci_cvb(ip_iconfs,noe*nconf)
      do jconf=nconf,iconf_add+2,-1
        call imove_cvb(iwork((jconf-2)*noe+ip_iconfs),iwork((jconf-1)*noe+ip_iconfs),noe)
      end do
      call izero(iwork(iconf_add*noe+ip_iconfs),noe)
      do i=1,min(nel_fr(ifrag),norb)
        iwork(i+iconf_add*noe+ip_iconfs-1) = 1
      end do
      do i=1,nel_fr(ifrag)-norb
        iwork(i+iconf_add*noe+ip_iconfs-1) = 2
      end do
    end if
    call cnfcheck_cvb(iwork(iconf_add*noe+ip_iconfs),nconf_fr(ifrag),nel_fr(ifrag))
    call cnfini_cvb(iwork(iconf_add*noe+ip_iconfs),nconf_fr(ifrag),nel_fr(ifrag),nS_fr(ifrag),i2s_fr(1,ifrag),nMs_fr(ifrag), &
                    nalf_fr(1,ifrag),nbet_fr(1,ifrag),nvbr_fr(ifrag),ndetvb_fr(ifrag),ndetvb2_fr(ifrag),mnion_fr(ifrag), &
                    mxion_fr(ifrag),nconfion_fr(0,ifrag),ifsc)
    iconf_add = iconf_add+nconf_fr(ifrag)
  end do
  ndetvb = 0
  ndetvb2 = 0
  nvbr = 0
  nelcheck = 0
  do i=1,nfrag
    ndetvb = ndetvb+ndetvb_fr(i)
    ndetvb2 = ndetvb2+ndetvb2_fr(i)
    nvbr = nvbr+nvbr_fr(i)
    nelcheck = nelcheck+nel_fr(i)
  end do
  if (nelcheck /= nel) then
    write(u6,*) ' Error: total number of electrons in fragment wavefunctions :',nelcheck,' not equal to number of electrons ',nel
    call abend_cvb()
  end if
  sc = (nfrag == 1) .and. (ifsc == 1)
  ! Set absym and use just lowest spin value if spinbas=determinants:
  absym(1) = (nalf == nbet)
  do ifrag=1,nfrag
    i2s_min = nel_fr(ifrag)
    do iS=1,nS_fr(ifrag)
      if (i2s_fr(iS,ifrag) /= 0) absym(1) = .false.
    end do
    if (kbasis == 6) then
      nS_fr(ifrag) = 1
      i2s_fr(1,ifrag) = i2s_min
    end if
  end do
  do i=2,5
    absym(i) = absym(1)
  end do
  nvb = nvb_cvb(kbasis)
  mnion = mnion_fr(1)
  mxion = mxion_fr(1)
  do i=2,nfrag
    mnion = min(mnion,mnion_fr(i))
    mxion = max(mxion,mxion_fr(i))
  end do
  ! ... Now remaining quantities that depend on NORB or NVB ...
  ! SYMELM
  ip_from = ip_symelm
  ip_to = ip_symelm
  do isyme=1,nsyme
    do iorb=1,norb
      if (ip_from /= ip_to) call fmove_cvb(work(ip_from),work(ip_to),norb)
      ip_from = ip_from+mxorb_cvb
      ip_to = ip_to+norb
    end do
    ip_from = ip_from+(mxorb_cvb-norb)*mxorb_cvb
  end do
  ! IORBREL
  ifrom = 1
  ito = 1
  do while (ifrom <= ndimrel)
    iorb = iorbrel(ifrom)
    jorb = iorbrel(ifrom+1)
    nmov = 3+iorbrel(ifrom+2)
    if ((iorb <= norb) .and. (jorb <= norb)) then
      if (ifrom /= ito) call imove_cvb(iorbrel(ifrom),iorbrel(ito),nmov)
      ito = ito+nmov
    end if
    ifrom = ifrom+nmov
  end do
  ndimrel = ito-1
  ! IFXSTR
  ito = 0
  do ifrom=1,nfxvb
    if (iwork(ifrom+ifxstr-1) <= nvb) then
      ito = ito+1
      iwork(ito+ifxstr-1) = iwork(ifrom+ifxstr-1)
    end if
  end do
  nfxvb = ito
  ! IDELSTR
  ito = 0
  do ifrom=1,nzrvb
    if (iwork(ifrom+idelstr-1) <= nvb) then
      ito = ito+1
      iwork(ito+idelstr-1) = iwork(ifrom+idelstr-1)
    end if
  end do
  nzrvb = ito
  ! IORTS
  ito = 0
  do ifrom=1,nort
    if ((iorts(1,ifrom) <= norb) .and. (iorts(2,ifrom) <= norb)) then
      ito = ito+1
      iorts(1,ito) = iorts(1,ifrom)
      iorts(2,ito) = iorts(2,ifrom)
    end if
  end do
  nort = ito
  ! IROTS
  ito = 0
  do ifrom=1,ndrot
    if ((irots(1,ifrom) <= norb) .and. (irots(2,ifrom) <= norb)) then
      ito = ito+1
      irots(1,ito) = irots(1,ifrom)
      irots(2,ito) = irots(2,ifrom)
    end if
  end do
  ndrot = ito
  ! Calling DEFS2 before INITOPT is required in order to set things such as ICRIT:
  call defs2_cvb(ifxorb)
  call initopt_cvb(icrit,lfxvb,nfxvb,iorts,nort,norb)

  call defs2_cvb(ifxorb)

  ! Try for new record
  call rdioff1_cvb(need)
  need = need+3*ihlf_cvb(1)+ihlf_cvb(noe*nconf)+mxaobf*norb+ihlf_cvb(norb)+nvbinp+nsyme*norb*norb+ihlf_cvb(ndimrel)+ &
         ihlf_cvb(norb)+ihlf_cvb(nfxvb)+ihlf_cvb(nzrvb)+ihlf_cvb(2*nort)+ihlf_cvb(2*ndrot)+ihlf_cvb(2*ndrot)+ihlf_cvb(nsyme)
  if (recinp == Zero) then
    recinp = recn_tmp01
  else if (recinp_old == Zero) then
    recinp_old = recn_tmp01
    recinp = recn_tmp02
  else
    swap = recinp_old
    recinp_old = recinp
    recinp = swap
  end if
  call reserv_cvb(need,recinp)
  call rdioff1_cvb(ioffs)
  call wrioff_cvb(1,recinp,ioffs)
  call wris_cvb([noe],1,recinp,ioffs)
  call wrioff_cvb(2,recinp,ioffs)
  call wris_cvb([nconf],1,recinp,ioffs)
  call wrioff_cvb(3,recinp,ioffs)
  call wris_cvb([kbasiscvb_inp],1,recinp,ioffs)
  call wrioff_cvb(4,recinp,ioffs)
  call wris_cvb(iwork(ip_iconfs),noe*nconf,recinp,ioffs)
  call wrioff_cvb(5,recinp,ioffs)
  call wrrs_cvb(orbs,mxaobf*norb,recinp,ioffs)
  call wrioff_cvb(6,recinp,ioffs)
  call wris_cvb(irdorbs,norb,recinp,ioffs)
  call wrioff_cvb(7,recinp,ioffs)
  call wrrs_cvb(work(ip_cvb),nvbinp,recinp,ioffs)
  call wrioff_cvb(8,recinp,ioffs)
  call wrrs_cvb(work(ip_symelm),nsyme*norb*norb,recinp,ioffs)
  call wrioff_cvb(9,recinp,ioffs)

  call dset_cvb(iorbrel,ifxorb,ifxstr,idelstr,iorts,irots,izeta)
else
  call maxdims_cvb()
end if
call mhpfreei_cvb(ip_iconfs)
call mhpfreer_cvb(ip_cvb)
call mhpfreer_cvb(ip_symelm)
call mhpfreei_cvb(ifxstr)
call mhpfreei_cvb(idelstr)
call hend_cvb()
call mfreei_cvb(ibase)

return

end subroutine input2_cvb
