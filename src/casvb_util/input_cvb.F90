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

subroutine input_cvb()

use Index_Functions, only: nTri_Elem
use casvb_global, only: absym, confsinp, gsinp, i2s_fr, iciweights, icrit, imethod, initial, inputmode, iorbprm, ipr, iprec, &
                        isaddle, ishstruc, isympr, isymv, ityp, ivbweights, iwidth, kbasis, lfxvb, lzrvb, mnion, mnion_fr, mxaobf, &
                        mxion, mxion_fr, mxirrep, mxiter, mxorb_cvb, mxops, mxsyme, nalf, nalf_fr, nbet, nbet_fr, nconf, nconf_fr, &
                        nconfion_fr, ndetvb, ndetvb_fr, ndetvb2_fr, ndimrel, ndrot, nel, nel_fr, nfrag, nfxorb, nfxvb, nijrel, &
                        nMs_fr, noe, norb, norbrel, nort, npcf, nS_fr, nspinb, nsyme, nvb, nvbinp, nvbr_fr, nzrvb, ploc, projcas, &
                        projsym, recinp, recinp_old, recn_tmp01, recn_tmp02, savvb, savvbci, sc, service, sij, spinbkw, strtci, &
                        strtint, strtmo, strtvb, tags
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoI

implicit none
integer(kind=iwp) :: i, i2s_min, iconf_add, idum(1), ifrag, ifrom, ifsc, igroup, io, ioffs, iorb, ip_from, ip_to, iS, istr, istr2, &
                     isyme, isymput, itag, ito, jorb, kbasis_old, kbasiscvb_inp, ltmp, mxalter, mxdimrel, mxortl, mxpair, mxread, &
                     nelcheck, nfrorb1, nfxorb1, nmov, nops, nread
real(kind=wp) :: swap
logical(kind=iwp) :: DoCholesky
character(len=50) :: inpstr
integer(kind=iwp), allocatable :: idelstr(:), ifxorb(:), ifxstr(:), iorbrel(:), iorts(:,:), irdorbs(:), irots(:,:), izeta(:), &
                                  itmp(:), itmp2(:), itmp3(:,:)
real(kind=wp), allocatable :: orbs(:,:), symelm(:), tmp(:)
integer(kind=iwp), parameter :: nbuf = 50, ncmp = 4, ncrit = 2, nendvb = 3, nglob = 5, nmeth = 12, nspec = 3, nstrin = 51, nwkw = 5
character(len=*), parameter :: crit(ncrit) = ['OVERLAP ','ENERGY  '], &
                               endvb(nendvb) = ['ENDVB   ','ENDCASVB','END     '], &
                               glbl(nglob) = ['XXXXxxxx','START   ','GUESS   ','PRINT   ','PREC    '], &
                               methkw(nmeth) = ['FLETCHER','TRIM    ','TRUSTOPT','DAVIDSON','STEEP   ','VB2CAS  ','AUGHESS ', &
                                                'AUG2    ','CHECK   ','DFLETCH ','NONE    ','SUPER   '], &
                               specl(nspec) = ['SERVICE ','MOSCOW  ','PERFLOC '], &
                               string(nstrin) = ['XXXXxxxx','XXXXxxxx','SAVE    ','XXXXxxxx','ORBPERM ','COUPLE  ','MAXITER ', &
                                                 'CRIT    ','CASPROJ ','PROJCAS ','NOCASPRO','NOPROJCA','XXXXxxxx','XXXXxxxx', &
                                                 'SYMELM  ','ORBREL  ','XXXXxxxx','SYMPROJ ','NOSYMPRO','FIXORB  ','FIXSTRUC', &
                                                 'DELSTRUC','FREORB  ','FRESTRUC','ORTHCON ','SADDLE  ','SHSTRUC ','VBWEIGHT', &
                                                 'CIWEIGHT','SCORR   ','NOSCORR ','METHOD  ','OPTIM   ','OPT     ','ENDOPTIM', &
                                                 'REPORT  ','ENDREPOR','XXXXxxxx','TUNE    ','XXXXxxxx','OPPOSITE','XXXXxxxx', &
                                                 'XXXXxxxx','STAT    ','INIT    ','NOINIT  ','TIDY    ','PLOC    ','NOPLOC  ', &
                                                 'ALTERNAT','ENDALTER'], &
                               weightkw(nwkw) = ['CHIRGWIN','LOWDIN  ','INVERSE ','NONE    ','ALL     ']
integer(kind=iwp), external :: nvb_cvb
logical(kind=iwp), external :: firsttime_cvb, valid_cvb ! ... Files/Hamiltonian available ...

call mma_allocate(ifxorb,mxorb_cvb,label='ifxorb')
call mma_allocate(izeta,mxsyme,label='izeta')

noe = 2*mxorb_cvb
mxpair = nTri_Elem(mxorb_cvb)
mxdimrel = mxpair*(3+mxops)

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
ifxorb(:) = 0
izeta(:) = 0
nfrag = 0

call mma_allocate(iorbrel,mxdimrel,label='iorbrel')
call mma_allocate(iorts,2,mxpair,label='iorts')
call mma_allocate(irots,2,mxpair,label='irots')
call mma_allocate(irdorbs,mxorb_cvb,label='irdorbs')
call mma_allocate(orbs,mxaobf,mxorb_cvb,label='orbs')
irdorbs(:) = 0
orbs(:,:) = Zero

call DecideOnCholesky(DoCholesky)
if (DoCholesky) then
  write(u6,*) '** Cholesky or RI/DF not yet implemented in CASVB **'
  call abend_cvb()
end if

call fstring_cvb(specl,nspec,istr,ncmp,2)
if (istr == 1) then
  ! 'SERVICE'
  service = .true.
  write(u6,'(1x,a,/)') '**** Service mode **** '
  !call service_cvb()
  write(u6,*) ' Casvb dummy routine called : SERV'
else if (istr == 2) then
  ! 'MOSCOW'
  service = .true.
  write(u6,'(1x,a,/)') '**** MOSCOW mode **** '
  !call moscow_cvb()
  write(u6,*) ' Casvb dummy routine called : MOSCOW'
else if (istr == 3) then
  service = .true.
  write(u6,'(1x,a,/)') '**** PERFLOC mode **** '
  !call perfloc_plc(3)
  write(u6,*) ' Molint dummy routine called : perfloc_plc'
else

  do

    ! CASSCF wavefunction information:
    call casinfoinp_cvb()
    ! VB wavefunction information:
    call fraginp_cvb()

    igroup = 0
    call fstring_cvb(glbl,nglob,istr,ncmp,2)
    if (istr /= 0) then
      igroup = 1
    else
      call fstring_cvb(string,nstrin,istr,ncmp,2)
      if (istr /= 0) then
        igroup = 2
      else
        call fstring_cvb(endvb,nendvb,istr,ncmp,2)
        if (istr /= 0) igroup = 3
      end if
    end if
    if (igroup == 3) then
      ! 'ENDVB', 'ENDCASVB' or 'END'
      istr = 0
    end if

    if (igroup /= 2) then
      if (istr == 2) then
        ! 'START'
        strtvb = Zero
        do
          call string_cvb(inpstr,1,nread,1)
          if (nread /= 1) exit
          if (inpstr(1:3) == 'CI=') then
            call setfn_cvb(strtci,inpstr(4:50))
          else if (inpstr(1:3) == 'VB=') then
            call setfn_cvb(strtvb,inpstr(4:50))
          else if (inpstr(1:3) == 'MO=') then
            call setfn_cvb(strtmo,inpstr(4:50))
          else if (inpstr(1:4) == 'INT=') then
            call setfn_cvb(strtint,inpstr(5:50))
          else
            exit
          end if
        end do
        if (valid_cvb(strtvb) .and. firsttime_cvb()) call touch_cvb('STRTGS')
      else if (istr == 3) then
        ! 'GUESS'
        call gsinp_cvb(orbs,irdorbs,nvbinp,kbasiscvb_inp,mxaobf,mxorb_cvb,kbasis,strtvb)
      else if (istr == 4) then
        ! 'PRINT'
        call int_cvb(ipr,10,nread,1)
      else if (istr == 5) then
        ! 'PREC'
        call int_cvb(idum,1,nread,1)
        iprec = idum(1)
        if (iprec < 0) then
          write(u6,*) ' Illegal precision :',iprec
          call abend_cvb()
        end if
        call int_cvb(idum,1,nread,1)
        iwidth = idum(1)
        call formats_cvb()
      end if

      ! 'ENDVB', 'ENDCASVB', 'END' or unrecognized keyword -- end of input:
      if (istr /= 0) cycle
    end if

    if (istr == 1) then
    else if (istr == 2) then
    else if (istr == 3) then
      ! 'SAVE'
      do
        call string_cvb(inpstr,1,nread,1)
        if (nread /= 1) exit
        if (inpstr(1:5) == 'VBCI=') then
          call setfn_cvb(savvbci,inpstr(6:50))
        else if (inpstr(1:3) == 'VB=') then
          call setfn_cvb(savvb,inpstr(4:50))
        else
          exit
        end if
      end do
    else if (istr == 4) then
    else if (istr == 5) then
      ! 'ORBPERM'
      if (firsttime_cvb()) call touch_cvb('ORBPERM')
      call int_cvb(iorbprm,mxorb_cvb,nread,0)
      if (nread > mxorb_cvb) then
        write(u6,*) ' Too many orbitals in ORBPERM keyword!'
        call abend_cvb()
      end if
      do iorb=1,nread
        if ((abs(iorbprm(iorb)) < 1) .or. (abs(iorbprm(iorb)) > mxorb_cvb)) then
          write(u6,'(a,40i3)') ' Illegal orbital label(s) in ORBPERM:',(iorbprm(io),io=1,nread)
          call abend_cvb()
        end if
      end do
    else if (istr == 6) then
      ! 'COUPLE'
      kbasis_old = kbasis
      call fstring_cvb(spinbkw,nspinb,kbasis,ncmp,1)
      if (kbasis == 0) kbasis = kbasis_old
      if (kbasis == 7) kbasis = 6
    else if (istr == 7) then
      ! 'MAXITER'
      call int_cvb(idum,1,nread,0)
      mxiter = idum(1)
    else if (istr == 8) then
      ! 'CRIT'
      call fstring_cvb(crit,ncrit,icrit,ncmp,1)
      if ((icrit /= 1) .and. (icrit /= 2)) then
        write(u6,*) ' Unrecognized CRIT keyword!'
        call abend_cvb()
      end if
    else if ((istr == 9) .or. (istr == 10)) then
      ! 'CASPROJ' or 'PROJCAS'
      projcas = .true.
    else if ((istr == 11) .or. (istr == 12)) then
      ! 'NOCASPROJ' or 'NOPROJCAS'
      projcas = .false.
    else if (istr == 15) then
      ! 'SYMELM'
      nsyme = nsyme+1
      if (nsyme > mxsyme) then
        write(u6,*) ' Too many symmetry elements found :',nsyme,mxsyme
        call abend_cvb()
      end if
      call mma_allocate(tmp,mxorb_cvb*mxorb_cvb*nsyme,label='symelm')
      if (allocated(symelm)) then
        tmp(1:mxorb_cvb*mxorb_cvb*(nsyme-1)) = symelm(:)
        if (allocated(symelm)) call mma_deallocate(symelm)
      end if
      call move_alloc(tmp,symelm)
      call symelminp_cvb(symelm,nsyme,tags,izeta,mxirrep,mxorb_cvb,mxsyme,ityp)
    else if (istr == 16) then
      ! 'ORBREL'
      iorb = 0
      jorb = 0
      call int_cvb(idum,1,nread,1)
      iorb = idum(1)
      call int_cvb(idum,1,nread,1)
      jorb = idum(1)
      if ((iorb < 1) .or. (iorb > mxorb_cvb) .or. (jorb < 1) .or. (jorb > mxorb_cvb)) then
        write(u6,*) ' Illegal orbital number(s) in ORBREL:',iorb,jorb
        call abend_cvb()
      end if
      iorbrel(1+ndimrel) = iorb
      iorbrel(2+ndimrel) = jorb
      nops = 0
      do
        call fstring_cvb(tags,nsyme,itag,3,1)
        if (itag == 0) exit
        nops = nops+1
        if (ndimrel+3+nops > mxdimrel) then
          write(u6,*) ' Too many symmetry elements in ORBREL keyword!'
          call abend_cvb()
        end if
        iorbrel(nops+3+ndimrel) = itag
      end do
      iorbrel(3+ndimrel) = nops
      norbrel = norbrel+1
      ndimrel = ndimrel+3+nops
    else if (istr == 18) then
      ! 'SYMPROJ'
      projsym = .true.
      isympr(:) = 0
      call int_cvb(idum,1,nread,1)
      isymput = idum(1)
      if (nread == 1) then
        isympr(isymput) = 1
        do
          call int_cvb(idum,1,nread,1)
          isymput = idum(1)
          if (nread /= 1) exit
          isympr(isymput) = 1
        end do
      else
        isympr(:) = isymv(:)
      end if
    else if (istr == 19) then
      ! 'NOSYMPROJ'
      projsym = .false.
    else if (istr == 20) then
      ! 'FIXORB'
      call mma_allocate(itmp,mxorb_cvb,label='itmp')
      ltmp = -1
      call intchk_cvb(itmp,mxorb_cvb,nfxorb,0,'FIXORB',ltmp)
      ifxorb(:) = 0
      do i=1,nfxorb
        ifxorb(itmp(i)) = 1
      end do
      call mma_deallocate(itmp)
    else if (istr == 21) then
      ! 'FIXSTRUC'
      lfxvb = 0
      call mma_maxINT(mxread)
      mxread = mxread/2
      call mma_allocate(itmp,mxread,label='itmp')
      call intchk_cvb(itmp,mxread,nfxvb,0,'FIXSTRUC',lfxvb)
      call mma_allocate(ifxstr,nfxvb,label='ifxstr')
      ifxstr(:) = itmp(1:nfxvb)
      call mma_deallocate(itmp)
    else if (istr == 22) then
      ! 'DELSTRUC'
      lzrvb = 0
      call mma_maxINT(mxread)
      mxread = mxread/2
      call mma_allocate(itmp,mxread,label='itmp')
      call intchk_cvb(itmp,mxread,nzrvb,0,'DELSTRUC',lzrvb)
      call mma_allocate(idelstr,nzrvb,label='idelstr')
      idelstr(:) = itmp(1:nzrvb)
      call mma_deallocate(itmp)
    else if (istr == 23) then
      ! 'FREORB' - not implemented
      call mma_allocate(itmp,mxorb_cvb,label='itmp')
      ltmp = -1
      call intchk_cvb(itmp,mxorb_cvb,nfrorb1,0,'FREORB',ltmp)
      call mma_allocate(itmp2,mxorb_cvb,label='itmp2')
      itmp2(:) = 0
      do i=1,nfrorb1
        itmp2(itmp(i)) = 1
      end do
      nfxorb1 = 0
      do i=1,mxorb_cvb
        if (itmp2(i) == 1) then
          nfxorb1 = nfxorb1+1
          itmp(nfxorb1) = i
        end if
      end do
      call mma_deallocate(itmp)
      nfxorb = max(nfxorb,nfxorb1)
    else if (istr == 24) then
      ! 'FRESTRUC' - not implemented (and code incomplete)
      lfxvb = 1
      call mma_maxINT(mxread)
      mxread = mxread/2
      call mma_allocate(itmp,mxread,label='itmp')
      call intchk_cvb(itmp,mxread,nfxvb,0,'FRESTRUC',lfxvb)
      call mma_allocate(ifxstr,nfxvb,label='ifxstr')
      ifxstr(:) = itmp(1:nfxvb)
      call mma_deallocate(itmp)
    else if (istr == 25) then
      ! 'ORTHCON'
      mxortl = 40
      call orthcon_cvb(iorts,mxortl,mxpair)
    else if (istr == 26) then
      ! 'SADDLE'
      call int_cvb(idum,1,nread,1)
      isaddle = idum(1)
    else if (istr == 27) then
      ! 'SHSTRUC'
      ishstruc = 1
    else if (istr == 28) then
      ! 'VBWEIGHT'
      ivbweights = 0
      do
        call fstring_cvb(weightkw,nwkw,istr2,ncmp,1)
        if (istr2 == 1) then
          if (mod(ivbweights,2) == 0) ivbweights = ivbweights+1
        else if (istr2 == 2) then
          if (mod(ivbweights,4) <= 1) ivbweights = ivbweights+2
        else if (istr2 == 3) then
          if (mod(ivbweights,8) <= 3) ivbweights = ivbweights+4
        else if (istr2 == 4) then
          ivbweights = 0
        else if (istr2 == 5) then
          ivbweights = 7
        end if
        if (istr2 <= 0) exit
      end do
    else if (istr == 29) then
      ! 'CIWEIGHT'
      npcf = 10
      iciweights = 0
      do
        call fstring_cvb(weightkw,nwkw,istr2,ncmp,1)
        if (istr2 == 1) then
          iciweights = ibset(iciweights,0)
        else if (istr2 == 2) then
          iciweights = ibset(iciweights,1)
        else if (istr2 == 3) then
          iciweights = ibset(iciweights,2)
        else if (istr2 == 4) then
          iciweights = 0
        else if (istr2 == 5) then
          iciweights = ibset(ibset(ibset(0,0),1),2)
        end if
        if (istr2 <= 0) exit
      end do
      call int_cvb(idum,1,nread,1)
      npcf = idum(1)
    else if (istr == 30) then
      ! 'SCORR'
      sij = .true.
    else if (istr == 31) then
      ! 'NOSCORR'
      sij = .false.
    else if (istr == 32) then
      ! 'METHOD'
      call fstring_cvb(methkw,nmeth,istr2,ncmp,1)
      if (istr2 /= 0) imethod = istr2
    else if ((istr == 33) .or. (istr == 34)) then
      ! 'OPTIM' or 'OPT'
      call maxdims_cvb()
      call loopcntr_cvb(1)
    else if (istr == 35) then
      ! 'ENDOPTIM'
      call maxdims_cvb()
      call loopcntr_cvb(2)
    else if (istr == 36) then
      ! 'REPORT '
      call maxdims_cvb()
      call loopcntr_cvb(3)
    else if (istr == 37) then
      ! 'ENDREPOR'
      call maxdims_cvb()
      call loopcntr_cvb(4)
    else if (istr == 39) then
      ! 'TUNE'
      call tuneinp_cvb()
    else if (istr == 41) then
      ! 'OPPOSITE'
      ! very good, but useless because this variable is unused
      !opposite = .true.
    else if (istr == 44) then
      ! 'STAT'
      if (firsttime_cvb()) call touch_cvb('STAT')
    else if (istr == 45) then
      ! 'INIT'
      initial = 1
    else if (istr == 46) then
      ! 'NOINIT'
      initial = 0
    else if (istr == 47) then
      ! 'TIDY'
    else if (istr == 48) then
      ! 'PLOC'
      ploc = .true.
    else if (istr == 49) then
      ! 'NOPLOC'
      ploc = .false.
    else if (istr == 50) then
      ! 'ALTERN'
      mxalter = 50
      call int_cvb(idum,1,nread,1)
      mxalter = idum(1)
      call loopcntr2_cvb(5,mxalter)
    else if (istr == 51) then
      ! 'ENDALTER'
      call loopcntr_cvb(6)
    end if

    ! 'ENDVB', 'ENDCASVB', 'END' or unrecognized keyword -- end of input:
    if (istr == 0) exit

  end do

end if

if (.not. allocated(ifxstr)) call mma_allocate(ifxstr,0,label='ifxstr')
if (.not. allocated(idelstr)) call mma_allocate(idelstr,0,label='ifxstr')
if (.not. allocated(symelm)) call mma_allocate(symelm,0,label='symelm')
if (.not. allocated(gsinp)) call mma_allocate(gsinp,0,label='gsinp')
if (.not. allocated(confsinp)) call mma_allocate(confsinp,0,0,label='confsinp')

if (inputmode == 2) then
  ! Input parsing complete for this step ...
  ! ... Work out NORB, NEL, S ...
  call casinfoset_cvb()
  ! ... Do ICONFS before others to get NVB and related info ...
  call mma_allocate(itmp3,noe,nconf,label='confsinp')
  itmp3(:,:) = confsinp(1:noe,:)
  call mma_deallocate(confsinp)
  call move_alloc(itmp3,confsinp)

  if (nfrag <= 1) then
    nMs_fr(1) = 1
    nalf_fr(1,1) = nalf
    nbet_fr(1,1) = nbet
  else
    nMs_fr(1:nfrag) = 1
    nalf_fr(1,1:nfrag) = (nel_fr(1:nfrag)+i2s_fr(1,1:nfrag))/2
    nbet_fr(1,1:nfrag) = nel_fr(1:nfrag)-nalf_fr(1,1:nfrag)
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
      call mma_allocate(itmp3,noe,nconf,label='confsinp')
      itmp3(:,1:nconf-1) = confsinp(:,:)
      call mma_deallocate(confsinp)
      call move_alloc(itmp3,confsinp)
      confsinp(:,iconf_add+2:) = confsinp(:,iconf_add+1:nconf-1)
      confsinp(:,iconf_add+1) = 0
      confsinp(1:min(nel_fr(ifrag),norb),iconf_add+1) = 1
      confsinp(1:nel_fr(ifrag)-norb,iconf_add+1) = 2
    end if
    call cnfcheck_cvb(confsinp(:,iconf_add+1),nconf_fr(ifrag),nel_fr(ifrag))
    call cnfini_cvb(confsinp(:,iconf_add+1),nconf_fr(ifrag),nel_fr(ifrag),nS_fr(ifrag),i2s_fr(1,ifrag),nMs_fr(ifrag), &
                    nalf_fr(1,ifrag),nvbr_fr(ifrag),ndetvb_fr(ifrag),ndetvb2_fr(ifrag),mnion_fr(ifrag),mxion_fr(ifrag), &
                    nconfion_fr(0,ifrag),ifsc)
    iconf_add = iconf_add+nconf_fr(ifrag)
  end do
  ndetvb = 0
  nelcheck = 0
  do i=1,nfrag
    ndetvb = ndetvb+ndetvb_fr(i)
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
      if (i2s_fr(iS,ifrag) /= 0) then
        absym(1) = .false.
        exit
      end if
    end do
    if (kbasis == 6) then
      nS_fr(ifrag) = 1
      i2s_fr(1,ifrag) = i2s_min
    end if
  end do
  absym(2:5) = absym(1)
  nvb = nvb_cvb(kbasis)
  mnion = mnion_fr(1)
  mxion = mxion_fr(1)
  do i=2,nfrag
    mnion = min(mnion,mnion_fr(i))
    mxion = max(mxion,mxion_fr(i))
  end do
  ! ... Now remaining quantities that depend on NORB or NVB ...
  ! SYMELM
  ip_from = 1
  ip_to = 1
  do isyme=1,nsyme
    do iorb=1,norb
      if (ip_from /= ip_to) symelm(ip_to:ip_to+norb-1) = symelm(ip_from:ip_from+norb-1)
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
      if (ifrom /= ito) iorbrel(ito:ito+nmov-1) = iorbrel(ifrom:ifrom+nmov-1)
      ito = ito+nmov
    end if
    ifrom = ifrom+nmov
  end do
  ndimrel = ito-1
  ! IFXSTR
  ito = 0
  do ifrom=1,nfxvb
    if (ifxstr(ifrom) <= nvb) then
      ito = ito+1
      ifxstr(ito) = ifxstr(ifrom)
    end if
  end do
  nfxvb = ito
  ! IDELSTR
  ito = 0
  do ifrom=1,nzrvb
    if (idelstr(ifrom) <= nvb) then
      ito = ito+1
      idelstr(ito) = idelstr(ifrom)
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
  ! (bug? irots is not initialized!)
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
  !need = ihlf_cvb(nbuf)
  !need = need+3*ihlf_cvb(1)+ihlf_cvb(noe*nconf)+mxaobf*norb+ihlf_cvb(norb)+nvbinp+nsyme*norb*norb+ihlf_cvb(ndimrel)+ &
  !       ihlf_cvb(norb)+ihlf_cvb(nfxvb)+ihlf_cvb(nzrvb)+ihlf_cvb(2*nort)+ihlf_cvb(2*ndrot)+ihlf_cvb(2*ndrot)+ihlf_cvb(nsyme)
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
  !call reserv_cvb(need,recinp)
  ioffs = (nbuf+RtoI-1)/RtoI
  call wrioff_cvb(1,recinp,ioffs)
  call wris_cvb([noe],1,recinp,ioffs)
  call wrioff_cvb(2,recinp,ioffs)
  call wris_cvb([nconf],1,recinp,ioffs)
  call wrioff_cvb(3,recinp,ioffs)
  call wris_cvb([kbasiscvb_inp],1,recinp,ioffs)
  call wrioff_cvb(4,recinp,ioffs)
  call wris_cvb(confsinp,noe*nconf,recinp,ioffs)
  call wrioff_cvb(5,recinp,ioffs)
  call wrrs_cvb(orbs,mxaobf*norb,recinp,ioffs)
  call wrioff_cvb(6,recinp,ioffs)
  call wris_cvb(irdorbs,norb,recinp,ioffs)
  call wrioff_cvb(7,recinp,ioffs)
  call wrrs_cvb(gsinp,nvbinp,recinp,ioffs)
  call wrioff_cvb(8,recinp,ioffs)
  call wrrs_cvb(symelm,nsyme*norb*norb,recinp,ioffs)
  call wrioff_cvb(9,recinp,ioffs)

  call dset_cvb(iorbrel,ifxorb,ifxstr,idelstr,iorts,irots,izeta)
else
  call maxdims_cvb()
end if
call bufio_end_cvb()

call mma_deallocate(iorbrel)
call mma_deallocate(ifxorb)
call mma_deallocate(iorts)
call mma_deallocate(irots)
call mma_deallocate(izeta)
call mma_deallocate(irdorbs)
call mma_deallocate(orbs)

call mma_deallocate(ifxstr)
call mma_deallocate(idelstr)
call mma_deallocate(symelm)
call mma_deallocate(gsinp)
call mma_deallocate(confsinp)

return

end subroutine input_cvb
