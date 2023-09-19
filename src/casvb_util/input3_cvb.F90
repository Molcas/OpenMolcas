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

subroutine input3_cvb(iorbrel,mxdimrel,ifxorb,ifxstr,izrstr,iorts,irots,izeta,ip_iconfs,orbs,irdorbs,ip_cvb,ip_symelm,kbasiscvb_inp)

implicit real*8(a-h,o-z)
! ... Files/Hamiltonian available ...
logical, external :: valid_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "inpmod_cvb.fh"
#include "spinb_cvb.fh"
#include "WrkSpc.fh"
parameter(nglob=5,nstrin=51,nendvb=3,nspec=3,ncrit=2,nmeth=12,nwkw=5,ncmp=4)
character*8 global, string, endvb, specl
character*8 crit, weightkw, methkw
character*50 inpstr
logical firsttime_cvb
logical DoCholesky
external firsttime_cvb
dimension global(nglob), string(nstrin), endvb(nendvb), specl(nspec)
dimension crit(ncrit)
dimension weightkw(nwkw), methkw(nmeth)
dimension iorbrel(mxdimrel), ifxorb(mxorb_cvb)
dimension iorts(*), irots(*), izeta(*)
dimension orbs(mxaobf,mxorb_cvb), irdorbs(mxorb_cvb)
dimension idum(1)
save global, string, endvb, specl
save crit, weightkw, methkw
data global/'XXXXxxxx','START   ','GUESS   ','PRINT   ','PREC    '/
data string/'XXXXxxxx','XXXXxxxx','SAVE    ','XXXXxxxx','ORBPERM ','COUPLE  ','MAXITER ','CRIT    ','CASPROJ ','PROJCAS ', &
            'NOCASPRO','NOPROJCA','XXXXxxxx','XXXXxxxx','SYMELM  ','ORBREL  ','XXXXxxxx','SYMPROJ ','NOSYMPRO','FIXORB  ', &
            'FIXSTRUC','DELSTRUC','FREORB  ','FRESTRUC','ORTHCON ','SADDLE  ','SHSTRUC ','VBWEIGHT','CIWEIGHT','SCORR   ', &
            'NOSCORR ','METHOD  ','OPTIM   ','OPT     ','ENDOPTIM','REPORT  ','ENDREPOR','XXXXxxxx','TUNE    ','XXXXxxxx', &
            'OPPOSITE','XXXXxxxx','XXXXxxxx','STAT    ','INIT    ','NOINIT  ','TIDY    ','PLOC    ','NOPLOC  ','ALTERNAT', &
            'ENDALTER'/
data endvb/'ENDVB   ','ENDCASVB','END     '/
data specl/'SERVICE ','MOSCOW  ','PERFLOC '/
data crit/'OVERLAP ','ENERGY  '/
data weightkw/'CHIRGWIN','LOWDIN  ','INVERSE ','NONE    ','ALL     '/
data methkw/'FLETCHER','TRIM    ','TRUSTOPT','DAVIDSON','STEEP   ','VB2CAS  ','AUGHESS ','AUG2    ','CHECK   ','DFLETCH ', &
            'NONE    ','SUPER   '/

call DecideOnCholesky(DoCholesky)
if (DoCholesky) then
  write(6,*) '** Cholesky or RI/DF not yet implemented in CASVB **'
  call abend_cvb()
end if

call fstring_cvb(specl,nspec,istr,ncmp,2)
if (istr == 1) then
  ! 'SERVICE'
  service = .true.
  write(6,'(1x,a,/)') '**** Service mode **** '
  !call service_cvb()
  write(6,*) ' Casvb dummy routine called : SERV'
  return
else if (istr == 2) then
  ! 'MOSCOW'
  service = .true.
  write(6,'(1x,a,/)') '**** MOSCOW mode **** '
  !call moscow_cvb()
  write(6,*) ' Casvb dummy routine called : MOSCOW'
  return
else if (istr == 3) then
  service = .true.
  write(6,'(1x,a,/)') '**** PERFLOC mode **** '
  !call perfloc_plc(3)
  write(6,*) ' Molint dummy routine called : perfloc_plc'
  return
end if

do

  ! CASSCF wavefunction information:
  call casinfoinp_cvb()
  ! VB wavefunction information:
  call fraginp_cvb(ip_iconfs)

  igroup = 0
  call fstring_cvb(global,nglob,istr,ncmp,2)
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
      strtvb = zero
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
      call gsinp_cvb(orbs,irdorbs,ip_cvb,nvbinp,kbasiscvb_inp,mxaobf,mxorb_cvb,kbasis,strtvb)
    else if (istr == 4) then
      ! 'PRINT'
      call int_cvb(ip,10,nread,1)
    else if (istr == 5) then
      ! 'PREC'
      call int_cvb(idum,1,nread,1)
      iprec = idum(1)
      if (iprec < 0) then
        write(6,*) ' Illegal precision :',iprec
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
      write(6,*) ' Too many orbitals in ORBPERM keyword!'
      call abend_cvb()
    end if
    do iorb=1,nread
      if ((abs(iorbprm(iorb)) < 1) .or. (abs(iorbprm(iorb)) > mxorb_cvb)) then
        write(6,'(a,40i3)') ' Illegal orbital label(s) in ORBPERM:',(iorbprm(ior),ior=1,nread)
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
      write(6,*) ' Unrecognized CRIT keyword!'
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
    call symelminp_cvb(ip_symelm,nsyme,tags,izeta,mxirrep,mxorb_cvb,mxsyme,ityp)
  else if (istr == 16) then
    ! 'ORBREL'
    iorb = 0
    jorb = 0
    call int_cvb(idum,1,nread,1)
    iorb = idum(1)
    call int_cvb(idum,1,nread,1)
    jorb = idum(1)
    if ((iorb < 1) .or. (iorb > mxorb_cvb) .or. (jorb < 1) .or. (jorb > mxorb_cvb)) then
      write(6,*) ' Illegal orbital number(s) in ORBREL:',iorb,jorb
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
        write(6,*) ' Too many symmetry elements in ORBREL keyword!'
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
    call izero(isympr,mxirrep)
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
      call imove_cvb(isymv,isympr,mxirrep)
    end if
  else if (istr == 19) then
    ! 'NOSYMPROJ'
    projsym = .false.
  else if (istr == 20) then
    ! 'FIXORB'
    itmp = mstacki_cvb(mxorb_cvb)
    call intchk_cvb(iwork(itmp),mxorb_cvb,nfxorb,0,'FIXORB',-1)
    call izero(ifxorb,mxorb_cvb)
    do i=1,nfxorb
      ifxorb(iwork(i+itmp-1)) = 1
    end do
    call mfreei_cvb(itmp)
  else if (istr == 21) then
    ! 'FIXSTRUC'
    lfxvb = 0
    call mhpfreei_cvb(ifxstr)
    mxread = mavaili_cvb()/2
    ifxstr = mheapi_cvb(mxread)
    call intchk_cvb(iwork(ifxstr),mxread,nfxvb,0,'FIXSTRUC',lfxvb)
    call mrealloci_cvb(ifxstr,nfxvb)
  else if (istr == 22) then
    ! 'DELSTRUC'
    lzrvb = 0
    call mhpfreei_cvb(izrstr)
    mxread = mavaili_cvb()/2
    izrstr = mheapi_cvb(mxread)
    call intchk_cvb(iwork(izrstr),mxread,nzrvb,0,'DELSTRUC',lzrvb)
    call mrealloci_cvb(izrstr,nzrvb)
  else if (istr == 23) then
    ! 'FREORB' - not implemented
    itmp = mstacki_cvb(mxorb_cvb)
    call intchk_cvb(iwork(itmp),mxorb_cvb,nfrorb1,0,'FREORB',-1)
    itmp2 = mstackiz_cvb(mxorb_cvb)
    do i=1,nfrorb1
      iwork(iwork(i+itmp-1)+itmp2-1) = 1
    end do
    nfxorb1 = 0
    do i=1,mxorb_cvb
      if (iwork(i+itmp2-1) == 1) then
        nfxorb1 = nfxorb1+1
        iwork(nfxorb1+itmp-1) = i
      end if
    end do
    call mfreei_cvb(itmp)
    nfxorb = max(nfxorb,nfxorb1)
  else if (istr == 24) then
    ! 'FRESTRUC' - not implemented (and code incomplete)
    lfxvb = 1
    call mhpfreei_cvb(ifxstr)
    mxread = mavaili_cvb()/2
    ifxstr = mheapi_cvb(mxread)
    call intchk_cvb(iwork(ifxstr),mxread,nfxvb,0,'FRESTRUC',lfxvb)
    call mrealloci_cvb(ifxstr,nfxvb)
  else if (istr == 25) then
    ! 'ORTHCON'
    mxgroup = 40
    mxortl = 40
    mxpair = mxorb_cvb*(mxorb_cvb-1)/2
    itmpa = mstacki_cvb(mxorb_cvb*mxorb_cvb)
    itmpb = mstacki_cvb(mxorb_cvb*mxgroup)
    itmpc = mstacki_cvb(mxgroup)
    itmpd = mstacki_cvb(mxortl)
    call orthcon_cvb(iorts,iwork(itmpa),iwork(itmpb),iwork(itmpc),iwork(itmpd),mxortl,mxpair)
    call mfreei_cvb(itmpa)
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
        if (mod(iciweights,2) == 0) iciweights = iciweights+1
      else if (istr2 == 2) then
        if (mod(iciweights,4) <= 1) iciweights = iciweights+2
      else if (istr2 == 3) then
        if (mod(iciweights,8) <= 3) iciweights = iciweights+4
      else if (istr2 == 4) then
        iciweights = 0
      else if (istr2 == 5) then
        iciweights = 7
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
    if (istr2 /= 0) then
      imethod = istr2
    end if
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
    opposite = .true.
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

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(irots)

end subroutine input3_cvb
