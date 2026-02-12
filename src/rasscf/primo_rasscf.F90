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

subroutine Primo_RasScf(VecTit,Ene,Occ,CMO)
!***********************************************************************
! purpose:
! Print MO-coefficients
!
!***********************************************************************

use rasscf_global, only: iPT2, OutFmt2, PreThr, ProThr, BName
use output_ras, only: LF
use general_data, only: NSYM, NASH, NBAS, NDEL, NFRO, NISH
use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
character(len=*) VecTit
real*8 CMO(*), Occ(*), Ene(*)
integer NSLCT(8)
logical PrOcc, PrEne
character(len=3) lIrrep(8)
character(len=8) Fmt1, Fmt2
character(len=132) Line, Blank
character(len=LenIn+8), external :: Clean_BName
integer NVSH(8)
integer, allocatable :: MrkIt(:), Slct(:)
real*8 CC
integer I, ib, iBas, iBOff, iCOff, iCol, IO, IORB, IS, ISEND, ISOFF, ISStart, IST, ISYM, left, lPaper, NB, NBTOT, nCols, ND, NFIA, &
        NO, NS, NSKIP, NSLCTT, lLine

call Get_cArray('Irreps',lIrrep,24)

! set print format
nCols = 10
lLine = 120
lPaper = 132
left = (lPaper-lLine)/2
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
Blank = ' '
! PAM Krapperup Nov 05: For the moment, selection of orbitals to be
! printed is ultimately determined here on the basis of PRETHR (and
! PROTHR) thresholds. These have either been set by the user, or
! determined in the CHKINP subroutine, possibly then on basis of
! user specification of OutFmt1 flags.
! Exception: OutFmt2='NOCORE  ' will inhibit printing of non-valece orbs.

! PAM Nov 09: Output marked as collapsible.
! print header
write(LF,*)
call CollapseOutput(1,'   Molecular orbitals:')
write(LF,'(6X,A)') '-------------------'
write(LF,*)
write(LF,Fmt2//'A)') trim(VecTit)
write(LF,*)

! Flag ipt2 in common in module rasscf_global.F90.
!   ipt2=0 means usual MO's, quasicanonical for
! inactives and virtuals, natural for active.
!   ipt2=1 means quasicanonical for actives also.
! PAM Apr 05: The rules for selecting orbitals for print
! appear a bit confused. I just follow the rules for now:
if (iPT2 == 0) then
  PrOcc = .true.
  PrEne = .true.
else
  PrOcc = .false.
  PrEne = .true.
end if

! Select orbitals to be printed.
! By default at least all occupied are printed.
! If iPT2==0 all orbitals with energy less than PrEThr
! and occupation bigger then PrOThr are also printed.
! If iPT2==1 all orbitals with occupations greater than
! PrOThr are also printed.

! Nr of orbitals=nr of basis functions
NBTOT = 0
do ISYM=1,NSYM
  NBTOT = NBTOT+NBAS(ISYM)
end do

! Put rasscf orbital energies on the runfile
call Put_darray('RASSCF OrbE',ENE,NBTOT)

! Initialize MARKIT.
call mma_allocate(MRKIT,NBTOT,Label='MrkIt')
MrkIt(:) = 0
if (OutFmt2 /= 'NOCORE  ') then
  ! Mark MARKIT as selected for occupied orbitals.
  IORB = 0
  do ISYM=1,NSYM
    NFIA = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
    do I=1,NFIA
      IORB = IORB+1
      MRKIT(IORB) = 1
    end do
    IORB = IORB+NBAS(ISYM)-NFIA
  end do
else
  ! Mark MARKIT as selected for valence or active orbitals.
  IORB = 0
  call Get_iArray('Non valence orbitals',NVSH,nSym)
  do ISYM=1,NSYM
    NFIA = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
    NSKIP = min(NFRO(ISYM)+NISH(ISYM),NVSH(ISYM))
    IORB = IORB+NSKIP
    do I=NSKIP+1,NFIA
      IORB = IORB+1
      MRKIT(IORB) = 1
    end do
    IORB = IORB+NBAS(ISYM)-NFIA
  end do
end if
! If PROCC, then only those orbitals that have occ no larger
! than or equal to PrOThr:
if (PROCC .and. (PrOThr >= 0.0d0)) then
  IORB = 0
  do ISYM=1,NSYM
    NFIA = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
    do I=1,NFIA
      IORB = IORB+1
      if (OCC(IORB) < PROTHR) MRKIT(IORB) = 0
    end do
    IORB = IORB+NBAS(ISYM)-NFIA
  end do
end if
! But if PRENE, then also those orbitals that have energy less
! than or equal to PrEThr, skipping deleted orbitals of course.
if (PRENE) then
  IORB = 0
  do ISYM=1,NSYM
    NFIA = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
    IORB = IORB+NFIA
    ND = NDEL(ISYM)
    do I=NFIA+1,NBAS(ISYM)-ND
      IORB = IORB+1
      if (ENE(IORB) <= PRETHR) MRKIT(IORB) = 1
    end do
    IORB = IORB+ND
  end do
end if
! Let ISELECT enumerate the orbitals to be printed, rather than
! just marking them:
NSLCTT = 0
NO = 0
do ISYM=1,NSYM
  NS = 0
  NB = NBAS(ISYM)
  do IO=NO+1,NO+NB
    if (MRKIT(IO) == 1) NS = NS+1
  end do
  NO = NO+NB
  NSLCT(ISYM) = NS
  NSLCTT = NSLCTT+NS
end do
call mma_allocate(SLCT,NSLCTT,Label='SLCT')
IS = 0
NO = 0
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  do IO=NO+1,NO+NB
    if (MRKIT(IO) == 1) then
      IS = IS+1
      SLCT(IS) = IO
    end if
  end do
  NO = NO+NB
end do
! Get rid of MARKIT.
call mma_deallocate(MRKIT)

! finally, print, the MOs
if (OutFmt2 == 'FULL    ') then

  ! print orbitals, using default format.
  ISOFF = 0
  IBOFF = 0
  ICOFF = 0
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    NS = NSLCT(ISYM)
    if (NS > 0) then
      write(LF,*)
      write(LF,*)
      write(LF,*)
      write(LF,Fmt2//'A,I2,A,A)') 'Molecular orbitals for symmetry species',iSym,': ',lIrrep(iSym)
      do ISSTART=1,NS,NCOLS
        ISEND = min(ISSTART+NCOLS-1,NS)
        write(LF,*)
        write(LF,*)
        write(LF,Fmt2//'A,7X,10I10)') 'Orbital ',(SLCT(ISOFF+I)-IBOFF,I=ISSTART,ISEND)
        if (PRENE) write(LF,Fmt2//'A,7X,10F10.4)') 'Energy  ',(ENE(SLCT(ISOFF+I)),I=ISSTART,ISEND)
        if (PROCC) write(LF,Fmt2//'A,7X,10F10.4)') 'Occ. No.',(OCC(SLCT(ISOFF+I)),I=ISSTART,ISEND)
        write(LF,*)
        do IB=1,NB
          write(LF,'(2X,I4,1X,A,10F10.4)') IB,Clean_BName(BName(IBOFF+IB),LENIN), &
                                           (CMO(ICOFF+(SLCT(ISOFF+I)-1-IBOFF)*NB+IB),I=ISSTART,ISEND)
        end do
      end do
    end if
    ISOFF = ISOFF+NS
    IBOFF = IBOFF+NB
    ICOFF = ICOFF+NB**2
  end do

else if (OutFmt2 == 'COMPACT ') then

  ! print orbitals, using compact format.

  ISOFF = 0
  IBOFF = 0
  ICOFF = 0
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NSLCT(ISYM) /= 0) then
      write(LF,*)
      write(LF,FMT2//'A,I2,A,A)') 'MOLECULAR ORBITALS FOR SYMMETRY SPECIES',ISYM,': ',LIRREP(ISYM)
      write(LF,*)
      if (PROCC .and. PRENE) then
        write(LF,FMT2//'A)') 'INDEX  ENERGY  OCCUPATION COEFFICIENTS ...'
      else if (PROCC) then
        write(LF,FMT2//'A)') 'INDEX  ENERGY  COEFFICIENTS ...'
      else if (PRENE) then
        write(LF,FMT2//'A)') 'INDEX  OCCUPATION  COEFFICIENTS ...'
      else
        write(LF,FMT2//'A)') 'INDEX  COEFFICIENTS ...'
      end if
      do IS=1,NSLCT(ISYM)
        IORB = SLCT(ISOFF+IS)
        ICOL = IORB-IBOFF
        LINE = BLANK
        IST = 1
        write(LINE(IST:132),'(I5)') ICOL
        IST = IST+5
        if (PRENE) then
          write(LINE(IST:132),'(F10.4)') ENE(IORB)
          IST = IST+10
        end if
        if (PROCC) then
          write(LINE(IST:132),'(F10.4)') OCC(IORB)
          IST = IST+10
        end if
        write(LF,FMT2//'A)') LINE
        LINE = BLANK
        IST = 9
        do IBAS=1,NB
          CC = CMO(ICOFF+(ICOL-1)*NB+IBAS)
          if (abs(CC) >= 0.1d0) then
            write(LINE(IST:132),'(I4,1X,A,A,F7.4,A)') IBAS,Clean_BName(BName(IBOFF+IBAS),LENIN),'(',CC,')'
            IST = IST+28
            if (IST > (132-LEFT-28)) then
              write(LF,FMT2//'A)') trim(LINE)
              LINE = BLANK
              IST = 9
            end if
          end if
        end do
        write(LF,FMT2//'A)') trim(LINE)
        LINE = BLANK
        IST = 9
      end do
    end if
    ISOFF = ISOFF+NSLCT(ISYM)
    IBOFF = IBOFF+NB
    ICOFF = ICOFF+NB**2
  end do

end if

call mma_deallocate(SLCT)
! PAM 09: Reset to non-collapsible output before return.
call CollapseOutput(0,'   Molecular orbitals:')
write(6,*)

end subroutine Primo_RasScf
