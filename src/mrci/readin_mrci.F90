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

subroutine READIN_MRCI()

use mrci_global, only: BNAME, CISEL, CSEL, CSPCK, CTRSH, ENP, ETHRE, GFAC, ICH, ICPF, IFIRST, INDX, INTSY, IORB, IPRINT, IRC, &
                       IREFCI, IREST, IROW, IROOT, ISAB, ITOC17, ITRANS, JJS, JREFX, KBUFF1, LN, LSYM, LUONE, LUSYMB, MAXIT, &
                       MEMPRM, MEMTOT, MEMWRK, MXREF, MXVEC, NASH, NBAS, NBAST, NBMAX, NBTRI, NCMO, NCSHT, NCSPCK, NCOMP, NCVAL, &
                       NDEL, NDMO, NELEC, NFMO, NFRO, NISH, NORB, NORBT, NREF, NRROOT, NSEL, NSM, NSYM, NVIR, NVIRP, NVIRT, &
                       POTNUC, SSEL, SQNLIM, THRORB
use guga_util_global, only: IAD10, nIOCR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp, u5, u6

implicit none
#include "Molcas.fh"
#include "warnings.h"
integer(kind=iwp) :: I, IADD10, iCmd, IDISK, IGFAC, IIN, ILIM, INTNUM, IO, IOM, iOpt, IORBS, IR, iRef, ISC(4), istatus, ISUM, &
                     ISYM, IT, IU, IV, IVA, IX1, IX2, IX3, IX4, IY1, IY2, IY3, IY4, J, jCmd, jEnd, JJ, jStart, LN1, LN2, LV, MXVC, &
                     NAMSIZ, NASHI, NASHT, NBASI, NC, NCSH(8), NCSHI, NDELI, NDELT, NDMOI, NDMOT, NFMOI, NFMOT, NFROI, NFROT, &
                     NINTSY, nIRC, NISHI, NISHT, NIWLK, nJJS, NORBI, NOTOT(8), NREFWR, NRF, NRLN1, nTit, NVAL(8), NVALI, NVALT, &
                     NVIRI, NVT, NVT2
real(kind=wp) :: SPIN
logical(kind=iwp) :: Skip
character(len=88) :: ModLine
character(len=72) :: Line, Title(10)
character(len=4) :: Command
integer(kind=iwp), allocatable :: IOCR(:)
character(len=4), parameter :: Cmd(19) = ['TITL','THRP','PRIN','FROZ','DELE','MAXI','ECON','REST','ROOT','ACPF','SDCI','GVAL', &
                                          'PROR','REFC','SELE','NRRO','MXVE','TRAN','END ']

! Initialize data and set defaults

IOM = MXORB
KBUFF1 = 2*9600
ETHRE = 1.0e-8_wp
SQNLIM = 1.0e-10_wp
CTRSH = 0.05_wp
THRORB = 1.0e-5_wp
ENP = One
NRROOT = 1
NSEL = 0
IPRINT = 1
MAXIT = 20
IREST = 0
ICPF = 0
IREFCI = 0
ITRANS = 0
IGFAC = 0
MXVC = 0
do I=1,8
  NFRO(I) = 0
  NDEL(I) = 0
  NBAS(I) = 0
  NORB(I) = 0
end do
do I=1,IOM+1
  IROW(I) = I*(I-1)/2
end do
do I=1,12
  IROOT(I) = I
end do
nTit = 0

! Read the header of the ONEINT file

NAMSIZ = LENIN8*MXORB
IDISK = 0
call WR_MOTRA_Info(LUONE,2,iDisk,ITOC17,64,POTNUC,NSYM,NBAS,NORB,NFMO,NDMO,8,BNAME,NAMSIZ)

!---  Read input from standard input ----------------------------------*
call RdNLst(u5,'MRCI')
Skip = .false.
jCmd = 0
do
  if (Skip) then
    Skip = .false.
  else
    read(u5,'(A)',iostat=istatus) Line
    if (istatus < 0) call Error(1)
    Command = Line(1:4)
    call UpCase(Command)
    if (Command(1:1) == '*') cycle
    if (Command == ' ') cycle
    jCmd = 0
    do iCmd=1,size(Cmd)
      if (Command == Cmd(iCmd)) jCmd = iCmd
    end do
  end if
  select case (jCmd)

    case default
      write(u6,*) 'READIN Error: Command not recognized.'
      write(u6,*) 'The command is:'//''''//Command//''''
      call QUIT(_RC_INPUT_ERROR_)

    case (1) !TITL
      !---  process TITL command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        Command = Line(1:4)
        call UpCase(Command)
        if (Command(1:1) == '*') cycle
        jCmd = 0
        do iCmd=1,size(Cmd)
          if (Command == Cmd(iCmd)) jCmd = iCmd
        end do
        if (jCmd /= 0) exit
        nTit = nTit+1
        if (nTit <= size(Title)) Title(nTit) = Line
      end do
      Skip = .true.

    case (2) !THRP
      !---  process THRP command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) CTRSH
      if (istatus > 0) call Error(2)

    case (3) !PRIN
      !---  process PRIN command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) IPRINT
      if (istatus > 0) call Error(2)

    case (4) !FROZ
      !---  process FROZ command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NFRO(I),I=1,8)
      if (istatus > 0) call Error(2)

    case (5) !DELE
      !---  process DELE command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NDEL(I),I=1,8)
      if (istatus > 0) call Error(2)

    case (6) !MAXI
      !---  process MAXI command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) MAXIT
      if (istatus > 0) call Error(2)

    case (7) !ECON
      !---  process ECON command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) ETHRE
      if (istatus > 0) call Error(2)

    case (8) !REST
      !---  process REST command --------------------------------------*
      IREST = 1

    case (9) !ROOT
      !---  process ROOT command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) (IROOT(I),I=1,NRROOT)
      if (istatus > 0) call Error(2)

    case (10) !ACPF
      !---  process ACPF command --------------------------------------*
      ICPF = 1

    case (11) !SDCI
      !---  process SDCI command --------------------------------------*
      ICPF = 0

    case (12) !GVAL
      !---  process GVAL command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) GFAC
      if (istatus > 0) call Error(2)
      IGFAC = 1

    case (13) !PROR
      !---  process PROR command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) THRORB
      if (istatus > 0) call Error(2)

    case (14) !REFC
      !---  process REFC command --------------------------------------*
      IREFCI = 1

    case (15) !SELE
      !---  process SELE command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) NSEL
      if (istatus > 0) call Error(2)
      JJ = 0
      do I=1,NSEL
        read(u5,*,iostat=istatus) NC,(CSEL(JJ+J),SSEL(JJ+J),J=1,NC)
        if (istatus < 0) call Error(1)
        if (istatus > 0) call Error(2)
        JJ = JJ+NC
        NCOMP(I) = NC
      end do

    case (16) !NRRO
      !---  process NRRO command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) NRROOT
      if (istatus > 0) call Error(2)
      if (NRROOT > MXVEC) then
        write(u6,1610) NRROOT,MXVEC
        call quit(_RC_INPUT_ERROR_)
      end if
      do I=1,NRROOT
        IROOT(I) = I
      end do

    case (17) !MXVE
      !---  process MXVE command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) MXVC
      if (istatus > 0) call Error(2)
      if (mxvc > mxvec) then
        write(u6,1710) mxvc,mxvec
        call quit(_RC_INPUT_ERROR_)
      end if

    case (18) !TRAN
      !---  process TRAN command --------------------------------------*
      ITRANS = 1

    case (19) !END
      exit

  end select
end do
!---  The end of the input is reached, print the title ----------------*
if (ntit == 0) then
  ntit = 1
  title(1) = ' ( No title was given )'
end if
write(u6,*)
write(u6,'(6X,A)') repeat('*',120)
write(u6,'(6X,A,A,A)') '*',repeat(' ',118),'*'
write(u6,'(6X,A,A,A6,A,A)') '*',repeat(' ',56),'Title:',repeat(' ',56),'*'
do i=1,nTit
  call Center_Text(Title(i))
  write(u6,'(6X,A,A,A72,A,A)') '*',repeat(' ',23),Title(i),repeat(' ',23),'*'
end do
write(u6,'(6X,A,A,A1)') '*',repeat(' ',118),'*'
write(u6,'(6X,A)') repeat('*',120)
write(u6,*)

!---  print the coordinates of the system -----------------------------*
call PrCoor()

!---  read the header of CIGUGA ---------------------------------------*

! Read the header of the CIGUGA file

IADD10 = 0
call iDAFILE(LUSYMB,2,IAD10,9,IADD10)
iOpt = 2
nJJS = 18
nIRC = 4
call mma_allocate(IOCR,nIOCR,label='IOCR')
call WR_GUGA(LUSYMB,iOpt,IADD10,NREF,SPIN,NELEC,LN,NSYM,NCSPCK,NINTSY,IFIRST,INTNUM,LSYM,NRF,LN1,NRLN1,NASH,NISH,8,IRC,nIRC,JJS, &
             nJJS,NVAL,IOCR,nIOCR)
if (ICPF == 1) then
  write(u6,*) '      THIS IS AN   A C P F   CALCULATION'
else
  write(u6,*) '      THIS IS AN   S D C I   CALCULATION'
  write(u6,*) '      (But an ACPF correction will be computed)'
end if
if (IGFAC == 0) then
  GFAC = Two/NELEC
  write(u6,*) '      USE THE DEFAULT ACPF G-VALUE GFAC=',GFAC
else
  write(u6,*) '      THE ACPF G-VALUE HAS BEEN SET TO GFAC=',GFAC
end if
write(u6,*)
if (IREST /= 0) write(u6,*) '      RESTARTED CALCULATION.'
write(u6,*) '      A SMALL CI IS PERFORMED INVOLVING ONLY THE REFERENCE STATES.'
write(u6,*) '      THIS REFERENCE CI WILL USE THE FOLLOWING ROOT SELECTION CRITERIA:'
write(u6,*)
write(u6,*)
if (MXVC == 0) MXVC = max(NRROOT,10)
if (NSEL == 0) then
  write(u6,*) '      ROOT SELECTION BY ENERGY ORDERING.'
  if (NRROOT == 1) then
    write(u6,'(A,I8)') '      ONE SINGLE ROOT, NUMBER ',IROOT(1)
  else
    write(u6,*) '      THE FOLLOWING ROOTS WILL BE SELECTED:'
    write(u6,'(1X,/(1x,12I3))') (IROOT(I),I=1,NRROOT)
  end if
else
  write(u6,*) '      ROOT SELECTION BY PROJECTION: THE EIGENVECTORS OF'
  write(u6,*) '      THE REFERENCE CI ARE ORDERED BY DECREASING SIZE OF'
  write(u6,*) '      THEIR PROJECTIONS ONTO A SELECTION SPACE.'
  if (NRROOT == 1) then
    write(u6,*) '     SELECT THE EIGENVECTOR WITH LARGEST PROJECTION.'
  else
    write(u6,'(A,I2,A)') '      SELECT THE ',NRROOT,' EIGENVECTORS WITH LARGEST PROJECTION.'
  end if
  write(u6,*) '      THE SELECTION SPACE IS SPANNED BY THE FOLLOWING VECTORS (NONZERO COMPONENTS ONLY):'
  JJ = 0
  do I=1,NSEL
    write(u6,'(6X,A,I2)') ' VECTOR NR. ',I
    NC = NCOMP(I)
    write(u6,'(11X,I2,5X,A20,F12.8)') (J,SSEL(JJ+J),CSEL(JJ+J),J=1,NC)
    JJ = JJ+NC
  end do
end if
write(u6,*)
if (IREFCI == 0) then
  if (IREST == 0) then
    write(u6,*) '      THE REFERENCE CI IS FOLLOWED BY THE FULL SPACE'
    write(u6,*) '      CALCULATION, WHERE THE SELECTION CRITERION'
    write(u6,*) '      IS MAXIMUM OVERLAP WITH THE ROOT(S) SELECTED IN'
    write(u6,*) '      THE REFERENCE CI.'
  else
    write(u6,*) '      THE REFERENCE CI IS FOLLOWED BY THE FULL SPACE'
    write(u6,*) '      CALCULATION, WITH ITERATIONS RESTARTED FROM'
    write(u6,*) '      CI VECTOR(S) READ FROM FILE. THE ROOT SELECTION'
    write(u6,*) '      CRITERION IS MAXIMUM OVERLAP WITH THE START'
    write(u6,*) '      VECTORS.'
  end if
else
  write(u6,*) '      ONLY THE REFERENCE CI WAS REQUESTED.'
end if
if (LN > IOM) then
  write(u6,*) 'READIN Error: Too many orbitals.'
  write(u6,'(1X,A,2I5)') 'actual,allowed:',LN,IOM
  call QUIT(_RC_INPUT_ERROR_)
end if
NISHT = 0
LV = 0
do I=1,NSYM
  NISHT = NISHT+NISH(I)
  LV = LV+NVAL(I)
end do
IIN = 0
IR = 0
IVA = 0
IU = NISHT+LV
IT = LV
IV = LN
NBAST = 0
NORBT = 0
NFMOT = 0
NFROT = 0
NASHT = 0
NVALT = 0
NVIRT = 0
NCSHT = 0
NDELT = 0
NDMOT = 0
do I=1,NSYM
  NORBI = NORB(I)
  NBASI = NBAS(I)
  NFMOI = NFMO(I)
  NFROI = NFRO(I)
  NISHI = NISH(I)
  NASHI = NASH(I)
  NVALI = NVAL(I)
  NDELI = NDEL(I)
  NDMOI = NDMO(I)
  NVIR(I) = NORBI-NFROI-NASHI-NISHI-NVALI-NDELI
  NVIRI = NVIR(I)
  NCSH(I) = NISHI+NASHI+NVALI+NVIRI
  NCSHI = NCSH(I)
  NBAST = NBAST+NBASI
  NORBT = NORBT+NORBI
  NFMOT = NFMOT+NFMOI
  NFROT = NFROT+NFROI
  NASHT = NASHT+NASHI
  NVALT = NVALT+NVALI
  NVIRT = NVIRT+NVIRI
  NCSHT = NCSHT+NCSHI
  NDELT = NDELT+NDELI
  NDMOT = NDMOT+NDMOI
  do J=1,NFROI
    IIN = IIN+1
    IR = IR-1
    ICH(IIN) = IR
  end do
  do J=1,NISHI
    IIN = IIN+1
    IT = IT+1
    ICH(IIN) = IT
    NSM(IT) = I
  end do
  do J=1,NASHI
    IIN = IIN+1
    IU = IU+1
    ICH(IIN) = IU
    NSM(IU) = I
  end do
  do J=1,NVALI
    IIN = IIN+1
    IVA = IVA+1
    ICH(IIN) = IVA
    NSM(IVA) = I
  end do
  do J=1,NVIRI
    IIN = IIN+1
    IV = IV+1
    ICH(IIN) = IV
    NSM(IV) = I
  end do
  do J=1,NDELI
    IIN = IIN+1
    ICH(IIN) = 0
  end do
end do
IORBS = 0
do ISYM=1,NSYM
  NOTOT(ISYM) = 0
end do
do ISYM=1,NSYM
  IO = NOTOT(ISYM)
  do I=1,NFMO(ISYM)+NFRO(ISYM)
    IO = IO+1
  end do
  NOTOT(ISYM) = IO
end do
do ISYM=1,NSYM
  IO = NOTOT(ISYM)
  do I=1,NISH(ISYM)
    IO = IO+1
    IORBS = IORBS+1
    IORB(IORBS) = IO
  end do
  NOTOT(ISYM) = IO
end do
do ISYM=1,NSYM
  IO = NOTOT(ISYM)
  do I=1,NASH(ISYM)
    IO = IO+1
    IORBS = IORBS+1
    IORB(IORBS) = IO
  end do
  NOTOT(ISYM) = IO
end do
do ISYM=1,NSYM
  IO = NOTOT(ISYM)
  do I=1,NVAL(ISYM)
    IO = IO+1
    IORBS = IORBS+1
    IORB(IORBS) = IO
  end do
  NOTOT(ISYM) = IO
end do
do ISYM=1,NSYM
  IO = NOTOT(ISYM)
  do I=1,NVIR(ISYM)
    IO = IO+1
    IORBS = IORBS+1
    IORB(IORBS) = IO
  end do
  NOTOT(ISYM) = IO
end do
! NR OF VIRTUALS IN PREVIOUS SYMMETRIES:
ISUM = 0
do I=1,NSYM
  NVIRP(I) = ISUM
  ISUM = ISUM+NVIR(I)
end do
NCMO = 0
NBMAX = 0
do I=1,NSYM
  if (NBAS(I) > NBMAX) NBMAX = NBAS(I)
  NCMO = NCMO+NBAS(I)**2
end do
NBTRI = (NBAST*(NBAST+1))/2
NVT = IROW(NVIRT+1)
NVT2 = IROW(NVIRT)
write(u6,*)
write(u6,'(A)') '      MALMQVIST DIAGONALIZATION'
write(u6,*)
write(u6,'(A,I8)') '      PRINT LEVEL                   ',IPRINT
write(u6,'(A,I12)') '      WORKSPACE WORDS, (Re(wp)) ',MEMTOT
write(u6,'(A,I8)') '      MAXIMUM NR OF ORBITALS        ',IOM
write(u6,'(A,I8)') '      MAX NR OF STORED CI/SGM ARR.  ',MXVC
write(u6,'(A,I8)') '      MAX NR OF ITERATIONS          ',MAXIT
write(u6,'(A,D9.2)') '      ENERGY CONVERGENCE THRESHOLD ',ETHRE
write(u6,'(A,F8.1)') '      SPIN QUANTUM NUMBER           ',SPIN
write(u6,'(A,I8)') '      CORRELATED ELECTRONS          ',NELEC
write(u6,'(A,I8)') '      WAVE FUNCTION SYMMETRY LABEL  ',LSYM
write(u6,'(A,I8)') '      POINT GROUP ORDER             ',NSYM
write(u6,*)
write(u6,101) 'SYMMETRY LABEL:',(I,I=1,NSYM)
write(u6,101) 'INACTIVE ORBITALS',(NISH(I),I=1,NSYM),NISHT
write(u6,101) 'ACTIVE ORBITALS',(NASH(I),I=1,NSYM),NASHT
write(u6,101) 'ADDED VALENCE ORB',(NVAL(I),I=1,NSYM),NVALT
write(u6,101) 'VIRTUAL ORBITALS',(NVIR(I),I=1,NSYM),NVIRT
write(u6,*)
write(u6,101) 'SUM:CORREL ORBITALS',(NCSH(I),I=1,NSYM),NCSHT
write(u6,*)
write(u6,101) 'FROZEN ORBITALS',(NFRO(I),I=1,NSYM),NFROT
write(u6,101) 'DELETED ORBITALS',(NDEL(I),I=1,NSYM),NDELT
write(u6,*)
write(u6,101) 'SUM:ORBITALS IN CI',(NORB(I),I=1,NSYM),NORBT
write(u6,*)
write(u6,101) 'PRE-FROZEN ORBITALS',(NFMO(I),I=1,NSYM),NFMOT
write(u6,101) 'PRE-DELETED ORBITALS',(NDMO(I),I=1,NSYM),NDMOT
write(u6,101) 'SUM:   TOTAL BASIS',(NBAS(I),I=1,NSYM),NBAST
write(u6,*)
if (LN1 == 0) then
  write(u6,*) '      ONE CLOSED SHELL REFERENCE STATE'
else
  write(u6,'(6X,I4,A)') NREF,' REFERENCE STATES'
  NREFWR = min(NREF,1000/LN1)
  LN2 = min(32,LN1)
  write(u6,'(6X,A,T47)') 'Occupation of the reference states'
  if (NREFWR < NREF) then
    write(u6,'(6X,A,I3,A)') '( Only the ',NREFWR,' first are listed)'
  end if
  write(u6,'(6X,A,T25,32I2)') 'Active orbital nr.',(I,I=1,LN2)
  jEnd = 0
  do iRef=1,NREFWR
    jStart = jEnd+1
    jEnd = jEnd+LN1
    write(u6,'(6X,A,I3,T25,32I2)') 'Ref nr',IREF,(IOCR(j),j=jStart,jStart-1+LN2)
  end do
end if
call mma_deallocate(IOCR)
write(u6,*)
if (INTNUM /= 0) write(u6,*) '      FIRST ORDER INTERACTING SPACE.'
IX1 = IRC(1)
IX2 = IRC(2)-IRC(1)
ISC(1) = IX1
ISC(2) = ISC(1)+IX2*NVIRT
IY1 = ISC(1)
IY2 = ISC(2)-ISC(1)
if (IFIRST == 0) then
  ILIM = 4
  IX3 = IRC(3)-IRC(2)
  IX4 = IRC(4)-IRC(3)
  ISC(3) = ISC(2)+IX3*NVT2
  ISC(4) = ISC(3)+IX4*NVT
  IY3 = ISC(3)-ISC(2)
  IY4 = ISC(4)-ISC(3)
  if (IPRINT >= 10) then
    write(u6,*)
    write(u6,*) '      INTERNAL WALKS:'
    write(u6,215) IX1,IX2,IX3,IX4
    write(u6,*)
    write(u6,*) '      FORMAL CONFIGURATIONS:'
    write(u6,215) IY1,IY2,IY3,IY4
    write(u6,'(6X,A,I7)') '                  TOTAL:',ISC(ILIM)
  end if
else
  ILIM = 2
  if (IPRINT >= 10) then
    write(u6,*)
    write(u6,*) '      INTERNAL WALKS:'
    write(u6,216) IX1,IX2
    write(u6,*)
    write(u6,*) '      FORMAL CONFIGURATIONS:'
    write(u6,216) IY1,IY2
    write(u6,'(6X,A,I7)') '                  TOTAL:',ISC(ILIM)
  end if
end if
NIWLK = IRC(ILIM)
NCVAL = IRC(1)
! ----------------------------------------------------------------------
if (NVIRT > 255) then
  write(u6,*)
  write(u6,*) ' Sorry -- The MRCI code uses internal integer codes'
  write(u6,*) ' where the index of virtual orbitals is kept in'
  write(u6,*) ' 8-bit fields. This cannot easily be increased'
  write(u6,*) ' and limits the number of virtual orbitals to '
  write(u6,*) ' 255. Your input asks for more virtuals than this.'
  write(u6,*) ' The program cannot run.'
  call Quit(_RC_INPUT_ERROR_)
end if
! ----------------------------------------------------------------------
! ALLOCATION OF DATA PERMANENTLY IN CORE
!

! CSPCK - ARRAY OF BIT-PACKED GUGA CASE NUMBERS OF INTERNAL WALKS.
! CONSISTS OF NCSPCK INTEGERS.

call mma_allocate(CSPCK,NCSPCK,label='CSPCK')
call iDAFILE(LUSYMB,2,CSPCK,NCSPCK,IADD10)

! INTSY - ARRAY OF BIT-PACKED SYMMETRY LABELS OF INTERNAL WALKS.
! CONSISTS OF NINTSY INTEGERS.

call mma_allocate(INTSY,NINTSY,label='INTSY')
call iDAFILE(LUSYMB,2,INTSY,NINTSY,IADD10)

! INDX - START POSITION IN CI ARRAY OF EACH INTERNAL-WALK-BLOCK

call mma_allocate(INDX,NIWLK,label='INDX')

! ISAB - ORDERING NR OF EACH VIRTUAL PAIR WITHIN ITS COMB-SYMM

call mma_allocate(ISAB,NVIRT,NVIRT,label='ISAB')

! JREFX - FOR EACH VALENCE CSF, EITHER 0 OR ITS REFERENCE NR.

call mma_allocate(JREFX,NCVAL,label='JREFX')
IADD10 = IAD10(2)
call iDAFILE(LUSYMB,2,JREFX,NCVAL,IADD10)

! PROJECTION SELECTION VECTORS

call mma_allocate(CISEL,NREF,NSEL,label='CISEL')

! START OF NON-PERMANENT AREA:

call INDMAT(CSPCK,INTSY,INDX,ISAB,JREFX,CISEL)
if (NREF > MXREF) then
  write(u6,*) 'READIN Error: Too many references.'
  write(u6,'(1X,A,2I6)') ' actual, allowed:',NREF,MXREF
  call QUIT(_RC_INPUT_ERROR_)
end if

! Total available memory (at start of program) is MEMTOT
! Available now is MEMWRK
! Already (permanently) allocated is MEMPRM

call mma_maxdble(MemWrk)
MEMPRM = MEMTOT-MEMWRK

call ALLOC_MRCI()

return

101 format(6X,A,T47,9I5)
215 format(/,6X,'                 VALENCE',I7,/,6X,' DOUBLET COUPLED SINGLES',I7,/,6X,' TRIPLET COUPLED DOUBLES',I7, &
           /,6X,' SINGLET COUPLED DOUBLES',I7)
216 format(/,6X,'                 VALENCE',I7,/,6X,' DOUBLET COUPLED SINGLES',I7)
1610 format('Too many roots,',i3,', max allowed is',i3)
1710 format('Too many vectors,',i3,', max allowed is',i3)

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      write(u6,*) 'READIN Error: Premature end of file while reading.'
    case (2)
      write(u6,*) 'READIN Error: I/O error during internal read.'
      write(u6,*) 'The line that could not be read is:'
      write(u6,*) Line
  end select
  call Quit(_RC_IO_ERROR_READ_)

end subroutine Error

end subroutine READIN_MRCI
