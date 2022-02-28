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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine ReadIn_CPF()
! Read input and allocate memory

use cpf_global, only: BNAME, CTRSH, ETHRE, ETOT, ICASE, ICH, ICONV, ICPF, IFIRST, ILIM, INCPF, INDX, IPRINT, IR1, IRC, IREST, &
                      IROW, ISAB, ISC, ISDCI, ISMAX, ITOC17, IV0, IV1, JJS, JSY, LN, LSYM, Lu_CIGuga, Lu_TraOne, LWSP, MAXIT, &
                      MAXITP, N, NASH, NBAS, NFRO, NISH, NORB, NORBT, NPFRO, NREF, NSM, NSYM, NVIR, NVIRT, POTNUC, WLEV
use guga_util_global, only: IAD10, nIOCR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u5, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), parameter :: mxTit = 10
integer(kind=iwp) :: I, IADD10, iCmd, IDISK, IIN, INTNUM, iOpt, IR, iRef, IRJ, istatus, iSym, IT, IU, IV, IVA, IX1, IX2, IX3, IX4, &
                     IY1, IY2, IY3, IY4, j, jCmd, jEnd, jStart, LN1, LN2, NAMSIZ, NASHI, NASHT, NBAST, NDEL(8), NDELI, NDELT, &
                     NFREF, NFROI, NFROT, nIRC, NISHI, NISHT, nJJS, NPDEL(8), NPDELT, NPFROT, NRLN1, nTit, NVAL(8), NVALI, NVALT, &
                     NVIR2, NVIRI, NVT, NVT2
real(kind=wp) :: S
logical(kind=iwp) :: Skip
character(len=88) :: ModLine
character(len=72) :: Line, Title(mxTit)
character(len=4) :: Command
integer(kind=iwp), allocatable :: IOCR(:), JREFX(:)
character(len=4), parameter :: Cmd(16) = ['TITL','MAXP','LEVS','THRP','PRIN','FROZ','DELE','MAXI','ECON','REST','MCPF','CPF ', &
                                          'SDCI','ACPF','LOW ','END ']

!---  Initialize arrays and variables ---------------------------------*
LWSP = .false.
ETHRE = 1.0e-6_wp
CTRSH = 5.0e-2_wp
IPRINT = 5
MAXIT = 20
IREST = 0
ICPF = 0
ISDCI = 0
INCPF = 0
ICONV = 0
MAXITP = 6
WLEV = 0.3_wp
ETOT = Zero
NPFRO(:) = 0
NFRO(:) = 0
NDEL(:) = 0
NPDEL(:) = 0
NISH(:) = 0
NASH(:) = 0
NVAL(:) = 0
NVIR(:) = 0
NORB(:) = 0
NBAS(:) = 0
NPFROT = 0
NFROT = 0
NDELT = 0
NPDELT = 0
NISHT = 0
NASHT = 0
NVALT = 0
NVIRT = 0
NORBT = 0
NBAST = 0
do I=1,size(IROW)
  IROW(I) = I*(I-1)/2
end do
nTit = 0

!---  read the header of TRAONE ---------------------------------------*
! Note: NORB(i)=NBAS(i)-NPFRO(i)-NPDEL(i)
NAMSIZ = LenIn8*MXORB
IDISK = 0
call WR_MOTRA_Info(Lu_TraOne,2,iDisk,ITOC17,64,POTNUC,NSYM,NBAS,NORB,NPFRO,NPDEL,8,BNAME,NAMSIZ)

!---  Read input from standard input ----------------------------------*
call RdNLst(u5,'CPF')
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
    jCmd = 0
    do iCmd=1,size(Cmd)
      if (Command == Cmd(iCmd)) jCmd = iCmd
    end do
  end if
  select case (jCmd)

    case default
      write(u6,*) 'READIN Error: Command not recognized.'
      write(u6,*) 'The command is:'//''''//Command//''''
      call QUIT_OnUserError()

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
        if (nTit <= mxTit) Title(nTit) = Line
      end do
      Skip = .true.

    case (2) !MAXP
      !---  process MAXP command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) MaxItP
      if (istatus > 0) call Error(2)

    case (3) !LEVS
      !---  process LEVS command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) WLev
      if (istatus > 0) call Error(2)

    case (4) !THRP
      !---  process THRP command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) CTrsh
      if (istatus > 0) call Error(2)

    case (5) !PRIN
      !---  process PRIN command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) iPrint
      if (istatus > 0) call Error(2)

    case (6) !FROZ
      !---  process FROZ command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (nFro(i),i=1,8)
      if (istatus > 0) call Error(2)

    case (7) !DELE
      !---  process DELE command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NDEL(i),i=1,8)
      if (istatus > 0) call Error(2)

    case (8) !MAXI
      !---  process MAXI command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) MaxIt
      if (istatus > 0) call Error(2)
      MaxIt = min(MaxIt,75)

    case (9) !ECON
      !---  process ECON command --------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) EThre
      if (istatus > 0) call Error(2)

    case (10) !REST
      !---  process REST command --------------------------------------*
      iRest = 1

    case (11) !MCPF
      !---  process MCPF command --------------------------------------*
      iCPF = 0
      iSDCI = 0
      iNCPF = 0

    case (12) !CPF
      !---  process CPF  command --------------------------------------*
      iCPF = 1
      iSDCI = 0
      iNCPF = 0

    case (13) !SDCI
      !---  process SDCI command --------------------------------------*
      iSDCI = 1
      iCPF = 0
      iNCPF = 0

    case (14) !ACPF
      !---  process ACPF command --------------------------------------*
      iNCPF = 1
      iCPF = 0
      iSDCI = 0

    case (15) !LOW
      !---  process LOW  command --------------------------------------*
      LWSP = .true.

    case (16) !END
      exit

  end select
end do
!---  The end of the input is reached, print the title ----------------*
if (ntit == 0) then
  ntit = 1
  title(1) = ' ( No title was given )'
end if
write(u6,*)
write(u6,'(6X,120A1)') ('*',i=1,120)
write(u6,'(6X,120A1)') '*',(' ',i=1,118),'*'
write(u6,'(6X,57A1,A6,57A1)') '*',(' ',i=1,56),'Title:',(' ',i=1,56),'*'
do i=1,nTit
  call Center_Text(Title(i))
  write(u6,'(6X,24A1,A72,24A1)') '*',(' ',j=1,23),Title(i),(' ',j=1,23),'*'
end do
write(u6,'(6X,120A1)') '*',(' ',i=1,118),'*'
write(u6,'(6X,120A1)') ('*',i=1,120)
write(u6,*)

!---  print the coordinates of the system -----------------------------*
call PrCoor()

!---  print the method used -------------------------------------------*
write(u6,*)
if (iSDCI == 1) then
  write(u6,'(6X,A)') 'This is an  S D C I  calculation'
else if (iCPF == 1) then
  write(u6,'(6X,A)') 'This is a  C P F  calculation'
else if (INCPF == 1) then
  write(u6,'(6X,A)') 'This is an  A C P F calculation'
else
  write(u6,'(6X,A)') 'This is an  M C P F  calculation'
end if
if (LWSP) write(u6,'(6X,A)') 'This is a LOW SPIN calculation'

!---  read the header of CIGUGA ---------------------------------------*
IADD10 = 0
call iDAFILE(Lu_CIGuga,2,IAD10,9,IADD10)
iOpt = 2
nJJS = 18
nIRC = 4
call mma_allocate(IOCR,nIOCR,label='IOCR')
call WR_GUGA(Lu_CIGuga,iOpt,IADD10,NFREF,S,N,LN,NSYM,IR1,IRJ,IFIRST,INTNUM,LSYM,NREF,LN1,NRLN1,NASH,NISH,8,IRC,nIRC,JJS,nJJS,NVAL, &
             IOCR,nIOCR)
if (LN >= MXORB) then
  write(u6,*) 'READIN Error: Too many orbitals.'
  write(u6,'(1X,A,2I5)') 'LN,MXORB:',LN,MXORB
  call QUIT_OnUserError()
end if
call mma_allocate(ICASE,IR1,label='ICASE')
call iDAFILE(Lu_CIGuga,2,ICASE,IR1,IADD10)
call mma_allocate(JSY,IRJ,label='JSY')
call iDAFILE(Lu_CIGuga,2,JSY,IRJ,IADD10)

!---  update orbital specifications -----------------------------------*
IV0 = 0
IV1 = 1
do I=1,NSYM
  NVIR(I) = NORB(I)-NFRO(I)-NISH(I)-NASH(I)-NVAL(I)-NDEL(I)
  NPFROT = NPFROT+NPFRO(I)
  NFROT = NFROT+NFRO(I)
  NISHT = NISHT+NISH(I)
  NASHT = NASHT+NASH(I)
  NVALT = NVALT+NVAL(I)
  NVIRT = NVIRT+NVIR(I)
  NDELT = NDELT+NDEL(I)
  NPDELT = NPDELT+NPDEL(I)
  NORBT = NORBT+NORB(I)
  NBAST = NBAST+NBAS(I)
end do
IIN = 0
IR = 0
IVA = 0
IU = NISHT+NVALT
IT = NVALT
IV = LN
NVIRT = 0
do I=1,NSYM
  NFROI = NFRO(I)
  NISHI = NISH(I)
  NASHI = NASH(I)
  NVALI = NVAL(I)
  NDELI = NDEL(I)
  !PAM97 NVIRDI = NORB(I)-NFROI-NASHI-NISHI-NVALI
  !PAM97 NVIR(I) = NVIRDI-NDEL(I)
  NVIRI = NVIR(I)
  NVIRT = NVIRT+NVIRI
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
NVT = IROW(NVIRT+1)
NVT2 = IROW(NVIRT)

!---  report input specifications -------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'ONE-ELECTRON BASIS:'
write(u6,'(6X,A)') '----------------------------'
write(u6,*)
write(u6,'(6X,A,T47,4X,4X,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
write(u6,*)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals pre-frozen in MOTRA',nPFroT,(nPFro(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals used by this program',norbT,(nOrb(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Pre-deleted in MOTRA',nPDelT,(nPDel(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Sum: No. of basis functions',nBasT,(nBas(iSym),iSym=1,nSym)
write(u6,*)
write(u6,'(6X,A)') 'ORBITAL SPECIFICATION:'
write(u6,'(6X,A)') '-------------------------------'
write(u6,*)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals frozen here',nFroT,(nFro(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Inactive orbitals',nIShT,(nISh(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Active orbitals',nAShT,(nASh(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Additional valence orbitals',nValT,(nVal(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Virtual orbitals',nVirT,(nVir(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals deleted here',nDelT,(nDel(iSym),iSym=1,nSym)
write(u6,'(6X,A,T47,I4,4x,8I4)') 'Sum: Total no. of orbitals',norbT,(nOrb(iSym),iSym=1,nSym)
write(u6,*)
write(u6,*)
write(u6,'(6X,A)') 'WAVE FUNCTION SPECIFICATION:'
write(u6,'(6X,A)') '----------------------------'
write(u6,*)
write(u6,'(6X,A,T47,I4)') 'Number of electrons in CI',N
write(u6,'(6X,A,T47,I4)') 'Internal orbitals in CI',LN
write(u6,'(6X,A,T47,I4)') 'External orbitals in CI',NVIRT
write(u6,'(6X,A,T47,I4)') 'Number of irreps',NSYM
write(u6,'(6X,A,T47,F4.1)') 'Spin quantum number',S
write(u6,'(6X,A,T47,I4)') 'State symmetry',LSYM
write(u6,*)
write(u6,'(6X,A)') 'REFERENCE STATE:'
write(u6,'(6X,A)') '------------------------------'
write(u6,*)
write(u6,'(6X,A,T47,I4)') 'Number of reference states',NREF
LN2 = min(16,LN1)
if (LN1 == 0) then
  write(u6,'(6X,A,T47)') 'One closed shell reference state'
else
  write(u6,'(6X,A,T47)') 'Occupation of active orbitals in the reference state:'
  write(u6,'(6X,A,T25,16I4)') 'Active orbital nr.',(I,I=1,LN2)
  jEnd = 0
  do iRef=1,nRef
    jStart = jEnd+1
    jEnd = jEnd+LN1
    write(u6,'(6X,A,I3,T25,16I4)') 'Ref nr',IREF,(IOCR(j),j=jStart,jEnd)
  end do
end if
call mma_deallocate(IOCR)
write(u6,*)
write(u6,'(6X,A)') 'OPTIONS:'
write(u6,'(6X,A)') '--------'
write(u6,*)
write(u6,'(6X,A,T47,I4)') 'Print parameter',iPrint
write(u6,'(6X,A,T47)') 'Pulay diagonalization'
if (INTNUM /= 0) write(u6,'(6X,A,T47)') 'First order interacting space'
IX1 = IRC(1)
IX2 = IRC(2)-IRC(1)
ISC(1) = IX1
ISC(2) = ISC(1)+IX2*NVIRT
IY1 = ISC(1)
IY2 = ISC(2)-ISC(1)
write(u6,214)
if (IFIRST == 0) then
  IX3 = IRC(3)-IRC(2)
  IX4 = IRC(4)-IRC(3)
  ISC(3) = ISC(2)+IX3*NVT2
  ISC(4) = ISC(3)+IX4*NVT
  IY3 = ISC(3)-ISC(2)
  IY4 = ISC(4)-ISC(3)
  write(u6,215) IX1,IX2,IX3,IX4
  write(u6,213)
  write(u6,215) IY1,IY2,IY3,IY4
else
  write(u6,216) IX1,IX2
  write(u6,213)
  write(u6,216) IY1,IY2
end if
ILIM = 4
if (IFIRST /= 0) ILIM = 2
! ERROR CONDITIONS:
!if (LN /= NISHT+NASHT+NVALT) then
!  write(u6,*) ' ERROR: Orbital specifications do not match'
!  write(u6,*) ' input to GUGA. The number of internal'
!  write(u6,*) ' orbitals must equal the number of inactive,'
!  write(u6,*) ' active, and additional valence orbitals.'
!  call QUIT(20)
!end if
! ALLOCATION FOR INDEX VECTORS
! THESE VECTORS ARE PERMANENTLY IN CORE
! -- INDEX
call mma_allocate(INDX,IRC(ILIM),label='INDX')
NVIR2 = NVIRT*NVIRT
call mma_allocate(ISAB,NVIR2,label='ISAB')
call mma_allocate(JREFX,ISC(1),label='JREFX')
IADD10 = IAD10(2)
call iDAFILE(Lu_CIGuga,2,JREFX,ISC(1),IADD10)
call INDMAT_CPF(JSY,INDX,ISAB,ISMAX,JREFX)
call mma_deallocate(JREFX)
call ALLOC_CPF()

return

214 format(//,6X,'INTERNAL CONFIGURATIONS')
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7, &
            /,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
213 format(//,6X,'FULL-SPACE CONFIGURATIONS (FORMAL)')
216 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)

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
  call Quit_OnUserError()

end subroutine Error

end subroutine ReadIn_CPF
