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

subroutine ReadIn_CPF(H,iH)

use iso_c_binding, only: c_f_pointer, c_loc
use Symmetry_Info, only: SMul => Mul

implicit real*8(A-H,O-Z)
! Read input and allocate memory
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
#include "niocr.fh"
dimension IOCR(nIOCR)
dimension H(*), iH(*)
#include "spin_cpf.fh"
parameter(nCmd=18)
parameter(mxTit=10)
character*4 Command, Cmd(nCmd)
character*72 Line, Title(mxTit)
character*88 ModLine
data Cmd/'TITL','MAXP','LEVS','THRP','PRIN','FROZ','DELE','MAXI','ECON','ETRS','REST','MCPF','CPF ','SDCI','ACPF','LOW ','EXTR', &
         'END '/
! Statement function
!---- convert a pointer in H to a pointer for iH
ipointer(i) = (i-1)*RtoI+1

!---  Initialize arrays and variables ---------------------------------*
Mul(:,:) = SMul(:,:)
KBUFF1 = 2*9600
D0 = 0.0d0
D1 = 1.0d0
D2 = 2.0d0
LWSP = .false.
SQ2 = sqrt(D2)
ETHRE = 1.0D-06
CTRSH = 5.0D-02
IPRINT = 5
MAXIT = 20
IREST = 0
!PAM97 IRHP = 0
ICPF = 0
ISDCI = 0
INCPF = 0
ICONV = 0
MAXITP = 6
WLEV = 0.3d0
ETOT = 0.0d0
do I=1,8
  NPFRO(I) = 0
  NFRO(I) = 0
  NDEL(I) = 0
  NPDEL(I) = 0
  NISH(I) = 0
  NASH(I) = 0
  NVAL(I) = 0
  NVIR(I) = 0
  NORB(I) = 0
  NBAS(I) = 0
end do
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
do I=1,MXORB+1
  IROW(I) = I*(I-1)/2
end do
do I=1,99
  LW(I) = 0
end do
nTit = 0

!---  read the header of TRAONE ---------------------------------------*
! Note: NORB(i)=NBAS(i)-NPFRO(i)-NPDEL(i)
NAMSIZ = LENIN8*MXORB
IDISK = 0
call WR_MOTRA_Info(Lu_TraOne,2,iDisk,ITOC17,64,POTNUC,NSYM,NBAS,NORB,NPFRO,NPDEL,8,NAME,NAMSIZ)

!---  Read input from standard input ----------------------------------*
call RdNLst(5,'CPF')
10 read(5,'(A)',end=991) Line
Command = Line(1:4)
call UpCase(Command)
if (Command(1:1) == '*') goto 10
jCmd = 0
do iCmd=1,nCmd
  if (Command == Cmd(iCmd)) jCmd = iCmd
end do
20 goto(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800) jCmd
write(6,*) 'READIN Error: Command not recognized.'
write(6,*) 'The command is:'//''''//Command//''''
call QUIT_OnUserError()

!---  process TITL command --------------------------------------------*
100 continue
read(5,'(A)',end=991) Line
Command = Line(1:4)
call UpCase(Command)
if (Command(1:1) == '*') goto 100
jCmd = 0
do iCmd=1,nCmd
  if (Command == Cmd(iCmd)) jCmd = iCmd
end do
if (jCmd /= 0) goto 20
nTit = nTit+1
if (nTit <= mxTit) Title(nTit) = Line
goto 100

!---  process MAXP command --------------------------------------------*
200 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 200
read(Line,*,Err=992) MaxItP
goto 10

!---  process LEVS command --------------------------------------------*
300 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 300
read(Line,*,Err=992) WLev
goto 10

!---  process THRP command --------------------------------------------*
400 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 400
read(Line,*,Err=992) CTrsh
goto 10

!---  process PRIN command --------------------------------------------*
500 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 500
read(Line,*,Err=992) iPrint
goto 10

!---  process FROZ command --------------------------------------------*
600 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 600
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (nFro(i),i=1,8)
goto 10

!---  process DELE command --------------------------------------------*
700 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 700
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (NDEL(i),i=1,8)
goto 10

!---  process MAXI command --------------------------------------------*
800 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 800
read(Line,*,Err=992) MaxIt
MaxIt = min(MaxIt,75)
goto 10

!---  process ECON command --------------------------------------------*
900 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 900
read(Line,*,Err=992) EThre
goto 10

!---  process ETRS command --------------------------------------------*
1000 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 1000
!PAM97 read(Line,*,Err=992) ETrsh
read(Line,*)
write(6,*) ' WARNING: The obsolete ETRS command is ignored.'
goto 10

!---  process REST command --------------------------------------------*
1100 continue
iRest = 1
goto 10

!---  process MCPF command --------------------------------------------*
1200 continue
iCPF = 0
iSDCI = 0
iNCPF = 0
goto 10

!---  process CPF  command --------------------------------------------*
1300 continue
iCPF = 1
iSDCI = 0
iNCPF = 0
goto 10

!---  process SDCI command --------------------------------------------*
1400 continue
iSDCI = 1
iCPF = 0
iNCPF = 0
goto 10

!---  process ACPF command --------------------------------------------*
1500 continue
iNCPF = 1
iCPF = 0
iSDCI = 0
goto 10

!---  process LOW  command --------------------------------------------*
1600 continue
LWSP = .true.
goto 10

!---  process EXTR command --------------------------------------------*
1700 write(6,*) 'The EXTRACT option is redundant and is ignored!'
goto 10

!---  The end of the input is reached, print the title ----------------*
1800 continue
if (ntit == 0) then
  ntit = 1
  title(1) = ' ( No title was given )'
end if
write(6,*)
write(6,'(6X,120A1)') ('*',i=1,120)
write(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
write(6,'(6X,57A1,A6,57A1)') '*',(' ',i=1,56),'Title:',(' ',i=1,56),'*'
do i=1,nTit
  call Center_Text(Title(i))
  write(6,'(6X,24A1,A72,24A1)') '*',(' ',j=1,23),Title(i),(' ',j=1,23),'*'
end do
write(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
write(6,'(6X,120A1)') ('*',i=1,120)
write(6,*)

!---  print the coordinates of the system -----------------------------*
call PrCoor()

!---  print the method used -------------------------------------------*
write(6,*)
if (iSDCI == 1) then
  write(6,'(6X,A)') 'This is an  S D C I  calculation'
else if (iCPF == 1) then
  write(6,'(6X,A)') 'This is a  C P F  calculation'
else if (INCPF == 1) then
  write(6,'(6X,A)') 'This is an  A C P F calculation'
else
  write(6,'(6X,A)') 'This is an  M C P F  calculation'
end if
if (LWSP) write(6,'(6X,A)') 'This is a LOW SPIN calculation'
call XFLUSH(6)

!---  read the header of CIGUGA ---------------------------------------*
IADD10 = 0
call iDAFILE(Lu_CIGuga,2,IAD10,9,IADD10)
iOpt = 2
nJJS = 18
nIRC = 4
call WR_GUGA(Lu_CIGuga,iOpt,IADD10,NFREF,S,N,LN,NSYM,IR1,IRJ,IFIRST,INTNUM,LSYM,NREF,LN1,NRLN1,NASH,NISH,8,IRC,nIRC,JJS,nJJS,NVAL, &
             IOCR,nIOCR)
if (LN >= MXORB) then
  write(6,*) 'READIN Error: Too many orbitals.'
  write(6,'(1X,A,2I5)') 'LN,MXORB:',LN,MXORB
  call QUIT_OnUserError()
end if
LW(1) = 1
call iDAFILE(Lu_CIGuga,2,iH(iPointer(LW(1))),IR1,IADD10)
LW(2) = LW(1)+(IR1+(RTOI-1))/RTOI
call iDAFILE(Lu_CIGuga,2,iH(iPointer(LW(2))),IRJ,IADD10)

!---  update orbital specifications -----------------------------------*
IV0 = 0
IV1 = 1
IV2 = 2
IV3 = 3
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
IN = 0
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
    IN = IN+1
    IR = IR-1
    ICH(IN) = IR
  end do
  do J=1,NISHI
    IN = IN+1
    IT = IT+1
    ICH(IN) = IT
    NSM(IT) = I
  end do
  do J=1,NASHI
    IN = IN+1
    IU = IU+1
    ICH(IN) = IU
    NSM(IU) = I
  end do
  do J=1,NVALI
    IN = IN+1
    IVA = IVA+1
    ICH(IN) = IVA
    NSM(IVA) = I
  end do
  do J=1,NVIRI
    IN = IN+1
    IV = IV+1
    ICH(IN) = IV
    NSM(IV) = I
  end do
  do J=1,NDELI
    IN = IN+1
    ICH(IN) = 0
  end do
end do
NVT = IROW(NVIRT+1)
NVT2 = IROW(NVIRT)

!---  report input specifications -------------------------------------*
write(6,*)
write(6,'(6X,A)') 'ONE-ELECTRON BASIS:'
write(6,'(6X,A)') '----------------------------'
write(6,*)
write(6,'(6X,A,T47,4X,4X,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
write(6,*)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals pre-frozen in MOTRA',nPFroT,(nPFro(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals used by this program',norbT,(nOrb(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Pre-deleted in MOTRA',nPDelT,(nPDel(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Sum: No. of basis functions',nBasT,(nBas(iSym),iSym=1,nSym)
call XFLUSH(6)
write(6,*)
write(6,'(6X,A)') 'ORBITAL SPECIFICATION:'
write(6,'(6X,A)') '-------------------------------'
write(6,*)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals frozen here',nFroT,(nFro(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Inactive orbitals',nIShT,(nISh(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Active orbitals',nAShT,(nASh(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Additional valence orbitals',nValT,(nVal(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Virtual orbitals',nVirT,(nVir(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals deleted here',nDelT,(nDel(iSym),iSym=1,nSym)
write(6,'(6X,A,T47,I4,4x,8I4)') 'Sum: Total no. of orbitals',norbT,(nOrb(iSym),iSym=1,nSym)
write(6,*)
call XFLUSH(6)
write(6,*)
write(6,'(6X,A)') 'WAVE FUNCTION SPECIFICATION:'
write(6,'(6X,A)') '----------------------------'
write(6,*)
write(6,'(6X,A,T47,I4)') 'Number of electrons in CI',N
write(6,'(6X,A,T47,I4)') 'Internal orbitals in CI',LN
write(6,'(6X,A,T47,I4)') 'External orbitals in CI',NVIRT
write(6,'(6X,A,T47,I4)') 'Number of irreps',NSYM
write(6,'(6X,A,T47,F4.1)') 'Spin quantum number',S
write(6,'(6X,A,T47,I4)') 'State symmetry',LSYM
call XFLUSH(6)
write(6,*)
write(6,'(6X,A)') 'REFERENCE STATE:'
write(6,'(6X,A)') '------------------------------'
write(6,*)
write(6,'(6X,A,T47,I4)') 'Number of reference states',NREF
LN2 = min(16,LN1)
if (LN1 == 0) then
  write(6,'(6X,A,T47)') 'One closed shell reference state'
else
  write(6,'(6X,A,T47)') 'Occupation of active orbitals in the reference state:'
  write(6,'(6X,A,T25,16I4)') 'Active orbital nr.',(I,I=1,LN2)
  jEnd = 0
  do iRef=1,nRef
    jStart = jEnd+1
    jEnd = jEnd+LN1
    write(6,'(6X,A,I3,T25,16I4)') 'Ref nr',IREF,(IOCR(j),j=jStart,jEnd)
  end do
end if
write(6,*)
call XFLUSH(6)
write(6,'(6X,A)') 'OPTIONS:'
write(6,'(6X,A)') '--------'
write(6,*)
write(6,'(6X,A,T47,I4)') 'Print parameter',iPrint
write(6,'(6X,A,T47)') 'Pulay diagonalization'
if (INTNUM /= 0) write(6,'(6X,A,T47)') 'First order interacting space'
!PAM97 if (IRHP /= 0) write(6,'(6X,A,T47)') 'Root homing'
IX1 = IRC(1)
IX2 = IRC(2)-IRC(1)
ISC(1) = IX1
ISC(2) = ISC(1)+IX2*NVIRT
IY1 = ISC(1)
IY2 = ISC(2)-ISC(1)
write(6,214)
call XFLUSH(6)
214 format(//,6X,'INTERNAL CONFIGURATIONS')
if (IFIRST /= 0) GO TO 205
IX3 = IRC(3)-IRC(2)
IX4 = IRC(4)-IRC(3)
ISC(3) = ISC(2)+IX3*NVT2
ISC(4) = ISC(3)+IX4*NVT
IY3 = ISC(3)-ISC(2)
IY4 = ISC(4)-ISC(3)
write(6,215) IX1,IX2,IX3,IX4
call XFLUSH(6)
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7, &
           /,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
write(6,213)
call XFLUSH(6)
213 format(//,6X,'FULL-SPACE CONFIGURATIONS (FORMAL)')
write(6,215) IY1,IY2,IY3,IY4
call XFLUSH(6)
GO TO 206
205 write(6,216) IX1,IX2
call XFLUSH(6)
216 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
write(6,213)
call XFLUSH(6)
write(6,216) IY1,IY2
call XFLUSH(6)
206 ILIM = 4
if (IFIRST /= 0) ILIM = 2
! ERROR CONDITIONS:
!if (LN /= NISHT+NASHT+NVALT) then
!  write(6,*) ' ERROR: Orbital specifications do not match'
!  write(6,*) ' input to GUGA. The number of internal'
!  write(6,*) ' orbitals must equal the number of inactive,'
!  write(6,*) ' active, and additional valence orbitals.'
!  call QUIT(20)
!end if
! ALLOCATION FOR INDEX VECTORS
! THESE VECTORS ARE PERMANENTLY IN CORE
! -- INDEX
!PAM97 LW(3) = LW(2)+IRJ
LW(3) = LW(2)+(IRJ+(RTOI-1))/RTOI
! ISAB
LW(4) = LW(3)+IRC(ILIM)
NVIR2 = NVIRT*NVIRT
! JREFX
LW(5) = LW(4)+NVIR2
IADD10 = IAD10(2)
call iDAFILE(Lu_CIGuga,2,iH(iPointer(LW(5))),ISC(1),IADD10)
! ADDRESSES NOT USED
LW(6) = LW(5)+ISC(1)
LW(7) = LW(6)
LW(8) = LW(7)
LW(9) = LW(8)
LW(10) = LW(9)
! LIMIT FOR PERMANENT VECTORS
LPERMA = LW(10)
call dINDMAT(H)
call ALLOC_CPF(ISMAX,LPERMA)

return

991 continue
write(6,*) 'READIN Error: Premature end of file while reading.'
call Quit_OnUserError()
992 continue
write(6,*) 'READIN Error: I/O error during internal read.'
write(6,*) 'The line that could not be read is:'
write(6,*) Line
call Quit_OnUserError()

! This is to allow type punning without an explicit interface
contains
subroutine dINDMAT(H)
  real*8, target :: H(*)
  integer, pointer :: iH2(:), iH3(:), iH4(:), iH5(:)
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(4))),iH4,[1])
  call c_f_pointer(c_loc(H(LW(5))),iH5,[1])
  call INDMAT_CPF(iH2,iH3,iH4,ISMAX,iH5)
  nullify(iH2,iH3,iH4,iH5)
end subroutine dINDMAT

end subroutine ReadIn_CPF
