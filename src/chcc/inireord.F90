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

subroutine IniReord(NaGrp,NaSGrp,NchBlk,LunAux,wrksize)
! nacitanie vsupu a inicializacia premnennych
! a tlac primitivnej hlavicky pre Reord procesz

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs
#endif
implicit none
#include "chcc1.fh"
#include "chcc_reord.fh"
!mp
#include "parcc.fh"
!mp

integer NaGrp, NaSGrp, NchBlk
integer LunAux, wrksize
!mp
integer nOrb(8), nOcc(8), nFro(8), nDel(8)
integer intkey1, intkey2
integer ndelvirt

integer LuSpool
character*80 LINE

#ifdef _MOLCAS_MPP_
integer jal1
#endif
integer NChLoc_min, NChLoc_max, NchBlk_tmp

character*3 msg
!mp

! setup defaults

call Get_iArray('nBas',nOrb,1) ! must read always nBas!!
call Get_iArray('nIsh',nOcc,1)
call Get_iArray('nFroPT',nFro,1) ! = 'nFro' in previous step
call Get_iArray('nDelPT',nDel,1) ! = 'nDel' in previous step

#ifdef _MOLCAS_MPP_
!mp get min/max
NChLoc_min = nc
NChLoc_max = 0
do jal1=0,Nprocs-1
  if (NChLoc(jal1) <= NChLoc_min) NChLoc_min = NChLoc(jal1)
  if (NChLoc(jal1) >= NChLoc_max) NChLoc_max = NChLoc(jal1)
end do

!mp calc reasonable starting value (200-300)
if (NChLoc_min /= NChLoc_max) then
  NChBlk = int(NChLoc_min/2)
else
  NChBlk = NChLoc_min
end if

if (NChBlk >= 300) NChBlk = min(200,int(NChBlk/2))

!mp fix num of ChV blocks to be less then 100
if (int(NChLoc_max/NChBlk) >= 100) NChBlk = int(NChLoc_max/100)-1
#else
NChLoc_max = nc
!mp calc reasonable starting value (200-300)
if (nc >= 300) then
  NChBlk = min(200,nc/2)
else
  NChBlk = nc
end if

NChLoc_min = NChBlk
!mp fix num of ChV blocks to be less then 100
if (int(nc/NChBlk) >= 100) NChBlk = int(nc/100)-1
#endif

nfr = nFro(1)
no = nOcc(1)-nFro(1)
ndelvirt = nDel(1)

nv = nOrb(1)-nDel(1)-nOcc(1) ! nOrb defined as = to nBas, right!

LunAux = 13
mhkey = 1
generkey = 1
intkey1 = 0
intkey2 = 0

NaGrp = 0
NaSGrp = 0
W34DistKey = 1
JoinLkey = 2 ! toto este nemam domyslene
restkey = 0
conv = 1.0d-6
printkey = 1
maxiter = 40

!mp read input file

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
5 continue
read(LuSpool,'(A80)') LINE
call UPCASE(LINE)
if (index(LINE,'&CHCC') == 0) goto 5
6 continue
read(LuSpool,'(A80)') LINE
if (LINE(1:1) == '*') goto 6
call UPCASE(LINE)

if (LINE(1:4) == 'TITL') then
  read(LuSpool,*)

else if (LINE(1:4) == 'FROZ') then
  read(LuSpool,*) nfr
  if ((nfr < 0) .or. (nfr >= no)) then
    write(6,*)
    write(6,*) 'Ilegal value for FROZen keyword : ',nfr
    call abend()
  end if
  no = no+nFro(1)-nfr

else if (LINE(1:4) == 'DELE') then
  read(LuSpool,*) ndelvirt
  if ((ndelvirt < 0) .or. (ndelvirt > nv)) then
    write(6,*)
    write(6,*) 'Ilegal value for DELETED keyword : ',ndelvirt
    call abend()
  end if
  nv = nv+nDel(1)-ndelvirt

else if (LINE(1:4) == 'LARG') then
  read(LuSpool,*) NaGrp
  if ((NaGrp < 0) .or. (NaGrp > maxGrp)) then
    write(6,*)
    write(6,*) 'Ilegal value for LARGE keyword : ',NaGrp
    write(6,*) 'Large segmentation must be -le 32'
    call abend()
  end if

else if (LINE(1:4) == 'SMAL') then
  read(LuSpool,*) NaSGrp
  if ((NaSGrp < 0) .or. (NaSGrp > 8)) then
    write(6,*)
    write(6,*) 'Ilegal value for SMALL keyword : ',NaSGrp
    write(6,*) 'Small segmentation must be -le 8'
    call abend()
  end if

  ! large == 0, small != 0 => quit
  if ((NaGrp == 0) .and. (NaSGrp /= 0)) then
    write(6,*)
    write(6,*) 'Small segmentation must be specified'
    write(6,*) 'with large segmentation, or both can'
    write(6,*) 'be left unspecified'
    call abend()
  end if

  if (NaGrp /= 0) then
    ! large != 0, small == 0 => small = 1
    if (NaSGrp == 0) NaSGrp = 1

    ! large * small <= 64
    if ((NaGrp*NaSGrp) > maxSGrp) then
      write(6,*)
      write(6,*) 'Product of Large and Small segmen-'
      write(6,*) 'tation must be less or equal to 64'
      call abend()
    end if
  end if

else if (LINE(1:4) == 'CHSE') then
  read(LuSpool,*) NchBlk_tmp
  if ((NchBlk_tmp < 1) .or. (NchBlk_tmp > NChLoc_min)) then
    write(6,*)
    write(6,*) 'Ilegal value for CHSegment keyword  : ',NchBlk_tmp
    write(6,*) 'Reseting to a reasonable value for    '
    write(6,*) 'this system :                         ',NchBlk
  else if (int(NChLoc_max/NchBlk_tmp) >= 100) then
    write(6,*) 'Number of block of the MO Cholesky vector'
    write(6,*) 'exceeded the limit. Increasing value of  '
    write(6,*) 'the CHSEgmentation keyword to : ',NchBlk
  else
    NchBlk = NchBlk_tmp
  end if

!mp else if (LINE(1:4) == 'LUNA') then  ... toto sa nikdy nevyuzivalo
!mp   read(LuSpool,*) LunAux

else if (LINE(1:4) == 'MHKE') then
  read(LuSpool,*) mhkey
  if ((mhkey < 0) .or. (mhkey > 2)) then
    mhkey = 1
    write(6,*)
    write(6,*) ' Warning!!!  Matrix handling key out of range'
    write(6,*) ' parameter mhkey changed to 1'
  end if

else if (LINE(1:4) == 'NOGE') then
  generkey = 0

else if (LINE(1:4) == 'ONTH') then
  intkey1 = 1

else if (LINE(1:4) == 'PREC') then
  intkey2 = 1

else if (LINE(1:4) == 'NODI') then
  W34DistKey = 0

else if (LINE(1:4) == 'JOIN') then
  read(LuSpool,*) JoinLkey
  if ((JoinLkey < 0) .or. (JoinLkey > 3)) then
    write(6,*)
    write(6,*) 'Ilegal value for Join keyword : ',JoinLkey
    write(6,*) 'Use one of 0, 1, 2, 3'
    write(6,*) 'For details, see the manual ...'
    call abend()
  end if

else if (LINE(1:4) == 'MAXI') then
  read(LuSpool,*) maxiter
  if (maxiter <= 0) then
    write(6,*)
    write(6,*) 'Ilegal value of the MAXITER keyword: ',maxiter
    write(6,*) 'Use integer > 0'
    call abend()
  end if

else if (LINE(1:4) == 'REST') then
  restkey = 1
  write(6,*)
  write(6,*) 'This option is temporary disabled'
  write(6,*) 'No Restart possible (... yet).'
  call abend()

else if (LINE(1:4) == 'THRE') then
  read(LuSpool,*) conv

else if (LINE(1:4) == 'PRIN') then
  read(LuSpool,*) printkey
  if (((printkey < 0) .or. (printkey > 10)) .or. ((printkey > 2) .and. (printkey < 10))) then

    write(6,*)
    write(6,*) 'Ilegal value of the PRINT keyword: ',printkey
    write(6,*) ' Use: 1  (Minimal) '
    write(6,*) '      2  (Minimal + Timings)'
    write(6,*) '      10 (Debug) '
    call abend()
  end if

else if (LINE(1:4) == 'END ') then
  goto 7
end if
goto 6
7 continue

call Close_LuSpool(LuSpool)

! take care of the algorithm keyword
if (intkey1 == intkey2) then
  if (intkey1 == 0) then
    write(6,*)
    write(6,*) 'None of OnTheFly/PreCalculate'
    write(6,*) 'algorithm was selected. Using'
    write(6,*) 'default: PreCalculate (1)'
    intkey = 1
  else
    write(6,*)
    write(6,*) 'OnTheFly and PreCalculate keywords'
    write(6,*) 'are mutually exclusive'
    call abend()
  end if
else
  if (intkey1 == 1) then
    intkey = 0
  else
    intkey = 1
  end if
end if

!2 tlac hlavicky
write(6,*)
write(6,*) '    Cholesky Based Closed-Shell CCSD code'
!mp write(6,*) ' Dedicated to the memory of Boris Jeltzin'
write(6,*)
write(6,*) '---------------------------------------------------'

write(6,'(A,i9)') ' Frozen Orbitals                   : ',nfr
write(6,'(A,i9)') ' Occupied Orbitals                 : ',no
write(6,'(A,i9)') ' Virtual Orbitals                  : ',nv
write(6,'(A,i9)') ' Total number of Cholesky Vectors  : ',nc

write(6,*) '---------------------------------------------------'

if (NaGrp /= 0) then
  write(6,'(A,i9)') ' Large Virtual Segmentation        : ',NaGrp
else
  write(6,'(A,A9)') ' Large Virtual Segmentation        :  auto'
end if

if (NaSGrp /= 0) then
  write(6,'(A,i9)') ' Small Virtual Segmentation        : ',NaSGrp
else
  write(6,'(A,A9)') ' Small Vectors Segmentation        :  auto'
end if

write(6,'(A,i9)') ' Cholesky Vectors Segmentation     : ',NchBlk

write(6,*) '---------------------------------------------------'

msg = 'No'
if (generkey == 1) msg = 'Yes'

write(6,'(A,A4)') ' Generate Scratch Files?                : ',msg
write(6,'(A,i4)') ' Precalculate (1) / On-the-Fly (0) Alg. : ',intkey
write(6,'(A,i4)') ' 3 and 4-ext. MO integrals distribute?  : ',W34DistKey
write(6,'(A,i4)') ' Parallel Join of varios MO integrals   : ',JoinLkey

write(6,*) '---------------------------------------------------'

write(6,'(A,E9.2)') ' Convergence Threshold             : ',conv
write(6,'(A,i9)') ' Maximum number of Iterations      : ',maxiter

write(6,*) '---------------------------------------------------'

write(6,'(A,i9)') ' Lun Number for Aux. Matrixes      : ',LunAux
write(6,'(A,i9)') ' BLAS/FTN Matrix Handling          : ',mhkey

msg = 'No'
if (restkey == 1) msg = 'Yes'

write(6,'(A,A10)') ' Start from RstFil ?               : ',msg
write(6,'(A,i9)') ' Print level                       : ',printkey

write(6,*) '---------------------------------------------------'
write(6,*)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(wrksize)

end subroutine IniReord
