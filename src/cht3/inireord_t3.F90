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
subroutine IniReord_t3(NaGrp,wrksize)
! nacitanie vsupu a inicializacia premnennych
! a tlac primitivnej hlavicky pre Reord procesz

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs
#endif
implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "cholesky.fh"
#include "ccsd_t3compat.fh"
integer NaGrp
integer wrksize
integer nOrb(8), nOcc(8)
integer ndelvirt
integer LuSpool
character*80 LINE
integer rc
real*8 FracMem
character*3 msg
#ifdef _MOLCAS_MPP_
integer jal1, jal2
#endif

! setup defaults

call Get_iArray('nOrb',nOrb,1)
call Get_iArray('nIsh',nOcc,1)

no = nOcc(1)
nv = nOrb(1)-nOcc(1)

FracMem = 0.0d0
call Cho_X_init(rc,FracMem) ! initialize cholesky info

! take local # of Cholesky Vectors on this node
#ifdef _MOLCAS_MPP_

do jal1=0,Nprocs-1
  NChLoc(jal1) = 0
end do

NChLoc(MyRank) = NumCho(1)

call gaigop(NChLoc(0),NProcs,'+')

jal2 = 0
do jal1=0,NProcs-1
  jal2 = jal2+NChLoc(jal1)
end do

nc = jal2
#else
nc = NumCho(1)
#endif

call Cho_X_final(rc)

ndelvirt = 0
LunAux = 13
mhkey = 1
generkey = 1
!mp !NaGrp = 1
call get_iScalar('CHCCLarge',NaGrp)
restkey = 0
printkey = 1

! t3 specific keywords

gen_files = .true.
run_triples = .true.
t3_starta = -1
t3_stopa = -1
t3_startb = -1
t3_stopb = -1

!mp read input file

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
5 read(LuSpool,'(A80)') LINE
call UPCASE(LINE)
if (index(LINE,'&CHT3') == 0) goto 5
6 read(LuSpool,'(A80)') LINE
if (LINE(1:1) == '*') goto 6
call UPCASE(LINE)
!
if (LINE(1:4) == 'TITL') then
  read(LuSpool,*)

else if (LINE(1:4) == 'FROZ') then ! FROZen
  read(LuSpool,*) nfr
  if ((nfr < 0) .or. (nfr >= no)) then
    write(6,*)
    write(6,*) 'Ilegal value for FROZen keyword : ',nfr
    call abend()
  end if
  no = no-nfr

else if (LINE(1:4) == 'DELE') then ! DELEted
  read(LuSpool,*) ndelvirt
  if ((ndelvirt < 0) .or. (ndelvirt >= nv)) then
    write(6,*)
    write(6,*) 'Ilegal value for DELEted keyword : ',ndelvirt
    call abend()
  end if
  nv = nv-ndelvirt

  !mp !else if (LINE(1:4) == 'LARG') then ! LARGegroup
  !mp !  read(LuSpool,*) NaGrp
  !mp !  if ((NaGrp < 1) .or. (NaGrp > 32)) then
  !mp !    write(6,*)
  !mp !    write(6,*) 'Ilegal value for LARGegroup keyword : ',NaGrp
  !mp !    write(6,*) 'Large segmentation must be <= 32'
  !mp !    call abend()
  !mp !  end if

  !mp !else if (LINE(1:4) == 'LUNA') then  !... toto sa nikdy nevyuzivalo
  !mp !  read(LuSpool,*) LunAux

else if (LINE(1:4) == 'MHKE') then ! MHKEy
  read(LuSpool,*) mhkey
  if ((mhkey < 0) .or. (mhkey > 2)) then
    mhkey = 1
    write(6,*)
    write(6,*) ' Warning!!! ',' MHKEy out of range, changed to 1'
  end if

else if (LINE(1:4) == 'REST') then ! RESTart
  restkey = 1
  write(6,*)
  write(6,*) 'RESTart option is temporary disabled'
  write(6,*) 'No Restart possible (... yet).'
  call abend()

else if (LINE(1:4) == 'PRIN') then ! PRINtkey
  read(LuSpool,*) printkey
  if (((printkey < 0) .or. (printkey > 10)) .or. ((printkey > 2) .and. (printkey < 10))) then

    write(6,*)
    write(6,*) 'Ilegal value of the PRINtkey keyword: ',printkey
    write(6,*) ' Use: 1  (Minimal) '
    write(6,*) '      2  (Minimal + Timings)'
    write(6,*) '      10 (Debug) '
    call abend()
  end if

else if (LINE(1:4) == 'NOGE') then ! NOGEnerate
  gen_files = .false.

else if (LINE(1:4) == 'NOTR') then ! NOTRiples
  run_triples = .false.

else if (LINE(1:4) == 'ALOO') then ! ALOOp
  read(LuSpool,*) t3_starta,t3_stopa
  if ((t3_starta < -1) .or. (t3_stopa < -1)) then
    write(6,*) 'ALOOp values can be either: '
    write(6,*) '-1 : indicating normal run, or'
    write(6,*) 'positive numbers!'
    call abend()
  end if

else if (LINE(1:4) == 'BLOO') then ! BLOOp
  read(LuSpool,*) t3_startb,t3_stopb
  if ((t3_startb < -1) .or. (t3_stopb < -1)) then
    write(6,*) 'BLOOp values can be either: '
    write(6,*) '-1 : indicating normal run, or'
    write(6,*) 'positive numbers!'
    call abend()
  end if

else if (LINE(1:4) == 'END ') then
  goto 7
end if
goto 6
7 continue

call Close_LuSpool(LuSpool)

!! take care of the cholesky vectors segmentation
!! to lead to < 100 blocks

!mp checks
if (t3_starta > t3_stopa) then
  write(6,*) 'Mismatch in input : '
  write(6,*) 'T3_STARTA = ',t3_starta
  write(6,*) 'T3_STOPA = ',t3_stopa
  call abend()
end if

if (t3_startb > t3_stopb) then
  write(6,*) 'Mismatch in input : '
  write(6,*) 'T3_STARTB = ',t3_startb
  write(6,*) 'T3_STOPB = ',t3_stopb
  call abend()
end if

if ((t3_starta < 0) .and. (t3_stopa > 0)) then
  write(6,*) 'Mismatch in input : '
  write(6,*) 'T3_STARTA = ',t3_starta
  write(6,*) 'T3_STOPA = ',t3_stopa
  call abend()
end if

if ((t3_startb < 0) .and. (t3_stopb > 0)) then
  write(6,*) 'Mismatch in input : '
  write(6,*) 'T3_STARTB = ',t3_startb
  write(6,*) 'T3_STOPB = ',t3_stopb
  call abend()
end if

!mp !if ((t3_starta < 0) .and. (t3_startb < 0)) then
!mp !  write(6,*) 'This restart combination not implemented'
!mp !  write(6,*) 'T3_STARTA = ',t3_starta
!mp !  write(6,*) 'T3_STARTB = ',t3_startb
!mp !  call abend()
!mp !end if

!2 tlac hlavicky
write(6,*)
write(6,*) '    Cholesky Based Closed-Shell (T) code'
write(6,*)
write(6,*) '--------------------------------------------------'

write(6,'(A,i9)') ' Frozen Orbitals                   : ',nfr
write(6,'(A,i9)') ' Occupied Orbitals                 : ',no
write(6,'(A,i9)') ' Virtual Orbitals                  : ',nv
write(6,'(A,i9)') ' Total number of Cholesky Vectors  : ',nc

write(6,*) '--------------------------------------------------'

write(6,'(A,i9)') ' Large Virtual Segmentation        : ',NaGrp

write(6,*) '--------------------------------------------------'

msg = 'No'
if (gen_files) msg = 'Yes'

write(6,'(A,A5)') ' Generate Triples Scratch Files?        : ',msg

msg = 'No'
if (.not. run_triples) msg = 'Yes'

write(6,'(A,A5)') ' Stop after Scratch Files generation?   : ',msg

write(6,*) '--------------------------------------------------'

if (t3_starta == -1) then
  write(6,'(A,i4)') ' Calculating full loop A                '
else
  write(6,'(A,i4)') ' VO index triplet to start with in loop A : ',t3_starta
  write(6,'(A,i4)') ' VO index triplet to stop  at   in loop A : ',t3_stopa
end if

if (t3_starta == -1) then
  write(6,'(A,i4)') ' Calculating full loop B                '
else
  write(6,'(A,i4)') ' VO index triplet to start with in loop B : ',t3_startb
  write(6,'(A,i4)') ' VO index triplet to stop  at   in loop B : ',t3_stopb
end if

write(6,*) '--------------------------------------------------'

write(6,'(A,i9)') ' Lun Number for Aux. Matrixes      : ',LunAux
write(6,'(A,i9)') ' BLAS/FTN Matrix Handling          : ',mhkey

msg = 'No'
if (restkey == 1) msg = 'Yes'

write(6,'(A,A10)') ' Start from RstFil ?               : ',msg
write(6,'(A,i9)') ' Print level                       : ',printkey

write(6,*) '--------------------------------------------------'
write(6,*)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(wrksize)

end subroutine IniReord_t3
