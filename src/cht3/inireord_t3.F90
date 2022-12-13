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

subroutine IniReord_t3(NaGrp)
! nacitanie vsupu a inicializacia premnennych
! a tlac primitivnej hlavicky pre Reord procesz

use ChT3_global, only: gen_files, LunAux, nc, nfr, no, nv, printkey, run_triples, t3_starta, t3_startb, t3_stopa, t3_stopb
#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: NaGrp
#include "cholesky.fh"
integer(kind=iwp) :: LuSpool, ndelvirt, nOcc(8), nOrb(8), rc
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: jal1, jal2
integer(kind=iwp), allocatable :: NChLoc(:)
#endif
real(kind=wp) :: FracMem
character(len=80) :: LINE
character(len=3) :: msg

! setup defaults

call Get_iArray('nOrb',nOrb,1)
call Get_iArray('nIsh',nOcc,1)

no = nOcc(1)
nv = nOrb(1)-nOcc(1)

FracMem = Zero
call Cho_X_init(rc,FracMem) ! initialize cholesky info

! take local # of Cholesky Vectors on this node
#ifdef _MOLCAS_MPP_

call mma_allocate(NChLoc,NProcs,label='NChLoc')
NChLoc(:) = 0

NChLoc(MyRank+1) = NumCho(1)

call gaigop(NChLoc,NProcs,'+')

jal2 = 0
do jal1=1,NProcs
  jal2 = jal2+NChLoc(jal1)
end do

call mma_deallocate(NChLoc)

nc = jal2
#else
nc = NumCho(1)
#endif

call Cho_X_final(rc)

ndelvirt = 0
LunAux = 13
!mhkey = 1
!generkey = 1
!mp !NaGrp = 1
call get_iScalar('CHCCLarge',NaGrp)
!restkey = 0
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
do
  read(LuSpool,'(A80)') LINE
  call UPCASE(LINE)
  if (index(LINE,'&CHT3') /= 0) exit
end do
do
  read(LuSpool,'(A80)') LINE
  if (LINE(1:1) == '*') cycle
  call UPCASE(LINE)

  select case (LINE(1:4))

    case ('TITL')
      read(LuSpool,*)

    case ('FROZ') ! FROZen
      read(LuSpool,*) nfr
      if ((nfr < 0) .or. (nfr >= no)) then
        write(u6,*)
        write(u6,*) 'Ilegal value for FROZen keyword : ',nfr
        call abend()
      end if
      no = no-nfr

    case ('DELE') ! DELEted
      read(LuSpool,*) ndelvirt
      if ((ndelvirt < 0) .or. (ndelvirt >= nv)) then
        write(u6,*)
        write(u6,*) 'Ilegal value for DELEted keyword : ',ndelvirt
        call abend()
      end if
      nv = nv-ndelvirt

    !mp !case ('LARG') ! LARGegroup
    !mp !  read(LuSpool,*) NaGrp
    !mp !  if ((NaGrp < 1) .or. (NaGrp > 32)) then
    !mp !    write(u6,*)
    !mp !    write(u6,*) 'Ilegal value for LARGegroup keyword : ',NaGrp
    !mp !    write(u6,*) 'Large segmentation must be <= 32'
    !mp !    call abend()
    !mp !  end if

    !mp !case ('LUNA') !... toto sa nikdy nevyuzivalo
    !mp !  read(LuSpool,*) LunAux

    !case ('MHKE') ! MHKEy
    !  read(LuSpool,*) mhkey
    !  if ((mhkey < 0) .or. (mhkey > 2)) then
    !    mhkey = 1
    !    write(u6,*)
    !    write(u6,*) ' Warning!!! ',' MHKEy out of range, changed to 1'
    !  end if

    case ('REST') ! RESTart
      !restkey = 1
      write(u6,*)
      write(u6,*) 'RESTart option is temporary disabled'
      write(u6,*) 'No Restart possible (... yet).'
      call abend()

    case ('PRIN') ! PRINtkey
      read(LuSpool,*) printkey
      if (((printkey < 0) .or. (printkey > 10)) .or. ((printkey > 2) .and. (printkey < 10))) then

        write(u6,*)
        write(u6,*) 'Ilegal value of the PRINtkey keyword: ',printkey
        write(u6,*) ' Use: 1  (Minimal) '
        write(u6,*) '      2  (Minimal + Timings)'
        write(u6,*) '      10 (Debug) '
        call abend()
      end if

    case ('NOGE') ! NOGEnerate
      gen_files = .false.

    case ('NOTR') ! NOTRiples
      run_triples = .false.

    case ('ALOO') ! ALOOp
      read(LuSpool,*) t3_starta,t3_stopa
      if ((t3_starta < -1) .or. (t3_stopa < -1)) then
        write(u6,*) 'ALOOp values can be either: '
        write(u6,*) '-1 : indicating normal run, or'
        write(u6,*) 'positive numbers!'
        call abend()
      end if

    case ('BLOO') ! BLOOp
      read(LuSpool,*) t3_startb,t3_stopb
      if ((t3_startb < -1) .or. (t3_stopb < -1)) then
        write(u6,*) 'BLOOp values can be either: '
        write(u6,*) '-1 : indicating normal run, or'
        write(u6,*) 'positive numbers!'
        call abend()
      end if

    case ('END ')
      exit

  end select
end do

call Close_LuSpool(LuSpool)

!! take care of the cholesky vectors segmentation
!! to lead to < 100 blocks

!mp checks
if (t3_starta > t3_stopa) then
  write(u6,*) 'Mismatch in input : '
  write(u6,*) 'T3_STARTA = ',t3_starta
  write(u6,*) 'T3_STOPA = ',t3_stopa
  call abend()
end if

if (t3_startb > t3_stopb) then
  write(u6,*) 'Mismatch in input : '
  write(u6,*) 'T3_STARTB = ',t3_startb
  write(u6,*) 'T3_STOPB = ',t3_stopb
  call abend()
end if

if ((t3_starta < 0) .and. (t3_stopa > 0)) then
  write(u6,*) 'Mismatch in input : '
  write(u6,*) 'T3_STARTA = ',t3_starta
  write(u6,*) 'T3_STOPA = ',t3_stopa
  call abend()
end if

if ((t3_startb < 0) .and. (t3_stopb > 0)) then
  write(u6,*) 'Mismatch in input : '
  write(u6,*) 'T3_STARTB = ',t3_startb
  write(u6,*) 'T3_STOPB = ',t3_stopb
  call abend()
end if

!mp !if ((t3_starta < 0) .and. (t3_startb < 0)) then
!mp !  write(u6,*) 'This restart combination not implemented'
!mp !  write(u6,*) 'T3_STARTA = ',t3_starta
!mp !  write(u6,*) 'T3_STARTB = ',t3_startb
!mp !  call abend()
!mp !end if

!2 tlac hlavicky
write(u6,*)
write(u6,*) '    Cholesky Based Closed-Shell (T) code'
write(u6,*)
write(u6,*) '--------------------------------------------------'

write(u6,'(A,i9)') ' Frozen Orbitals                   : ',nfr
write(u6,'(A,i9)') ' Occupied Orbitals                 : ',no
write(u6,'(A,i9)') ' Virtual Orbitals                  : ',nv
write(u6,'(A,i9)') ' Total number of Cholesky Vectors  : ',nc

write(u6,*) '--------------------------------------------------'

write(u6,'(A,i9)') ' Large Virtual Segmentation        : ',NaGrp

write(u6,*) '--------------------------------------------------'

msg = 'No'
if (gen_files) msg = 'Yes'

write(u6,'(A,A5)') ' Generate Triples Scratch Files?        : ',msg

msg = 'No'
if (.not. run_triples) msg = 'Yes'

write(u6,'(A,A5)') ' Stop after Scratch Files generation?   : ',msg

write(u6,*) '--------------------------------------------------'

if (t3_starta == -1) then
  write(u6,'(A,i4)') ' Calculating full loop A                '
else
  write(u6,'(A,i4)') ' VO index triplet to start with in loop A : ',t3_starta
  write(u6,'(A,i4)') ' VO index triplet to stop  at   in loop A : ',t3_stopa
end if

if (t3_starta == -1) then
  write(u6,'(A,i4)') ' Calculating full loop B                '
else
  write(u6,'(A,i4)') ' VO index triplet to start with in loop B : ',t3_startb
  write(u6,'(A,i4)') ' VO index triplet to stop  at   in loop B : ',t3_stopb
end if

write(u6,*) '--------------------------------------------------'

write(u6,'(A,i9)') ' Lun Number for Aux. Matrixes      : ',LunAux
!write(u6,'(A,i9)') ' BLAS/FTN Matrix Handling          : ',mhkey

!msg = 'No'
!if (restkey == 1) msg = 'Yes'

!write(u6,'(A,A10)') ' Start from RstFil ?               : ',msg
write(u6,'(A,i9)') ' Print level                       : ',printkey

write(u6,*) '--------------------------------------------------'
write(u6,*)

return

end subroutine IniReord_t3
