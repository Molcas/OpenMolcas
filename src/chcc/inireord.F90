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

subroutine IniReord(NaGrp,NaSGrp,NchBlk,LunAux)
! nacitanie vsupu a inicializacia premnennych
! a tlac primitivnej hlavicky pre Reord procesz

use chcc_global, only: conv, generkey, intkey, JoinLkey, maxiter, maxGrp, maxSGrp, mhkey, nc, nfr, no, nv, printkey, restkey, &
                       W34DistKey
#ifdef _MOLCAS_MPP_
use chcc_global, only: NChLoc
use Para_Info, only: nProcs
#endif
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: NaGrp, NaSGrp, NchBlk, LunAux
integer(kind=iwp) :: intkey1, intkey2, LuSpool, NchBlk_tmp, NChLoc_max, NChLoc_min, nDel(8), ndelvirt, nFro(8), nOcc(8), nOrb(8)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: jal1
#endif
character(len=80) :: LINE
character(len=3) :: msg

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
conv = 1.0e-6_wp
printkey = 1
maxiter = 40

!mp read input file

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
do
  read(LuSpool,'(A80)') LINE
  call UPCASE(LINE)
  if (index(LINE,'&CHCC') /= 0) exit
end do
do
  read(LuSpool,'(A80)') LINE
  if (LINE(1:1) == '*') cycle
  call UPCASE(LINE)

  select case (LINE(1:4))

    case ('TITL')
      read(LuSpool,*)

    case ('FROZ')
      read(LuSpool,*) nfr
      if ((nfr < 0) .or. (nfr >= no)) then
        write(u6,*)
        write(u6,*) 'Illegal value for FROZen keyword : ',nfr
        call abend()
      end if
      no = no+nFro(1)-nfr

    case ('DELE')
      read(LuSpool,*) ndelvirt
      if ((ndelvirt < 0) .or. (ndelvirt > nv)) then
        write(u6,*)
        write(u6,*) 'Illegal value for DELETED keyword : ',ndelvirt
        call abend()
      end if
      nv = nv+nDel(1)-ndelvirt

    case ('LARG')
      read(LuSpool,*) NaGrp
      if ((NaGrp < 0) .or. (NaGrp > maxGrp)) then
        write(u6,*)
        write(u6,*) 'Illegal value for LARGE keyword : ',NaGrp
        write(u6,*) 'Large segmentation must be <= ',maxGrp
        call abend()
      end if

    case ('SMAL')
      read(LuSpool,*) NaSGrp
      if ((NaSGrp < 0) .or. (NaSGrp > 8)) then
        write(u6,*)
        write(u6,*) 'Illegal value for SMALL keyword : ',NaSGrp
        write(u6,*) 'Small segmentation must be <= 8'
        call abend()
      end if

      ! large == 0, small != 0 => quit
      if ((NaGrp == 0) .and. (NaSGrp /= 0)) then
        write(u6,*)
        write(u6,*) 'Small segmentation must be specified'
        write(u6,*) 'with large segmentation, or both can'
        write(u6,*) 'be left unspecified'
        call abend()
      end if

      if (NaGrp /= 0) then
        ! large != 0, small == 0 => small = 1
        if (NaSGrp == 0) NaSGrp = 1

        ! large * small <= 64
        if ((NaGrp*NaSGrp) > maxSGrp) then
          write(u6,*)
          write(u6,*) 'Product of Large and Small segmen-'
          write(u6,*) 'tation must be less or equal to ',maxSGrp
          call abend()
        end if
      end if

    case ('CHSE')
      read(LuSpool,*) NchBlk_tmp
      if ((NchBlk_tmp < 1) .or. (NchBlk_tmp > NChLoc_min)) then
        write(u6,*)
        write(u6,*) 'Illegal value for CHSegment keyword  : ',NchBlk_tmp
        write(u6,*) 'Reseting to a reasonable value for    '
        write(u6,*) 'this system :                         ',NchBlk
      else if (int(NChLoc_max/NchBlk_tmp) >= 100) then
        write(u6,*) 'Number of block of the MO Cholesky vector'
        write(u6,*) 'exceeded the limit. Increasing value of  '
        write(u6,*) 'the CHSEgmentation keyword to : ',NchBlk
      else
        NchBlk = NchBlk_tmp
      end if

    !mp case ('LUNA')  ... toto sa nikdy nevyuzivalo
    !mp   read(LuSpool,*) LunAux

    case ('MHKE')
      read(LuSpool,*) mhkey
      if ((mhkey < 0) .or. (mhkey > 2)) then
        mhkey = 1
        write(u6,*)
        write(u6,*) ' Warning!!!  Matrix handling key out of range'
        write(u6,*) ' parameter mhkey changed to 1'
      end if

    case ('NOGE')
      generkey = 0

    case ('ONTH')
      intkey1 = 1

    case ('PREC')
      intkey2 = 1

    case ('NODI')
      W34DistKey = 0

    case ('JOIN')
      read(LuSpool,*) JoinLkey
      if ((JoinLkey < 0) .or. (JoinLkey > 3)) then
        write(u6,*)
        write(u6,*) 'Illegal value for Join keyword : ',JoinLkey
        write(u6,*) 'Use one of 0, 1, 2, 3'
        write(u6,*) 'For details, see the manual ...'
        call abend()
      end if

    case ('MAXI')
      read(LuSpool,*) maxiter
      if (maxiter <= 0) then
        write(u6,*)
        write(u6,*) 'Illegal value of the MAXITER keyword: ',maxiter
        write(u6,*) 'Use integer > 0'
        call abend()
      end if

    case ('REST')
      restkey = 1
      write(u6,*)
      write(u6,*) 'This option is temporary disabled'
      write(u6,*) 'No Restart possible (... yet).'
      call abend()

    case ('THRE')
      read(LuSpool,*) conv

    case ('PRIN')
      read(LuSpool,*) printkey
      if (((printkey < 0) .or. (printkey > 10)) .or. ((printkey > 2) .and. (printkey < 10))) then

        write(u6,*)
        write(u6,*) 'Illegal value of the PRINT keyword: ',printkey
        write(u6,*) ' Use: 1  (Minimal) '
        write(u6,*) '      2  (Minimal + Timings)'
        write(u6,*) '      10 (Debug) '
        call abend()
      end if

    case ('END ')
      exit

  end select
end do

call Close_LuSpool(LuSpool)

! take care of the algorithm keyword
if (intkey1 == intkey2) then
  if (intkey1 == 0) then
    write(u6,*)
    write(u6,*) 'None of OnTheFly/PreCalculate'
    write(u6,*) 'algorithm was selected. Using'
    write(u6,*) 'default: PreCalculate (1)'
    intkey = 1
  else
    write(u6,*)
    write(u6,*) 'OnTheFly and PreCalculate keywords'
    write(u6,*) 'are mutually exclusive'
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
write(u6,*)
write(u6,*) '    Cholesky Based Closed-Shell CCSD code'
!mp write(u6,*) ' Dedicated to the memory of Boris Jeltzin'
write(u6,*)
write(u6,*) '---------------------------------------------------'

write(u6,'(A,i9)') ' Frozen Orbitals                   : ',nfr
write(u6,'(A,i9)') ' Occupied Orbitals                 : ',no
write(u6,'(A,i9)') ' Virtual Orbitals                  : ',nv
write(u6,'(A,i9)') ' Total number of Cholesky Vectors  : ',nc

write(u6,*) '---------------------------------------------------'

if (NaGrp /= 0) then
  write(u6,'(A,i9)') ' Large Virtual Segmentation        : ',NaGrp
else
  write(u6,'(A,A9)') ' Large Virtual Segmentation        :  auto'
end if

if (NaSGrp /= 0) then
  write(u6,'(A,i9)') ' Small Virtual Segmentation        : ',NaSGrp
else
  write(u6,'(A,A9)') ' Small Vectors Segmentation        :  auto'
end if

write(u6,'(A,i9)') ' Cholesky Vectors Segmentation     : ',NchBlk

write(u6,*) '---------------------------------------------------'

msg = 'No'
if (generkey == 1) msg = 'Yes'

write(u6,'(A,A4)') ' Generate Scratch Files?                : ',msg
write(u6,'(A,i4)') ' Precalculate (1) / On-the-Fly (0) Alg. : ',intkey
write(u6,'(A,i4)') ' 3 and 4-ext. MO integrals distribute?  : ',W34DistKey
write(u6,'(A,i4)') ' Parallel Join of varios MO integrals   : ',JoinLkey

write(u6,*) '---------------------------------------------------'

write(u6,'(A,E9.2)') ' Convergence Threshold             : ',conv
write(u6,'(A,i9)') ' Maximum number of Iterations      : ',maxiter

write(u6,*) '---------------------------------------------------'

write(u6,'(A,i9)') ' Lun Number for Aux. Matrixes      : ',LunAux
write(u6,'(A,i9)') ' BLAS/FTN Matrix Handling          : ',mhkey

msg = 'No'
if (restkey == 1) msg = 'Yes'

write(u6,'(A,A10)') ' Start from RstFil ?               : ',msg
write(u6,'(A,i9)') ' Print level                       : ',printkey

write(u6,*) '---------------------------------------------------'
write(u6,*)

return

end subroutine IniReord
