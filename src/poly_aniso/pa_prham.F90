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

subroutine pa_prham(exch,npair,i_pair,nneq,neq,nexch,nmax,lmax,eso,HLIN1,HLIN3,HLIN9,HDIP,HKEX,HDMO,HITO,Dipol,AnisoLines1, &
                    AnisoLines3,AnisoLines9,KE,DM_Exchange,JITO_exchange)
! this function prints the exchange Hamiltonian
! it does not compute any new infromation

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, iwp, u6, CtoB

implicit none
integer(kind=iwp), intent(in) :: exch, npair, i_pair(npair,2), nneq, neq(nneq), nexch(nneq), nmax, lmax
real(kind=wp), intent(in) :: eso(nneq,nmax)
complex(kind=wp), intent(in) :: HLIN1(npair,nmax,nmax,nmax,nmax), HLIN3(npair,nmax,nmax,nmax,nmax), &
                                HLIN9(npair,nmax,nmax,nmax,nmax), HDIP(npair,nmax,nmax,nmax,nmax), &
                                HKEX(npair,nmax,nmax,nmax,nmax), HDMO(npair,nmax,nmax,nmax,nmax), HITO(npair,nmax,nmax,nmax,nmax)
logical(kind=iwp), intent(in) :: Dipol, AnisoLines1, AnisoLines3, AnisoLines9, KE, DM_exchange, JITO_exchange
integer(kind=iwp) :: i, i1, i2, ibas(exch,lmax), ibuf, icoord(lmax), intc(lmax), is1, is2, j, js1, js2, l, lb, lb1, lb2, lp, lpr, &
                     mem_local, nb, nb1, nb2, nind(lmax,2)
complex(kind=wp), allocatable :: H1(:), H2(:), H3(:), H4(:), HTOT(:)
integer(kind=iwp), external :: norder
real(kind=wp), external :: dznrm2_

!=======================================================================
if (npair == 0) then
  call WarningMessage(2,'PA_PRHAM: npair = 0')
  return
end if
if (nmax == 0) then
  call WarningMessage(2,'PA_PRHAM:  nmax = 0')
  return
end if
if (exch == 0) then
  call WarningMessage(2,'PA_PRHAM:  exch = 0')
  return
end if
if (lmax == 0) then
  call WarningMessage(2,'PA_PRHAM:  lmax = 0')
  return
end if
ibuf = npair*nmax*nmax*nmax*nmax
if (ibuf > 0) then
  if ((dznrm2_(ibuf,HLIN1,1) == Zero) .and. AnisoLines1) call WarningMessage(2,'PA_PRHAM:  HLIN1 is empty')
  if ((dznrm2_(ibuf,HLIN3,1) == Zero) .and. AnisoLines3) call WarningMessage(2,'PA_PRHAM:  HLIN3 is empty')
  if ((dznrm2_(ibuf,HLIN9,1) == Zero) .and. AnisoLines9) call WarningMessage(2,'PA_PRHAM:  HLIN9 is empty')
  if ((dznrm2_(ibuf,HDIP,1) == Zero) .and. Dipol) call WarningMessage(2,'PA_PRHAM:  HDIP is empty')
  if ((dznrm2_(ibuf,HKEX,1) == Zero) .and. KE) call WarningMessage(2,'PA_PRHAM:  HKEX is empty')
  if ((dznrm2_(ibuf,HDMO,1) == Zero) .and. DM_exchange) call WarningMessage(2,'PA_PRHAM:  HDMO is empty')
  if ((dznrm2_(ibuf,HITO,1) == Zero) .and. JITO_exchange) call WarningMessage(2,'PA_PRHAM:  HITO is empty')
end if !ibuf
!=======================================================================
! allocate memory
mem_local = 0
if (exch >= 0) then
  call mma_allocate(htot,exch,'htot')
  call mma_allocate(h1,exch,'h1')
  call mma_allocate(h2,exch,'h2')
  call mma_allocate(h3,exch,'h3')
  call mma_allocate(h4,exch,'h4')
  mem_local = mem_local+size(htot)*CtoB
  mem_local = mem_local+size(h1)*CtoB
  mem_local = mem_local+size(h2)*CtoB
  mem_local = mem_local+size(h3)*CtoB
  mem_local = mem_local+size(h4)*CtoB
end if

!do lp=1,nPair
!  write(u6,'(A,i2,A,i3)') 'i_Pair(',lp,',1)=',i_pair(lp,1)
!  write(u6,'(A,i2,A,i3)') 'i_Pair(',lp,',2)=',i_pair(lp,2)
!end do
! generate the tables:
lpr = 0
l = 0
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do
nind(l+1:,:) = 0
intc(1) = 1
if (lmax > 1) then
  do i=2,lmax
    i1 = nind(i-1,1)
    intc(i) = intc(i-1)*nexch(i1)
  end do
end if
do nb=1,exch
  nb1 = nb-1
  do i=1,lmax
    ibas(nb,lmax-i+1) = nb1/intc(lmax-i+1)
    nb1 = nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
  end do
end do

write(u6,*)
write(u6,'(a)') repeat('%',100)
write(u6,'(30x,a)') 'Hamiltonian of the Total Magnetic Interaction'
write(u6,'(a)') repeat('%',100)

write(u6,'(A)') 'Only matrix elements with non-zero absolute value are listed below'

if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|      Dipolar Magnetic Interaction   | ', &
                   '|      Kinetic Exchange Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 4
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. KE) .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|      Dipolar Magnetic Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 3
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|      Kinetic Exchange Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 3
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|      Dipolar Magnetic Interaction   | ','|      Kinetic Exchange Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 3
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. KE) .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 2
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. KE) .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|      Dipolar Magnetic Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 2
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. KE .and. (.not. JITO_exchange)) then
  write(u6,'(6A)') '                 ','|      Kinetic Exchange Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 2
! JITO is active below:
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|          ITO Exchange Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 2
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|      Dipolar Magnetic Interaction   | ','|          ITO Exchange Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 3
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. KE .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|      Kinetic Exchange Interaction   | ','|          ITO Exchange Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 3
else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. KE .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|      Dipolar Magnetic Interaction   | ','|      Kinetic Exchange Interaction   | ', &
                   '|          ITO Exchange Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 4
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|          ITO Exchange Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 3
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. KE) .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|      Dipolar Magnetic Interaction   | ', &
                   '|          ITO Exchange Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 4
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. KE .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|      Kinetic Exchange Interaction   | ', &
                   '|          ITO Exchange Interaction   | ','|       TOTAL  Magnetic Interaction   | '
  lpr = 4
else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. KE .and. JITO_exchange) then
  write(u6,'(6A)') '                 ','|       Lines  Model of Interaction   | ','|      Dipolar Magnetic Interaction   | ', &
                   '|      Kinetic Exchange Interaction   | ','|          ITO Exchange Interaction   | ', &
                   '|       TOTAL  Magnetic Interaction   | '
  lpr = 5

end if
write(u6,'(110A)') '-----------------',('|-----  Real  ---------  Imaginary  --| ',i=1,lpr)

do nb1=1,exch
  H1(:) = cZero
  H2(:) = cZero
  H3(:) = cZero
  H4(:) = cZero
  HTOT(:) = cZero
  do lp=1,npair
    icoord(:) = ibas(nb1,:)

    lb1 = i_pair(lp,1)
    lb2 = i_pair(lp,2)
    i1 = nind(lb1,1)
    i2 = nind(lb2,1)
    is1 = icoord(lb1)+1
    is2 = icoord(lb2)+1

    do js1=1,nexch(i1)
      icoord(lb1) = js1-1
      do js2=1,nexch(i2)
        icoord(lb2) = js2-1
        nb2 = norder(icoord,intc,lmax)

        H1(nb2) = H1(nb2)+HLIN1(lp,is1,js1,is2,js2)+HLIN3(lp,is1,js1,is2,js2)+HLIN9(lp,is1,js1,is2,js2)

        H2(nb2) = H2(nb2)+HDIP(lp,is1,js1,is2,js2)

        H3(nb2) = H3(nb2)+HKEX(lp,is1,js1,is2,js2)

        H4(nb2) = H4(nb2)+HITO(lp,is1,js1,is2,js2)

        HTOT(nb2) = HTOT(nb2)+HLIN1(lp,is1,js1,is2,js2)+HLIN3(lp,is1,js1,is2,js2)+HLIN9(lp,is1,js1,is2,js2)+ &
                    HDIP(lp,is1,js1,is2,js2)+HKEX(lp,is1,js1,is2,js2)+HDMO(lp,is1,js1,is2,js2)+HITO(lp,is1,js1,is2,js2)

      end do ! js2
    end do ! js1
  end do ! lp

  if (.not. KE) then
    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          HTOT(nb1) = HTOT(nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end if !.not.KE
  ! all matrices for HEXCH( Nb1, xxx) are known:
  ! proceed to print
  if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. KE) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H2(nb2)) > Zero) .or. (abs(H3(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H2(nb2),H3(nb2),HTOT(nb2)
    end do

  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. KE)) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H2(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H2(nb2),HTOT(nb2)
    end do

  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. KE) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H3(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H3(nb2),HTOT(nb2)
    end do

  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. KE) then
    do nb2=nb1,exch
      if ((abs(H2(nb2)) > Zero) .or. (abs(H3(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H2(nb2),H3(nb2),HTOT(nb2)
    end do

  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. KE)) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),HTOT(nb2)
    end do

  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. KE)) then
    do nb2=nb1,exch
      if ((abs(H2(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H2(nb2),HTOT(nb2)
    end do

  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. KE) then
    do nb2=nb1,exch
      if ((abs(H3(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H3(nb2),HTOT(nb2)
    end do
  ! JITO is active below:
  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. (.not. KE) .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H4(nb2),HTOT(nb2)
    end do
  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. KE) .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H2(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H2(nb2),H4(nb2),HTOT(nb2)
    end do
  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. (.not. Dipol) .and. KE .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H3(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H3(nb2),H4(nb2),HTOT(nb2)
    end do
  else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. KE .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H2(nb2)) > Zero) .or. (abs(H3(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H2(nb2),H3(nb2),H4(nb2),HTOT(nb2)
    end do
  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. KE) .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H4(nb2),HTOT(nb2)
    end do
  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. KE) .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H2(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H2(nb2),H4(nb2),HTOT(nb2)
    end do
  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. KE .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H3(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,4(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H3(nb2),H4(nb2),HTOT(nb2)
    end do
  else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. KE .and. JITO_exchange) then
    do nb2=nb1,exch
      if ((abs(H1(nb2)) > Zero) .or. (abs(H2(nb2)) > Zero) .or. (abs(H3(nb2)) > Zero) .or. (abs(H4(nb2)) > Zero) .or. &
          (abs(HTOT(nb2)) > Zero)) &
        write(u6,'(A,i4,A,i4,A,5(2ES19.11,2x))') '<',nb1,'| H |',nb2,'> =',H1(nb2),H2(nb2),H3(nb2),H4(nb2),HTOT(nb2)
    end do
  end if

end do !nb1

write(u6,'(10A)') '-----------------',('|-------------------------------------| ',i=1,lpr)

!=======================================================================
! deallocate memory
if (exch >= 0) then
  call mma_deallocate(htot)
  call mma_deallocate(h1)
  call mma_deallocate(h2)
  call mma_deallocate(h3)
  call mma_deallocate(h4)
end if

return

end subroutine pa_prham
