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

subroutine action_ccsort(foka,fokb,fi,eps)

use ccsort_global, only: clopkey, daddr, Escf, fullprint, iokey, ISPIN, LSYM, luna1, luna2, luna3, luna4, lunab, lunda1, lunda2, &
                         lunt3, mapdri, mapiri, maxspace, mbas, NACTEL, noa, nob, NORB, nsize, NSYM, nva, nvb, posri0, reclen, &
                         t3key, typ
use CCT3_global, only: T3IntPos, T3Off
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoB

implicit none
real(kind=wp), intent(out) :: foka(mbas*(mbas+1)/2), fokb(mbas*(mbas+1)/2)
real(kind=wp), intent(in) :: fi(*), eps(mbas)
integer(kind=iwp) :: a, freespace, ickey, keyred, ndimv1, ndimv2, ndimv3, ndimvi, p, post, rc, symp, sympq, sympqr, symq, symr, &
                     syms, t3help1, t3help2, t3help3, t3help4, vsize, wrksize
integer(kind=iwp), allocatable :: AMMAP(:,:,:), ABMAP(:,:,:), JN(:,:), KN(:,:), LN(:,:), PQIND(:,:)
real(kind=wp), allocatable :: CCSORT(:), CCSORT2(:), VALN(:,:)

! distribute memory

!.1 calc. work space requirements
call initwrk(wrksize)

!.2 test if allocation of required memory is possible

!.2.1 allocate work space
call mma_maxdble(maxspace)
maxspace = maxspace-4
if (maxspace < wrksize) then
  write(u6,*) ' Allocation of work space failed!'
  write(u6,*) ' Increase the size of the variable MOLCAS_MEM'
  call Abend()
end if
call mma_allocate(CCSORT,wrksize,label='CCSORT')

!.3 set wrk = 0
CCSORT(:) = Zero

! def foka,fokb

ndimv1 = 0
do symp=1,nsym
  ndimv1 = ndimv1+(norb(symp)+1)*norb(symp)/2
end do

foka(1:ndimv1) = fi(1:ndimv1)
fokb(1:ndimv1) = fi(1:ndimv1)

! make names

call mktempanam()

! if T3 are requited, define T3IntPos and calc T3Off
if (t3key == 1) call DefT3par(noa,nsym)

! open files INTA1,INTA2,INTA3,INTA4 and INTAB

if (iokey == 1) then
  ! Fortran IO
  call molcas_binaryopen_vanilla(luna1,'INTA1')
  call molcas_binaryopen_vanilla(luna2,'INTA2')
  call molcas_binaryopen_vanilla(luna3,'INTA3')
  call molcas_binaryopen_vanilla(luna4,'INTA4')
  call molcas_binaryopen_vanilla(lunab,'INTAB')
  !open(unit=luna1,file='INTA1',form='unformatted')
  !open(unit=luna2,file='INTA2',form='unformatted')
  !open(unit=luna3,file='INTA3',form='unformatted')
  !open(unit=luna4,file='INTA4',form='unformatted')
  !open(unit=lunab,file='INTAB',form='unformatted')

else
  ! MOLCAS IO
  call daname(luna1,'INTA1')
  call daname(luna2,'INTA2')
  call daname(luna3,'INTA3')
  call daname(luna4,'INTA4')
  call daname(lunab,'INTAB')
  daddr(luna1) = 0
  daddr(luna2) = 0
  daddr(luna3) = 0
  daddr(luna4) = 0
  daddr(lunab) = 0
end if

! define #V1 (for <p,q,i,j>
call mkmappqij()

! make head of INTAB file for nonsymmetrical (C1) case
if (nsym == 1) call initintabc1()

! allocate space for ammap,abmap
call mma_allocate(AMMAP,mbas,8,8,Label='AMMAP')
call mma_allocate(ABMAP,mbas,mbas,8,Label='ABMAP')

do symp=1,nsym

  if (fullprint > 0) write(u6,'(6X,A,2X,I2)') 'Symmetry of the pivot index',symp

  ! define #V2 (for <_a,m,p,q>)
  call mkmapampq(symp)

  ! if T3 are required, make maps for R_i
  ! get mapd and mapi for  R_i(a,bc)
  if ((t3key == 1) .and. (noa(symp) > 0)) call ccsort_t3grc0(3,8,4,4,4,0,symp,posri0,post,mapdri,mapiri)

  ! open TEMPDA2 fils for <am|rs> integrals, if there are some virtiuals
  ! in symp symmetry
  if (nvb(symp) > 0) then
    call mkampqmap(AMMAP,symp,rc)
    call daopen('TEMPDA2 ',lunda2,reclen)
  end if

  do symq=1,nsym
    sympq = mul(symp,symq)

    ! open direct access file here, to enable exact specification of
    ! the number of records (only for symmetrical cases; syma>=symb)
    ! N.B. nrec is not needed now
    if ((nsym > 1) .and. (symp >= symq)) call daopen('TEMPDA1 ',lunda1,reclen)

    ! make abmap for syma>=symb
    if (symp >= symq) call mkabpqmap(ABMAP,symp,symq,rc)

    do symr=1,nsym
      sympqr = mul(sympq,symr)
      syms = sympqr

      ! calc size of the integral file
      if (symp == symr) then
        if (symp == symq) then
          vsize = (norb(symp)*(norb(symp)+1))/2
          vsize = vsize*(vsize+1)/2
          ickey = 3
        else
          vsize = (norb(symp)*(norb(symp)+1)*norb(symq)*(norb(symq)+1))/4
          ickey = 2
        end if
      else
        vsize = norb(symp)*norb(symq)*norb(symr)*norb(syms)
        ickey = 1
      end if

      if (vsize == 0) cycle

      if (fullprint > 1) write(u6,'(6X,A,I4,4X,4I2)') 'Block',typ(symp,symq,symr),symp,symq,symr,syms

      ! test for incore expansion
      call mma_maxdble(freespace)
      freespace = freespace-4
      if (fullprint >= 2) then
        write(u6,*)
        write(u6,'(6X,A,I10)') 'Available freespace   ',freespace
        write(u6,'(6X,A,I10)') 'Available freespace/MB',freespace*RtoB/1024**2
        write(u6,'(6X,A,I10)') 'Incore expansion      ',vsize+mbas*mbas
        write(u6,'(6X,A,I10)') 'Out of core expansion ',4*nsize*mbas
      end if

      if (freespace >= (vsize+mbas*mbas)) then
        ! INCORE EXPANSION
        if (fullprint >= 1) then
          write(u6,*)
          write(u6,'(6X,A)') 'Incore expansion      '
        end if
        call mma_allocate(CCSORT2,vsize,label='CCSORT2')

        if (ickey == 1) then
          ! case V(p,q,r,s)
          call esb_ic_1(symp,symq,symr,syms,CCSORT2,norb(symp),norb(symq),norb(symr),norb(syms))

        else if (ickey == 2) then
          ! case V(pr,qs)
          call mma_allocate(PQIND,mbas,mbas,label='PQIND')
          call esb_ic_2(symp,symq,CCSORT2,norb(symp),norb(symq),PQIND)
          call mma_deallocate(PQIND)

        else
          ! case V(prqs)
          call mma_allocate(PQIND,mbas,mbas,label='PQIND')
          call esb_ic_3(symp,CCSORT2,norb(symp),PQIND)
          call mma_deallocate(PQIND)

        end if

      else
        ! OUT OF CORE EXPANSION
        ! init temp files and realize expansion of this block
        ickey = 0
        if (fullprint >= 1) write(u6,'(6X,A)') 'Out of core expansion '
        call inittemp(norb(symp))
        if (freespace < (4*nsize*mbas)) then
          write(u6,*) ' Allocation of work space for Out-of-core failed!'
          write(u6,*) ' Increase the size of the variable MOLCAS_MEM'
          call Abend()
        end if

        ! allocate space for valn,jn,kn,ln
        call mma_allocate(VALN,nsize,mbas,label='VALN')
        call mma_allocate(JN,nsize,mbas,label='JN')
        call mma_allocate(KN,nsize,mbas,label='KN')
        call mma_allocate(LN,nsize,mbas,label='LN')

        call exppsb(symp,symq,symr,syms,VALN,JN,KN,LN)

        ! release space for valn,jn,kn,ln
        call mma_deallocate(VALN)
        call mma_deallocate(JN)
        call mma_deallocate(KN)
        call mma_deallocate(LN)

      end if

      ! run over all pivot indices

      do p=1,norb(symp)

        ! def dimensions of vint and read this block of integrals into vint
        ndimvi = norb(symp)
        ndimv1 = norb(symq)
        ndimv2 = norb(symr)
        ndimv3 = norb(syms)

        if (ickey == 0) then
          ! Out of core expansion

          if (symq == syms) then
            keyred = 1
          else
            keyred = 0
          end if
          call unpackk(p,CCSORT,ndimv1,ndimv2,ndimv3,keyred)

        else if (ickey == 1) then
          ! else Incore expansions

          ! case V(p,q,r,s)
          call unpackk_ic_1(p,CCSORT,ndimv1,ndimv2,ndimv3,CCSORT2,ndimvi)
        else if (ickey == 2) then
          ! case V(pr,qs)
          call unpackk_ic_2(p,CCSORT,ndimvi,ndimv1,CCSORT2)
        else
          ! case V(prqs)
          call unpackk_ic_3(p,CCSORT,ndimvi,CCSORT2)
        end if

        ! cycle

        ! add integrals to T3nam if needed (v nacechranej forme)

        if (t3key == 1) then
          if (p <= noa(symp)) then
            if (symq > syms) then

              ! calc proper address in t3nam file
              t3help4 = 0
              do t3help1=1,symp-1
                t3help4 = t3help4+noa(t3help1)
              end do
              t3help4 = t3help4+p
              t3help1 = mapiri(symr,symq,1)
              daddr(lunt3) = T3IntPos(t3help4)+T3Off(t3help1,symp)

              ! def required parameters
              t3help1 = nvb(symr)
              t3help2 = nvb(symq)
              t3help3 = nvb(syms)

              ! do packing
              call t3intpck2(CCSORT,CCSORT(posri0),ndimv1,ndimv2,ndimv3,t3help1,t3help2,t3help3,symq,symr,syms,nob,nvb)

            else if (symq == syms) then

              ! calc proper address in t3nam file
              t3help4 = 0
              do t3help1=1,symp-1
                t3help4 = t3help4+noa(t3help1)
              end do
              t3help4 = t3help4+p
              t3help1 = mapiri(symr,symq,1)
              daddr(lunt3) = T3IntPos(t3help4)+T3Off(t3help1,symp)

              ! def required parameters
              t3help1 = nvb(symr)
              t3help2 = nvb(symq)*(nvb(symq)+1)/2

              ! do packing
              call t3intpck1(CCSORT,CCSORT(posri0),ndimv1,ndimv2,ndimv3,t3help1,t3help2,symq,symr,syms,nob,nvb)

            end if
          end if
        end if

        ! add integrals to #1 <pq|ij> if needed

        ! contributions only for symi(r)>=symj(s)
        if (symr >= syms) call addpqij(CCSORT,wrksize,symp,symq,symr,syms,p,CCSORT,ndimv1,ndimv2,ndimv3)

        ! updete fok if necessary (only for open shell case)

        if (clopkey == 1) then

          if ((symp == symr) .and. (symq == syms) .and. (p > nob(symp)) .and. (p <= (noa(symp)))) then
            call fokupdate1(foka,fokb,symq,p,CCSORT,ndimv1,ndimv2,ndimv3)
          end if

          if ((symp == syms) .and. (symq == symr) .and. (p > nob(symp)) .and. (p <= (noa(symp)))) then
            call fokupdate2(foka,symq,p,CCSORT,ndimv1,ndimv2,ndimv3)
          end if

        end if

        ! add corresponding <am|pq> integrals to TEMPDA2
        ! and pack _a_brs to direct access file TEMPDA1 if needed and symm in not C1
        if (p > nob(symp)) then
          a = p-nob(symp)
          call ampack(CCSORT,wrksize,symp,symq,symr,syms,a,CCSORT,ndimv1,ndimv2,ndimv3,AMMAP)
          if ((nsym > 1) .and. (symp >= symq)) then
            call abpack(CCSORT,wrksize,symp,symq,symr,syms,a,CCSORT,ndimv1,ndimv2,ndimv3,ABMAP)
          end if
        end if

        ! add INTAB file (for nonsymmetrical (C1) state)

        if ((nsym == 1) .and. (p > nob(1))) call addintabc1(CCSORT,wrksize,p-nob(1),CCSORT,ndimv1)

      end do

      ! cycle

      ! close temp files
      !call closetemp(norb(symp))

      if (ickey >= 1) call mma_deallocate(CCSORT2)

    end do

    ! cycle

    if ((nsym > 1) .and. (symp >= symq)) then
      ! add contributions to INTAB comming from symp,sumq and close TEMPDA1 file
      ! only for symmetrical cases; only for syma>=symb
      call addintab(CCSORT,wrksize,symp,symq,ABMAP)
      close(lunda1)
      call vf('TEMPDA1 ',lunda1)
    end if

  end do

  ! cycle

  ! add contributions to INTA1-4 if there are some virtuals in symp symmetry
  ! and close TEMPDA2 files

  !if (nvb(symp) > 0) then
  call addinta(CCSORT,wrksize,symp,AMMAP)
  close(lunda2)
  call vf('TEMPDA2 ',lunda2)
  !end if

end do

! if T3 are required, reorganize T3nam file
if (t3key == 1) call t3reorg(CCSORT,wrksize,noa,nsym)

! release space for ammap,abmap
call mma_deallocate(AMMAP)
call mma_deallocate(ABMAP)

! close files INTA1,INTA2,INTA3 and INTA4, INTAB1

if (iokey == 1) then
  ! Fortran IO
  close(luna1)
  close(luna2)
  close(luna3)
  close(luna4)
  close(lunab)

else
  ! MOLCAS IO
  call daclos(luna1)
  call daclos(luna2)
  call daclos(luna3)
  call daclos(luna4)
  call daclos(lunab)
end if

!return

! def static integrals (file INTSTA)

call mkintsta(CCSORT,wrksize,foka,fokb)

! write general informations to INPDAT

call molcas_binaryopen_vanilla(1,'INPDAT')
!open(unit=1,file='INPDAT',form='unformatted')
write(1) NACTEL,ISPIN,NSYM,LSYM,mul,noa,nob,nva,nvb,norb,eps,Escf
close(1)

! Release the memory
call mma_deallocate(CCSORT)

return

end subroutine action_ccsort
