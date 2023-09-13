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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine permvb2_cvb(v1,iperm,vb,iapr,ixapr,xalf,xbet,mingrph,maxgrph,nk,locc,lunocc,inewocc,inocc2,negs,inda,phsa,indb,phsb,v2, &
                       ialg)

implicit real*8(a-h,o-w,y-z),integer(x)
logical vb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension iperm(norb)
dimension iapr(ndetvb), ixapr(nda+1)
dimension xalf(0:norb,0:nalf), xbet(0:norb,0:nbet)
dimension mingrph(0:norb), maxgrph(0:norb)
dimension nk(0:norb), locc(norb+1), lunocc(norb+1)
dimension inewocc(norb), inocc2(norb), negs(norb)
dimension inda(nda), phsa(nda), indb(ndb), phsb(ndb)
! V1 is dimensioned either NDET or NDETVB according to CI/VB
! V2 is dimensioned NDET/NDA or NDETVB according to CI/VB
dimension v1(*), v2(*)

! Some tests of permutation
! Valid?
call izero(negs,norb)
do i=1,norb
  iprm = abs(iperm(i))
  if ((iprm < 1) .or. (iprm > norb)) then
    write(6,*) ' Illegal orbital permutation!'
    call abend_cvb()
  end if
  negs(iprm) = negs(iprm)+1
end do
do iorb=1,norb
  if (negs(iorb) /= 1) then
    write(6,*) ' Illegal orbital permutation!'
    call abend_cvb()
  end if
end do
! Return if identity
do iorb=1,norb
  if (iperm(iorb) /= iorb) goto 35
end do
return
35 continue
! Use IALG=2 if only phase changes
do iorb=1,norb
  if (abs(iperm(iorb)) /= iorb) goto 45
end do
ialg = 2
45 continue
call izero(negs,norb)
do i=1,norb
  if (iperm(i) < 0) negs(abs(iperm(i))) = 1
end do
! Alpha loop:
call izero(inocc2,norb)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf,0)
  maxgrph(iorb) = min(iorb,nalf)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
index = 1
200 continue
call izero(inewocc,norb)
do ialf=1,nalf
  inewocc(abs(iperm(locc(ialf)))) = ialf
end do
ineg = 0
ia = 0
do iorb=1,norb
  if (inewocc(iorb) /= 0) then
    ia = ia+1
    inocc2(ia) = inewocc(iorb)
    inewocc(iorb) = 1
    if (negs(iorb) == 1) ineg = ineg+1
  end if
end do
if (mod(ineg,2) == 0) then
  phsa(index) = party_cvb(inocc2,nalf)
else
  phsa(index) = -party_cvb(inocc2,nalf)
end if
inda(index) = indget_cvb(inewocc,nalf,norb,xalf)

call loind_cvb(norb,nalf,nk,mingrph,maxgrph,locc,lunocc,index,xalf,*200)
! Beta loop:
call izero(inocc2,norb)
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet,0)
  maxgrph(iorb) = min(iorb,nbet)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)
call imove_cvb(maxgrph,nk,norb+1)
call occupy_cvb(nk,norb,locc,lunocc)
index = 1
500 continue
call izero(inewocc,norb)
do ibet=1,nbet
  inewocc(abs(iperm(locc(ibet)))) = ibet
end do
ineg = 0
ib = 0
do iorb=1,norb
  if (inewocc(iorb) /= 0) then
    ib = ib+1
    inocc2(ib) = inewocc(iorb)
    inewocc(iorb) = 1
    if (negs(iorb) == 1) ineg = ineg+1
  end if
end do
if (mod(ineg,2) == 0) then
  phsb(index) = party_cvb(inocc2,nbet)
else
  phsb(index) = -party_cvb(inocc2,nbet)
end if
indb(index) = indget_cvb(inewocc,nbet,norb,xbet)

call loind_cvb(norb,nbet,nk,mingrph,maxgrph,locc,lunocc,index,xbet,*500)

if (vb) then
  call fzero(v2,ndetvb)
  do ia=1,nda
    iato = inda(ia)
    do ixa=ixapr(ia),ixapr(ia+1)-1
      ib = iapr(ixa)
      ibto = indb(ib)
      do ixato=ixapr(iato),ixapr(iato+1)-1
        if (iapr(ixato) == ibto) goto 1200
      end do
      ! Shouldn't get here ...
      write(6,'(a,100i3)') ' Error, VB determinants not closed under permutation :',iperm
      call abend_cvb()
1200  continue
      v2(ixa) = phsa(ia)*phsb(ib)*v1(ixato)
    end do
  end do
  call fmove_cvb(v2,v1,ndetvb)
else if (ialg == 1) then
  ! Brute force strategy if enough memory (x1.5 faster):
  do ib=1,ndb
    iboff = (ib-1)*nda
    inboff = (indb(ib)-1)*nda
    do ia=1,nda
      v2(ia+iboff) = phsa(ia)*phsb(ib)*v1(inda(ia)+inboff)
    end do
  end do
  call fmove_cvb(v2,v1,ndet)
else if (ialg == 2) then
  ! More-or-less in-place update of V1:
  do ia=1,nda
    if (ia == inda(ia)) then
      if (phsa(ia) == -one) then
        ioffs = ia-nda
        do ib=1,ndb
          v1(ib*nda+ioffs) = -v1(ib*nda+ioffs)
        end do
      end if
    else if (inda(ia) /= 0) then
      ! Cyclic permutation involving IA:
      ioffs = ia-nda
      do ib=1,ndb
        v2(ib) = v1(ib*nda+ioffs)
      end do
      iat = ia
3400  continue
      if (phsa(iat) == one) then
        ioffs1 = iat-nda
        ioffs2 = inda(iat)-nda
        do ib=1,ndb
          v1(ib*nda+ioffs1) = v1(ib*nda+ioffs2)
        end do
      else
        ioffs1 = iat-nda
        ioffs2 = inda(iat)-nda
        do ib=1,ndb
          v1(ib*nda+ioffs1) = -v1(ib*nda+ioffs2)
        end do
      end if
      iatold = iat
      iat = inda(iat)
      inda(iatold) = 0
      if (inda(iat) /= ia) goto 3400
      if (phsa(iat) == one) then
        ioffs = iat-nda
        do ib=1,ndb
          v1(ib*nda+ioffs) = v2(ib)
        end do
      else
        ioffs = iat-nda
        do ib=1,ndb
          v1(ib*nda+ioffs) = -v2(ib)
        end do
      end if
      inda(iat) = 0
    end if
  end do
  do ib=1,ndb
    if (ib == indb(ib)) then
      if (phsb(ib) == -one) then
        ioffs = (ib-1)*nda
        do ia=1,nda
          v1(ia+ioffs) = -v1(ia+ioffs)
        end do
      end if
    else if (indb(ib) /= 0) then
      ! Cyclic permutation involving IB:
      ioffs = (ib-1)*nda
      do ia=1,nda
        v2(ia) = v1(ia+ioffs)
      end do
      ibt = ib
4400  continue
      if (phsb(ibt) == one) then
        ioffs1 = (ibt-1)*nda
        ioffs2 = (indb(ibt)-1)*nda
        do ia=1,nda
          v1(ia+ioffs1) = v1(ia+ioffs2)
        end do
      else
        ioffs1 = (ibt-1)*nda
        ioffs2 = (indb(ibt)-1)*nda
        do ia=1,nda
          v1(ia+ioffs1) = -v1(ia+ioffs2)
        end do
      end if
      ibtold = ibt
      ibt = indb(ibt)
      indb(ibtold) = 0
      if (indb(ibt) /= ib) goto 4400
      if (phsb(ibt) == one) then
        ioffs = (ibt-1)*nda
        do ia=1,nda
          v1(ia+ioffs) = v2(ia)
        end do
      else
        ioffs = (ibt-1)*nda
        do ia=1,nda
          v1(ia+ioffs) = -v2(ia)
        end do
      end if
      indb(ibt) = 0
    end if
  end do
end if

return

end subroutine permvb2_cvb
!**********************************
!** Routines involving CI and VB **
!**********************************
