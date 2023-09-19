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

subroutine symelminp_cvb(ip_symelm,nsyme,tags,izeta,mxirrep,mxorb,mxsyme,ityp)

implicit real*8(a-h,o-z)
! ... Matrices (orthogonal/determinant) ...
logical, external :: mxorth_cvb
#include "WrkSpc.fh"
parameter(nsymelm=5,nsign=2,ncmp=4)
character*8 symelm(nsymelm), sign(nsign)
character*3 tags(mxsyme)
dimension izeta(*)
dimension ityp(mxorb)
dimension iaux(1), daux(1)
save symelm, sign
data symelm/'IRREPS  ','COEFFS  ','TRANS   ','END     ','ENDSYMEL'/
data sign/'+       ','-       '/
save zero, one
data zero/0d0/,one/1d0/

nsyme = nsyme+1
if (nsyme > mxsyme) then
  write(6,*) ' Too many symmetry elements found :',nsyme,mxsyme
  call abend_cvb()
end if
tags(nsyme) = ' '
call string_cvb(tags(nsyme),1,nread,1)
call fstring_cvb(sign,nsign,isign,ncmp,1)
if (isign == 1) then
  izeta(nsyme) = 1
else if (isign == 2) then
  izeta(nsyme) = -1
else
  izeta(nsyme) = 0
end if
call mreallocr_cvb(ip_symelm,mxorb*mxorb*nsyme)
ishft = mxorb*mxorb*(nsyme-1)
call mxunit_cvb(work(ishft+ip_symelm),mxorb)

do
  call fstring_cvb(symelm,nsymelm,istr2,ncmp,2)
  if (istr2 == 1) then
    ! 'IRREPS'
    do i=1,mxirrep
      iaux = 0
      call int_cvb(iaux,1,nread,0)
      irrep = iaux(1)
      if (irrep /= 0) then
        do iorb=1,mxorb
          if (irrep == ityp(iorb)) work(iorb+(iorb-1)*mxorb+ishft+ip_symelm-1) = -one
        end do
      end if
    end do
  else if (istr2 == 2) then
    ! 'COEFFS'
    do i=1,mxorb
      iaux = 0
      call int_cvb(iaux,1,nread,0)
      iorb = iaux(1)
      if (iorb /= 0) then
        work(iorb+(iorb-1)*mxorb+ishft+ip_symelm-1) = -one
      else
        exit
      end if
    end do
  else if (istr2 == 3) then
    ! 'TRANS'
    iaux = 0
    call int_cvb(iaux,1,nread,0)
    idim = iaux(1)
    if ((idim < 1) .or. (idim > mxorb)) then
      write(6,*) ' Illegal dimension in TRANS:',idim,mxorb
      call abend_cvb()
    end if
    itmp = mstacki_cvb(idim)
    do i=1,idim
      call int_cvb(iaux,1,nread,0)
      iorb = iaux(1)
      if ((iorb < 1) .or. (iorb > mxorb)) then
        write(6,*) ' Illegal orbital number in TRANS:',iorb
        call abend_cvb()
      end if
      iwork(i+itmp-1) = iorb
    end do
    do ior=1,idim
      iorb = iwork(ior+itmp-1)
      do jor=1,idim
        jorb = iwork(jor+itmp-1)
        daux = zero
        call real_cvb(daux(1),1,nread,0)
        work(iorb+(jorb-1)*mxorb+ishft+ip_symelm-1) = daux(1)
      end do
    end do
    call mfreei_cvb(itmp)
  end if
  ! 'END' , 'ENDSYMEL' or unrecognized keyword -- end SYMELM input:
  if ((istr2 == 4) .or. (istr2 == 5) .or. (istr2 == 0)) exit
end do
if (.not. mxorth_cvb(work(ishft+ip_symelm),mxorb)) then
  write(6,*) ' Symmetry element ',tags(nsyme),' not orthogonal!'
  write(6,*) ' Check usage of TRANS keyword.'
  call abend_cvb()
end if

return

end subroutine symelminp_cvb
