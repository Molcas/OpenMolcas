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

subroutine zasun_pck(i1,length,valn,jn,kn,ln)
! this routine writes one block of 3-indices and appropriate
! values of integrals into an open TEMP-file
!
! i1 - number of pivot index (I)
! length - number of valid integrals in block (I)
! this routine has also jn,kn,ln,valn
! and stattemp and tmpnam as inputs, but they are
! transported through commons  in reorg.fh

implicit real*8(a-h,o-z)
integer length, i1
#include "reorg.fh"
#include "SysDef.fh"
real*8 valn(1:nsize,1:mbas)
integer jn(1:nsize,1:mbas)
integer kn(1:nsize,1:mbas)
integer ln(1:nsize,1:mbas)
! help variable
integer m2, iRec
integer jkl(1:nsize)
integer constj
parameter(constj=1048576)
integer constk
parameter(constk=1024)
character*(RtoB+ItoB) pp(1:nsize), pphelp
real*8 rhelp
integer ihelp

! pack indices and integral values

do m2=1,length
  jkl(m2) = ln(m2,i1)+constj*jn(m2,i1)
end do
do m2=1,length
  jkl(m2) = jkl(m2)+constk*kn(m2,i1)
end do

do m2=1,length
  rhelp = valn(m2,i1)
  ihelp = jkl(m2)
  pphelp(1:RtoB) = transfer(rhelp,pphelp(1:RtoB))
  pphelp(RtoB+1:) = transfer(ihelp,pphelp(RtoB+1:))
  pp(m2) = pphelp
end do

! open corresponding TEMP file in corresponding form

if (iokey == 1) then

  ! Fortran IO

  if (stattemp(i1) == 0) then
    ! file will be opened first time, it must be opened
    ! with the pointer at then first position
    call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
    !open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown')
    stattemp(i1) = 1

  else
    ! file was already used in expansion of this block, it must
    ! be opened with the pointer at the end of the file
#   ifdef _DECAXP_
    call molcas_open_ext2(lunpublic,tmpnam(i1),'append','unformatted',f_iostat,.false.,1,'unknown',is_error)

    !open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown',access='append')
#   else
    call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
    !open(unit=lunpublic,file=tmpnam(i1),form='unformatted',status='unknown')
    do iRec=1,nrectemp(i1)
      read(lunpublic) m2
    end do
#   endif

  end if

  call zashlp1(lunpublic,pp,length)
  close(lunpublic)

else

  ! MOLCAS IO

  call daname(lunpublic,tmpnam(i1))
  call cdafile(lunpublic,1,pp,(RtoB+ItoB)*length,stattemp(i1))
  call daclos(lunpublic)

end if

nrectemp(i1) = nrectemp(i1)+1
lrectemp(i1) = length

return

end subroutine zasun_pck
