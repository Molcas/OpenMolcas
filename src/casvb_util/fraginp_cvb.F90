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

subroutine fraginp_cvb(ip_iconfs)

implicit real*8(a-h,o-z)
parameter(nstrin=2,ncmp=4)
character*8 string(nstrin)
dimension dum(1)
save string
data string/'WAVE    ','CON    '/
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "frag_cvb.fh"
#include "WrkSpc.fh"

1000 call fstring_cvb(string,nstrin,istr,ncmp,2)
if (istr == 1) then
  ! 'WAVE'
  nfrag = nfrag+1
  nel_fr(nfrag) = 0
  call int_cvb(nel_fr(nfrag),1,nread,1)
  nMs_fr(nfrag) = 0
  nS_fr(nfrag) = 0
  nalf_fr(1,nfrag) = 0
  nbet_fr(1,nfrag) = 0
  i2s_fr(1,nfrag) = -1
100 Scurr = -one
  call real_cvb(dum,1,nread,1)
  Scurr = dum(1)
  if (Scurr /= -one) then
    nS_fr(nfrag) = nS_fr(nfrag)+1
    i2s_fr(nS_fr(nfrag),nfrag) = nint(2d0*Scurr)
    goto 100
  end if
else if (istr == 2) then
  ! 'CON'
  if (nfrag == 0) then
    nfrag = 1
    nel_fr(nfrag) = 0
    nMs_fr(nfrag) = 0
    nS_fr(nfrag) = 0
    nalf_fr(1,nfrag) = 0
    nbet_fr(1,nfrag) = 0
    i2s_fr(1,nfrag) = -1
  end if

  mxconf = max(mavaili_cvb()-1000,0)/noe
  call mrealloci_cvb(ip_iconfs,noe*mxconf)
  nconf_fr(nfrag) = 1
  nconf = nconf+1
200 continue
  if (mxconf < nconf) then
    write(6,*) ' Insufficient memory for configuration read',mavaili_cvb(),mxconf,nconf
    call abend_cvb()
  end if
  call izero(iwork(noe*(nconf-1)+ip_iconfs),noe)
  call int_cvb(iwork(noe*(nconf-1)+ip_iconfs),noe,nread,1)
  call fstring_cvb('CON',1,istr2,3,2)
  if (istr2 /= 0) then
    nconf_fr(nfrag) = nconf_fr(nfrag)+1
    nconf = nconf+1
    goto 200
  end if
  call mrealloci_cvb(ip_iconfs,noe*nconf)
end if
if (istr /= 0) goto 1000

return

end subroutine fraginp_cvb
