************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine inisewm(prgnam,ndiff)
      Implicit Real*8 (A-H,O-Z)
      character*(*) prgnam
      character*16 pgnm_local
#include "unixinfo.fh"
#include "para_info.fh"
      Logical DoRys

c      ProgName=prgnam
      pgnm_local=prgnam
      call locase(pgnm_local)
*
      info=0
      If(pgnm_local.eq.'seward') then
      else if(pgnm_local.eq.'scf') then
        DoRys=.true.
        call inisew(info,DoRys,ndiff)
      else if(pgnm_local.eq.'dtraf') then
        DoRys=.true.
        call inisew(info,DoRys,ndiff)
      else if(pgnm_local.eq.'dkext') then
        DoRys=.true.
        call inisew(info,DoRys,ndiff)
      else if(pgnm_local.eq.'mltpl') then
        DoRys=.true.  ! for Schwarz prescreening
        call inisew(info,DoRys,ndiff)
      else if(pgnm_local.eq.'alaska') then
        DoRys=.true.
        call inisew(info,DoRys,ndiff)
      else if(pgnm_local.eq.'mckinley') then
        DoRys=.true.
        call inisew(info,DoRys,ndiff)
      else if(pgnm_local.eq.'slapaf') then
      else if(pgnm_local.eq.'espf') Then
        DoRys=.true.
        call inisew(info,DoRys,ndiff)
      else
        DoRys=.false.
        call inisew(info,DoRys,ndiff)
      EndIf
*
      Return
      End
