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
      subroutine gugadrt_njexcit(indjk,ljk,iextbit,nextbit,ivalid,      &
     &                           jstep,kttmp,k0)
#include "gendrt.fh"
#include "Sysdrt.fh"
#include "refstate.fh"
#include "ref.fh"
      dimension indjk(ljk),itexcit(n_ref)

      kp=k0
      do idxref=1,n_ref
        call upacknod(indjk,idxref,ival,nextbit,iextbit,ljk)
        if(jstep.eq.1.or.jstep.eq.2) then
          if(iref_occ(kp+1,idxref).eq.0) ival=ival+1
        endif
        if(jstep.eq.3) then
          if(iref_occ(kp+1,idxref).eq.0) ival=ival+2
          if(iref_occ(kp+1,idxref).eq.1) ival=ival+1
        endif
        if(ival.gt.2) ival=3
        itexcit(idxref)=ival
      enddo
      inm=itexcit(1)
      do idxref=2,n_ref
        inm=min(inm,itexcit(idxref))
      enddo

      if(inm.gt.2) then
        ivalid=0
      else
        kttmp=inm
        ivalid=1
        if(jstep.ne.0) then
          do i=1,n_ref
            ival=itexcit(i)
            call packnod(indjk,i,ival,nextbit,iextbit,ljk)
          enddo
        endif
      endif

      return
      end
