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
        subroutine defcommon (nfr,no,nv)
!
! this routine do :
!
! define commons needed in DIRCC routines
!
        implicit none
        integer nfr,nv,no
!mp!        integer me,nprocs
!
#include "uhf.fh"
#include "param_cht3.fh"
!mp!        common /my_mpi_world_com/ me, nprocs
!
!mp!        include 'task_info_inc'
!mp!        include 'ws_conn_inc'
!
!       logical llmpi
!
#include "ioind.fh"
!
! ----  UHF -----
!
        noab(1)=no
        noab(2)=no
!
        nnoab(1)=noab(1)*(noab(1)-1)/2
        nnoab(2)=noab(2)*(noab(2)-1)/2
        nnoab(3)=noab(1)*noab(2)
!
        nuab(1)=nv
        nuab(2)=nv
!
        nnuab(1)=(nuab(1)*(nuab(1)-1))/2
        nnuab(2)=(nuab(2)*(nuab(2)-1))/2
        nnuab(3)=nuab(1)*nuab(2)
!
        ich(1)="A"
        ich(2)="B"
        ich(3)="C"
!
! ----  PARAM ----
!
        nso=0
!??????????????????
        it=1
        itlast=1
        nbf=nfr+no+nv
        nomx=nbf
        nu=nv
        mx2=(nbf*(nbf+1))/2
        nno=(no*(no+1))/2
        nnu=(nu*(nu+1))/2
        nuo=no*nu
!
! ------ my_mpi_world_com --------
!
! zatial pre sekvencny chod
!
!mp!      me=0
!mp!      nprocs=1
!mp!      llmpi=.false.
!mp!      nws=1
!mp!      iws(1)=1
!mp!      lws(1)=.true.
!
! ------ pre Get3DM -----------
!
! ------    IOPT    -----------
!
        IOPT(14)=6
        IOPT(30)=0
        IOPT(76)=0
        IOPT(93)=0
        IOPT(93)=64
        IOPT(95)=0
!
!  sets IOPT(27) to an extreme number to force one file - temporary !!!  preskumat !!!
!  in prder to be compatible with RHF i=2^31-1
      iopt(27)=2147483647
!
!
!
! zatial tolko
!
        return
        end
