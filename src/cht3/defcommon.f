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
        subroutine defcommon (nfr,no,nv)
c
c this routine do :
c
c define commons needed in DIRCC routines
c
        implicit none
        integer nfr,nv,no
cmp!        integer me,nprocs
c
#include "uhf.fh"
#include "param_cht3.fh"
cmp!        common /my_mpi_world_com/ me, nprocs
c
cmp!        include 'task_info_inc'
cmp!        include 'ws_conn_inc'
c
c       logical llmpi
c
#include "ioind.fh"
c
c ----  UHF -----
c
        noab(1)=no
        noab(2)=no
c
        nnoab(1)=noab(1)*(noab(1)-1)/2
        nnoab(2)=noab(2)*(noab(2)-1)/2
        nnoab(3)=noab(1)*noab(2)
c
        nuab(1)=nv
        nuab(2)=nv
c
        nnuab(1)=(nuab(1)*(nuab(1)-1))/2
        nnuab(2)=(nuab(2)*(nuab(2)-1))/2
        nnuab(3)=nuab(1)*nuab(2)
c
        ich(1)="A"
        ich(2)="B"
        ich(3)="C"
c
c ----  PARAM ----
c
        nso=0
c??????????????????
        it=1
        itlast=1
        nbf=nfr+no+nv
        nomx=nbf
        nu=nv
        mx2=(nbf*(nbf+1))/2
        nno=(no*(no+1))/2
        nnu=(nu*(nu+1))/2
        nuo=no*nu
c
c ------ my_mpi_world_com --------
c
c zatial pre sekvencny chod
c
cmp!      me=0
cmp!      nprocs=1
cmp!      llmpi=.false.
cmp!      nws=1
cmp!      iws(1)=1
cmp!      lws(1)=.true.
c
c ------ pre Get3DM -----------
c
c ------    IOPT    -----------
c
        IOPT(14)=6
        IOPT(30)=0
        IOPT(76)=0
        IOPT(93)=0
        IOPT(93)=64
        IOPT(95)=0
c
C  sets IOPT(27) to an extreme number to force one file - temporary !!!  preskumat !!!
C  in prder to be compatible with RHF i=2^31-1
      iopt(27)=2147483647
c
c
c
c zatial tolko
c
        return
        end
