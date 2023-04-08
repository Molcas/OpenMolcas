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
        subroutine checkMem(NvGrp, NvSGrp, NchBlk,                      &
     & Jal1, Jal2, wrksize, maxdim)
!
! Check memory consumption according to the specified
! orbital segmentation
!
        implicit none
#include "chcc1.fh"
!
        integer NvGrp,NvSGrp,NchBlk
        integer wrksize
!
        integer maxdim,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        integer PossV1,PossV2,PossV3,PossV4
        integer PossL11,PossL12
        integer PossH1,PossH2,PossH3,PossH4,PossH5
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
        integer PossT,PossMax
!        jalove
        integer Jal1,Jal2,Jal3
!
!2      Distribute memory (and redefine wrksize)
!
!*.1        Distribute Permanent arrays
!
        PossT=1
        call DistMemPerm (PossT)
        PossMax=PossT
!
!*.2    Distribute Work arrays for Reord
!
        call DefParReord (NvGrp,maxdim)
!        also
        call DefParo2v4 (NvGrp,NvGrp,NvSGrp,NvSGrp,                     &
     &                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
        PossT=PossFree
        call DistMemReord (NvGrp,maxdim,mdSGrpa,NchBlk,                 &
     &       PossV1,PossV2,PossV3,PossV4,PossM1,PossM2,                 &
     &       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
!
!*.3    Distribute Work arrays for o3v3
!
        call DefParo3v3 (NvGrp,maxdim)
!
        PossT=PossFree
        call DistMemo3v3jk (NvGrp,maxdim,                               &
     &       PossV1,PossV2,PossV3,PossV4,                               &
     &       PossH1,PossH2,PossH3,PossH4,PossH5,                        &
     &       PossK,PossQ,                                               &
     &       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
!
        PossT=PossFree
        call DistMemo3v3t2 (NvGrp,maxdim,                               &
     &       PossV1,PossV2,PossV3,PossV4,                               &
     &       PossH1,PossH2,PossH3,PossH4,                               &
     &       PossK,PossQ,                                               &
     &       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
!
        PossT=PossFree
        call DistMemo3v3chol (NvGrp,maxdim,                             &
     &       PossV1,PossV2,PossV3,PossV4,                               &
     &       PossH1,PossH2,PossH3,PossH4,                               &
     &       PossM1,PossM2,PossM3,PossM4,PossM5,                        &
     &       PossK,PossQ,                                               &
     &       PossT)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
!
!*.4    Distribute Work arrays for o2v4
!
        call DefParo2v4 (NvGrp,NvGrp,NvSGrp,NvSGrp,                     &
     &                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
        PossT=PossFree
!        Nazvy premennych tu davam ine, lebo ide len o vypocet narokov na Mem
        call DistMemo2v4 (NvGrp,NvGrp,NvSGrp,NvSGrp,                    &
     &                 mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,                 &
     &                 PossV1,PossV2,PossV3,PossV4,                     &
     &                 PossL11,PossL12,                                 &
     &                 PossH1,PossH2,PossH3,PossH4,PossQ,               &
     &                 Jal1,Jal2,                                       &
     &                 PossM1,PossM2,PossM3,PossM4,PossM5,Jal3,         &
     &                 PossT,PossK)
        if (PossT.gt.PossMax) then
          PossMax=PossT
        end if
!
!*.5        Memory requirements for Summary step are upperbounded by o3v3
!        and also o2v4, this it is skipped
!
        wrksize=PossMax
!
        return
        end
