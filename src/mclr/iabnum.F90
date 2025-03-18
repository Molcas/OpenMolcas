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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      Integer FUNCTION IABNUM(IASTR,IBSTR,IAGRP,IBGRP,IGENSG,           &
     &                ISGNA,ISGNB,ISGNAB,IOOS,NORB,IPSFAC,PSSIGN,       &
     &                IPRNT)
      Use Str_info, only: STR,nElec,NoCTyp
!
! Encapsulation routine for IABNUS
!
      IMPLICIT None
      INTEGER IASTR(*),IBSTR(*)
      INTEGER IAGRP,IBGRP,IGENSG
      INTEGER ISGNA(*),ISGNB(*)
      INTEGER ISGNAB
      INTEGER IOOS(NOCTYP(IAGRP),NOCTYP(IBGRP),*)
      INTEGER NORB,IPSFAC
      REAL*8 PSSIGN
      INTEGER IPRNT

      INTEGER, EXTERNAL:: IABNUS
!
      IABNUM = IABNUS(IASTR,NELEC(IAGRP),Str(IAGRP)%STREO,              &
     &                Str(IAGRP)%STCL,Str(IAGRP)%STSM,                  &
     &                NOCTYP(IAGRP),                                    &
     &                Str(IAGRP)%Z,Str(IAGRP)%ISTSO,                    &
     &                Str(IAGRP)%NSTSO,                                 &
     &                IBSTR,NELEC(IBGRP),Str(IBGRP)%STREO,              &
     &                Str(IBGRP)%STCL,Str(IBGRP)%STSM,                  &
     &                NOCTYP(IBGRP),                                    &
     &                Str(IBGRP)%Z,Str(IBGRP)%ISTSO,                    &
     &                Str(IBGRP)%NSTSO,                                 &
     &                IOOS,NORB,IGENSG,ISGNA,ISGNB,ISGNAB,PSSIGN,       &
     &                IPSFAC,IPRNT)
      END FUNCTION IABNUM
