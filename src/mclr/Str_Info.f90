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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
!        OCSTR        :        Offsets for occupation of strings
!        STREO        :        reordering array
!        STSM         :        Symmetry of each string
!        STCL         :        Class of each string
!        NSTSO        :        Number of strings per symmetry and occupation
!        ISTSO        :        Offset of strings per symmetry and occupation
!        EL1          :        Number of electrons in RAS1 per sub type
!        EL3          :        Number of electrons in RAS3 per sub type
!        ACTP         :        is sub-type active
!        Z            :         Lexical adressing of arrays
!        EL123        :        -"-    But array
!        STSTMI       :        Explicit offsets and lengths
!        STSTMN       :                  "







Module Str_Info
      Implicit None
      Type String_Info
           Sequence
           Integer, Pointer:: OCSTR(:)=>Null()
           Integer, Allocatable:: OCSTR_hidden(:)
           Integer, Pointer:: STREO(:)=>Null()
           Integer, Allocatable:: STREO_hidden(:)
           Integer, Pointer:: STSM(:)=>Null()
           Integer, Allocatable:: STSM_hidden(:)
           Integer, Pointer:: STCL(:)=>Null()
           Integer, Allocatable:: STCL_hidden(:)
           Integer, Pointer:: NSTSO(:)=>Null()
           Integer, Allocatable:: NSTSO_hidden(:)
           Integer, Pointer:: ISTSO(:)=>Null()
           Integer, Allocatable:: ISTSO_hidden(:)
           Integer, Pointer:: EL1(:)=>Null()
           Integer, Allocatable:: EL1_hidden(:)
           Integer, Pointer:: EL3(:)=>Null()
           Integer, Allocatable:: EL3_hidden(:)
           Integer, Pointer:: ACTP(:)=>Null()
           Integer, Allocatable:: ACTP_hidden(:)
           Integer, Pointer:: Z(:)=>Null()
           Integer, Allocatable:: Z_hidden(:)
           Integer, Pointer:: EL123(:)=>Null()
           Integer, Allocatable:: EL123_hidden(:)
           Integer, Allocatable:: STSTMI(:)
           Integer, Allocatable:: STSTMN(:)
      End Type String_Info

      Type (String_Info), Allocatable, Target:: Str(:)
End Module Str_Info
