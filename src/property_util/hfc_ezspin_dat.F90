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
! Copyright (C) 2026, Dong Q. Le                                       *
!***********************************************************************
! Data derived from EasySpin (https://easyspin.org/)                   *
! Copyright (c) 2026                                                   *
! Licensed under the MIT License.                                      *
! Please see the file readme_easyspin_license.md in the                *
! Tools/create_easy_spin_data directory for licence details.           *
!***********************************************************************


module hfc_data
  use stdalloc,      only: mma_allocate, mma_deallocate
  use Definitions,   only: iwp, wp
  implicit none

  integer, parameter :: numb_records = 352

  type :: easyspin_data

    integer(kind=iwp)             :: AtNumb
    integer(kind=iwp)             :: MassNum
    real(kind=wp)                 :: abundance
    real(kind=wp)                 :: nucspin
    real(kind=wp)                 :: gfactor
  end type easyspin_data

  type :: element_map
    integer(kind=iwp)             :: first_rec
    integer(kind=iwp)             :: last_rec
  end type element_map

  type(element_map),allocatable   :: elem_map_db(:)
  type(easyspin_data),allocatable :: ezspin_db(:)

    public :: init_isotope_data, gfac_spin_by_mass, get_first_nonzero_gfactor, &
              gfactor_by_nucspin, nucspin_by_gfactor, free_isotope_data

  contains

!======================================================================
  subroutine init_isotope_data()
    integer :: istat
    allocate(ezspin_db(352), stat=istat)

    if (istat > 0) then
      write(6,*) ' Error in allocating EasySpin ezspin_db.'
    end if

    allocate(elem_map_db(117), stat=istat)

    if (istat > 0) then
      write(6,*) ' Error in allocating EasySpin elem_map_db.'
    end if


!--> Element H
      elem_map_db(1)%first_rec = 1
      elem_map_db(1)%last_rec  = 3
!----------------------------
!--> Element He
      elem_map_db(2)%first_rec = 4
      elem_map_db(2)%last_rec  = 5
!----------------------------
!--> Element Li
      elem_map_db(3)%first_rec = 6
      elem_map_db(3)%last_rec  = 7
!----------------------------
!--> Element Be
      elem_map_db(4)%first_rec = 8
      elem_map_db(4)%last_rec  = 8
!----------------------------
!--> Element B
      elem_map_db(5)%first_rec = 9
      elem_map_db(5)%last_rec  = 10
!----------------------------
!--> Element C
      elem_map_db(6)%first_rec = 11
      elem_map_db(6)%last_rec  = 13
!----------------------------
!--> Element N
      elem_map_db(7)%first_rec = 14
      elem_map_db(7)%last_rec  = 15
!----------------------------
!--> Element O
      elem_map_db(8)%first_rec = 16
      elem_map_db(8)%last_rec  = 18
!----------------------------
!--> Element F
      elem_map_db(9)%first_rec = 19
      elem_map_db(9)%last_rec  = 19
!----------------------------
!--> Element Ne
      elem_map_db(10)%first_rec = 20
      elem_map_db(10)%last_rec  = 22
!----------------------------
!--> Element Na
      elem_map_db(11)%first_rec = 23
      elem_map_db(11)%last_rec  = 24
!----------------------------
!--> Element Mg
      elem_map_db(12)%first_rec = 25
      elem_map_db(12)%last_rec  = 27
!----------------------------
!--> Element Al
      elem_map_db(13)%first_rec = 28
      elem_map_db(13)%last_rec  = 28
!----------------------------
!--> Element Si
      elem_map_db(14)%first_rec = 29
      elem_map_db(14)%last_rec  = 31
!----------------------------
!--> Element P
      elem_map_db(15)%first_rec = 32
      elem_map_db(15)%last_rec  = 32
!----------------------------
!--> Element S
      elem_map_db(16)%first_rec = 33
      elem_map_db(16)%last_rec  = 36
!----------------------------
!--> Element Cl
      elem_map_db(17)%first_rec = 37
      elem_map_db(17)%last_rec  = 39
!----------------------------
!--> Element Ar
      elem_map_db(18)%first_rec = 40
      elem_map_db(18)%last_rec  = 43
!----------------------------
!--> Element K
      elem_map_db(19)%first_rec = 44
      elem_map_db(19)%last_rec  = 46
!----------------------------
!--> Element Ca
      elem_map_db(20)%first_rec = 47
      elem_map_db(20)%last_rec  = 53
!----------------------------
!--> Element Sc
      elem_map_db(21)%first_rec = 54
      elem_map_db(21)%last_rec  = 54
!----------------------------
!--> Element Ti
      elem_map_db(22)%first_rec = 55
      elem_map_db(22)%last_rec  = 59
!----------------------------
!--> Element V
      elem_map_db(23)%first_rec = 60
      elem_map_db(23)%last_rec  = 61
!----------------------------
!--> Element Cr
      elem_map_db(24)%first_rec = 62
      elem_map_db(24)%last_rec  = 65
!----------------------------
!--> Element Mn
      elem_map_db(25)%first_rec = 66
      elem_map_db(25)%last_rec  = 67
!----------------------------
!--> Element Fe
      elem_map_db(26)%first_rec = 68
      elem_map_db(26)%last_rec  = 71
!----------------------------
!--> Element Co
      elem_map_db(27)%first_rec = 72
      elem_map_db(27)%last_rec  = 73
!----------------------------
!--> Element Ni
      elem_map_db(28)%first_rec = 74
      elem_map_db(28)%last_rec  = 78
!----------------------------
!--> Element Cu
      elem_map_db(29)%first_rec = 79
      elem_map_db(29)%last_rec  = 80
!----------------------------
!--> Element Zn
      elem_map_db(30)%first_rec = 81
      elem_map_db(30)%last_rec  = 85
!----------------------------
!--> Element Ga
      elem_map_db(31)%first_rec = 86
      elem_map_db(31)%last_rec  = 87
!----------------------------
!--> Element Ge
      elem_map_db(32)%first_rec = 88
      elem_map_db(32)%last_rec  = 92
!----------------------------
!--> Element As
      elem_map_db(33)%first_rec = 93
      elem_map_db(33)%last_rec  = 93
!----------------------------
!--> Element Se
      elem_map_db(34)%first_rec = 94
      elem_map_db(34)%last_rec  = 100
!----------------------------
!--> Element Br
      elem_map_db(35)%first_rec = 101
      elem_map_db(35)%last_rec  = 102
!----------------------------
!--> Element Kr
      elem_map_db(36)%first_rec = 103
      elem_map_db(36)%last_rec  = 109
!----------------------------
!--> Element Rb
      elem_map_db(37)%first_rec = 110
      elem_map_db(37)%last_rec  = 111
!----------------------------
!--> Element Sr
      elem_map_db(38)%first_rec = 112
      elem_map_db(38)%last_rec  = 115
!----------------------------
!--> Element Y
      elem_map_db(39)%first_rec = 116
      elem_map_db(39)%last_rec  = 116
!----------------------------
!--> Element Zr
      elem_map_db(40)%first_rec = 117
      elem_map_db(40)%last_rec  = 121
!----------------------------
!--> Element Nb
      elem_map_db(41)%first_rec = 122
      elem_map_db(41)%last_rec  = 122
!----------------------------
!--> Element Mo
      elem_map_db(42)%first_rec = 123
      elem_map_db(42)%last_rec  = 129
!----------------------------
!--> Element Tc
      elem_map_db(43)%first_rec = 130
      elem_map_db(43)%last_rec  = 130
!----------------------------
!--> Element Ru
      elem_map_db(44)%first_rec = 131
      elem_map_db(44)%last_rec  = 137
!----------------------------
!--> Element Rh
      elem_map_db(45)%first_rec = 138
      elem_map_db(45)%last_rec  = 138
!----------------------------
!--> Element Pd
      elem_map_db(46)%first_rec = 139
      elem_map_db(46)%last_rec  = 144
!----------------------------
!--> Element Ag
      elem_map_db(47)%first_rec = 145
      elem_map_db(47)%last_rec  = 146
!----------------------------
!--> Element Cd
      elem_map_db(48)%first_rec = 147
      elem_map_db(48)%last_rec  = 154
!----------------------------
!--> Element In
      elem_map_db(49)%first_rec = 155
      elem_map_db(49)%last_rec  = 156
!----------------------------
!--> Element Sn
      elem_map_db(50)%first_rec = 157
      elem_map_db(50)%last_rec  = 166
!----------------------------
!--> Element Sb
      elem_map_db(51)%first_rec = 167
      elem_map_db(51)%last_rec  = 169
!----------------------------
!--> Element Te
      elem_map_db(52)%first_rec = 170
      elem_map_db(52)%last_rec  = 177
!----------------------------
!--> Element I
      elem_map_db(53)%first_rec = 178
      elem_map_db(53)%last_rec  = 179
!----------------------------
!--> Element Xe
      elem_map_db(54)%first_rec = 180
      elem_map_db(54)%last_rec  = 188
!----------------------------
!--> Element Cs
      elem_map_db(55)%first_rec = 189
      elem_map_db(55)%last_rec  = 192
!----------------------------
!--> Element Ba
      elem_map_db(56)%first_rec = 193
      elem_map_db(56)%last_rec  = 200
!----------------------------
!--> Element La
      elem_map_db(57)%first_rec = 201
      elem_map_db(57)%last_rec  = 203
!----------------------------
!--> Element Ce
      elem_map_db(58)%first_rec = 204
      elem_map_db(58)%last_rec  = 208
!----------------------------
!--> Element Pr
      elem_map_db(59)%first_rec = 209
      elem_map_db(59)%last_rec  = 209
!----------------------------
!--> Element Nd
      elem_map_db(60)%first_rec = 210
      elem_map_db(60)%last_rec  = 216
!----------------------------
!--> Element Pm
      elem_map_db(61)%first_rec = 217
      elem_map_db(61)%last_rec  = 217
!----------------------------
!--> Element Sm
      elem_map_db(62)%first_rec = 218
      elem_map_db(62)%last_rec  = 225
!----------------------------
!--> Element Eu
      elem_map_db(63)%first_rec = 226
      elem_map_db(63)%last_rec  = 230
!----------------------------
!--> Element Gd
      elem_map_db(64)%first_rec = 231
      elem_map_db(64)%last_rec  = 237
!----------------------------
!--> Element Tb
      elem_map_db(65)%first_rec = 238
      elem_map_db(65)%last_rec  = 240
!----------------------------
!--> Element Dy
      elem_map_db(66)%first_rec = 241
      elem_map_db(66)%last_rec  = 247
!----------------------------
!--> Element Ho
      elem_map_db(67)%first_rec = 248
      elem_map_db(67)%last_rec  = 248
!----------------------------
!--> Element Er
      elem_map_db(68)%first_rec = 249
      elem_map_db(68)%last_rec  = 254
!----------------------------
!--> Element Tm
      elem_map_db(69)%first_rec = 255
      elem_map_db(69)%last_rec  = 256
!----------------------------
!--> Element Yb
      elem_map_db(70)%first_rec = 257
      elem_map_db(70)%last_rec  = 263
!----------------------------
!--> Element Lu
      elem_map_db(71)%first_rec = 264
      elem_map_db(71)%last_rec  = 267
!----------------------------
!--> Element Hf
      elem_map_db(72)%first_rec = 268
      elem_map_db(72)%last_rec  = 273
!----------------------------
!--> Element Ta
      elem_map_db(73)%first_rec = 274
      elem_map_db(73)%last_rec  = 275
!----------------------------
!--> Element W
      elem_map_db(74)%first_rec = 276
      elem_map_db(74)%last_rec  = 280
!----------------------------
!--> Element Re
      elem_map_db(75)%first_rec = 281
      elem_map_db(75)%last_rec  = 282
!----------------------------
!--> Element Os
      elem_map_db(76)%first_rec = 283
      elem_map_db(76)%last_rec  = 289
!----------------------------
!--> Element Ir
      elem_map_db(77)%first_rec = 290
      elem_map_db(77)%last_rec  = 291
!----------------------------
!--> Element Pt
      elem_map_db(78)%first_rec = 292
      elem_map_db(78)%last_rec  = 297
!----------------------------
!--> Element Au
      elem_map_db(79)%first_rec = 298
      elem_map_db(79)%last_rec  = 298
!----------------------------
!--> Element Hg
      elem_map_db(80)%first_rec = 299
      elem_map_db(80)%last_rec  = 305
!----------------------------
!--> Element Tl
      elem_map_db(81)%first_rec = 306
      elem_map_db(81)%last_rec  = 308
!----------------------------
!--> Element Pb
      elem_map_db(82)%first_rec = 309
      elem_map_db(82)%last_rec  = 312
!----------------------------
!--> Element Bi
      elem_map_db(83)%first_rec = 313
      elem_map_db(83)%last_rec  = 314
!----------------------------
!--> Element Po
      elem_map_db(84)%first_rec = 315
      elem_map_db(84)%last_rec  = 315
!----------------------------
!--> Element At
      elem_map_db(85)%first_rec = 316
      elem_map_db(85)%last_rec  = 316
!----------------------------
!--> Element Rn
      elem_map_db(86)%first_rec = 317
      elem_map_db(86)%last_rec  = 317
!----------------------------
!--> Element Fr
      elem_map_db(87)%first_rec = 318
      elem_map_db(87)%last_rec  = 318
!----------------------------
!--> Element Ra
      elem_map_db(88)%first_rec = 319
      elem_map_db(88)%last_rec  = 319
!----------------------------
!--> Element Ac
      elem_map_db(89)%first_rec = 320
      elem_map_db(89)%last_rec  = 320
!----------------------------
!--> Element Th
      elem_map_db(90)%first_rec = 321
      elem_map_db(90)%last_rec  = 322
!----------------------------
!--> Element Pa
      elem_map_db(91)%first_rec = 323
      elem_map_db(91)%last_rec  = 323
!----------------------------
!--> Element U
      elem_map_db(92)%first_rec = 324
      elem_map_db(92)%last_rec  = 326
!----------------------------
!--> Element Np
      elem_map_db(93)%first_rec = 327
      elem_map_db(93)%last_rec  = 327
!----------------------------
!--> Element Pu
      elem_map_db(94)%first_rec = 328
      elem_map_db(94)%last_rec  = 328
!----------------------------
!--> Element Am
      elem_map_db(95)%first_rec = 329
      elem_map_db(95)%last_rec  = 329
!----------------------------
!--> Element Cm
      elem_map_db(96)%first_rec = 330
      elem_map_db(96)%last_rec  = 330
!----------------------------
!--> Element Bk
      elem_map_db(97)%first_rec = 331
      elem_map_db(97)%last_rec  = 331
!----------------------------
!--> Element Cf
      elem_map_db(98)%first_rec = 332
      elem_map_db(98)%last_rec  = 332
!----------------------------
!--> Element Es
      elem_map_db(99)%first_rec = 333
      elem_map_db(99)%last_rec  = 333
!----------------------------
!--> Element Fm
      elem_map_db(100)%first_rec = 334
      elem_map_db(100)%last_rec  = 334
!----------------------------
!--> Element Md
      elem_map_db(101)%first_rec = 335
      elem_map_db(101)%last_rec  = 335
!----------------------------
!--> Element No
      elem_map_db(102)%first_rec = 336
      elem_map_db(102)%last_rec  = 336
!----------------------------
!--> Element Lr
      elem_map_db(103)%first_rec = 337
      elem_map_db(103)%last_rec  = 337
!----------------------------
!--> Element Rf
      elem_map_db(104)%first_rec = 338
      elem_map_db(104)%last_rec  = 338
!----------------------------
!--> Element Db
      elem_map_db(105)%first_rec = 339
      elem_map_db(105)%last_rec  = 339
!----------------------------
!--> Element Sg
      elem_map_db(106)%first_rec = 340
      elem_map_db(106)%last_rec  = 340
!----------------------------
!--> Element Bh
      elem_map_db(107)%first_rec = 341
      elem_map_db(107)%last_rec  = 341
!----------------------------
!--> Element Hs
      elem_map_db(108)%first_rec = 342
      elem_map_db(108)%last_rec  = 342
!----------------------------
!--> Element Mt
      elem_map_db(109)%first_rec = 343
      elem_map_db(109)%last_rec  = 343
!----------------------------
!--> Element Ds
      elem_map_db(110)%first_rec = 344
      elem_map_db(110)%last_rec  = 344
!----------------------------
!--> Element Rg
      elem_map_db(111)%first_rec = 345
      elem_map_db(111)%last_rec  = 345
!----------------------------
!--> Element Cn
      elem_map_db(112)%first_rec = 346
      elem_map_db(112)%last_rec  = 346
!----------------------------
!--> Element Nh
      elem_map_db(113)%first_rec = 347
      elem_map_db(113)%last_rec  = 347
!----------------------------
!--> Element Fl
      elem_map_db(114)%first_rec = 348
      elem_map_db(114)%last_rec  = 348
!----------------------------
!--> Element Mc
      elem_map_db(115)%first_rec = 349
      elem_map_db(115)%last_rec  = 349
!----------------------------
!--> Element Lv
      elem_map_db(116)%first_rec = 350
      elem_map_db(116)%last_rec  = 350
!----------------------------
!--> Element Ts
      elem_map_db(117)%first_rec = 351
      elem_map_db(117)%last_rec  = 352
!----------------------------
!--> Record 001    = 1H
    ezspin_db(1)%AtNumb      = 1_iwp
    ezspin_db(1)%MassNum     = 1_iwp
    ezspin_db(1)%abundance    = 99.9885_wp
    ezspin_db(1)%nucspin      = 0.5_wp
    ezspin_db(1)%gfactor      = +5.58569468_wp
!----------------------------
!--> Record 002    = 2H
    ezspin_db(2)%AtNumb      = 1_iwp
    ezspin_db(2)%MassNum     = 2_iwp
    ezspin_db(2)%abundance    = 0.0115_wp
    ezspin_db(2)%nucspin      = 1.0_wp
    ezspin_db(2)%gfactor      = +0.8574382_wp
!----------------------------
!--> Record 003    = 3H
    ezspin_db(3)%AtNumb      = 1_iwp
    ezspin_db(3)%MassNum     = 3_iwp
    ezspin_db(3)%abundance    = 0.0_wp
    ezspin_db(3)%nucspin      = 0.5_wp
    ezspin_db(3)%gfactor      = +5.95799369_wp
!----------------------------
!--> Record 004    = 4He
    ezspin_db(4)%AtNumb      = 2_iwp
    ezspin_db(4)%MassNum     = 4_iwp
    ezspin_db(4)%abundance    = 99.999863_wp
    ezspin_db(4)%nucspin      = 0.0_wp
    ezspin_db(4)%gfactor      = 0.0_wp
!----------------------------
!--> Record 005    = 3He
    ezspin_db(5)%AtNumb      = 2_iwp
    ezspin_db(5)%MassNum     = 3_iwp
    ezspin_db(5)%abundance    = 0.000137_wp
    ezspin_db(5)%nucspin      = 0.5_wp
    ezspin_db(5)%gfactor      = -4.25499544_wp
!----------------------------
!--> Record 006    = 7Li
    ezspin_db(6)%AtNumb      = 3_iwp
    ezspin_db(6)%MassNum     = 7_iwp
    ezspin_db(6)%abundance    = 92.41_wp
    ezspin_db(6)%nucspin      = 1.5_wp
    ezspin_db(6)%gfactor      = +2.170951_wp
!----------------------------
!--> Record 007    = 6Li
    ezspin_db(7)%AtNumb      = 3_iwp
    ezspin_db(7)%MassNum     = 6_iwp
    ezspin_db(7)%abundance    = 7.59_wp
    ezspin_db(7)%nucspin      = 1.0_wp
    ezspin_db(7)%gfactor      = +0.8220473_wp
!----------------------------
!--> Record 008    = 9Be
    ezspin_db(8)%AtNumb      = 4_iwp
    ezspin_db(8)%MassNum     = 9_iwp
    ezspin_db(8)%abundance    = 100.0_wp
    ezspin_db(8)%nucspin      = 1.5_wp
    ezspin_db(8)%gfactor      = -0.78495_wp
!----------------------------
!--> Record 009    = 11B
    ezspin_db(9)%AtNumb      = 5_iwp
    ezspin_db(9)%MassNum     = 11_iwp
    ezspin_db(9)%abundance    = 80.1_wp
    ezspin_db(9)%nucspin      = 1.5_wp
    ezspin_db(9)%gfactor      = +1.7924326_wp
!----------------------------
!--> Record 010    = 10B
    ezspin_db(10)%AtNumb      = 5_iwp
    ezspin_db(10)%MassNum     = 10_iwp
    ezspin_db(10)%abundance    = 19.9_wp
    ezspin_db(10)%nucspin      = 3.0_wp
    ezspin_db(10)%gfactor      = +0.600215_wp
!----------------------------
!--> Record 011    = 12C
    ezspin_db(11)%AtNumb      = 6_iwp
    ezspin_db(11)%MassNum     = 12_iwp
    ezspin_db(11)%abundance    = 98.93_wp
    ezspin_db(11)%nucspin      = 0.0_wp
    ezspin_db(11)%gfactor      = 0.0_wp
!----------------------------
!--> Record 012    = 13C
    ezspin_db(12)%AtNumb      = 6_iwp
    ezspin_db(12)%MassNum     = 13_iwp
    ezspin_db(12)%abundance    = 1.07_wp
    ezspin_db(12)%nucspin      = 0.5_wp
    ezspin_db(12)%gfactor      = +1.4048236_wp
!----------------------------
!--> Record 013    = 14C
    ezspin_db(13)%AtNumb      = 6_iwp
    ezspin_db(13)%MassNum     = 14_iwp
    ezspin_db(13)%abundance    = 0.0_wp
    ezspin_db(13)%nucspin      = 0.0_wp
    ezspin_db(13)%gfactor      = 0.0_wp
!----------------------------
!--> Record 014    = 14N
    ezspin_db(14)%AtNumb      = 7_iwp
    ezspin_db(14)%MassNum     = 14_iwp
    ezspin_db(14)%abundance    = 99.632_wp
    ezspin_db(14)%nucspin      = 1.0_wp
    ezspin_db(14)%gfactor      = +0.40376100_wp
!----------------------------
!--> Record 015    = 15N
    ezspin_db(15)%AtNumb      = 7_iwp
    ezspin_db(15)%MassNum     = 15_iwp
    ezspin_db(15)%abundance    = 0.368_wp
    ezspin_db(15)%nucspin      = 0.5_wp
    ezspin_db(15)%gfactor      = -0.56637768_wp
!----------------------------
!--> Record 016    = 16O
    ezspin_db(16)%AtNumb      = 8_iwp
    ezspin_db(16)%MassNum     = 16_iwp
    ezspin_db(16)%abundance    = 99.757_wp
    ezspin_db(16)%nucspin      = 0.0_wp
    ezspin_db(16)%gfactor      = 0.0_wp
!----------------------------
!--> Record 017    = 18O
    ezspin_db(17)%AtNumb      = 8_iwp
    ezspin_db(17)%MassNum     = 18_iwp
    ezspin_db(17)%abundance    = 0.205_wp
    ezspin_db(17)%nucspin      = 0.0_wp
    ezspin_db(17)%gfactor      = 0.0_wp
!----------------------------
!--> Record 018    = 17O
    ezspin_db(18)%AtNumb      = 8_iwp
    ezspin_db(18)%MassNum     = 17_iwp
    ezspin_db(18)%abundance    = 0.038_wp
    ezspin_db(18)%nucspin      = 2.5_wp
    ezspin_db(18)%gfactor      = -0.757516_wp
!----------------------------
!--> Record 019    = 19F
    ezspin_db(19)%AtNumb      = 9_iwp
    ezspin_db(19)%MassNum     = 19_iwp
    ezspin_db(19)%abundance    = 100.0_wp
    ezspin_db(19)%nucspin      = 0.5_wp
    ezspin_db(19)%gfactor      = +5.257736_wp
!----------------------------
!--> Record 020    = 20Ne
    ezspin_db(20)%AtNumb      = 10_iwp
    ezspin_db(20)%MassNum     = 20_iwp
    ezspin_db(20)%abundance    = 90.48_wp
    ezspin_db(20)%nucspin      = 0.0_wp
    ezspin_db(20)%gfactor      = 0.0_wp
!----------------------------
!--> Record 021    = 22Ne
    ezspin_db(21)%AtNumb      = 10_iwp
    ezspin_db(21)%MassNum     = 22_iwp
    ezspin_db(21)%abundance    = 9.25_wp
    ezspin_db(21)%nucspin      = 0.0_wp
    ezspin_db(21)%gfactor      = 0.0_wp
!----------------------------
!--> Record 022    = 21Ne
    ezspin_db(22)%AtNumb      = 10_iwp
    ezspin_db(22)%MassNum     = 21_iwp
    ezspin_db(22)%abundance    = 0.27_wp
    ezspin_db(22)%nucspin      = 1.5_wp
    ezspin_db(22)%gfactor      = -0.441198_wp
!----------------------------
!--> Record 023    = 23Na
    ezspin_db(23)%AtNumb      = 11_iwp
    ezspin_db(23)%MassNum     = 23_iwp
    ezspin_db(23)%abundance    = 100.0_wp
    ezspin_db(23)%nucspin      = 1.5_wp
    ezspin_db(23)%gfactor      = +1.478348_wp
!----------------------------
!--> Record 024    = 22Na
    ezspin_db(24)%AtNumb      = 11_iwp
    ezspin_db(24)%MassNum     = 22_iwp
    ezspin_db(24)%abundance    = 0.0_wp
    ezspin_db(24)%nucspin      = 3.0_wp
    ezspin_db(24)%gfactor      = +0.582_wp
!----------------------------
!--> Record 025    = 24Mg
    ezspin_db(25)%AtNumb      = 12_iwp
    ezspin_db(25)%MassNum     = 24_iwp
    ezspin_db(25)%abundance    = 78.99_wp
    ezspin_db(25)%nucspin      = 0.0_wp
    ezspin_db(25)%gfactor      = 0.0_wp
!----------------------------
!--> Record 026    = 26Mg
    ezspin_db(26)%AtNumb      = 12_iwp
    ezspin_db(26)%MassNum     = 26_iwp
    ezspin_db(26)%abundance    = 11.01_wp
    ezspin_db(26)%nucspin      = 0.0_wp
    ezspin_db(26)%gfactor      = 0.0_wp
!----------------------------
!--> Record 027    = 25Mg
    ezspin_db(27)%AtNumb      = 12_iwp
    ezspin_db(27)%MassNum     = 25_iwp
    ezspin_db(27)%abundance    = 10.00_wp
    ezspin_db(27)%nucspin      = 2.5_wp
    ezspin_db(27)%gfactor      = -0.34218_wp
!----------------------------
!--> Record 028    = 27Al
    ezspin_db(28)%AtNumb      = 13_iwp
    ezspin_db(28)%MassNum     = 27_iwp
    ezspin_db(28)%abundance    = 100.0_wp
    ezspin_db(28)%nucspin      = 2.5_wp
    ezspin_db(28)%gfactor      = +1.4566028_wp
!----------------------------
!--> Record 029    = 28Si
    ezspin_db(29)%AtNumb      = 14_iwp
    ezspin_db(29)%MassNum     = 28_iwp
    ezspin_db(29)%abundance    = 92.2297_wp
    ezspin_db(29)%nucspin      = 0.0_wp
    ezspin_db(29)%gfactor      = 0.0_wp
!----------------------------
!--> Record 030    = 29Si
    ezspin_db(30)%AtNumb      = 14_iwp
    ezspin_db(30)%MassNum     = 29_iwp
    ezspin_db(30)%abundance    = 4.6832_wp
    ezspin_db(30)%nucspin      = 0.5_wp
    ezspin_db(30)%gfactor      = -1.11058_wp
!----------------------------
!--> Record 031    = 30Si
    ezspin_db(31)%AtNumb      = 14_iwp
    ezspin_db(31)%MassNum     = 30_iwp
    ezspin_db(31)%abundance    = 3.0872_wp
    ezspin_db(31)%nucspin      = 0.0_wp
    ezspin_db(31)%gfactor      = 0.0_wp
!----------------------------
!--> Record 032    = 31P
    ezspin_db(32)%AtNumb      = 15_iwp
    ezspin_db(32)%MassNum     = 31_iwp
    ezspin_db(32)%abundance    = 100.0_wp
    ezspin_db(32)%nucspin      = 0.5_wp
    ezspin_db(32)%gfactor      = +2.26320_wp
!----------------------------
!--> Record 033    = 32S
    ezspin_db(33)%AtNumb      = 16_iwp
    ezspin_db(33)%MassNum     = 32_iwp
    ezspin_db(33)%abundance    = 94.93_wp
    ezspin_db(33)%nucspin      = 0.0_wp
    ezspin_db(33)%gfactor      = 0.0_wp
!----------------------------
!--> Record 034    = 34S
    ezspin_db(34)%AtNumb      = 16_iwp
    ezspin_db(34)%MassNum     = 34_iwp
    ezspin_db(34)%abundance    = 4.29_wp
    ezspin_db(34)%nucspin      = 0.0_wp
    ezspin_db(34)%gfactor      = 0.0_wp
!----------------------------
!--> Record 035    = 33S
    ezspin_db(35)%AtNumb      = 16_iwp
    ezspin_db(35)%MassNum     = 33_iwp
    ezspin_db(35)%abundance    = 0.76_wp
    ezspin_db(35)%nucspin      = 1.5_wp
    ezspin_db(35)%gfactor      = +0.429214_wp
!----------------------------
!--> Record 036    = 36S
    ezspin_db(36)%AtNumb      = 16_iwp
    ezspin_db(36)%MassNum     = 36_iwp
    ezspin_db(36)%abundance    = 0.02_wp
    ezspin_db(36)%nucspin      = 0.0_wp
    ezspin_db(36)%gfactor      = 0.0_wp
!----------------------------
!--> Record 037    = 35Cl
    ezspin_db(37)%AtNumb      = 17_iwp
    ezspin_db(37)%MassNum     = 35_iwp
    ezspin_db(37)%abundance    = 75.78_wp
    ezspin_db(37)%nucspin      = 1.5_wp
    ezspin_db(37)%gfactor      = +0.5479162_wp
!----------------------------
!--> Record 038    = 37Cl
    ezspin_db(38)%AtNumb      = 17_iwp
    ezspin_db(38)%MassNum     = 37_iwp
    ezspin_db(38)%abundance    = 24.22_wp
    ezspin_db(38)%nucspin      = 1.5_wp
    ezspin_db(38)%gfactor      = +0.4560824_wp
!----------------------------
!--> Record 039    = 36Cl
    ezspin_db(39)%AtNumb      = 17_iwp
    ezspin_db(39)%MassNum     = 36_iwp
    ezspin_db(39)%abundance    = 0.0_wp
    ezspin_db(39)%nucspin      = 2.0_wp
    ezspin_db(39)%gfactor      = +0.642735_wp
!----------------------------
!--> Record 040    = 40Ar
    ezspin_db(40)%AtNumb      = 18_iwp
    ezspin_db(40)%MassNum     = 40_iwp
    ezspin_db(40)%abundance    = 99.6003_wp
    ezspin_db(40)%nucspin      = 0.0_wp
    ezspin_db(40)%gfactor      = 0.0_wp
!----------------------------
!--> Record 041    = 36Ar
    ezspin_db(41)%AtNumb      = 18_iwp
    ezspin_db(41)%MassNum     = 36_iwp
    ezspin_db(41)%abundance    = 0.3365_wp
    ezspin_db(41)%nucspin      = 0.0_wp
    ezspin_db(41)%gfactor      = 0.0_wp
!----------------------------
!--> Record 042    = 38Ar
    ezspin_db(42)%AtNumb      = 18_iwp
    ezspin_db(42)%MassNum     = 38_iwp
    ezspin_db(42)%abundance    = 0.0632_wp
    ezspin_db(42)%nucspin      = 0.0_wp
    ezspin_db(42)%gfactor      = 0.0_wp
!----------------------------
!--> Record 043    = 39Ar
    ezspin_db(43)%AtNumb      = 18_iwp
    ezspin_db(43)%MassNum     = 39_iwp
    ezspin_db(43)%abundance    = 0.0_wp
    ezspin_db(43)%nucspin      = 3.5_wp
    ezspin_db(43)%gfactor      = -0.4537_wp
!----------------------------
!--> Record 044    = 39K
    ezspin_db(44)%AtNumb      = 19_iwp
    ezspin_db(44)%MassNum     = 39_iwp
    ezspin_db(44)%abundance    = 93.2581_wp
    ezspin_db(44)%nucspin      = 1.5_wp
    ezspin_db(44)%gfactor      = +0.26098_wp
!----------------------------
!--> Record 045    = 41K
    ezspin_db(45)%AtNumb      = 19_iwp
    ezspin_db(45)%MassNum     = 41_iwp
    ezspin_db(45)%abundance    = 6.7302_wp
    ezspin_db(45)%nucspin      = 1.5_wp
    ezspin_db(45)%gfactor      = +0.1432467_wp
!----------------------------
!--> Record 046    = 40K
    ezspin_db(46)%AtNumb      = 19_iwp
    ezspin_db(46)%MassNum     = 40_iwp
    ezspin_db(46)%abundance    = 0.0117_wp
    ezspin_db(46)%nucspin      = 4.0_wp
    ezspin_db(46)%gfactor      = -0.324525_wp
!----------------------------
!--> Record 047    = 40Ca
    ezspin_db(47)%AtNumb      = 20_iwp
    ezspin_db(47)%MassNum     = 40_iwp
    ezspin_db(47)%abundance    = 96.941_wp
    ezspin_db(47)%nucspin      = 0.0_wp
    ezspin_db(47)%gfactor      = 0.0_wp
!----------------------------
!--> Record 048    = 44Ca
    ezspin_db(48)%AtNumb      = 20_iwp
    ezspin_db(48)%MassNum     = 44_iwp
    ezspin_db(48)%abundance    = 2.086_wp
    ezspin_db(48)%nucspin      = 0.0_wp
    ezspin_db(48)%gfactor      = 0.0_wp
!----------------------------
!--> Record 049    = 42Ca
    ezspin_db(49)%AtNumb      = 20_iwp
    ezspin_db(49)%MassNum     = 42_iwp
    ezspin_db(49)%abundance    = 0.647_wp
    ezspin_db(49)%nucspin      = 0.0_wp
    ezspin_db(49)%gfactor      = 0.0_wp
!----------------------------
!--> Record 050    = 48Ca
    ezspin_db(50)%AtNumb      = 20_iwp
    ezspin_db(50)%MassNum     = 48_iwp
    ezspin_db(50)%abundance    = 0.187_wp
    ezspin_db(50)%nucspin      = 0.0_wp
    ezspin_db(50)%gfactor      = 0.0_wp
!----------------------------
!--> Record 051    = 43Ca
    ezspin_db(51)%AtNumb      = 20_iwp
    ezspin_db(51)%MassNum     = 43_iwp
    ezspin_db(51)%abundance    = 0.135_wp
    ezspin_db(51)%nucspin      = 3.5_wp
    ezspin_db(51)%gfactor      = -0.37637_wp
!----------------------------
!--> Record 052    = 46Ca
    ezspin_db(52)%AtNumb      = 20_iwp
    ezspin_db(52)%MassNum     = 46_iwp
    ezspin_db(52)%abundance    = 0.004_wp
    ezspin_db(52)%nucspin      = 0.0_wp
    ezspin_db(52)%gfactor      = 0.0_wp
!----------------------------
!--> Record 053    = 41Ca
    ezspin_db(53)%AtNumb      = 20_iwp
    ezspin_db(53)%MassNum     = 41_iwp
    ezspin_db(53)%abundance    = 0.0_wp
    ezspin_db(53)%nucspin      = 3.5_wp
    ezspin_db(53)%gfactor      = -0.4556517_wp
!----------------------------
!--> Record 054    = 45Sc
    ezspin_db(54)%AtNumb      = 21_iwp
    ezspin_db(54)%MassNum     = 45_iwp
    ezspin_db(54)%abundance    = 100.0_wp
    ezspin_db(54)%nucspin      = 3.5_wp
    ezspin_db(54)%gfactor      = +1.35899_wp
!----------------------------
!--> Record 055    = 48Ti
    ezspin_db(55)%AtNumb      = 22_iwp
    ezspin_db(55)%MassNum     = 48_iwp
    ezspin_db(55)%abundance    = 73.72_wp
    ezspin_db(55)%nucspin      = 0.0_wp
    ezspin_db(55)%gfactor      = 0.0_wp
!----------------------------
!--> Record 056    = 46Ti
    ezspin_db(56)%AtNumb      = 22_iwp
    ezspin_db(56)%MassNum     = 46_iwp
    ezspin_db(56)%abundance    = 8.25_wp
    ezspin_db(56)%nucspin      = 0.0_wp
    ezspin_db(56)%gfactor      = 0.0_wp
!----------------------------
!--> Record 057    = 47Ti
    ezspin_db(57)%AtNumb      = 22_iwp
    ezspin_db(57)%MassNum     = 47_iwp
    ezspin_db(57)%abundance    = 7.44_wp
    ezspin_db(57)%nucspin      = 2.5_wp
    ezspin_db(57)%gfactor      = -0.31539_wp
!----------------------------
!--> Record 058    = 49Ti
    ezspin_db(58)%AtNumb      = 22_iwp
    ezspin_db(58)%MassNum     = 49_iwp
    ezspin_db(58)%abundance    = 5.41_wp
    ezspin_db(58)%nucspin      = 3.5_wp
    ezspin_db(58)%gfactor      = -0.315477_wp
!----------------------------
!--> Record 059    = 50Ti
    ezspin_db(59)%AtNumb      = 22_iwp
    ezspin_db(59)%MassNum     = 50_iwp
    ezspin_db(59)%abundance    = 5.18_wp
    ezspin_db(59)%nucspin      = 0.0_wp
    ezspin_db(59)%gfactor      = 0.0_wp
!----------------------------
!--> Record 060    = 51V
    ezspin_db(60)%AtNumb      = 23_iwp
    ezspin_db(60)%MassNum     = 51_iwp
    ezspin_db(60)%abundance    = 99.75_wp
    ezspin_db(60)%nucspin      = 3.5_wp
    ezspin_db(60)%gfactor      = +1.47106_wp
!----------------------------
!--> Record 061    = 50V
    ezspin_db(61)%AtNumb      = 23_iwp
    ezspin_db(61)%MassNum     = 50_iwp
    ezspin_db(61)%abundance    = 0.25_wp
    ezspin_db(61)%nucspin      = 6.0_wp
    ezspin_db(61)%gfactor      = +0.5576148_wp
!----------------------------
!--> Record 062    = 52Cr
    ezspin_db(62)%AtNumb      = 24_iwp
    ezspin_db(62)%MassNum     = 52_iwp
    ezspin_db(62)%abundance    = 83.789_wp
    ezspin_db(62)%nucspin      = 0.0_wp
    ezspin_db(62)%gfactor      = 0.0_wp
!----------------------------
!--> Record 063    = 53Cr
    ezspin_db(63)%AtNumb      = 24_iwp
    ezspin_db(63)%MassNum     = 53_iwp
    ezspin_db(63)%abundance    = 9.501_wp
    ezspin_db(63)%nucspin      = 1.5_wp
    ezspin_db(63)%gfactor      = -0.31636_wp
!----------------------------
!--> Record 064    = 50Cr
    ezspin_db(64)%AtNumb      = 24_iwp
    ezspin_db(64)%MassNum     = 50_iwp
    ezspin_db(64)%abundance    = 4.345_wp
    ezspin_db(64)%nucspin      = 0.0_wp
    ezspin_db(64)%gfactor      = 0.0_wp
!----------------------------
!--> Record 065    = 54Cr
    ezspin_db(65)%AtNumb      = 24_iwp
    ezspin_db(65)%MassNum     = 54_iwp
    ezspin_db(65)%abundance    = 2.365_wp
    ezspin_db(65)%nucspin      = 0.0_wp
    ezspin_db(65)%gfactor      = 0.0_wp
!----------------------------
!--> Record 066    = 55Mn
    ezspin_db(66)%AtNumb      = 25_iwp
    ezspin_db(66)%MassNum     = 55_iwp
    ezspin_db(66)%abundance    = 100.0_wp
    ezspin_db(66)%nucspin      = 2.5_wp
    ezspin_db(66)%gfactor      = +1.3813_wp
!----------------------------
!--> Record 067    = 53Mn
    ezspin_db(67)%AtNumb      = 25_iwp
    ezspin_db(67)%MassNum     = 53_iwp
    ezspin_db(67)%abundance    = 0.0_wp
    ezspin_db(67)%nucspin      = 3.5_wp
    ezspin_db(67)%gfactor      = +1.439_wp
!----------------------------
!--> Record 068    = 56Fe
    ezspin_db(68)%AtNumb      = 26_iwp
    ezspin_db(68)%MassNum     = 56_iwp
    ezspin_db(68)%abundance    = 91.754_wp
    ezspin_db(68)%nucspin      = 0.0_wp
    ezspin_db(68)%gfactor      = 0.0_wp
!----------------------------
!--> Record 069    = 54Fe
    ezspin_db(69)%AtNumb      = 26_iwp
    ezspin_db(69)%MassNum     = 54_iwp
    ezspin_db(69)%abundance    = 5.845_wp
    ezspin_db(69)%nucspin      = 0.0_wp
    ezspin_db(69)%gfactor      = 0.0_wp
!----------------------------
!--> Record 070    = 57Fe
    ezspin_db(70)%AtNumb      = 26_iwp
    ezspin_db(70)%MassNum     = 57_iwp
    ezspin_db(70)%abundance    = 2.119_wp
    ezspin_db(70)%nucspin      = 0.5_wp
    ezspin_db(70)%gfactor      = +0.1809_wp
!----------------------------
!--> Record 071    = 58Fe
    ezspin_db(71)%AtNumb      = 26_iwp
    ezspin_db(71)%MassNum     = 58_iwp
    ezspin_db(71)%abundance    = 0.282_wp
    ezspin_db(71)%nucspin      = 0.0_wp
    ezspin_db(71)%gfactor      = 0.0_wp
!----------------------------
!--> Record 072    = 59Co
    ezspin_db(72)%AtNumb      = 27_iwp
    ezspin_db(72)%MassNum     = 59_iwp
    ezspin_db(72)%abundance    = 100.0_wp
    ezspin_db(72)%nucspin      = 3.5_wp
    ezspin_db(72)%gfactor      = +1.322_wp
!----------------------------
!--> Record 073    = 60Co
    ezspin_db(73)%AtNumb      = 27_iwp
    ezspin_db(73)%MassNum     = 60_iwp
    ezspin_db(73)%abundance    = 0.0_wp
    ezspin_db(73)%nucspin      = 5.0_wp
    ezspin_db(73)%gfactor      = +0.7598_wp
!----------------------------
!--> Record 074    = 58Ni
    ezspin_db(74)%AtNumb      = 28_iwp
    ezspin_db(74)%MassNum     = 58_iwp
    ezspin_db(74)%abundance    = 68.0769_wp
    ezspin_db(74)%nucspin      = 0.0_wp
    ezspin_db(74)%gfactor      = 0.0_wp
!----------------------------
!--> Record 075    = 60Ni
    ezspin_db(75)%AtNumb      = 28_iwp
    ezspin_db(75)%MassNum     = 60_iwp
    ezspin_db(75)%abundance    = 26.2231_wp
    ezspin_db(75)%nucspin      = 0.0_wp
    ezspin_db(75)%gfactor      = 0.0_wp
!----------------------------
!--> Record 076    = 62Ni
    ezspin_db(76)%AtNumb      = 28_iwp
    ezspin_db(76)%MassNum     = 62_iwp
    ezspin_db(76)%abundance    = 3.6345_wp
    ezspin_db(76)%nucspin      = 0.0_wp
    ezspin_db(76)%gfactor      = 0.0_wp
!----------------------------
!--> Record 077    = 61Ni
    ezspin_db(77)%AtNumb      = 28_iwp
    ezspin_db(77)%MassNum     = 61_iwp
    ezspin_db(77)%abundance    = 1.1399_wp
    ezspin_db(77)%nucspin      = 1.5_wp
    ezspin_db(77)%gfactor      = -0.50001_wp
!----------------------------
!--> Record 078    = 64Ni
    ezspin_db(78)%AtNumb      = 28_iwp
    ezspin_db(78)%MassNum     = 64_iwp
    ezspin_db(78)%abundance    = 0.9256_wp
    ezspin_db(78)%nucspin      = 0.0_wp
    ezspin_db(78)%gfactor      = 0.0_wp
!----------------------------
!--> Record 079    = 63Cu
    ezspin_db(79)%AtNumb      = 29_iwp
    ezspin_db(79)%MassNum     = 63_iwp
    ezspin_db(79)%abundance    = 69.17_wp
    ezspin_db(79)%nucspin      = 1.5_wp
    ezspin_db(79)%gfactor      = +1.4824_wp
!----------------------------
!--> Record 080    = 65Cu
    ezspin_db(80)%AtNumb      = 29_iwp
    ezspin_db(80)%MassNum     = 65_iwp
    ezspin_db(80)%abundance    = 30.83_wp
    ezspin_db(80)%nucspin      = 1.5_wp
    ezspin_db(80)%gfactor      = +1.5878_wp
!----------------------------
!--> Record 081    = 64Zn
    ezspin_db(81)%AtNumb      = 30_iwp
    ezspin_db(81)%MassNum     = 64_iwp
    ezspin_db(81)%abundance    = 48.63_wp
    ezspin_db(81)%nucspin      = 0.0_wp
    ezspin_db(81)%gfactor      = 0.0_wp
!----------------------------
!--> Record 082    = 66Zn
    ezspin_db(82)%AtNumb      = 30_iwp
    ezspin_db(82)%MassNum     = 66_iwp
    ezspin_db(82)%abundance    = 27.90_wp
    ezspin_db(82)%nucspin      = 0.0_wp
    ezspin_db(82)%gfactor      = 0.0_wp
!----------------------------
!--> Record 083    = 68Zn
    ezspin_db(83)%AtNumb      = 30_iwp
    ezspin_db(83)%MassNum     = 68_iwp
    ezspin_db(83)%abundance    = 18.75_wp
    ezspin_db(83)%nucspin      = 0.0_wp
    ezspin_db(83)%gfactor      = 0.0_wp
!----------------------------
!--> Record 084    = 67Zn
    ezspin_db(84)%AtNumb      = 30_iwp
    ezspin_db(84)%MassNum     = 67_iwp
    ezspin_db(84)%abundance    = 4.10_wp
    ezspin_db(84)%nucspin      = 2.5_wp
    ezspin_db(84)%gfactor      = +0.350192_wp
!----------------------------
!--> Record 085    = 70Zn
    ezspin_db(85)%AtNumb      = 30_iwp
    ezspin_db(85)%MassNum     = 70_iwp
    ezspin_db(85)%abundance    = 0.62_wp
    ezspin_db(85)%nucspin      = 0.0_wp
    ezspin_db(85)%gfactor      = 0.0_wp
!----------------------------
!--> Record 086    = 69Ga
    ezspin_db(86)%AtNumb      = 31_iwp
    ezspin_db(86)%MassNum     = 69_iwp
    ezspin_db(86)%abundance    = 60.108_wp
    ezspin_db(86)%nucspin      = 1.5_wp
    ezspin_db(86)%gfactor      = +1.34439_wp
!----------------------------
!--> Record 087    = 71Ga
    ezspin_db(87)%AtNumb      = 31_iwp
    ezspin_db(87)%MassNum     = 71_iwp
    ezspin_db(87)%abundance    = 39.892_wp
    ezspin_db(87)%nucspin      = 1.5_wp
    ezspin_db(87)%gfactor      = +1.70818_wp
!----------------------------
!--> Record 088    = 74Ge
    ezspin_db(88)%AtNumb      = 32_iwp
    ezspin_db(88)%MassNum     = 74_iwp
    ezspin_db(88)%abundance    = 36.28_wp
    ezspin_db(88)%nucspin      = 0.0_wp
    ezspin_db(88)%gfactor      = 0.0_wp
!----------------------------
!--> Record 089    = 72Ge
    ezspin_db(89)%AtNumb      = 32_iwp
    ezspin_db(89)%MassNum     = 72_iwp
    ezspin_db(89)%abundance    = 27.54_wp
    ezspin_db(89)%nucspin      = 0.0_wp
    ezspin_db(89)%gfactor      = 0.0_wp
!----------------------------
!--> Record 090    = 70Ge
    ezspin_db(90)%AtNumb      = 32_iwp
    ezspin_db(90)%MassNum     = 70_iwp
    ezspin_db(90)%abundance    = 20.84_wp
    ezspin_db(90)%nucspin      = 0.0_wp
    ezspin_db(90)%gfactor      = 0.0_wp
!----------------------------
!--> Record 091    = 73Ge
    ezspin_db(91)%AtNumb      = 32_iwp
    ezspin_db(91)%MassNum     = 73_iwp
    ezspin_db(91)%abundance    = 7.73_wp
    ezspin_db(91)%nucspin      = 4.5_wp
    ezspin_db(91)%gfactor      = -0.1954373_wp
!----------------------------
!--> Record 092    = 76Ge
    ezspin_db(92)%AtNumb      = 32_iwp
    ezspin_db(92)%MassNum     = 76_iwp
    ezspin_db(92)%abundance    = 7.61_wp
    ezspin_db(92)%nucspin      = 0.0_wp
    ezspin_db(92)%gfactor      = 0.0_wp
!----------------------------
!--> Record 093    = 75As
    ezspin_db(93)%AtNumb      = 33_iwp
    ezspin_db(93)%MassNum     = 75_iwp
    ezspin_db(93)%abundance    = 100.0_wp
    ezspin_db(93)%nucspin      = 1.5_wp
    ezspin_db(93)%gfactor      = +0.95965_wp
!----------------------------
!--> Record 094    = 80Se
    ezspin_db(94)%AtNumb      = 34_iwp
    ezspin_db(94)%MassNum     = 80_iwp
    ezspin_db(94)%abundance    = 49.61_wp
    ezspin_db(94)%nucspin      = 0.0_wp
    ezspin_db(94)%gfactor      = 0.0_wp
!----------------------------
!--> Record 095    = 78Se
    ezspin_db(95)%AtNumb      = 34_iwp
    ezspin_db(95)%MassNum     = 78_iwp
    ezspin_db(95)%abundance    = 23.77_wp
    ezspin_db(95)%nucspin      = 0.0_wp
    ezspin_db(95)%gfactor      = 0.0_wp
!----------------------------
!--> Record 096    = 76Se
    ezspin_db(96)%AtNumb      = 34_iwp
    ezspin_db(96)%MassNum     = 76_iwp
    ezspin_db(96)%abundance    = 9.37_wp
    ezspin_db(96)%nucspin      = 0.0_wp
    ezspin_db(96)%gfactor      = 0.0_wp
!----------------------------
!--> Record 097    = 82Se
    ezspin_db(97)%AtNumb      = 34_iwp
    ezspin_db(97)%MassNum     = 82_iwp
    ezspin_db(97)%abundance    = 8.73_wp
    ezspin_db(97)%nucspin      = 0.0_wp
    ezspin_db(97)%gfactor      = 0.0_wp
!----------------------------
!--> Record 098    = 77Se
    ezspin_db(98)%AtNumb      = 34_iwp
    ezspin_db(98)%MassNum     = 77_iwp
    ezspin_db(98)%abundance    = 7.63_wp
    ezspin_db(98)%nucspin      = 0.5_wp
    ezspin_db(98)%gfactor      = +1.07008_wp
!----------------------------
!--> Record 099    = 74Se
    ezspin_db(99)%AtNumb      = 34_iwp
    ezspin_db(99)%MassNum     = 74_iwp
    ezspin_db(99)%abundance    = 0.89_wp
    ezspin_db(99)%nucspin      = 0.0_wp
    ezspin_db(99)%gfactor      = 0.0_wp
!----------------------------
!--> Record 100    = 79Se
    ezspin_db(100)%AtNumb      = 34_iwp
    ezspin_db(100)%MassNum     = 79_iwp
    ezspin_db(100)%abundance    = 0.0_wp
    ezspin_db(100)%nucspin      = 3.5_wp
    ezspin_db(100)%gfactor      = -0.29_wp
!----------------------------
!--> Record 101    = 79Br
    ezspin_db(101)%AtNumb      = 35_iwp
    ezspin_db(101)%MassNum     = 79_iwp
    ezspin_db(101)%abundance    = 50.69_wp
    ezspin_db(101)%nucspin      = 1.5_wp
    ezspin_db(101)%gfactor      = +1.404267_wp
!----------------------------
!--> Record 102    = 81Br
    ezspin_db(102)%AtNumb      = 35_iwp
    ezspin_db(102)%MassNum     = 81_iwp
    ezspin_db(102)%abundance    = 49.31_wp
    ezspin_db(102)%nucspin      = 1.5_wp
    ezspin_db(102)%gfactor      = +1.513708_wp
!----------------------------
!--> Record 103    = 84Kr
    ezspin_db(103)%AtNumb      = 36_iwp
    ezspin_db(103)%MassNum     = 84_iwp
    ezspin_db(103)%abundance    = 57.00_wp
    ezspin_db(103)%nucspin      = 0.0_wp
    ezspin_db(103)%gfactor      = 0.0_wp
!----------------------------
!--> Record 104    = 86Kr
    ezspin_db(104)%AtNumb      = 36_iwp
    ezspin_db(104)%MassNum     = 86_iwp
    ezspin_db(104)%abundance    = 17.30_wp
    ezspin_db(104)%nucspin      = 0.0_wp
    ezspin_db(104)%gfactor      = 0.0_wp
!----------------------------
!--> Record 105    = 82Kr
    ezspin_db(105)%AtNumb      = 36_iwp
    ezspin_db(105)%MassNum     = 82_iwp
    ezspin_db(105)%abundance    = 11.58_wp
    ezspin_db(105)%nucspin      = 0.0_wp
    ezspin_db(105)%gfactor      = 0.0_wp
!----------------------------
!--> Record 106    = 83Kr
    ezspin_db(106)%AtNumb      = 36_iwp
    ezspin_db(106)%MassNum     = 83_iwp
    ezspin_db(106)%abundance    = 11.49_wp
    ezspin_db(106)%nucspin      = 4.5_wp
    ezspin_db(106)%gfactor      = -0.215704_wp
!----------------------------
!--> Record 107    = 80Kr
    ezspin_db(107)%AtNumb      = 36_iwp
    ezspin_db(107)%MassNum     = 80_iwp
    ezspin_db(107)%abundance    = 2.28_wp
    ezspin_db(107)%nucspin      = 0.0_wp
    ezspin_db(107)%gfactor      = 0.0_wp
!----------------------------
!--> Record 108    = 78Kr
    ezspin_db(108)%AtNumb      = 36_iwp
    ezspin_db(108)%MassNum     = 78_iwp
    ezspin_db(108)%abundance    = 0.35_wp
    ezspin_db(108)%nucspin      = 0.0_wp
    ezspin_db(108)%gfactor      = 0.0_wp
!----------------------------
!--> Record 109    = 85Kr
    ezspin_db(109)%AtNumb      = 36_iwp
    ezspin_db(109)%MassNum     = 85_iwp
    ezspin_db(109)%abundance    = 0.0_wp
    ezspin_db(109)%nucspin      = 4.5_wp
    ezspin_db(109)%gfactor      = -0.2233_wp
!----------------------------
!--> Record 110    = 85Rb
    ezspin_db(110)%AtNumb      = 37_iwp
    ezspin_db(110)%MassNum     = 85_iwp
    ezspin_db(110)%abundance    = 72.17_wp
    ezspin_db(110)%nucspin      = 2.5_wp
    ezspin_db(110)%gfactor      = +0.541192_wp
!----------------------------
!--> Record 111    = 87Rb
    ezspin_db(111)%AtNumb      = 37_iwp
    ezspin_db(111)%MassNum     = 87_iwp
    ezspin_db(111)%abundance    = 27.83_wp
    ezspin_db(111)%nucspin      = 1.5_wp
    ezspin_db(111)%gfactor      = +1.83421_wp
!----------------------------
!--> Record 112    = 88Sr
    ezspin_db(112)%AtNumb      = 38_iwp
    ezspin_db(112)%MassNum     = 88_iwp
    ezspin_db(112)%abundance    = 82.58_wp
    ezspin_db(112)%nucspin      = 0.0_wp
    ezspin_db(112)%gfactor      = 0.0_wp
!----------------------------
!--> Record 113    = 86Sr
    ezspin_db(113)%AtNumb      = 38_iwp
    ezspin_db(113)%MassNum     = 86_iwp
    ezspin_db(113)%abundance    = 9.86_wp
    ezspin_db(113)%nucspin      = 0.0_wp
    ezspin_db(113)%gfactor      = 0.0_wp
!----------------------------
!--> Record 114    = 87Sr
    ezspin_db(114)%AtNumb      = 38_iwp
    ezspin_db(114)%MassNum     = 87_iwp
    ezspin_db(114)%abundance    = 7.00_wp
    ezspin_db(114)%nucspin      = 4.5_wp
    ezspin_db(114)%gfactor      = -0.24284_wp
!----------------------------
!--> Record 115    = 84Sr
    ezspin_db(115)%AtNumb      = 38_iwp
    ezspin_db(115)%MassNum     = 84_iwp
    ezspin_db(115)%abundance    = 0.56_wp
    ezspin_db(115)%nucspin      = 0.0_wp
    ezspin_db(115)%gfactor      = 0.0_wp
!----------------------------
!--> Record 116    = 89Y
    ezspin_db(116)%AtNumb      = 39_iwp
    ezspin_db(116)%MassNum     = 89_iwp
    ezspin_db(116)%abundance    = 100.0_wp
    ezspin_db(116)%nucspin      = 0.5_wp
    ezspin_db(116)%gfactor      = -0.2748308_wp
!----------------------------
!--> Record 117    = 90Zr
    ezspin_db(117)%AtNumb      = 40_iwp
    ezspin_db(117)%MassNum     = 90_iwp
    ezspin_db(117)%abundance    = 51.45_wp
    ezspin_db(117)%nucspin      = 0.0_wp
    ezspin_db(117)%gfactor      = 0.0_wp
!----------------------------
!--> Record 118    = 94Zr
    ezspin_db(118)%AtNumb      = 40_iwp
    ezspin_db(118)%MassNum     = 94_iwp
    ezspin_db(118)%abundance    = 17.38_wp
    ezspin_db(118)%nucspin      = 0.0_wp
    ezspin_db(118)%gfactor      = 0.0_wp
!----------------------------
!--> Record 119    = 92Zr
    ezspin_db(119)%AtNumb      = 40_iwp
    ezspin_db(119)%MassNum     = 92_iwp
    ezspin_db(119)%abundance    = 17.15_wp
    ezspin_db(119)%nucspin      = 0.0_wp
    ezspin_db(119)%gfactor      = 0.0_wp
!----------------------------
!--> Record 120    = 91Zr
    ezspin_db(120)%AtNumb      = 40_iwp
    ezspin_db(120)%MassNum     = 91_iwp
    ezspin_db(120)%abundance    = 11.22_wp
    ezspin_db(120)%nucspin      = 2.5_wp
    ezspin_db(120)%gfactor      = -0.521448_wp
!----------------------------
!--> Record 121    = 96Zr
    ezspin_db(121)%AtNumb      = 40_iwp
    ezspin_db(121)%MassNum     = 96_iwp
    ezspin_db(121)%abundance    = 2.80_wp
    ezspin_db(121)%nucspin      = 0.0_wp
    ezspin_db(121)%gfactor      = 0.0_wp
!----------------------------
!--> Record 122    = 93Nb
    ezspin_db(122)%AtNumb      = 41_iwp
    ezspin_db(122)%MassNum     = 93_iwp
    ezspin_db(122)%abundance    = 100.0_wp
    ezspin_db(122)%nucspin      = 4.5_wp
    ezspin_db(122)%gfactor      = +1.3712_wp
!----------------------------
!--> Record 123    = 98Mo
    ezspin_db(123)%AtNumb      = 42_iwp
    ezspin_db(123)%MassNum     = 98_iwp
    ezspin_db(123)%abundance    = 24.13_wp
    ezspin_db(123)%nucspin      = 0.0_wp
    ezspin_db(123)%gfactor      = 0.0_wp
!----------------------------
!--> Record 124    = 96Mo
    ezspin_db(124)%AtNumb      = 42_iwp
    ezspin_db(124)%MassNum     = 96_iwp
    ezspin_db(124)%abundance    = 16.68_wp
    ezspin_db(124)%nucspin      = 0.0_wp
    ezspin_db(124)%gfactor      = 0.0_wp
!----------------------------
!--> Record 125    = 95Mo
    ezspin_db(125)%AtNumb      = 42_iwp
    ezspin_db(125)%MassNum     = 95_iwp
    ezspin_db(125)%abundance    = 15.92_wp
    ezspin_db(125)%nucspin      = 2.5_wp
    ezspin_db(125)%gfactor      = -0.3657_wp
!----------------------------
!--> Record 126    = 92Mo
    ezspin_db(126)%AtNumb      = 42_iwp
    ezspin_db(126)%MassNum     = 92_iwp
    ezspin_db(126)%abundance    = 14.84_wp
    ezspin_db(126)%nucspin      = 0.0_wp
    ezspin_db(126)%gfactor      = 0.0_wp
!----------------------------
!--> Record 127    = 100Mo
    ezspin_db(127)%AtNumb      = 42_iwp
    ezspin_db(127)%MassNum     = 100_iwp
    ezspin_db(127)%abundance    = 9.63_wp
    ezspin_db(127)%nucspin      = 0.0_wp
    ezspin_db(127)%gfactor      = 0.0_wp
!----------------------------
!--> Record 128    = 97Mo
    ezspin_db(128)%AtNumb      = 42_iwp
    ezspin_db(128)%MassNum     = 97_iwp
    ezspin_db(128)%abundance    = 9.55_wp
    ezspin_db(128)%nucspin      = 2.5_wp
    ezspin_db(128)%gfactor      = -0.3734_wp
!----------------------------
!--> Record 129    = 94Mo
    ezspin_db(129)%AtNumb      = 42_iwp
    ezspin_db(129)%MassNum     = 94_iwp
    ezspin_db(129)%abundance    = 9.25_wp
    ezspin_db(129)%nucspin      = 0.0_wp
    ezspin_db(129)%gfactor      = 0.0_wp
!----------------------------
!--> Record 130    = 99Tc
    ezspin_db(130)%AtNumb      = 43_iwp
    ezspin_db(130)%MassNum     = 99_iwp
    ezspin_db(130)%abundance    = 0.0_wp
    ezspin_db(130)%nucspin      = 4.5_wp
    ezspin_db(130)%gfactor      = +1.2632_wp
!----------------------------
!--> Record 131    = 102Ru
    ezspin_db(131)%AtNumb      = 44_iwp
    ezspin_db(131)%MassNum     = 102_iwp
    ezspin_db(131)%abundance    = 31.55_wp
    ezspin_db(131)%nucspin      = 0.0_wp
    ezspin_db(131)%gfactor      = 0.0_wp
!----------------------------
!--> Record 132    = 104Ru
    ezspin_db(132)%AtNumb      = 44_iwp
    ezspin_db(132)%MassNum     = 104_iwp
    ezspin_db(132)%abundance    = 18.62_wp
    ezspin_db(132)%nucspin      = 0.0_wp
    ezspin_db(132)%gfactor      = 0.0_wp
!----------------------------
!--> Record 133    = 101Ru
    ezspin_db(133)%AtNumb      = 44_iwp
    ezspin_db(133)%MassNum     = 101_iwp
    ezspin_db(133)%abundance    = 17.06_wp
    ezspin_db(133)%nucspin      = 2.5_wp
    ezspin_db(133)%gfactor      = -0.288_wp
!----------------------------
!--> Record 134    = 99Ru
    ezspin_db(134)%AtNumb      = 44_iwp
    ezspin_db(134)%MassNum     = 99_iwp
    ezspin_db(134)%abundance    = 12.76_wp
    ezspin_db(134)%nucspin      = 2.5_wp
    ezspin_db(134)%gfactor      = -0.256_wp
!----------------------------
!--> Record 135    = 100Ru
    ezspin_db(135)%AtNumb      = 44_iwp
    ezspin_db(135)%MassNum     = 100_iwp
    ezspin_db(135)%abundance    = 12.60_wp
    ezspin_db(135)%nucspin      = 0.0_wp
    ezspin_db(135)%gfactor      = 0.0_wp
!----------------------------
!--> Record 136    = 96Ru
    ezspin_db(136)%AtNumb      = 44_iwp
    ezspin_db(136)%MassNum     = 96_iwp
    ezspin_db(136)%abundance    = 5.54_wp
    ezspin_db(136)%nucspin      = 0.0_wp
    ezspin_db(136)%gfactor      = 0.0_wp
!----------------------------
!--> Record 137    = 98Ru
    ezspin_db(137)%AtNumb      = 44_iwp
    ezspin_db(137)%MassNum     = 98_iwp
    ezspin_db(137)%abundance    = 1.87_wp
    ezspin_db(137)%nucspin      = 0.0_wp
    ezspin_db(137)%gfactor      = 0.0_wp
!----------------------------
!--> Record 138    = 103Rh
    ezspin_db(138)%AtNumb      = 45_iwp
    ezspin_db(138)%MassNum     = 103_iwp
    ezspin_db(138)%abundance    = 100.0_wp
    ezspin_db(138)%nucspin      = 0.5_wp
    ezspin_db(138)%gfactor      = -0.1768_wp
!----------------------------
!--> Record 139    = 106Pd
    ezspin_db(139)%AtNumb      = 46_iwp
    ezspin_db(139)%MassNum     = 106_iwp
    ezspin_db(139)%abundance    = 27.33_wp
    ezspin_db(139)%nucspin      = 0.0_wp
    ezspin_db(139)%gfactor      = 0.0_wp
!----------------------------
!--> Record 140    = 108Pd
    ezspin_db(140)%AtNumb      = 46_iwp
    ezspin_db(140)%MassNum     = 108_iwp
    ezspin_db(140)%abundance    = 26.46_wp
    ezspin_db(140)%nucspin      = 0.0_wp
    ezspin_db(140)%gfactor      = 0.0_wp
!----------------------------
!--> Record 141    = 105Pd
    ezspin_db(141)%AtNumb      = 46_iwp
    ezspin_db(141)%MassNum     = 105_iwp
    ezspin_db(141)%abundance    = 22.33_wp
    ezspin_db(141)%nucspin      = 2.5_wp
    ezspin_db(141)%gfactor      = -0.257_wp
!----------------------------
!--> Record 142    = 110Pd
    ezspin_db(142)%AtNumb      = 46_iwp
    ezspin_db(142)%MassNum     = 110_iwp
    ezspin_db(142)%abundance    = 11.72_wp
    ezspin_db(142)%nucspin      = 0.0_wp
    ezspin_db(142)%gfactor      = 0.0_wp
!----------------------------
!--> Record 143    = 104Pd
    ezspin_db(143)%AtNumb      = 46_iwp
    ezspin_db(143)%MassNum     = 104_iwp
    ezspin_db(143)%abundance    = 11.14_wp
    ezspin_db(143)%nucspin      = 0.0_wp
    ezspin_db(143)%gfactor      = 0.0_wp
!----------------------------
!--> Record 144    = 102Pd
    ezspin_db(144)%AtNumb      = 46_iwp
    ezspin_db(144)%MassNum     = 102_iwp
    ezspin_db(144)%abundance    = 1.02_wp
    ezspin_db(144)%nucspin      = 0.0_wp
    ezspin_db(144)%gfactor      = 0.0_wp
!----------------------------
!--> Record 145    = 107Ag
    ezspin_db(145)%AtNumb      = 47_iwp
    ezspin_db(145)%MassNum     = 107_iwp
    ezspin_db(145)%abundance    = 51.839_wp
    ezspin_db(145)%nucspin      = 0.5_wp
    ezspin_db(145)%gfactor      = -0.22714_wp
!----------------------------
!--> Record 146    = 109Ag
    ezspin_db(146)%AtNumb      = 47_iwp
    ezspin_db(146)%MassNum     = 109_iwp
    ezspin_db(146)%abundance    = 48.161_wp
    ezspin_db(146)%nucspin      = 0.5_wp
    ezspin_db(146)%gfactor      = -0.26112_wp
!----------------------------
!--> Record 147    = 114Cd
    ezspin_db(147)%AtNumb      = 48_iwp
    ezspin_db(147)%MassNum     = 114_iwp
    ezspin_db(147)%abundance    = 28.73_wp
    ezspin_db(147)%nucspin      = 0.0_wp
    ezspin_db(147)%gfactor      = 0.0_wp
!----------------------------
!--> Record 148    = 112Cd
    ezspin_db(148)%AtNumb      = 48_iwp
    ezspin_db(148)%MassNum     = 112_iwp
    ezspin_db(148)%abundance    = 24.13_wp
    ezspin_db(148)%nucspin      = 0.0_wp
    ezspin_db(148)%gfactor      = 0.0_wp
!----------------------------
!--> Record 149    = 111Cd
    ezspin_db(149)%AtNumb      = 48_iwp
    ezspin_db(149)%MassNum     = 111_iwp
    ezspin_db(149)%abundance    = 12.80_wp
    ezspin_db(149)%nucspin      = 0.5_wp
    ezspin_db(149)%gfactor      = -1.18977_wp
!----------------------------
!--> Record 150    = 110Cd
    ezspin_db(150)%AtNumb      = 48_iwp
    ezspin_db(150)%MassNum     = 110_iwp
    ezspin_db(150)%abundance    = 12.49_wp
    ezspin_db(150)%nucspin      = 0.0_wp
    ezspin_db(150)%gfactor      = 0.0_wp
!----------------------------
!--> Record 151    = 113Cd
    ezspin_db(151)%AtNumb      = 48_iwp
    ezspin_db(151)%MassNum     = 113_iwp
    ezspin_db(151)%abundance    = 12.22_wp
    ezspin_db(151)%nucspin      = 0.5_wp
    ezspin_db(151)%gfactor      = -1.244602_wp
!----------------------------
!--> Record 152    = 116Cd
    ezspin_db(152)%AtNumb      = 48_iwp
    ezspin_db(152)%MassNum     = 116_iwp
    ezspin_db(152)%abundance    = 7.49_wp
    ezspin_db(152)%nucspin      = 0.0_wp
    ezspin_db(152)%gfactor      = 0.0_wp
!----------------------------
!--> Record 153    = 106Cd
    ezspin_db(153)%AtNumb      = 48_iwp
    ezspin_db(153)%MassNum     = 106_iwp
    ezspin_db(153)%abundance    = 1.25_wp
    ezspin_db(153)%nucspin      = 0.0_wp
    ezspin_db(153)%gfactor      = 0.0_wp
!----------------------------
!--> Record 154    = 108Cd
    ezspin_db(154)%AtNumb      = 48_iwp
    ezspin_db(154)%MassNum     = 108_iwp
    ezspin_db(154)%abundance    = 0.89_wp
    ezspin_db(154)%nucspin      = 0.0_wp
    ezspin_db(154)%gfactor      = 0.0_wp
!----------------------------
!--> Record 155    = 115In
    ezspin_db(155)%AtNumb      = 49_iwp
    ezspin_db(155)%MassNum     = 115_iwp
    ezspin_db(155)%abundance    = 95.71_wp
    ezspin_db(155)%nucspin      = 4.5_wp
    ezspin_db(155)%gfactor      = +1.2313_wp
!----------------------------
!--> Record 156    = 113In
    ezspin_db(156)%AtNumb      = 49_iwp
    ezspin_db(156)%MassNum     = 113_iwp
    ezspin_db(156)%abundance    = 4.29_wp
    ezspin_db(156)%nucspin      = 4.5_wp
    ezspin_db(156)%gfactor      = +1.2286_wp
!----------------------------
!--> Record 157    = 120Sn
    ezspin_db(157)%AtNumb      = 50_iwp
    ezspin_db(157)%MassNum     = 120_iwp
    ezspin_db(157)%abundance    = 32.58_wp
    ezspin_db(157)%nucspin      = 0.0_wp
    ezspin_db(157)%gfactor      = 0.0_wp
!----------------------------
!--> Record 158    = 118Sn
    ezspin_db(158)%AtNumb      = 50_iwp
    ezspin_db(158)%MassNum     = 118_iwp
    ezspin_db(158)%abundance    = 24.22_wp
    ezspin_db(158)%nucspin      = 0.0_wp
    ezspin_db(158)%gfactor      = 0.0_wp
!----------------------------
!--> Record 159    = 116Sn
    ezspin_db(159)%AtNumb      = 50_iwp
    ezspin_db(159)%MassNum     = 116_iwp
    ezspin_db(159)%abundance    = 14.54_wp
    ezspin_db(159)%nucspin      = 0.0_wp
    ezspin_db(159)%gfactor      = 0.0_wp
!----------------------------
!--> Record 160    = 119Sn
    ezspin_db(160)%AtNumb      = 50_iwp
    ezspin_db(160)%MassNum     = 119_iwp
    ezspin_db(160)%abundance    = 8.59_wp
    ezspin_db(160)%nucspin      = 0.5_wp
    ezspin_db(160)%gfactor      = -2.09456_wp
!----------------------------
!--> Record 161    = 117Sn
    ezspin_db(161)%AtNumb      = 50_iwp
    ezspin_db(161)%MassNum     = 117_iwp
    ezspin_db(161)%abundance    = 7.68_wp
    ezspin_db(161)%nucspin      = 0.5_wp
    ezspin_db(161)%gfactor      = -2.00208_wp
!----------------------------
!--> Record 162    = 124Sn
    ezspin_db(162)%AtNumb      = 50_iwp
    ezspin_db(162)%MassNum     = 124_iwp
    ezspin_db(162)%abundance    = 5.79_wp
    ezspin_db(162)%nucspin      = 0.0_wp
    ezspin_db(162)%gfactor      = 0.0_wp
!----------------------------
!--> Record 163    = 122Sn
    ezspin_db(163)%AtNumb      = 50_iwp
    ezspin_db(163)%MassNum     = 122_iwp
    ezspin_db(163)%abundance    = 4.63_wp
    ezspin_db(163)%nucspin      = 0.0_wp
    ezspin_db(163)%gfactor      = 0.0_wp
!----------------------------
!--> Record 164    = 112Sn
    ezspin_db(164)%AtNumb      = 50_iwp
    ezspin_db(164)%MassNum     = 112_iwp
    ezspin_db(164)%abundance    = 0.97_wp
    ezspin_db(164)%nucspin      = 0.0_wp
    ezspin_db(164)%gfactor      = 0.0_wp
!----------------------------
!--> Record 165    = 114Sn
    ezspin_db(165)%AtNumb      = 50_iwp
    ezspin_db(165)%MassNum     = 114_iwp
    ezspin_db(165)%abundance    = 0.66_wp
    ezspin_db(165)%nucspin      = 0.0_wp
    ezspin_db(165)%gfactor      = 0.0_wp
!----------------------------
!--> Record 166    = 115Sn
    ezspin_db(166)%AtNumb      = 50_iwp
    ezspin_db(166)%MassNum     = 115_iwp
    ezspin_db(166)%abundance    = 0.34_wp
    ezspin_db(166)%nucspin      = 0.5_wp
    ezspin_db(166)%gfactor      = -1.8377_wp
!----------------------------
!--> Record 167    = 121Sb
    ezspin_db(167)%AtNumb      = 51_iwp
    ezspin_db(167)%MassNum     = 121_iwp
    ezspin_db(167)%abundance    = 57.21_wp
    ezspin_db(167)%nucspin      = 2.5_wp
    ezspin_db(167)%gfactor      = +1.3454_wp
!----------------------------
!--> Record 168    = 123Sb
    ezspin_db(168)%AtNumb      = 51_iwp
    ezspin_db(168)%MassNum     = 123_iwp
    ezspin_db(168)%abundance    = 42.79_wp
    ezspin_db(168)%nucspin      = 3.5_wp
    ezspin_db(168)%gfactor      = +0.72851_wp
!----------------------------
!--> Record 169    = 125Sb
    ezspin_db(169)%AtNumb      = 51_iwp
    ezspin_db(169)%MassNum     = 125_iwp
    ezspin_db(169)%abundance    = 0.0_wp
    ezspin_db(169)%nucspin      = 3.5_wp
    ezspin_db(169)%gfactor      = +0.751_wp
!----------------------------
!--> Record 170    = 130Te
    ezspin_db(170)%AtNumb      = 52_iwp
    ezspin_db(170)%MassNum     = 130_iwp
    ezspin_db(170)%abundance    = 34.08_wp
    ezspin_db(170)%nucspin      = 0.0_wp
    ezspin_db(170)%gfactor      = 0.0_wp
!----------------------------
!--> Record 171    = 128Te
    ezspin_db(171)%AtNumb      = 52_iwp
    ezspin_db(171)%MassNum     = 128_iwp
    ezspin_db(171)%abundance    = 31.74_wp
    ezspin_db(171)%nucspin      = 0.0_wp
    ezspin_db(171)%gfactor      = 0.0_wp
!----------------------------
!--> Record 172    = 126Te
    ezspin_db(172)%AtNumb      = 52_iwp
    ezspin_db(172)%MassNum     = 126_iwp
    ezspin_db(172)%abundance    = 18.84_wp
    ezspin_db(172)%nucspin      = 0.0_wp
    ezspin_db(172)%gfactor      = 0.0_wp
!----------------------------
!--> Record 173    = 125Te
    ezspin_db(173)%AtNumb      = 52_iwp
    ezspin_db(173)%MassNum     = 125_iwp
    ezspin_db(173)%abundance    = 7.07_wp
    ezspin_db(173)%nucspin      = 0.5_wp
    ezspin_db(173)%gfactor      = -1.7770102_wp
!----------------------------
!--> Record 174    = 124Te
    ezspin_db(174)%AtNumb      = 52_iwp
    ezspin_db(174)%MassNum     = 124_iwp
    ezspin_db(174)%abundance    = 4.74_wp
    ezspin_db(174)%nucspin      = 0.0_wp
    ezspin_db(174)%gfactor      = 0.0_wp
!----------------------------
!--> Record 175    = 122Te
    ezspin_db(175)%AtNumb      = 52_iwp
    ezspin_db(175)%MassNum     = 122_iwp
    ezspin_db(175)%abundance    = 2.55_wp
    ezspin_db(175)%nucspin      = 0.0_wp
    ezspin_db(175)%gfactor      = 0.0_wp
!----------------------------
!--> Record 176    = 123Te
    ezspin_db(176)%AtNumb      = 52_iwp
    ezspin_db(176)%MassNum     = 123_iwp
    ezspin_db(176)%abundance    = 0.89_wp
    ezspin_db(176)%nucspin      = 0.5_wp
    ezspin_db(176)%gfactor      = -1.473896_wp
!----------------------------
!--> Record 177    = 120Te
    ezspin_db(177)%AtNumb      = 52_iwp
    ezspin_db(177)%MassNum     = 120_iwp
    ezspin_db(177)%abundance    = 0.09_wp
    ezspin_db(177)%nucspin      = 0.0_wp
    ezspin_db(177)%gfactor      = 0.0_wp
!----------------------------
!--> Record 178    = 127I
    ezspin_db(178)%AtNumb      = 53_iwp
    ezspin_db(178)%MassNum     = 127_iwp
    ezspin_db(178)%abundance    = 100.0_wp
    ezspin_db(178)%nucspin      = 2.5_wp
    ezspin_db(178)%gfactor      = +1.12531_wp
!----------------------------
!--> Record 179    = 129I
    ezspin_db(179)%AtNumb      = 53_iwp
    ezspin_db(179)%MassNum     = 129_iwp
    ezspin_db(179)%abundance    = 0.0_wp
    ezspin_db(179)%nucspin      = 3.5_wp
    ezspin_db(179)%gfactor      = +0.74886_wp
!----------------------------
!--> Record 180    = 132Xe
    ezspin_db(180)%AtNumb      = 54_iwp
    ezspin_db(180)%MassNum     = 132_iwp
    ezspin_db(180)%abundance    = 26.89_wp
    ezspin_db(180)%nucspin      = 0.0_wp
    ezspin_db(180)%gfactor      = 0.0_wp
!----------------------------
!--> Record 181    = 129Xe
    ezspin_db(181)%AtNumb      = 54_iwp
    ezspin_db(181)%MassNum     = 129_iwp
    ezspin_db(181)%abundance    = 26.44_wp
    ezspin_db(181)%nucspin      = 0.5_wp
    ezspin_db(181)%gfactor      = -1.55595_wp
!----------------------------
!--> Record 182    = 131Xe
    ezspin_db(182)%AtNumb      = 54_iwp
    ezspin_db(182)%MassNum     = 131_iwp
    ezspin_db(182)%abundance    = 21.18_wp
    ezspin_db(182)%nucspin      = 1.5_wp
    ezspin_db(182)%gfactor      = +0.461_wp
!----------------------------
!--> Record 183    = 134Xe
    ezspin_db(183)%AtNumb      = 54_iwp
    ezspin_db(183)%MassNum     = 134_iwp
    ezspin_db(183)%abundance    = 10.44_wp
    ezspin_db(183)%nucspin      = 0.0_wp
    ezspin_db(183)%gfactor      = 0.0_wp
!----------------------------
!--> Record 184    = 136Xe
    ezspin_db(184)%AtNumb      = 54_iwp
    ezspin_db(184)%MassNum     = 136_iwp
    ezspin_db(184)%abundance    = 8.87_wp
    ezspin_db(184)%nucspin      = 0.0_wp
    ezspin_db(184)%gfactor      = 0.0_wp
!----------------------------
!--> Record 185    = 130Xe
    ezspin_db(185)%AtNumb      = 54_iwp
    ezspin_db(185)%MassNum     = 130_iwp
    ezspin_db(185)%abundance    = 4.08_wp
    ezspin_db(185)%nucspin      = 0.0_wp
    ezspin_db(185)%gfactor      = 0.0_wp
!----------------------------
!--> Record 186    = 128Xe
    ezspin_db(186)%AtNumb      = 54_iwp
    ezspin_db(186)%MassNum     = 128_iwp
    ezspin_db(186)%abundance    = 1.92_wp
    ezspin_db(186)%nucspin      = 0.0_wp
    ezspin_db(186)%gfactor      = 0.0_wp
!----------------------------
!--> Record 187    = 124Xe
    ezspin_db(187)%AtNumb      = 54_iwp
    ezspin_db(187)%MassNum     = 124_iwp
    ezspin_db(187)%abundance    = 0.09_wp
    ezspin_db(187)%nucspin      = 0.0_wp
    ezspin_db(187)%gfactor      = 0.0_wp
!----------------------------
!--> Record 188    = 126Xe
    ezspin_db(188)%AtNumb      = 54_iwp
    ezspin_db(188)%MassNum     = 126_iwp
    ezspin_db(188)%abundance    = 0.09_wp
    ezspin_db(188)%nucspin      = 0.0_wp
    ezspin_db(188)%gfactor      = 0.0_wp
!----------------------------
!--> Record 189    = 133Cs
    ezspin_db(189)%AtNumb      = 55_iwp
    ezspin_db(189)%MassNum     = 133_iwp
    ezspin_db(189)%abundance    = 100.0_wp
    ezspin_db(189)%nucspin      = 3.5_wp
    ezspin_db(189)%gfactor      = +0.7377214_wp
!----------------------------
!--> Record 190    = 134Cs
    ezspin_db(190)%AtNumb      = 55_iwp
    ezspin_db(190)%MassNum     = 134_iwp
    ezspin_db(190)%abundance    = 0.0_wp
    ezspin_db(190)%nucspin      = 4.0_wp
    ezspin_db(190)%gfactor      = +0.74843_wp
!----------------------------
!--> Record 191    = 135Cs
    ezspin_db(191)%AtNumb      = 55_iwp
    ezspin_db(191)%MassNum     = 135_iwp
    ezspin_db(191)%abundance    = 0.0_wp
    ezspin_db(191)%nucspin      = 3.5_wp
    ezspin_db(191)%gfactor      = +0.78069_wp
!----------------------------
!--> Record 192    = 137Cs
    ezspin_db(192)%AtNumb      = 55_iwp
    ezspin_db(192)%MassNum     = 137_iwp
    ezspin_db(192)%abundance    = 0.0_wp
    ezspin_db(192)%nucspin      = 3.5_wp
    ezspin_db(192)%gfactor      = +0.81466_wp
!----------------------------
!--> Record 193    = 138Ba
    ezspin_db(193)%AtNumb      = 56_iwp
    ezspin_db(193)%MassNum     = 138_iwp
    ezspin_db(193)%abundance    = 71.698_wp
    ezspin_db(193)%nucspin      = 0.0_wp
    ezspin_db(193)%gfactor      = 0.0_wp
!----------------------------
!--> Record 194    = 137Ba
    ezspin_db(194)%AtNumb      = 56_iwp
    ezspin_db(194)%MassNum     = 137_iwp
    ezspin_db(194)%abundance    = 11.232_wp
    ezspin_db(194)%nucspin      = 1.5_wp
    ezspin_db(194)%gfactor      = +0.62491_wp
!----------------------------
!--> Record 195    = 136Ba
    ezspin_db(195)%AtNumb      = 56_iwp
    ezspin_db(195)%MassNum     = 136_iwp
    ezspin_db(195)%abundance    = 7.854_wp
    ezspin_db(195)%nucspin      = 0.0_wp
    ezspin_db(195)%gfactor      = 0.0_wp
!----------------------------
!--> Record 196    = 135Ba
    ezspin_db(196)%AtNumb      = 56_iwp
    ezspin_db(196)%MassNum     = 135_iwp
    ezspin_db(196)%abundance    = 6.592_wp
    ezspin_db(196)%nucspin      = 1.5_wp
    ezspin_db(196)%gfactor      = +0.55863_wp
!----------------------------
!--> Record 197    = 134Ba
    ezspin_db(197)%AtNumb      = 56_iwp
    ezspin_db(197)%MassNum     = 134_iwp
    ezspin_db(197)%abundance    = 2.417_wp
    ezspin_db(197)%nucspin      = 0.0_wp
    ezspin_db(197)%gfactor      = 0.0_wp
!----------------------------
!--> Record 198    = 130Ba
    ezspin_db(198)%AtNumb      = 56_iwp
    ezspin_db(198)%MassNum     = 130_iwp
    ezspin_db(198)%abundance    = 0.106_wp
    ezspin_db(198)%nucspin      = 0.0_wp
    ezspin_db(198)%gfactor      = 0.0_wp
!----------------------------
!--> Record 199    = 132Ba
    ezspin_db(199)%AtNumb      = 56_iwp
    ezspin_db(199)%MassNum     = 132_iwp
    ezspin_db(199)%abundance    = 0.101_wp
    ezspin_db(199)%nucspin      = 0.0_wp
    ezspin_db(199)%gfactor      = 0.0_wp
!----------------------------
!--> Record 200    = 133Ba
    ezspin_db(200)%AtNumb      = 56_iwp
    ezspin_db(200)%MassNum     = 133_iwp
    ezspin_db(200)%abundance    = 0.0_wp
    ezspin_db(200)%nucspin      = 0.5_wp
    ezspin_db(200)%gfactor      = -1.5433_wp
!----------------------------
!--> Record 201    = 139La
    ezspin_db(201)%AtNumb      = 57_iwp
    ezspin_db(201)%MassNum     = 139_iwp
    ezspin_db(201)%abundance    = 99.910_wp
    ezspin_db(201)%nucspin      = 3.5_wp
    ezspin_db(201)%gfactor      = +0.795156_wp
!----------------------------
!--> Record 202    = 138La
    ezspin_db(202)%AtNumb      = 57_iwp
    ezspin_db(202)%MassNum     = 138_iwp
    ezspin_db(202)%abundance    = 0.090_wp
    ezspin_db(202)%nucspin      = 5.0_wp
    ezspin_db(202)%gfactor      = +0.742729_wp
!----------------------------
!--> Record 203    = 137La
    ezspin_db(203)%AtNumb      = 57_iwp
    ezspin_db(203)%MassNum     = 137_iwp
    ezspin_db(203)%abundance    = 0.0_wp
    ezspin_db(203)%nucspin      = 3.5_wp
    ezspin_db(203)%gfactor      = +0.7714_wp
!----------------------------
!--> Record 204    = 140Ce
    ezspin_db(204)%AtNumb      = 58_iwp
    ezspin_db(204)%MassNum     = 140_iwp
    ezspin_db(204)%abundance    = 88.450_wp
    ezspin_db(204)%nucspin      = 0.0_wp
    ezspin_db(204)%gfactor      = 0.0_wp
!----------------------------
!--> Record 205    = 142Ce
    ezspin_db(205)%AtNumb      = 58_iwp
    ezspin_db(205)%MassNum     = 142_iwp
    ezspin_db(205)%abundance    = 11.114_wp
    ezspin_db(205)%nucspin      = 0.0_wp
    ezspin_db(205)%gfactor      = 0.0_wp
!----------------------------
!--> Record 206    = 138Ce
    ezspin_db(206)%AtNumb      = 58_iwp
    ezspin_db(206)%MassNum     = 138_iwp
    ezspin_db(206)%abundance    = 0.251_wp
    ezspin_db(206)%nucspin      = 0.0_wp
    ezspin_db(206)%gfactor      = 0.0_wp
!----------------------------
!--> Record 207    = 136Ce
    ezspin_db(207)%AtNumb      = 58_iwp
    ezspin_db(207)%MassNum     = 136_iwp
    ezspin_db(207)%abundance    = 0.185_wp
    ezspin_db(207)%nucspin      = 0.0_wp
    ezspin_db(207)%gfactor      = 0.0_wp
!----------------------------
!--> Record 208    = 141Ce
    ezspin_db(208)%AtNumb      = 58_iwp
    ezspin_db(208)%MassNum     = 141_iwp
    ezspin_db(208)%abundance    = 0.00_wp
    ezspin_db(208)%nucspin      = 3.5_wp
    ezspin_db(208)%gfactor      = 0.3143_wp
!----------------------------
!--> Record 209    = 141Pr
    ezspin_db(209)%AtNumb      = 59_iwp
    ezspin_db(209)%MassNum     = 141_iwp
    ezspin_db(209)%abundance    = 100.0_wp
    ezspin_db(209)%nucspin      = 2.5_wp
    ezspin_db(209)%gfactor      = +1.7102_wp
!----------------------------
!--> Record 210    = 142Nd
    ezspin_db(210)%AtNumb      = 60_iwp
    ezspin_db(210)%MassNum     = 142_iwp
    ezspin_db(210)%abundance    = 27.2_wp
    ezspin_db(210)%nucspin      = 0.0_wp
    ezspin_db(210)%gfactor      = 0.0_wp
!----------------------------
!--> Record 211    = 144Nd
    ezspin_db(211)%AtNumb      = 60_iwp
    ezspin_db(211)%MassNum     = 144_iwp
    ezspin_db(211)%abundance    = 23.8_wp
    ezspin_db(211)%nucspin      = 0.0_wp
    ezspin_db(211)%gfactor      = 0.0_wp
!----------------------------
!--> Record 212    = 146Nd
    ezspin_db(212)%AtNumb      = 60_iwp
    ezspin_db(212)%MassNum     = 146_iwp
    ezspin_db(212)%abundance    = 17.2_wp
    ezspin_db(212)%nucspin      = 0.0_wp
    ezspin_db(212)%gfactor      = 0.0_wp
!----------------------------
!--> Record 213    = 143Nd
    ezspin_db(213)%AtNumb      = 60_iwp
    ezspin_db(213)%MassNum     = 143_iwp
    ezspin_db(213)%abundance    = 12.2_wp
    ezspin_db(213)%nucspin      = 3.5_wp
    ezspin_db(213)%gfactor      = -0.3043_wp
!----------------------------
!--> Record 214    = 145Nd
    ezspin_db(214)%AtNumb      = 60_iwp
    ezspin_db(214)%MassNum     = 145_iwp
    ezspin_db(214)%abundance    = 8.3_wp
    ezspin_db(214)%nucspin      = 3.5_wp
    ezspin_db(214)%gfactor      = -0.187_wp
!----------------------------
!--> Record 215    = 148Nd
    ezspin_db(215)%AtNumb      = 60_iwp
    ezspin_db(215)%MassNum     = 148_iwp
    ezspin_db(215)%abundance    = 5.7_wp
    ezspin_db(215)%nucspin      = 0.0_wp
    ezspin_db(215)%gfactor      = 0.0_wp
!----------------------------
!--> Record 216    = 150Nd
    ezspin_db(216)%AtNumb      = 60_iwp
    ezspin_db(216)%MassNum     = 150_iwp
    ezspin_db(216)%abundance    = 5.6_wp
    ezspin_db(216)%nucspin      = 0.0_wp
    ezspin_db(216)%gfactor      = 0.0_wp
!----------------------------
!--> Record 217    = 147Pm
    ezspin_db(217)%AtNumb      = 61_iwp
    ezspin_db(217)%MassNum     = 147_iwp
    ezspin_db(217)%abundance    = 0.0_wp
    ezspin_db(217)%nucspin      = 3.5_wp
    ezspin_db(217)%gfactor      = +0.737_wp
!----------------------------
!--> Record 218    = 152Sm
    ezspin_db(218)%AtNumb      = 62_iwp
    ezspin_db(218)%MassNum     = 152_iwp
    ezspin_db(218)%abundance    = 26.75_wp
    ezspin_db(218)%nucspin      = 0.0_wp
    ezspin_db(218)%gfactor      = 0.0_wp
!----------------------------
!--> Record 219    = 154Sm
    ezspin_db(219)%AtNumb      = 62_iwp
    ezspin_db(219)%MassNum     = 154_iwp
    ezspin_db(219)%abundance    = 22.75_wp
    ezspin_db(219)%nucspin      = 0.0_wp
    ezspin_db(219)%gfactor      = 0.0_wp
!----------------------------
!--> Record 220    = 147Sm
    ezspin_db(220)%AtNumb      = 62_iwp
    ezspin_db(220)%MassNum     = 147_iwp
    ezspin_db(220)%abundance    = 14.99_wp
    ezspin_db(220)%nucspin      = 3.5_wp
    ezspin_db(220)%gfactor      = -0.232_wp
!----------------------------
!--> Record 221    = 149Sm
    ezspin_db(221)%AtNumb      = 62_iwp
    ezspin_db(221)%MassNum     = 149_iwp
    ezspin_db(221)%abundance    = 13.82_wp
    ezspin_db(221)%nucspin      = 3.5_wp
    ezspin_db(221)%gfactor      = -0.1908_wp
!----------------------------
!--> Record 222    = 148Sm
    ezspin_db(222)%AtNumb      = 62_iwp
    ezspin_db(222)%MassNum     = 148_iwp
    ezspin_db(222)%abundance    = 11.24_wp
    ezspin_db(222)%nucspin      = 0.0_wp
    ezspin_db(222)%gfactor      = 0.0_wp
!----------------------------
!--> Record 223    = 150Sm
    ezspin_db(223)%AtNumb      = 62_iwp
    ezspin_db(223)%MassNum     = 150_iwp
    ezspin_db(223)%abundance    = 7.38_wp
    ezspin_db(223)%nucspin      = 0.0_wp
    ezspin_db(223)%gfactor      = 0.0_wp
!----------------------------
!--> Record 224    = 144Sm
    ezspin_db(224)%AtNumb      = 62_iwp
    ezspin_db(224)%MassNum     = 144_iwp
    ezspin_db(224)%abundance    = 3.07_wp
    ezspin_db(224)%nucspin      = 0.0_wp
    ezspin_db(224)%gfactor      = 0.0_wp
!----------------------------
!--> Record 225    = 151Sm
    ezspin_db(225)%AtNumb      = 62_iwp
    ezspin_db(225)%MassNum     = 151_iwp
    ezspin_db(225)%abundance    = 0.0_wp
    ezspin_db(225)%nucspin      = 2.5_wp
    ezspin_db(225)%gfactor      = +0.1444_wp
!----------------------------
!--> Record 226    = 153Eu
    ezspin_db(226)%AtNumb      = 63_iwp
    ezspin_db(226)%MassNum     = 153_iwp
    ezspin_db(226)%abundance    = 52.19_wp
    ezspin_db(226)%nucspin      = 2.5_wp
    ezspin_db(226)%gfactor      = +0.6134_wp
!----------------------------
!--> Record 227    = 151Eu
    ezspin_db(227)%AtNumb      = 63_iwp
    ezspin_db(227)%MassNum     = 151_iwp
    ezspin_db(227)%abundance    = 47.81_wp
    ezspin_db(227)%nucspin      = 2.5_wp
    ezspin_db(227)%gfactor      = +1.3887_wp
!----------------------------
!--> Record 228    = 152Eu
    ezspin_db(228)%AtNumb      = 63_iwp
    ezspin_db(228)%MassNum     = 152_iwp
    ezspin_db(228)%abundance    = 0.0_wp
    ezspin_db(228)%nucspin      = 3.0_wp
    ezspin_db(228)%gfactor      = -0.6467_wp
!----------------------------
!--> Record 229    = 154Eu
    ezspin_db(229)%AtNumb      = 63_iwp
    ezspin_db(229)%MassNum     = 154_iwp
    ezspin_db(229)%abundance    = 0.0_wp
    ezspin_db(229)%nucspin      = 3.0_wp
    ezspin_db(229)%gfactor      = -0.6683_wp
!----------------------------
!--> Record 230    = 155Eu
    ezspin_db(230)%AtNumb      = 63_iwp
    ezspin_db(230)%MassNum     = 155_iwp
    ezspin_db(230)%abundance    = 0.0_wp
    ezspin_db(230)%nucspin      = 2.5_wp
    ezspin_db(230)%gfactor      = +0.608_wp
!----------------------------
!--> Record 231    = 158Gd
    ezspin_db(231)%AtNumb      = 64_iwp
    ezspin_db(231)%MassNum     = 158_iwp
    ezspin_db(231)%abundance    = 24.84_wp
    ezspin_db(231)%nucspin      = 0.0_wp
    ezspin_db(231)%gfactor      = 0.0_wp
!----------------------------
!--> Record 232    = 160Gd
    ezspin_db(232)%AtNumb      = 64_iwp
    ezspin_db(232)%MassNum     = 160_iwp
    ezspin_db(232)%abundance    = 21.86_wp
    ezspin_db(232)%nucspin      = 0.0_wp
    ezspin_db(232)%gfactor      = 0.0_wp
!----------------------------
!--> Record 233    = 156Gd
    ezspin_db(233)%AtNumb      = 64_iwp
    ezspin_db(233)%MassNum     = 156_iwp
    ezspin_db(233)%abundance    = 20.47_wp
    ezspin_db(233)%nucspin      = 0.0_wp
    ezspin_db(233)%gfactor      = 0.0_wp
!----------------------------
!--> Record 234    = 157Gd
    ezspin_db(234)%AtNumb      = 64_iwp
    ezspin_db(234)%MassNum     = 157_iwp
    ezspin_db(234)%abundance    = 15.65_wp
    ezspin_db(234)%nucspin      = 1.5_wp
    ezspin_db(234)%gfactor      = -0.2265_wp
!----------------------------
!--> Record 235    = 155Gd
    ezspin_db(235)%AtNumb      = 64_iwp
    ezspin_db(235)%MassNum     = 155_iwp
    ezspin_db(235)%abundance    = 14.80_wp
    ezspin_db(235)%nucspin      = 1.5_wp
    ezspin_db(235)%gfactor      = -0.1715_wp
!----------------------------
!--> Record 236    = 154Gd
    ezspin_db(236)%AtNumb      = 64_iwp
    ezspin_db(236)%MassNum     = 154_iwp
    ezspin_db(236)%abundance    = 2.18_wp
    ezspin_db(236)%nucspin      = 0.0_wp
    ezspin_db(236)%gfactor      = 0.0_wp
!----------------------------
!--> Record 237    = 152Gd
    ezspin_db(237)%AtNumb      = 64_iwp
    ezspin_db(237)%MassNum     = 152_iwp
    ezspin_db(237)%abundance    = 0.20_wp
    ezspin_db(237)%nucspin      = 0.0_wp
    ezspin_db(237)%gfactor      = 0.0_wp
!----------------------------
!--> Record 238    = 159Tb
    ezspin_db(238)%AtNumb      = 65_iwp
    ezspin_db(238)%MassNum     = 159_iwp
    ezspin_db(238)%abundance    = 100.0_wp
    ezspin_db(238)%nucspin      = 1.5_wp
    ezspin_db(238)%gfactor      = +1.343_wp
!----------------------------
!--> Record 239    = 157Tb
    ezspin_db(239)%AtNumb      = 65_iwp
    ezspin_db(239)%MassNum     = 157_iwp
    ezspin_db(239)%abundance    = 0.0_wp
    ezspin_db(239)%nucspin      = 1.5_wp
    ezspin_db(239)%gfactor      = +1.34_wp
!----------------------------
!--> Record 240    = 160Tb
    ezspin_db(240)%AtNumb      = 65_iwp
    ezspin_db(240)%MassNum     = 160_iwp
    ezspin_db(240)%abundance    = 0.0_wp
    ezspin_db(240)%nucspin      = 3.0_wp
    ezspin_db(240)%gfactor      = +0.5967_wp
!----------------------------
!--> Record 241    = 164Dy
    ezspin_db(241)%AtNumb      = 66_iwp
    ezspin_db(241)%MassNum     = 164_iwp
    ezspin_db(241)%abundance    = 28.18_wp
    ezspin_db(241)%nucspin      = 0.0_wp
    ezspin_db(241)%gfactor      = 0.0_wp
!----------------------------
!--> Record 242    = 162Dy
    ezspin_db(242)%AtNumb      = 66_iwp
    ezspin_db(242)%MassNum     = 162_iwp
    ezspin_db(242)%abundance    = 25.51_wp
    ezspin_db(242)%nucspin      = 0.0_wp
    ezspin_db(242)%gfactor      = 0.0_wp
!----------------------------
!--> Record 243    = 163Dy
    ezspin_db(243)%AtNumb      = 66_iwp
    ezspin_db(243)%MassNum     = 163_iwp
    ezspin_db(243)%abundance    = 24.90_wp
    ezspin_db(243)%nucspin      = 2.5_wp
    ezspin_db(243)%gfactor      = +0.269_wp
!----------------------------
!--> Record 244    = 161Dy
    ezspin_db(244)%AtNumb      = 66_iwp
    ezspin_db(244)%MassNum     = 161_iwp
    ezspin_db(244)%abundance    = 18.91_wp
    ezspin_db(244)%nucspin      = 2.5_wp
    ezspin_db(244)%gfactor      = -0.192_wp
!----------------------------
!--> Record 245    = 160Dy
    ezspin_db(245)%AtNumb      = 66_iwp
    ezspin_db(245)%MassNum     = 160_iwp
    ezspin_db(245)%abundance    = 2.34_wp
    ezspin_db(245)%nucspin      = 0.0_wp
    ezspin_db(245)%gfactor      = 0.0_wp
!----------------------------
!--> Record 246    = 158Dy
    ezspin_db(246)%AtNumb      = 66_iwp
    ezspin_db(246)%MassNum     = 158_iwp
    ezspin_db(246)%abundance    = 0.10_wp
    ezspin_db(246)%nucspin      = 0.0_wp
    ezspin_db(246)%gfactor      = 0.0_wp
!----------------------------
!--> Record 247    = 156Dy
    ezspin_db(247)%AtNumb      = 66_iwp
    ezspin_db(247)%MassNum     = 156_iwp
    ezspin_db(247)%abundance    = 0.06_wp
    ezspin_db(247)%nucspin      = 0.0_wp
    ezspin_db(247)%gfactor      = 0.0_wp
!----------------------------
!--> Record 248    = 165Ho
    ezspin_db(248)%AtNumb      = 67_iwp
    ezspin_db(248)%MassNum     = 165_iwp
    ezspin_db(248)%abundance    = 100.0_wp
    ezspin_db(248)%nucspin      = 3.5_wp
    ezspin_db(248)%gfactor      = +1.668_wp
!----------------------------
!--> Record 249    = 166Er
    ezspin_db(249)%AtNumb      = 68_iwp
    ezspin_db(249)%MassNum     = 166_iwp
    ezspin_db(249)%abundance    = 33.61_wp
    ezspin_db(249)%nucspin      = 0.0_wp
    ezspin_db(249)%gfactor      = 0.0_wp
!----------------------------
!--> Record 250    = 168Er
    ezspin_db(250)%AtNumb      = 68_iwp
    ezspin_db(250)%MassNum     = 168_iwp
    ezspin_db(250)%abundance    = 26.78_wp
    ezspin_db(250)%nucspin      = 0.0_wp
    ezspin_db(250)%gfactor      = 0.0_wp
!----------------------------
!--> Record 251    = 167Er
    ezspin_db(251)%AtNumb      = 68_iwp
    ezspin_db(251)%MassNum     = 167_iwp
    ezspin_db(251)%abundance    = 22.93_wp
    ezspin_db(251)%nucspin      = 3.5_wp
    ezspin_db(251)%gfactor      = -0.1611_wp
!----------------------------
!--> Record 252    = 170Er
    ezspin_db(252)%AtNumb      = 68_iwp
    ezspin_db(252)%MassNum     = 170_iwp
    ezspin_db(252)%abundance    = 14.93_wp
    ezspin_db(252)%nucspin      = 0.0_wp
    ezspin_db(252)%gfactor      = 0.0_wp
!----------------------------
!--> Record 253    = 164Er
    ezspin_db(253)%AtNumb      = 68_iwp
    ezspin_db(253)%MassNum     = 164_iwp
    ezspin_db(253)%abundance    = 1.61_wp
    ezspin_db(253)%nucspin      = 0.0_wp
    ezspin_db(253)%gfactor      = 0.0_wp
!----------------------------
!--> Record 254    = 162Er
    ezspin_db(254)%AtNumb      = 68_iwp
    ezspin_db(254)%MassNum     = 162_iwp
    ezspin_db(254)%abundance    = 0.14_wp
    ezspin_db(254)%nucspin      = 0.0_wp
    ezspin_db(254)%gfactor      = 0.0_wp
!----------------------------
!--> Record 255    = 169Tm
    ezspin_db(255)%AtNumb      = 69_iwp
    ezspin_db(255)%MassNum     = 169_iwp
    ezspin_db(255)%abundance    = 100.0_wp
    ezspin_db(255)%nucspin      = 0.5_wp
    ezspin_db(255)%gfactor      = -0.462_wp
!----------------------------
!--> Record 256    = 171Tm
    ezspin_db(256)%AtNumb      = 69_iwp
    ezspin_db(256)%MassNum     = 171_iwp
    ezspin_db(256)%abundance    = 0.0_wp
    ezspin_db(256)%nucspin      = 0.5_wp
    ezspin_db(256)%gfactor      = -0.456_wp
!----------------------------
!--> Record 257    = 174Yb
    ezspin_db(257)%AtNumb      = 70_iwp
    ezspin_db(257)%MassNum     = 174_iwp
    ezspin_db(257)%abundance    = 31.83_wp
    ezspin_db(257)%nucspin      = 0.0_wp
    ezspin_db(257)%gfactor      = 0.0_wp
!----------------------------
!--> Record 258    = 172Yb
    ezspin_db(258)%AtNumb      = 70_iwp
    ezspin_db(258)%MassNum     = 172_iwp
    ezspin_db(258)%abundance    = 21.83_wp
    ezspin_db(258)%nucspin      = 0.0_wp
    ezspin_db(258)%gfactor      = 0.0_wp
!----------------------------
!--> Record 259    = 173Yb
    ezspin_db(259)%AtNumb      = 70_iwp
    ezspin_db(259)%MassNum     = 173_iwp
    ezspin_db(259)%abundance    = 16.13_wp
    ezspin_db(259)%nucspin      = 2.5_wp
    ezspin_db(259)%gfactor      = -0.2592_wp
!----------------------------
!--> Record 260    = 171Yb
    ezspin_db(260)%AtNumb      = 70_iwp
    ezspin_db(260)%MassNum     = 171_iwp
    ezspin_db(260)%abundance    = 14.28_wp
    ezspin_db(260)%nucspin      = 0.5_wp
    ezspin_db(260)%gfactor      = +0.98734_wp
!----------------------------
!--> Record 261    = 176Yb
    ezspin_db(261)%AtNumb      = 70_iwp
    ezspin_db(261)%MassNum     = 176_iwp
    ezspin_db(261)%abundance    = 12.76_wp
    ezspin_db(261)%nucspin      = 0.0_wp
    ezspin_db(261)%gfactor      = 0.0_wp
!----------------------------
!--> Record 262    = 170Yb
    ezspin_db(262)%AtNumb      = 70_iwp
    ezspin_db(262)%MassNum     = 170_iwp
    ezspin_db(262)%abundance    = 3.04_wp
    ezspin_db(262)%nucspin      = 0.0_wp
    ezspin_db(262)%gfactor      = 0.0_wp
!----------------------------
!--> Record 263    = 168Yb
    ezspin_db(263)%AtNumb      = 70_iwp
    ezspin_db(263)%MassNum     = 168_iwp
    ezspin_db(263)%abundance    = 0.13_wp
    ezspin_db(263)%nucspin      = 0.0_wp
    ezspin_db(263)%gfactor      = 0.0_wp
!----------------------------
!--> Record 264    = 175Lu
    ezspin_db(264)%AtNumb      = 71_iwp
    ezspin_db(264)%MassNum     = 175_iwp
    ezspin_db(264)%abundance    = 97.41_wp
    ezspin_db(264)%nucspin      = 3.5_wp
    ezspin_db(264)%gfactor      = +0.6378_wp
!----------------------------
!--> Record 265    = 176Lu
    ezspin_db(265)%AtNumb      = 71_iwp
    ezspin_db(265)%MassNum     = 176_iwp
    ezspin_db(265)%abundance    = 2.59_wp
    ezspin_db(265)%nucspin      = 7.0_wp
    ezspin_db(265)%gfactor      = +0.4517_wp
!----------------------------
!--> Record 266    = 173Lu
    ezspin_db(266)%AtNumb      = 71_iwp
    ezspin_db(266)%MassNum     = 173_iwp
    ezspin_db(266)%abundance    = 0.0_wp
    ezspin_db(266)%nucspin      = 3.5_wp
    ezspin_db(266)%gfactor      = +0.6517_wp
!----------------------------
!--> Record 267    = 174Lu
    ezspin_db(267)%AtNumb      = 71_iwp
    ezspin_db(267)%MassNum     = 174_iwp
    ezspin_db(267)%abundance    = 0.0_wp
    ezspin_db(267)%nucspin      = 1.0_wp
    ezspin_db(267)%gfactor      = +1.988_wp
!----------------------------
!--> Record 268    = 180Hf
    ezspin_db(268)%AtNumb      = 72_iwp
    ezspin_db(268)%MassNum     = 180_iwp
    ezspin_db(268)%abundance    = 35.08_wp
    ezspin_db(268)%nucspin      = 0.0_wp
    ezspin_db(268)%gfactor      = 0.0_wp
!----------------------------
!--> Record 269    = 178Hf
    ezspin_db(269)%AtNumb      = 72_iwp
    ezspin_db(269)%MassNum     = 178_iwp
    ezspin_db(269)%abundance    = 27.28_wp
    ezspin_db(269)%nucspin      = 0.0_wp
    ezspin_db(269)%gfactor      = 0.0_wp
!----------------------------
!--> Record 270    = 177Hf
    ezspin_db(270)%AtNumb      = 72_iwp
    ezspin_db(270)%MassNum     = 177_iwp
    ezspin_db(270)%abundance    = 18.60_wp
    ezspin_db(270)%nucspin      = 3.5_wp
    ezspin_db(270)%gfactor      = +0.2267_wp
!----------------------------
!--> Record 271    = 179Hf
    ezspin_db(271)%AtNumb      = 72_iwp
    ezspin_db(271)%MassNum     = 179_iwp
    ezspin_db(271)%abundance    = 13.62_wp
    ezspin_db(271)%nucspin      = 4.5_wp
    ezspin_db(271)%gfactor      = -0.1424_wp
!----------------------------
!--> Record 272    = 176Hf
    ezspin_db(272)%AtNumb      = 72_iwp
    ezspin_db(272)%MassNum     = 176_iwp
    ezspin_db(272)%abundance    = 5.26_wp
    ezspin_db(272)%nucspin      = 0.0_wp
    ezspin_db(272)%gfactor      = 0.0_wp
!----------------------------
!--> Record 273    = 174Hf
    ezspin_db(273)%AtNumb      = 72_iwp
    ezspin_db(273)%MassNum     = 174_iwp
    ezspin_db(273)%abundance    = 0.16_wp
    ezspin_db(273)%nucspin      = 0.0_wp
    ezspin_db(273)%gfactor      = 0.0_wp
!----------------------------
!--> Record 274    = 181Ta
    ezspin_db(274)%AtNumb      = 73_iwp
    ezspin_db(274)%MassNum     = 181_iwp
    ezspin_db(274)%abundance    = 99.988_wp
    ezspin_db(274)%nucspin      = 3.5_wp
    ezspin_db(274)%gfactor      = +0.67729_wp
!----------------------------
!--> Record 275    = 180Ta
    ezspin_db(275)%AtNumb      = 73_iwp
    ezspin_db(275)%MassNum     = 180_iwp
    ezspin_db(275)%abundance    = 0.012_wp
    ezspin_db(275)%nucspin      = 9.0_wp
    ezspin_db(275)%gfactor      = 0.5361_wp
!----------------------------
!--> Record 276    = 184W
    ezspin_db(276)%AtNumb      = 74_iwp
    ezspin_db(276)%MassNum     = 184_iwp
    ezspin_db(276)%abundance    = 30.64_wp
    ezspin_db(276)%nucspin      = 0.0_wp
    ezspin_db(276)%gfactor      = 0.0_wp
!----------------------------
!--> Record 277    = 186W
    ezspin_db(277)%AtNumb      = 74_iwp
    ezspin_db(277)%MassNum     = 186_iwp
    ezspin_db(277)%abundance    = 28.43_wp
    ezspin_db(277)%nucspin      = 0.0_wp
    ezspin_db(277)%gfactor      = 0.0_wp
!----------------------------
!--> Record 278    = 182W
    ezspin_db(278)%AtNumb      = 74_iwp
    ezspin_db(278)%MassNum     = 182_iwp
    ezspin_db(278)%abundance    = 26.50_wp
    ezspin_db(278)%nucspin      = 0.0_wp
    ezspin_db(278)%gfactor      = 0.0_wp
!----------------------------
!--> Record 279    = 183W
    ezspin_db(279)%AtNumb      = 74_iwp
    ezspin_db(279)%MassNum     = 183_iwp
    ezspin_db(279)%abundance    = 14.31_wp
    ezspin_db(279)%nucspin      = 0.5_wp
    ezspin_db(279)%gfactor      = +0.2355695_wp
!----------------------------
!--> Record 280    = 180W
    ezspin_db(280)%AtNumb      = 74_iwp
    ezspin_db(280)%MassNum     = 180_iwp
    ezspin_db(280)%abundance    = 0.12_wp
    ezspin_db(280)%nucspin      = 0.0_wp
    ezspin_db(280)%gfactor      = 0.0_wp
!----------------------------
!--> Record 281    = 187Re
    ezspin_db(281)%AtNumb      = 75_iwp
    ezspin_db(281)%MassNum     = 187_iwp
    ezspin_db(281)%abundance    = 62.60_wp
    ezspin_db(281)%nucspin      = 2.5_wp
    ezspin_db(281)%gfactor      = +1.2879_wp
!----------------------------
!--> Record 282    = 185Re
    ezspin_db(282)%AtNumb      = 75_iwp
    ezspin_db(282)%MassNum     = 185_iwp
    ezspin_db(282)%abundance    = 37.40_wp
    ezspin_db(282)%nucspin      = 2.5_wp
    ezspin_db(282)%gfactor      = +1.2748_wp
!----------------------------
!--> Record 283    = 192Os
    ezspin_db(283)%AtNumb      = 76_iwp
    ezspin_db(283)%MassNum     = 192_iwp
    ezspin_db(283)%abundance    = 40.78_wp
    ezspin_db(283)%nucspin      = 0.0_wp
    ezspin_db(283)%gfactor      = 0.0_wp
!----------------------------
!--> Record 284    = 190Os
    ezspin_db(284)%AtNumb      = 76_iwp
    ezspin_db(284)%MassNum     = 190_iwp
    ezspin_db(284)%abundance    = 26.26_wp
    ezspin_db(284)%nucspin      = 0.0_wp
    ezspin_db(284)%gfactor      = 0.0_wp
!----------------------------
!--> Record 285    = 189Os
    ezspin_db(285)%AtNumb      = 76_iwp
    ezspin_db(285)%MassNum     = 189_iwp
    ezspin_db(285)%abundance    = 16.15_wp
    ezspin_db(285)%nucspin      = 1.5_wp
    ezspin_db(285)%gfactor      = +0.439956_wp
!----------------------------
!--> Record 286    = 188Os
    ezspin_db(286)%AtNumb      = 76_iwp
    ezspin_db(286)%MassNum     = 188_iwp
    ezspin_db(286)%abundance    = 13.24_wp
    ezspin_db(286)%nucspin      = 0.0_wp
    ezspin_db(286)%gfactor      = 0.0_wp
!----------------------------
!--> Record 287    = 187Os
    ezspin_db(287)%AtNumb      = 76_iwp
    ezspin_db(287)%MassNum     = 187_iwp
    ezspin_db(287)%abundance    = 1.96_wp
    ezspin_db(287)%nucspin      = 0.5_wp
    ezspin_db(287)%gfactor      = +0.1293038_wp
!----------------------------
!--> Record 288    = 186Os
    ezspin_db(288)%AtNumb      = 76_iwp
    ezspin_db(288)%MassNum     = 186_iwp
    ezspin_db(288)%abundance    = 1.59_wp
    ezspin_db(288)%nucspin      = 0.0_wp
    ezspin_db(288)%gfactor      = 0.0_wp
!----------------------------
!--> Record 289    = 184Os
    ezspin_db(289)%AtNumb      = 76_iwp
    ezspin_db(289)%MassNum     = 184_iwp
    ezspin_db(289)%abundance    = 0.02_wp
    ezspin_db(289)%nucspin      = 0.0_wp
    ezspin_db(289)%gfactor      = 0.0_wp
!----------------------------
!--> Record 290    = 193Ir
    ezspin_db(290)%AtNumb      = 77_iwp
    ezspin_db(290)%MassNum     = 193_iwp
    ezspin_db(290)%abundance    = 62.7_wp
    ezspin_db(290)%nucspin      = 1.5_wp
    ezspin_db(290)%gfactor      = +0.1091_wp
!----------------------------
!--> Record 291    = 191Ir
    ezspin_db(291)%AtNumb      = 77_iwp
    ezspin_db(291)%MassNum     = 191_iwp
    ezspin_db(291)%abundance    = 37.3_wp
    ezspin_db(291)%nucspin      = 1.5_wp
    ezspin_db(291)%gfactor      = +0.1005_wp
!----------------------------
!--> Record 292    = 195Pt
    ezspin_db(292)%AtNumb      = 78_iwp
    ezspin_db(292)%MassNum     = 195_iwp
    ezspin_db(292)%abundance    = 33.832_wp
    ezspin_db(292)%nucspin      = 0.5_wp
    ezspin_db(292)%gfactor      = +1.2190_wp
!----------------------------
!--> Record 293    = 194Pt
    ezspin_db(293)%AtNumb      = 78_iwp
    ezspin_db(293)%MassNum     = 194_iwp
    ezspin_db(293)%abundance    = 32.967_wp
    ezspin_db(293)%nucspin      = 0.0_wp
    ezspin_db(293)%gfactor      = 0.0_wp
!----------------------------
!--> Record 294    = 196Pt
    ezspin_db(294)%AtNumb      = 78_iwp
    ezspin_db(294)%MassNum     = 196_iwp
    ezspin_db(294)%abundance    = 25.242_wp
    ezspin_db(294)%nucspin      = 0.0_wp
    ezspin_db(294)%gfactor      = 0.0_wp
!----------------------------
!--> Record 295    = 198Pt
    ezspin_db(295)%AtNumb      = 78_iwp
    ezspin_db(295)%MassNum     = 198_iwp
    ezspin_db(295)%abundance    = 7.163_wp
    ezspin_db(295)%nucspin      = 0.0_wp
    ezspin_db(295)%gfactor      = 0.0_wp
!----------------------------
!--> Record 296    = 192Pt
    ezspin_db(296)%AtNumb      = 78_iwp
    ezspin_db(296)%MassNum     = 192_iwp
    ezspin_db(296)%abundance    = 0.784_wp
    ezspin_db(296)%nucspin      = 0.0_wp
    ezspin_db(296)%gfactor      = 0.0_wp
!----------------------------
!--> Record 297    = 190Pt
    ezspin_db(297)%AtNumb      = 78_iwp
    ezspin_db(297)%MassNum     = 190_iwp
    ezspin_db(297)%abundance    = 0.014_wp
    ezspin_db(297)%nucspin      = 0.0_wp
    ezspin_db(297)%gfactor      = 0.0_wp
!----------------------------
!--> Record 298    = 197Au
    ezspin_db(298)%AtNumb      = 79_iwp
    ezspin_db(298)%MassNum     = 197_iwp
    ezspin_db(298)%abundance    = 100.0_wp
    ezspin_db(298)%nucspin      = 1.5_wp
    ezspin_db(298)%gfactor      = +0.097164_wp
!----------------------------
!--> Record 299    = 202Hg
    ezspin_db(299)%AtNumb      = 80_iwp
    ezspin_db(299)%MassNum     = 202_iwp
    ezspin_db(299)%abundance    = 29.86_wp
    ezspin_db(299)%nucspin      = 0.0_wp
    ezspin_db(299)%gfactor      = 0.0_wp
!----------------------------
!--> Record 300    = 200Hg
    ezspin_db(300)%AtNumb      = 80_iwp
    ezspin_db(300)%MassNum     = 200_iwp
    ezspin_db(300)%abundance    = 23.10_wp
    ezspin_db(300)%nucspin      = 0.0_wp
    ezspin_db(300)%gfactor      = 0.0_wp
!----------------------------
!--> Record 301    = 199Hg
    ezspin_db(301)%AtNumb      = 80_iwp
    ezspin_db(301)%MassNum     = 199_iwp
    ezspin_db(301)%abundance    = 16.87_wp
    ezspin_db(301)%nucspin      = 0.5_wp
    ezspin_db(301)%gfactor      = +1.011771_wp
!----------------------------
!--> Record 302    = 201Hg
    ezspin_db(302)%AtNumb      = 80_iwp
    ezspin_db(302)%MassNum     = 201_iwp
    ezspin_db(302)%abundance    = 13.18_wp
    ezspin_db(302)%nucspin      = 1.5_wp
    ezspin_db(302)%gfactor      = -0.373484_wp
!----------------------------
!--> Record 303    = 198Hg
    ezspin_db(303)%AtNumb      = 80_iwp
    ezspin_db(303)%MassNum     = 198_iwp
    ezspin_db(303)%abundance    = 9.97_wp
    ezspin_db(303)%nucspin      = 0.0_wp
    ezspin_db(303)%gfactor      = 0.0_wp
!----------------------------
!--> Record 304    = 204Hg
    ezspin_db(304)%AtNumb      = 80_iwp
    ezspin_db(304)%MassNum     = 204_iwp
    ezspin_db(304)%abundance    = 6.87_wp
    ezspin_db(304)%nucspin      = 0.0_wp
    ezspin_db(304)%gfactor      = 0.0_wp
!----------------------------
!--> Record 305    = 196Hg
    ezspin_db(305)%AtNumb      = 80_iwp
    ezspin_db(305)%MassNum     = 196_iwp
    ezspin_db(305)%abundance    = 0.15_wp
    ezspin_db(305)%nucspin      = 0.0_wp
    ezspin_db(305)%gfactor      = 0.0_wp
!----------------------------
!--> Record 306    = 205Tl
    ezspin_db(306)%AtNumb      = 81_iwp
    ezspin_db(306)%MassNum     = 205_iwp
    ezspin_db(306)%abundance    = 70.476_wp
    ezspin_db(306)%nucspin      = 0.5_wp
    ezspin_db(306)%gfactor      = +3.2764292_wp
!----------------------------
!--> Record 307    = 203Tl
    ezspin_db(307)%AtNumb      = 81_iwp
    ezspin_db(307)%MassNum     = 203_iwp
    ezspin_db(307)%abundance    = 29.524_wp
    ezspin_db(307)%nucspin      = 0.5_wp
    ezspin_db(307)%gfactor      = +3.24451574_wp
!----------------------------
!--> Record 308    = 204Tl
    ezspin_db(308)%AtNumb      = 81_iwp
    ezspin_db(308)%MassNum     = 204_iwp
    ezspin_db(308)%abundance    = 0.0_wp
    ezspin_db(308)%nucspin      = 2.0_wp
    ezspin_db(308)%gfactor      = +0.045_wp
!----------------------------
!--> Record 309    = 208Pb
    ezspin_db(309)%AtNumb      = 82_iwp
    ezspin_db(309)%MassNum     = 208_iwp
    ezspin_db(309)%abundance    = 52.4_wp
    ezspin_db(309)%nucspin      = 0.0_wp
    ezspin_db(309)%gfactor      = 0.0_wp
!----------------------------
!--> Record 310    = 206Pb
    ezspin_db(310)%AtNumb      = 82_iwp
    ezspin_db(310)%MassNum     = 206_iwp
    ezspin_db(310)%abundance    = 24.1_wp
    ezspin_db(310)%nucspin      = 0.0_wp
    ezspin_db(310)%gfactor      = 0.0_wp
!----------------------------
!--> Record 311    = 207Pb
    ezspin_db(311)%AtNumb      = 82_iwp
    ezspin_db(311)%MassNum     = 207_iwp
    ezspin_db(311)%abundance    = 22.1_wp
    ezspin_db(311)%nucspin      = 0.5_wp
    ezspin_db(311)%gfactor      = +1.18512_wp
!----------------------------
!--> Record 312    = 204Pb
    ezspin_db(312)%AtNumb      = 82_iwp
    ezspin_db(312)%MassNum     = 204_iwp
    ezspin_db(312)%abundance    = 1.4_wp
    ezspin_db(312)%nucspin      = 0.0_wp
    ezspin_db(312)%gfactor      = 0.0_wp
!----------------------------
!--> Record 313    = 209Bi
    ezspin_db(313)%AtNumb      = 83_iwp
    ezspin_db(313)%MassNum     = 209_iwp
    ezspin_db(313)%abundance    = 100.0_wp
    ezspin_db(313)%nucspin      = 4.5_wp
    ezspin_db(313)%gfactor      = +0.9134_wp
!----------------------------
!--> Record 314    = 207Bi
    ezspin_db(314)%AtNumb      = 83_iwp
    ezspin_db(314)%MassNum     = 207_iwp
    ezspin_db(314)%abundance    = 0.0_wp
    ezspin_db(314)%nucspin      = 4.5_wp
    ezspin_db(314)%gfactor      = +0.9092_wp
!----------------------------
!--> Record 315    = 209Po
    ezspin_db(315)%AtNumb      = 84_iwp
    ezspin_db(315)%MassNum     = 209_iwp
    ezspin_db(315)%abundance    = 0.0_wp
    ezspin_db(315)%nucspin      = 0.5_wp
    ezspin_db(315)%gfactor      = +1.5_wp
!----------------------------
!--> Record 316    = 0At
    ezspin_db(316)%AtNumb      = 85_iwp
    ezspin_db(316)%MassNum     = 0_iwp
    ezspin_db(316)%abundance    = 0.0_wp
    ezspin_db(316)%nucspin      = -1.0_wp
    ezspin_db(316)%gfactor      = 0.0_wp
!----------------------------
!--> Record 317    = 0Rn
    ezspin_db(317)%AtNumb      = 86_iwp
    ezspin_db(317)%MassNum     = 0_iwp
    ezspin_db(317)%abundance    = 0.0_wp
    ezspin_db(317)%nucspin      = -1.0_wp
    ezspin_db(317)%gfactor      = 0.0_wp
!----------------------------
!--> Record 318    = 0Fr
    ezspin_db(318)%AtNumb      = 87_iwp
    ezspin_db(318)%MassNum     = 0_iwp
    ezspin_db(318)%abundance    = 0.0_wp
    ezspin_db(318)%nucspin      = -1.0_wp
    ezspin_db(318)%gfactor      = 0.0_wp
!----------------------------
!--> Record 319    = 0Ra
    ezspin_db(319)%AtNumb      = 88_iwp
    ezspin_db(319)%MassNum     = 0_iwp
    ezspin_db(319)%abundance    = 0.0_wp
    ezspin_db(319)%nucspin      = -1.0_wp
    ezspin_db(319)%gfactor      = 0.0_wp
!----------------------------
!--> Record 320    = 227Ac
    ezspin_db(320)%AtNumb      = 89_iwp
    ezspin_db(320)%MassNum     = 227_iwp
    ezspin_db(320)%abundance    = 0.0_wp
    ezspin_db(320)%nucspin      = 1.5_wp
    ezspin_db(320)%gfactor      = +0.73_wp
!----------------------------
!--> Record 321    = 232Th
    ezspin_db(321)%AtNumb      = 90_iwp
    ezspin_db(321)%MassNum     = 232_iwp
    ezspin_db(321)%abundance    = 100.0_wp
    ezspin_db(321)%nucspin      = 0.0_wp
    ezspin_db(321)%gfactor      = 0.0_wp
!----------------------------
!--> Record 322    = 229Th
    ezspin_db(322)%AtNumb      = 90_iwp
    ezspin_db(322)%MassNum     = 229_iwp
    ezspin_db(322)%abundance    = 0.0_wp
    ezspin_db(322)%nucspin      = 2.5_wp
    ezspin_db(322)%gfactor      = +0.18_wp
!----------------------------
!--> Record 323    = 231Pa
    ezspin_db(323)%AtNumb      = 91_iwp
    ezspin_db(323)%MassNum     = 231_iwp
    ezspin_db(323)%abundance    = 100.0_wp
    ezspin_db(323)%nucspin      = 1.5_wp
    ezspin_db(323)%gfactor      = 1.34_wp
!----------------------------
!--> Record 324    = 238U
    ezspin_db(324)%AtNumb      = 92_iwp
    ezspin_db(324)%MassNum     = 238_iwp
    ezspin_db(324)%abundance    = 99.2745_wp
    ezspin_db(324)%nucspin      = 0.0_wp
    ezspin_db(324)%gfactor      = 0.0_wp
!----------------------------
!--> Record 325    = 235U
    ezspin_db(325)%AtNumb      = 92_iwp
    ezspin_db(325)%MassNum     = 235_iwp
    ezspin_db(325)%abundance    = 0.7200_wp
    ezspin_db(325)%nucspin      = 3.5_wp
    ezspin_db(325)%gfactor      = -0.109_wp
!----------------------------
!--> Record 326    = 234U
    ezspin_db(326)%AtNumb      = 92_iwp
    ezspin_db(326)%MassNum     = 234_iwp
    ezspin_db(326)%abundance    = 0.0055_wp
    ezspin_db(326)%nucspin      = 0.0_wp
    ezspin_db(326)%gfactor      = 0.0_wp
!----------------------------
!--> Record 327    = 237Np
    ezspin_db(327)%AtNumb      = 93_iwp
    ezspin_db(327)%MassNum     = 237_iwp
    ezspin_db(327)%abundance    = 0.0_wp
    ezspin_db(327)%nucspin      = 2.5_wp
    ezspin_db(327)%gfactor      = +1.256_wp
!----------------------------
!--> Record 328    = 239Pu
    ezspin_db(328)%AtNumb      = 94_iwp
    ezspin_db(328)%MassNum     = 239_iwp
    ezspin_db(328)%abundance    = 0.0_wp
    ezspin_db(328)%nucspin      = 0.5_wp
    ezspin_db(328)%gfactor      = +0.406_wp
!----------------------------
!--> Record 329    = 243Am
    ezspin_db(329)%AtNumb      = 95_iwp
    ezspin_db(329)%MassNum     = 243_iwp
    ezspin_db(329)%abundance    = 0.0_wp
    ezspin_db(329)%nucspin      = 2.5_wp
    ezspin_db(329)%gfactor      = +0.6_wp
!----------------------------
!--> Record 330    = 0Cm
    ezspin_db(330)%AtNumb      = 96_iwp
    ezspin_db(330)%MassNum     = 0_iwp
    ezspin_db(330)%abundance    = 0.0_wp
    ezspin_db(330)%nucspin      = -1.0_wp
    ezspin_db(330)%gfactor      = 0.0_wp
!----------------------------
!--> Record 331    = 0Bk
    ezspin_db(331)%AtNumb      = 97_iwp
    ezspin_db(331)%MassNum     = 0_iwp
    ezspin_db(331)%abundance    = 0.0_wp
    ezspin_db(331)%nucspin      = -1.0_wp
    ezspin_db(331)%gfactor      = 0.0_wp
!----------------------------
!--> Record 332    = 0Cf
    ezspin_db(332)%AtNumb      = 98_iwp
    ezspin_db(332)%MassNum     = 0_iwp
    ezspin_db(332)%abundance    = 0.0_wp
    ezspin_db(332)%nucspin      = -1.0_wp
    ezspin_db(332)%gfactor      = 0.0_wp
!----------------------------
!--> Record 333    = 0Es
    ezspin_db(333)%AtNumb      = 99_iwp
    ezspin_db(333)%MassNum     = 0_iwp
    ezspin_db(333)%abundance    = 0.0_wp
    ezspin_db(333)%nucspin      = -1.0_wp
    ezspin_db(333)%gfactor      = 0.0_wp
!----------------------------
!--> Record 334    = 0Fm
    ezspin_db(334)%AtNumb      = 100_iwp
    ezspin_db(334)%MassNum     = 0_iwp
    ezspin_db(334)%abundance    = 0.0_wp
    ezspin_db(334)%nucspin      = -1.0_wp
    ezspin_db(334)%gfactor      = 0.0_wp
!----------------------------
!--> Record 335    = 0Md
    ezspin_db(335)%AtNumb      = 101_iwp
    ezspin_db(335)%MassNum     = 0_iwp
    ezspin_db(335)%abundance    = 0.0_wp
    ezspin_db(335)%nucspin      = -1.0_wp
    ezspin_db(335)%gfactor      = 0.0_wp
!----------------------------
!--> Record 336    = 0No
    ezspin_db(336)%AtNumb      = 102_iwp
    ezspin_db(336)%MassNum     = 0_iwp
    ezspin_db(336)%abundance    = 0.0_wp
    ezspin_db(336)%nucspin      = -1.0_wp
    ezspin_db(336)%gfactor      = 0.0_wp
!----------------------------
!--> Record 337    = 0Lr
    ezspin_db(337)%AtNumb      = 103_iwp
    ezspin_db(337)%MassNum     = 0_iwp
    ezspin_db(337)%abundance    = 0.0_wp
    ezspin_db(337)%nucspin      = -1.0_wp
    ezspin_db(337)%gfactor      = 0.0_wp
!----------------------------
!--> Record 338    = 0Rf
    ezspin_db(338)%AtNumb      = 104_iwp
    ezspin_db(338)%MassNum     = 0_iwp
    ezspin_db(338)%abundance    = 0.0_wp
    ezspin_db(338)%nucspin      = -1.0_wp
    ezspin_db(338)%gfactor      = 0.0_wp
!----------------------------
!--> Record 339    = 0Db
    ezspin_db(339)%AtNumb      = 105_iwp
    ezspin_db(339)%MassNum     = 0_iwp
    ezspin_db(339)%abundance    = 0.0_wp
    ezspin_db(339)%nucspin      = -1.0_wp
    ezspin_db(339)%gfactor      = 0.0_wp
!----------------------------
!--> Record 340    = 0Sg
    ezspin_db(340)%AtNumb      = 106_iwp
    ezspin_db(340)%MassNum     = 0_iwp
    ezspin_db(340)%abundance    = 0.0_wp
    ezspin_db(340)%nucspin      = -1.0_wp
    ezspin_db(340)%gfactor      = 0.0_wp
!----------------------------
!--> Record 341    = 0Bh
    ezspin_db(341)%AtNumb      = 107_iwp
    ezspin_db(341)%MassNum     = 0_iwp
    ezspin_db(341)%abundance    = 0.0_wp
    ezspin_db(341)%nucspin      = -1.0_wp
    ezspin_db(341)%gfactor      = 0.0_wp
!----------------------------
!--> Record 342    = 0Hs
    ezspin_db(342)%AtNumb      = 108_iwp
    ezspin_db(342)%MassNum     = 0_iwp
    ezspin_db(342)%abundance    = 0.0_wp
    ezspin_db(342)%nucspin      = -1.0_wp
    ezspin_db(342)%gfactor      = 0.0_wp
!----------------------------
!--> Record 343    = 0Mt
    ezspin_db(343)%AtNumb      = 109_iwp
    ezspin_db(343)%MassNum     = 0_iwp
    ezspin_db(343)%abundance    = 0.0_wp
    ezspin_db(343)%nucspin      = -1.0_wp
    ezspin_db(343)%gfactor      = 0.0_wp
!----------------------------
!--> Record 344    = 0Ds
    ezspin_db(344)%AtNumb      = 110_iwp
    ezspin_db(344)%MassNum     = 0_iwp
    ezspin_db(344)%abundance    = 0.0_wp
    ezspin_db(344)%nucspin      = -1.0_wp
    ezspin_db(344)%gfactor      = 0.0_wp
!----------------------------
!--> Record 345    = 0Rg
    ezspin_db(345)%AtNumb      = 111_iwp
    ezspin_db(345)%MassNum     = 0_iwp
    ezspin_db(345)%abundance    = 0.0_wp
    ezspin_db(345)%nucspin      = -1.0_wp
    ezspin_db(345)%gfactor      = 0.0_wp
!----------------------------
!--> Record 346    = 0Cn
    ezspin_db(346)%AtNumb      = 112_iwp
    ezspin_db(346)%MassNum     = 0_iwp
    ezspin_db(346)%abundance    = 0.0_wp
    ezspin_db(346)%nucspin      = -1.0_wp
    ezspin_db(346)%gfactor      = 0.0_wp
!----------------------------
!--> Record 347    = 0Nh
    ezspin_db(347)%AtNumb      = 113_iwp
    ezspin_db(347)%MassNum     = 0_iwp
    ezspin_db(347)%abundance    = 0.0_wp
    ezspin_db(347)%nucspin      = -1.0_wp
    ezspin_db(347)%gfactor      = 0.0_wp
!----------------------------
!--> Record 348    = 0Fl
    ezspin_db(348)%AtNumb      = 114_iwp
    ezspin_db(348)%MassNum     = 0_iwp
    ezspin_db(348)%abundance    = 0.0_wp
    ezspin_db(348)%nucspin      = -1.0_wp
    ezspin_db(348)%gfactor      = 0.0_wp
!----------------------------
!--> Record 349    = 0Mc
    ezspin_db(349)%AtNumb      = 115_iwp
    ezspin_db(349)%MassNum     = 0_iwp
    ezspin_db(349)%abundance    = 0.0_wp
    ezspin_db(349)%nucspin      = -1.0_wp
    ezspin_db(349)%gfactor      = 0.0_wp
!----------------------------
!--> Record 350    = 0Lv
    ezspin_db(350)%AtNumb      = 116_iwp
    ezspin_db(350)%MassNum     = 0_iwp
    ezspin_db(350)%abundance    = 0.0_wp
    ezspin_db(350)%nucspin      = -1.0_wp
    ezspin_db(350)%gfactor      = 0.0_wp
!----------------------------
!--> Record 351    = 0Ts
    ezspin_db(351)%AtNumb      = 117_iwp
    ezspin_db(351)%MassNum     = 0_iwp
    ezspin_db(351)%abundance    = 0.0_wp
    ezspin_db(351)%nucspin      = -1.0_wp
    ezspin_db(351)%gfactor      = 0.0_wp
!----------------------------
!--> Record 352    = 0Og
    ezspin_db(352)%AtNumb      = 118_iwp
    ezspin_db(352)%MassNum     = 0_iwp
    ezspin_db(352)%abundance    = 0.0_wp
    ezspin_db(352)%nucspin      = -1.0_wp
    ezspin_db(352)%gfactor      = 0.0_wp
!----------------------------

  end subroutine init_isotope_data
!======================================================================
  function gfactor_by_nucspin(AtNumb,nucspin) result(gfactor)
    integer(kind=iwp),intent(in) :: AtNumb
    real(kind=wp), intent(in)    :: nucspin
    real(kind=wp)                :: gfactor
    integer(kind=iwp)            :: iRec, first_rec, last_rec
    logical(kind=iwp)            :: not_found_gfac


    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    gfactor = 0.0_wp
    not_found_gfac = .true.
    do iRec=first_rec,last_rec
      if (abs(ezspin_db(iRec)%nucspin - nucspin) < 1e-5_wp) then
        gfactor = ezspin_db(iRec)%gfactor
        exit
      endif
    enddo

    if (not_found_gfac) call warning_not_found(1)
    call isot_info(AtNumb, gfactor, nucspin)

  end function gfactor_by_nucspin
!======================================================================


!======================================================================
  function nucspin_by_gfactor(AtNumb,gfactor) result(nucspin)
    integer(kind=iwp),intent(in) :: AtNumb
    real(kind=wp), intent(in)    :: gfactor
    real(kind=wp)                :: nucspin
    integer(kind=iwp)            :: iRec, first_rec, last_rec
    logical(kind=iwp)            :: not_found_spin

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    nucspin = 0.5_wp
    not_found_spin = .true.
    do iRec=first_rec,last_rec
      if (abs(ezspin_db(iRec)%gfactor - gfactor) < 0.001_wp) then
        not_found_spin = .false.
        nucspin = ezspin_db(iRec)%nucspin
        exit
      endif
    enddo

    if (not_found_spin) call warning_not_found(0)
    call isot_info(AtNumb, gfactor, nucspin)

  end function nucspin_by_gfactor
!======================================================================


!======================================================================
  subroutine gfac_spin_by_mass(AtNumb,MassNum,gfactor,nucspin)
    integer(kind=iwp),intent(in) :: AtNumb, MassNum
    real(kind=wp), intent(out)   :: nucspin, gfactor
    integer(kind=iwp)            :: iRec, first_rec, last_rec
    logical(kind=iwp)            :: not_found_gfac , not_found_spin

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    nucspin = 0.5_wp
    gfactor = 0.0_wp
    not_found_gfac = .true.
    not_found_spin = .true.
    do iRec=first_rec,last_rec
      if (ezspin_db(iRec)%MassNum == MassNum) then
        not_found_gfac = .false.
        not_found_spin = .false.
        nucspin = ezspin_db(iRec)%nucspin
        gfactor = ezspin_db(iRec)%gfactor
        exit
      endif
    enddo


    if (not_found_spin) call warning_not_found(0)
    if (not_found_gfac) call warning_not_found(1)
    call isot_info(AtNumb, gfactor, nucspin, MassNum)

  end subroutine gfac_spin_by_mass
!======================================================================


!======================================================================
  subroutine get_first_nonzero_gfactor(AtNumb, gfactor,nucspin)
    integer(kind=iwp),intent(in)           :: AtNumb
    real(kind=wp), intent(out)             :: gfactor, nucspin
    integer(kind=iwp)                      :: iRec, first_rec, last_rec
    logical(kind=iwp)                      :: not_found_gfac, not_found_spin

    first_rec = elem_map_db(AtNumb)%first_rec
    last_rec  = elem_map_db(AtNumb)%last_rec

    nucspin = 0.5_wp
    gfactor = 0.0_wp
    not_found_gfac = .true.
    not_found_spin = .true.
    do iRec=first_rec,last_rec
      if (ezspin_db(iRec)%gfactor /= 0.0_wp) then
        not_found_gfac = .false.
        not_found_spin = .false.
        gfactor = ezspin_db(iRec)%gfactor
        nucspin = ezspin_db(iRec)%nucspin
        exit
      endif
    enddo

    if (not_found_spin) call warning_not_found(0)
    if (not_found_gfac) call warning_not_found(1)
    call isot_info(AtNumb, gfactor, nucspin)

  end subroutine get_first_nonzero_gfactor
!======================================================================


!======================================================================
  subroutine isot_info(AtNumb, gfactor, nucspin,MassNum)
    integer(kind=iwp), intent(in) :: AtNumb
    integer(kind=iwp), intent(in), optional :: MassNum
    real(kind=wp), intent(in)     :: gfactor, nucspin

    write(6,'(11X,A18)')    repeat('_', 18)
    write(6,'(11X,A18,I0)')    'INFO          Z = ', AtNumb
    if (present(MassNum)) then
      write(6,'(11X,A18,I0)')    '              A = ', MassNum
    endif
    write(6,'(11X,A18,F12.8)') '       g-factor = ', gfactor
    write(6,'(11X,A18,F5.2)')  '        nucspin = ', nucspin
    write(6,*) ''
  end subroutine isot_info
!======================================================================


!======================================================================
  subroutine warning_not_found(option)
    integer(kind=iwp), intent(in) :: option

    if(option == 0) then
      write(6,'(11X,A37)')  'WARNING: Not found a matching record!'
      write(6,'(20X,A61)')  'Defaulting to nucspin = 0.5. This value might be fictitious. '
    else if(option == 1) then
      write(6,'(11X,A37)')  'WARNING: Not found a matching record!'
      write(6,'(20X,A61)')  'Defaulting to g-factor = 0.0. This value might be fictitious.'
    endif
  end subroutine warning_not_found
!======================================================================


!======================================================================
  subroutine free_isotope_data()
    integer :: istat

    deallocate(ezspin_db, stat=istat)
    if (istat > 0) then
      write(6,*) ' Error in deallocating EasySpin ezspin_db.'
    end if

    deallocate(elem_map_db, stat=istat)
    if (istat > 0) then
      write(6,*) ' Error in deallocating EasySpin elem_map_db.'
    end if

  end subroutine free_isotope_data

end module hfc_data
