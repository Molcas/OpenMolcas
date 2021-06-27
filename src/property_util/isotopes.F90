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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************
!
! Isotope numbers and masses
!
! Taken from NIST: "Atomic Weights and Isotopic Compositions with
! Relative Atomic Masses" (https://www.nist.gov/pml/data/comp.cfm)
! Last updated: March 2017 (Version 4.1, July 2015)
!
! Each item in the the array ElementList contains data for an element:
!  - %Symbol: symbol
!  - %Natural: number of natural ocurring isotopes
!  - %Isotopes: array with all the isotopes for the element, sorted
!               by order of abundance (most to least, artificial
!               sorted by increasing mass number); elements with no
!               natural isotopes have the most stable first. Each item
!               in this array contains:
!    - %A: mass number (protons + neutrons)
!    - %m: isotopic mass in Da
!
! The "default" isotope for each element is simply the first item in
! the %Isotopes member.
!
! Manual changes from NIST data:
!  - Definitive symbols for all elements
!  - Most stable isotope (from Wikipedia) selected for Z > 94
!  - 3H and 14C included as natural

module Isotopes

use stdalloc, only: mma_Allocate, mma_Deallocate
use Definitions, only: wp, iwp, u6

implicit none
private
type Iso_t
  integer(kind=iwp) :: A
  real(kind=wp) :: m
end type Iso_t
type Element_t
  character(len=2) :: Symbol
  integer(kind=iwp) :: Z, Natural
  type(Iso_t), allocatable :: Isotopes(:)
end type Element_t
integer(kind=iwp), parameter :: MaxAtomNum = 118
type(Element_t), allocatable :: ElementList(:)
character(len=2), parameter :: PTab(MaxAtomNum) = [' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
                                                   'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
                                                   'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                                                   'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr', &
                                                   'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                                                   'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                                                   'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                                                   'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
                                                   'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                                                   'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                                                   'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
                                                   'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' &
                                                  ]
#include "constants2.fh"

interface Isotope
  module procedure Isotope_sym, Isotope_num
end interface Isotope

protected :: ElementList
public :: MaxAtomNum, Isotope, ElementList, Initialize_Isotopes, Free_Isotopes, NuclideMass, PTab

! Private extensions to mma interfaces

interface cptr2loff
  module procedure elm_cptr2loff
  module procedure iso_cptr2loff
end interface
interface mma_Allocate
  module procedure element_mma_allo_1D, element_mma_allo_1D_lim
  module procedure isotope_mma_allo_1D, isotope_mma_allo_1D_lim
end interface
interface mma_Deallocate
  module procedure element_mma_free_1D
  module procedure isotope_mma_free_1D
end interface

contains

! This subroutine allocates and fills the data in ElementList
! Since each array has a different size, it has to be done dynamically

subroutine Initialize_Isotopes()
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A)
      import :: Iso_t
      type(Iso_t), allocatable :: A(:)
    end subroutine c_null_alloc
  end interface
  integer(kind=iwp) :: i
# endif
  if (allocated(ElementList)) return
  call mma_Allocate(ElementList,MaxAtomNum,'ElmList')
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  do i=1,size(ElementList,1)
    call c_null_alloc(ElementList(i)%Isotopes)
  end do
# endif

# include "macros.fh"
  unused_proc(mma_allocate(ElementList,[0,0]))
  unused_proc(mma_allocate(ElementList(0)%Isotopes,[0,0]))

  ElementList(1)%Symbol = adjustl(PTab(1)) ! H
  ElementList(1)%Natural = 3
  call mma_Allocate(ElementList(1)%Isotopes,7)
  ElementList(1)%Isotopes(:) = [ &
                                Iso_t(1,1.00782503223_wp), &
                                Iso_t(2,2.01410177812_wp), &
                                Iso_t(3,3.0160492779_wp), &
                                Iso_t(4,4.02643_wp), &
                                Iso_t(5,5.035311_wp), &
                                Iso_t(6,6.04496_wp), &
                                Iso_t(7,7.0527_wp)]

  ElementList(2)%Symbol = adjustl(PTab(2)) ! He
  ElementList(2)%Natural = 2
  call mma_Allocate(ElementList(2)%Isotopes,8)
  ElementList(2)%Isotopes(:) = [ &
                                Iso_t(4,4.00260325413_wp), &
                                Iso_t(3,3.0160293201_wp), &
                                Iso_t(5,5.012057_wp), &
                                Iso_t(6,6.018885891_wp), &
                                Iso_t(7,7.0279907_wp), &
                                Iso_t(8,8.03393439_wp), &
                                Iso_t(9,9.043946_wp), &
                                Iso_t(10,10.05279_wp)]

  ElementList(3)%Symbol = adjustl(PTab(3)) ! Li
  ElementList(3)%Natural = 2
  call mma_Allocate(ElementList(3)%Isotopes,11)
  ElementList(3)%Isotopes(:) = [ &
                                Iso_t(7,7.0160034366_wp), &
                                Iso_t(6,6.0151228874_wp), &
                                Iso_t(3,3.0308_wp), &
                                Iso_t(4,4.02719_wp), &
                                Iso_t(5,5.012538_wp), &
                                Iso_t(8,8.022486246_wp), &
                                Iso_t(9,9.02679019_wp), &
                                Iso_t(10,10.035483_wp), &
                                Iso_t(11,11.04372358_wp), &
                                Iso_t(12,12.052517_wp), &
                                Iso_t(13,13.06263_wp)]

  ElementList(4)%Symbol = adjustl(PTab(4)) ! Be
  ElementList(4)%Natural = 1
  call mma_Allocate(ElementList(4)%Isotopes,12)
  ElementList(4)%Isotopes(:) = [ &
                                Iso_t(9,9.012183065_wp), &
                                Iso_t(5,5.0399_wp), &
                                Iso_t(6,6.0197264_wp), &
                                Iso_t(7,7.016928717_wp), &
                                Iso_t(8,8.005305102_wp), &
                                Iso_t(10,10.013534695_wp), &
                                Iso_t(11,11.02166108_wp), &
                                Iso_t(12,12.0269221_wp), &
                                Iso_t(13,13.036135_wp), &
                                Iso_t(14,14.04289_wp), &
                                Iso_t(15,15.05342_wp), &
                                Iso_t(16,16.06167_wp)]

  ElementList(5)%Symbol = adjustl(PTab(5)) ! B
  ElementList(5)%Natural = 2
  call mma_Allocate(ElementList(5)%Isotopes,16)
  ElementList(5)%Isotopes(:) = [ &
                                Iso_t(11,11.00930536_wp), &
                                Iso_t(10,10.01293695_wp), &
                                Iso_t(6,6.0508_wp), &
                                Iso_t(7,7.029712_wp), &
                                Iso_t(8,8.0246073_wp), &
                                Iso_t(9,9.01332965_wp), &
                                Iso_t(12,12.0143527_wp), &
                                Iso_t(13,13.0177802_wp), &
                                Iso_t(14,14.025404_wp), &
                                Iso_t(15,15.031088_wp), &
                                Iso_t(16,16.039842_wp), &
                                Iso_t(17,17.04699_wp), &
                                Iso_t(18,18.05566_wp), &
                                Iso_t(19,19.0631_wp), &
                                Iso_t(20,20.07207_wp), &
                                Iso_t(21,21.08129_wp)]

  ElementList(6)%Symbol = adjustl(PTab(6)) ! C
  ElementList(6)%Natural = 3
  call mma_Allocate(ElementList(6)%Isotopes,16)
  ElementList(6)%Isotopes(:) = [ &
                                Iso_t(12,12.0_wp), &
                                Iso_t(13,13.00335483507_wp), &
                                Iso_t(14,14.0032419884_wp), &
                                Iso_t(8,8.037643_wp), &
                                Iso_t(9,9.0310372_wp), &
                                Iso_t(10,10.01685331_wp), &
                                Iso_t(11,11.0114336_wp), &
                                Iso_t(15,15.01059926_wp), &
                                Iso_t(16,16.0147013_wp), &
                                Iso_t(17,17.022577_wp), &
                                Iso_t(18,18.026751_wp), &
                                Iso_t(19,19.0348_wp), &
                                Iso_t(20,20.04032_wp), &
                                Iso_t(21,21.049_wp), &
                                Iso_t(22,22.05753_wp), &
                                Iso_t(23,23.0689_wp)]

  ElementList(7)%Symbol = adjustl(PTab(7)) ! N
  ElementList(7)%Natural = 2
  call mma_Allocate(ElementList(7)%Isotopes,16)
  ElementList(7)%Isotopes(:) = [ &
                                Iso_t(14,14.00307400443_wp), &
                                Iso_t(15,15.00010889888_wp), &
                                Iso_t(10,10.04165_wp), &
                                Iso_t(11,11.026091_wp), &
                                Iso_t(12,12.0186132_wp), &
                                Iso_t(13,13.00573861_wp), &
                                Iso_t(16,16.0061019_wp), &
                                Iso_t(17,17.008449_wp), &
                                Iso_t(18,18.014078_wp), &
                                Iso_t(19,19.017022_wp), &
                                Iso_t(20,20.023366_wp), &
                                Iso_t(21,21.02711_wp), &
                                Iso_t(22,22.03439_wp), &
                                Iso_t(23,23.04114_wp), &
                                Iso_t(24,24.05039_wp), &
                                Iso_t(25,25.0601_wp)]

  ElementList(8)%Symbol = adjustl(PTab(8)) ! O
  ElementList(8)%Natural = 3
  call mma_Allocate(ElementList(8)%Isotopes,17)
  ElementList(8)%Isotopes(:) = [ &
                                Iso_t(16,15.99491461957_wp), &
                                Iso_t(18,17.99915961286_wp), &
                                Iso_t(17,16.9991317565_wp), &
                                Iso_t(12,12.034262_wp), &
                                Iso_t(13,13.024815_wp), &
                                Iso_t(14,14.00859636_wp), &
                                Iso_t(15,15.00306562_wp), &
                                Iso_t(19,19.003578_wp), &
                                Iso_t(20,20.00407535_wp), &
                                Iso_t(21,21.008655_wp), &
                                Iso_t(22,22.009966_wp), &
                                Iso_t(23,23.015696_wp), &
                                Iso_t(24,24.01986_wp), &
                                Iso_t(25,25.02936_wp), &
                                Iso_t(26,26.03729_wp), &
                                Iso_t(27,27.04772_wp), &
                                Iso_t(28,28.05591_wp)]

  ElementList(9)%Symbol = adjustl(PTab(9)) ! F
  ElementList(9)%Natural = 1
  call mma_Allocate(ElementList(9)%Isotopes,18)
  ElementList(9)%Isotopes(:) = [ &
                                Iso_t(19,18.99840316273_wp), &
                                Iso_t(14,14.034315_wp), &
                                Iso_t(15,15.018043_wp), &
                                Iso_t(16,16.0114657_wp), &
                                Iso_t(17,17.00209524_wp), &
                                Iso_t(18,18.00093733_wp), &
                                Iso_t(20,19.999981252_wp), &
                                Iso_t(21,20.9999489_wp), &
                                Iso_t(22,22.002999_wp), &
                                Iso_t(23,23.003557_wp), &
                                Iso_t(24,24.008115_wp), &
                                Iso_t(25,25.012199_wp), &
                                Iso_t(26,26.020038_wp), &
                                Iso_t(27,27.02644_wp), &
                                Iso_t(28,28.03534_wp), &
                                Iso_t(29,29.04254_wp), &
                                Iso_t(30,30.05165_wp), &
                                Iso_t(31,31.05971_wp)]

  ElementList(10)%Symbol = adjustl(PTab(10)) ! Ne
  ElementList(10)%Natural = 3
  call mma_Allocate(ElementList(10)%Isotopes,19)
  ElementList(10)%Isotopes(:) = [ &
                                 Iso_t(20,19.9924401762_wp), &
                                 Iso_t(22,21.991385114_wp), &
                                 Iso_t(21,20.993846685_wp), &
                                 Iso_t(16,16.02575_wp), &
                                 Iso_t(17,17.01771396_wp), &
                                 Iso_t(18,18.0057087_wp), &
                                 Iso_t(19,19.00188091_wp), &
                                 Iso_t(23,22.99446691_wp), &
                                 Iso_t(24,23.99361065_wp), &
                                 Iso_t(25,24.997789_wp), &
                                 Iso_t(26,26.000515_wp), &
                                 Iso_t(27,27.007553_wp), &
                                 Iso_t(28,28.01212_wp), &
                                 Iso_t(29,29.01975_wp), &
                                 Iso_t(30,30.02473_wp), &
                                 Iso_t(31,31.0331_wp), &
                                 Iso_t(32,32.03972_wp), &
                                 Iso_t(33,33.04938_wp), &
                                 Iso_t(34,34.05673_wp)]

  ElementList(11)%Symbol = adjustl(PTab(11)) ! Na
  ElementList(11)%Natural = 1
  call mma_Allocate(ElementList(11)%Isotopes,20)
  ElementList(11)%Isotopes(:) = [ &
                                 Iso_t(23,22.989769282_wp), &
                                 Iso_t(18,18.02688_wp), &
                                 Iso_t(19,19.01388_wp), &
                                 Iso_t(20,20.0073544_wp), &
                                 Iso_t(21,20.99765469_wp), &
                                 Iso_t(22,21.99443741_wp), &
                                 Iso_t(24,23.99096295_wp), &
                                 Iso_t(25,24.989954_wp), &
                                 Iso_t(26,25.9926346_wp), &
                                 Iso_t(27,26.9940765_wp), &
                                 Iso_t(28,27.998939_wp), &
                                 Iso_t(29,29.0028771_wp), &
                                 Iso_t(30,30.0090979_wp), &
                                 Iso_t(31,31.013163_wp), &
                                 Iso_t(32,32.02019_wp), &
                                 Iso_t(33,33.02573_wp), &
                                 Iso_t(34,34.03359_wp), &
                                 Iso_t(35,35.04062_wp), &
                                 Iso_t(36,36.04929_wp), &
                                 Iso_t(37,37.05705_wp)]

  ElementList(12)%Symbol = adjustl(PTab(12)) ! Mg
  ElementList(12)%Natural = 3
  call mma_Allocate(ElementList(12)%Isotopes,22)
  ElementList(12)%Isotopes(:) = [ &
                                 Iso_t(24,23.985041697_wp), &
                                 Iso_t(26,25.982592968_wp), &
                                 Iso_t(25,24.985836976_wp), &
                                 Iso_t(19,19.034169_wp), &
                                 Iso_t(20,20.01885_wp), &
                                 Iso_t(21,21.011716_wp), &
                                 Iso_t(22,21.99957065_wp), &
                                 Iso_t(23,22.99412421_wp), &
                                 Iso_t(27,26.984340624_wp), &
                                 Iso_t(28,27.9838767_wp), &
                                 Iso_t(29,28.988617_wp), &
                                 Iso_t(30,29.9904629_wp), &
                                 Iso_t(31,30.996648_wp), &
                                 Iso_t(32,31.9991102_wp), &
                                 Iso_t(33,33.0053271_wp), &
                                 Iso_t(34,34.008935_wp), &
                                 Iso_t(35,35.01679_wp), &
                                 Iso_t(36,36.02188_wp), &
                                 Iso_t(37,37.03037_wp), &
                                 Iso_t(38,38.03658_wp), &
                                 Iso_t(39,39.04538_wp), &
                                 Iso_t(40,40.05218_wp)]

  ElementList(13)%Symbol = adjustl(PTab(13)) ! Al
  ElementList(13)%Natural = 1
  call mma_Allocate(ElementList(13)%Isotopes,23)
  ElementList(13)%Isotopes(:) = [ &
                                 Iso_t(27,26.98153853_wp), &
                                 Iso_t(21,21.02897_wp), &
                                 Iso_t(22,22.01954_wp), &
                                 Iso_t(23,23.00724435_wp), &
                                 Iso_t(24,23.9999489_wp), &
                                 Iso_t(25,24.9904281_wp), &
                                 Iso_t(26,25.986891904_wp), &
                                 Iso_t(28,27.98191021_wp), &
                                 Iso_t(29,28.9804565_wp), &
                                 Iso_t(30,29.98296_wp), &
                                 Iso_t(31,30.983945_wp), &
                                 Iso_t(32,31.988085_wp), &
                                 Iso_t(33,32.990909_wp), &
                                 Iso_t(34,33.996705_wp), &
                                 Iso_t(35,34.999764_wp), &
                                 Iso_t(36,36.00639_wp), &
                                 Iso_t(37,37.01053_wp), &
                                 Iso_t(38,38.0174_wp), &
                                 Iso_t(39,39.02254_wp), &
                                 Iso_t(40,40.03003_wp), &
                                 Iso_t(41,41.03638_wp), &
                                 Iso_t(42,42.04384_wp), &
                                 Iso_t(43,43.05147_wp)]

  ElementList(14)%Symbol = adjustl(PTab(14)) ! Si
  ElementList(14)%Natural = 3
  call mma_Allocate(ElementList(14)%Isotopes,24)
  ElementList(14)%Isotopes(:) = [ &
                                 Iso_t(28,27.97692653465_wp), &
                                 Iso_t(29,28.9764946649_wp), &
                                 Iso_t(30,29.973770136_wp), &
                                 Iso_t(22,22.03579_wp), &
                                 Iso_t(23,23.02544_wp), &
                                 Iso_t(24,24.011535_wp), &
                                 Iso_t(25,25.004109_wp), &
                                 Iso_t(26,25.99233384_wp), &
                                 Iso_t(27,26.98670481_wp), &
                                 Iso_t(31,30.975363194_wp), &
                                 Iso_t(32,31.97415154_wp), &
                                 Iso_t(33,32.97797696_wp), &
                                 Iso_t(34,33.978576_wp), &
                                 Iso_t(35,34.984583_wp), &
                                 Iso_t(36,35.986695_wp), &
                                 Iso_t(37,36.992921_wp), &
                                 Iso_t(38,37.995523_wp), &
                                 Iso_t(39,39.002491_wp), &
                                 Iso_t(40,40.00583_wp), &
                                 Iso_t(41,41.01301_wp), &
                                 Iso_t(42,42.01778_wp), &
                                 Iso_t(43,43.0248_wp), &
                                 Iso_t(44,44.03061_wp), &
                                 Iso_t(45,45.03995_wp)]

  ElementList(15)%Symbol = adjustl(PTab(15)) ! P
  ElementList(15)%Natural = 1
  call mma_Allocate(ElementList(15)%Isotopes,24)
  ElementList(15)%Isotopes(:) = [ &
                                 Iso_t(31,30.97376199842_wp), &
                                 Iso_t(24,24.03577_wp), &
                                 Iso_t(25,25.02119_wp), &
                                 Iso_t(26,26.01178_wp), &
                                 Iso_t(27,26.999224_wp), &
                                 Iso_t(28,27.9923266_wp), &
                                 Iso_t(29,28.98180079_wp), &
                                 Iso_t(30,29.97831375_wp), &
                                 Iso_t(32,31.973907643_wp), &
                                 Iso_t(33,32.9717257_wp), &
                                 Iso_t(34,33.97364589_wp), &
                                 Iso_t(35,34.9733141_wp), &
                                 Iso_t(36,35.97826_wp), &
                                 Iso_t(37,36.979607_wp), &
                                 Iso_t(38,37.984252_wp), &
                                 Iso_t(39,38.986227_wp), &
                                 Iso_t(40,39.99133_wp), &
                                 Iso_t(41,40.994654_wp), &
                                 Iso_t(42,42.00108_wp), &
                                 Iso_t(43,43.00502_wp), &
                                 Iso_t(44,44.01121_wp), &
                                 Iso_t(45,45.01645_wp), &
                                 Iso_t(46,46.02446_wp), &
                                 Iso_t(47,47.03139_wp)]

  ElementList(16)%Symbol = adjustl(PTab(16)) ! S
  ElementList(16)%Natural = 4
  call mma_Allocate(ElementList(16)%Isotopes,24)
  ElementList(16)%Isotopes(:) = [ &
                                 Iso_t(32,31.9720711744_wp), &
                                 Iso_t(34,33.967867004_wp), &
                                 Iso_t(33,32.9714589098_wp), &
                                 Iso_t(36,35.96708071_wp), &
                                 Iso_t(26,26.02907_wp), &
                                 Iso_t(27,27.01828_wp), &
                                 Iso_t(28,28.00437_wp), &
                                 Iso_t(29,28.996611_wp), &
                                 Iso_t(30,29.98490703_wp), &
                                 Iso_t(31,30.97955701_wp), &
                                 Iso_t(35,34.96903231_wp), &
                                 Iso_t(37,36.97112551_wp), &
                                 Iso_t(38,37.9711633_wp), &
                                 Iso_t(39,38.975134_wp), &
                                 Iso_t(40,39.9754826_wp), &
                                 Iso_t(41,40.9795935_wp), &
                                 Iso_t(42,41.9810651_wp), &
                                 Iso_t(43,42.9869076_wp), &
                                 Iso_t(44,43.9901188_wp), &
                                 Iso_t(45,44.99572_wp), &
                                 Iso_t(46,46.00004_wp), &
                                 Iso_t(47,47.00795_wp), &
                                 Iso_t(48,48.0137_wp), &
                                 Iso_t(49,49.02276_wp)]

  ElementList(17)%Symbol = adjustl(PTab(17)) ! Cl
  ElementList(17)%Natural = 2
  call mma_Allocate(ElementList(17)%Isotopes,24)
  ElementList(17)%Isotopes(:) = [ &
                                 Iso_t(35,34.968852682_wp), &
                                 Iso_t(37,36.965902602_wp), &
                                 Iso_t(28,28.02954_wp), &
                                 Iso_t(29,29.01478_wp), &
                                 Iso_t(30,30.00477_wp), &
                                 Iso_t(31,30.992414_wp), &
                                 Iso_t(32,31.98568464_wp), &
                                 Iso_t(33,32.97745199_wp), &
                                 Iso_t(34,33.973762485_wp), &
                                 Iso_t(36,35.968306809_wp), &
                                 Iso_t(38,37.96801044_wp), &
                                 Iso_t(39,38.9680082_wp), &
                                 Iso_t(40,39.970415_wp), &
                                 Iso_t(41,40.970685_wp), &
                                 Iso_t(42,41.97325_wp), &
                                 Iso_t(43,42.97389_wp), &
                                 Iso_t(44,43.97787_wp), &
                                 Iso_t(45,44.98029_wp), &
                                 Iso_t(46,45.98517_wp), &
                                 Iso_t(47,46.98916_wp), &
                                 Iso_t(48,47.99564_wp), &
                                 Iso_t(49,49.00123_wp), &
                                 Iso_t(50,50.00905_wp), &
                                 Iso_t(51,51.01554_wp)]

  ElementList(18)%Symbol = adjustl(PTab(18)) ! Ar
  ElementList(18)%Natural = 3
  call mma_Allocate(ElementList(18)%Isotopes,24)
  ElementList(18)%Isotopes(:) = [ &
                                 Iso_t(40,39.9623831237_wp), &
                                 Iso_t(36,35.967545105_wp), &
                                 Iso_t(38,37.96273211_wp), &
                                 Iso_t(30,30.02307_wp), &
                                 Iso_t(31,31.01212_wp), &
                                 Iso_t(32,31.9976378_wp), &
                                 Iso_t(33,32.98992555_wp), &
                                 Iso_t(34,33.98027009_wp), &
                                 Iso_t(35,34.97525759_wp), &
                                 Iso_t(37,36.96677633_wp), &
                                 Iso_t(39,38.964313_wp), &
                                 Iso_t(41,40.96450057_wp), &
                                 Iso_t(42,41.9630457_wp), &
                                 Iso_t(43,42.9656361_wp), &
                                 Iso_t(44,43.9649238_wp), &
                                 Iso_t(45,44.96803973_wp), &
                                 Iso_t(46,45.968083_wp), &
                                 Iso_t(47,46.972935_wp), &
                                 Iso_t(48,47.97591_wp), &
                                 Iso_t(49,48.9819_wp), &
                                 Iso_t(50,49.98613_wp), &
                                 Iso_t(51,50.9937_wp), &
                                 Iso_t(52,51.99896_wp), &
                                 Iso_t(53,53.00729_wp)]

  ElementList(19)%Symbol = adjustl(PTab(19)) ! K
  ElementList(19)%Natural = 3
  call mma_Allocate(ElementList(19)%Isotopes,25)
  ElementList(19)%Isotopes(:) = [ &
                                 Iso_t(39,38.9637064864_wp), &
                                 Iso_t(41,40.9618252579_wp), &
                                 Iso_t(40,39.963998166_wp), &
                                 Iso_t(32,32.02265_wp), &
                                 Iso_t(33,33.00756_wp), &
                                 Iso_t(34,33.99869_wp), &
                                 Iso_t(35,34.98800541_wp), &
                                 Iso_t(36,35.98130201_wp), &
                                 Iso_t(37,36.97337589_wp), &
                                 Iso_t(38,37.96908112_wp), &
                                 Iso_t(42,41.96240231_wp), &
                                 Iso_t(43,42.9607347_wp), &
                                 Iso_t(44,43.96158699_wp), &
                                 Iso_t(45,44.96069149_wp), &
                                 Iso_t(46,45.96198159_wp), &
                                 Iso_t(47,46.9616616_wp), &
                                 Iso_t(48,47.96534119_wp), &
                                 Iso_t(49,48.96821075_wp), &
                                 Iso_t(50,49.97238_wp), &
                                 Iso_t(51,50.975828_wp), &
                                 Iso_t(52,51.98224_wp), &
                                 Iso_t(53,52.98746_wp), &
                                 Iso_t(54,53.99463_wp), &
                                 Iso_t(55,55.00076_wp), &
                                 Iso_t(56,56.00851_wp)]

  ElementList(20)%Symbol = adjustl(PTab(20)) ! Ca
  ElementList(20)%Natural = 6
  call mma_Allocate(ElementList(20)%Isotopes,25)
  ElementList(20)%Isotopes(:) = [ &
                                 Iso_t(40,39.962590863_wp), &
                                 Iso_t(44,43.95548156_wp), &
                                 Iso_t(42,41.95861783_wp), &
                                 Iso_t(48,47.95252276_wp), &
                                 Iso_t(43,42.95876644_wp), &
                                 Iso_t(46,45.953689_wp), &
                                 Iso_t(34,34.01487_wp), &
                                 Iso_t(35,35.00514_wp), &
                                 Iso_t(36,35.993074_wp), &
                                 Iso_t(37,36.98589785_wp), &
                                 Iso_t(38,37.97631922_wp), &
                                 Iso_t(39,38.97071081_wp), &
                                 Iso_t(41,40.96227792_wp), &
                                 Iso_t(45,44.95618635_wp), &
                                 Iso_t(47,46.9545424_wp), &
                                 Iso_t(49,48.95566274_wp), &
                                 Iso_t(50,49.9574992_wp), &
                                 Iso_t(51,50.960989_wp), &
                                 Iso_t(52,51.963217_wp), &
                                 Iso_t(53,52.96945_wp), &
                                 Iso_t(54,53.9734_wp), &
                                 Iso_t(55,54.9803_wp), &
                                 Iso_t(56,55.98508_wp), &
                                 Iso_t(57,56.99262_wp), &
                                 Iso_t(58,57.99794_wp)]

  ElementList(21)%Symbol = adjustl(PTab(21)) ! Sc
  ElementList(21)%Natural = 1
  call mma_Allocate(ElementList(21)%Isotopes,26)
  ElementList(21)%Isotopes(:) = [ &
                                 Iso_t(45,44.95590828_wp), &
                                 Iso_t(36,36.01648_wp), &
                                 Iso_t(37,37.00374_wp), &
                                 Iso_t(38,37.99512_wp), &
                                 Iso_t(39,38.984785_wp), &
                                 Iso_t(40,39.9779673_wp), &
                                 Iso_t(41,40.969251105_wp), &
                                 Iso_t(42,41.96551653_wp), &
                                 Iso_t(43,42.9611505_wp), &
                                 Iso_t(44,43.9594029_wp), &
                                 Iso_t(46,45.95516826_wp), &
                                 Iso_t(47,46.9524037_wp), &
                                 Iso_t(48,47.9522236_wp), &
                                 Iso_t(49,48.9500146_wp), &
                                 Iso_t(50,49.952176_wp), &
                                 Iso_t(51,50.953592_wp), &
                                 Iso_t(52,51.95688_wp), &
                                 Iso_t(53,52.95909_wp), &
                                 Iso_t(54,53.96393_wp), &
                                 Iso_t(55,54.96782_wp), &
                                 Iso_t(56,55.97345_wp), &
                                 Iso_t(57,56.97777_wp), &
                                 Iso_t(58,57.98403_wp), &
                                 Iso_t(59,58.98894_wp), &
                                 Iso_t(60,59.99565_wp), &
                                 Iso_t(61,61.001_wp)]

  ElementList(22)%Symbol = adjustl(PTab(22)) ! Ti
  ElementList(22)%Natural = 5
  call mma_Allocate(ElementList(22)%Isotopes,26)
  ElementList(22)%Isotopes(:) = [ &
                                 Iso_t(48,47.94794198_wp), &
                                 Iso_t(46,45.95262772_wp), &
                                 Iso_t(47,46.95175879_wp), &
                                 Iso_t(49,48.94786568_wp), &
                                 Iso_t(50,49.94478689_wp), &
                                 Iso_t(38,38.01145_wp), &
                                 Iso_t(39,39.00236_wp), &
                                 Iso_t(40,39.9905_wp), &
                                 Iso_t(41,40.983148_wp), &
                                 Iso_t(42,41.97304903_wp), &
                                 Iso_t(43,42.9685225_wp), &
                                 Iso_t(44,43.95968995_wp), &
                                 Iso_t(45,44.95812198_wp), &
                                 Iso_t(51,50.94661065_wp), &
                                 Iso_t(52,51.946893_wp), &
                                 Iso_t(53,52.94973_wp), &
                                 Iso_t(54,53.95105_wp), &
                                 Iso_t(55,54.95527_wp), &
                                 Iso_t(56,55.95791_wp), &
                                 Iso_t(57,56.96364_wp), &
                                 Iso_t(58,57.9666_wp), &
                                 Iso_t(59,58.97247_wp), &
                                 Iso_t(60,59.97603_wp), &
                                 Iso_t(61,60.98245_wp), &
                                 Iso_t(62,61.98651_wp), &
                                 Iso_t(63,62.99375_wp)]

  ElementList(23)%Symbol = adjustl(PTab(23)) ! V
  ElementList(23)%Natural = 2
  call mma_Allocate(ElementList(23)%Isotopes,27)
  ElementList(23)%Isotopes(:) = [ &
                                 Iso_t(51,50.94395704_wp), &
                                 Iso_t(50,49.94715601_wp), &
                                 Iso_t(40,40.01276_wp), &
                                 Iso_t(41,41.00021_wp), &
                                 Iso_t(42,41.99182_wp), &
                                 Iso_t(43,42.980766_wp), &
                                 Iso_t(44,43.97411_wp), &
                                 Iso_t(45,44.9657748_wp), &
                                 Iso_t(46,45.96019878_wp), &
                                 Iso_t(47,46.95490491_wp), &
                                 Iso_t(48,47.9522522_wp), &
                                 Iso_t(49,48.9485118_wp), &
                                 Iso_t(52,51.94477301_wp), &
                                 Iso_t(53,52.9443367_wp), &
                                 Iso_t(54,53.946439_wp), &
                                 Iso_t(55,54.94724_wp), &
                                 Iso_t(56,55.95048_wp), &
                                 Iso_t(57,56.95252_wp), &
                                 Iso_t(58,57.95672_wp), &
                                 Iso_t(59,58.95939_wp), &
                                 Iso_t(60,59.96431_wp), &
                                 Iso_t(61,60.96725_wp), &
                                 Iso_t(62,61.97265_wp), &
                                 Iso_t(63,62.97639_wp), &
                                 Iso_t(64,63.98264_wp), &
                                 Iso_t(65,64.9875_wp), &
                                 Iso_t(66,65.99398_wp)]

  ElementList(24)%Symbol = adjustl(PTab(24)) ! Cr
  ElementList(24)%Natural = 4
  call mma_Allocate(ElementList(24)%Isotopes,27)
  ElementList(24)%Isotopes(:) = [ &
                                 Iso_t(52,51.94050623_wp), &
                                 Iso_t(53,52.94064815_wp), &
                                 Iso_t(50,49.94604183_wp), &
                                 Iso_t(54,53.93887916_wp), &
                                 Iso_t(42,42.0067_wp), &
                                 Iso_t(43,42.99753_wp), &
                                 Iso_t(44,43.98536_wp), &
                                 Iso_t(45,44.97905_wp), &
                                 Iso_t(46,45.968359_wp), &
                                 Iso_t(47,46.9628974_wp), &
                                 Iso_t(48,47.9540291_wp), &
                                 Iso_t(49,48.9513333_wp), &
                                 Iso_t(51,50.94476502_wp), &
                                 Iso_t(55,54.94083843_wp), &
                                 Iso_t(56,55.9406531_wp), &
                                 Iso_t(57,56.943613_wp), &
                                 Iso_t(58,57.94435_wp), &
                                 Iso_t(59,58.94859_wp), &
                                 Iso_t(60,59.95008_wp), &
                                 Iso_t(61,60.95442_wp), &
                                 Iso_t(62,61.9561_wp), &
                                 Iso_t(63,62.96165_wp), &
                                 Iso_t(64,63.96408_wp), &
                                 Iso_t(65,64.96996_wp), &
                                 Iso_t(66,65.97366_wp), &
                                 Iso_t(67,66.98016_wp), &
                                 Iso_t(68,67.98403_wp)]

  ElementList(25)%Symbol = adjustl(PTab(25)) ! Mn
  ElementList(25)%Natural = 1
  call mma_Allocate(ElementList(25)%Isotopes,28)
  ElementList(25)%Isotopes(:) = [ &
                                 Iso_t(55,54.93804391_wp), &
                                 Iso_t(44,44.00715_wp), &
                                 Iso_t(45,44.99449_wp), &
                                 Iso_t(46,45.98609_wp), &
                                 Iso_t(47,46.975775_wp), &
                                 Iso_t(48,47.96852_wp), &
                                 Iso_t(49,48.959595_wp), &
                                 Iso_t(50,49.95423778_wp), &
                                 Iso_t(51,50.94820847_wp), &
                                 Iso_t(52,51.9455639_wp), &
                                 Iso_t(53,52.94128889_wp), &
                                 Iso_t(54,53.9403576_wp), &
                                 Iso_t(56,55.93890369_wp), &
                                 Iso_t(57,56.9382861_wp), &
                                 Iso_t(58,57.9400666_wp), &
                                 Iso_t(59,58.9403911_wp), &
                                 Iso_t(60,59.9431366_wp), &
                                 Iso_t(61,60.9444525_wp), &
                                 Iso_t(62,61.94795_wp), &
                                 Iso_t(63,62.9496647_wp), &
                                 Iso_t(64,63.9538494_wp), &
                                 Iso_t(65,64.9560198_wp), &
                                 Iso_t(66,65.960547_wp), &
                                 Iso_t(67,66.96424_wp), &
                                 Iso_t(68,67.96962_wp), &
                                 Iso_t(69,68.97366_wp), &
                                 Iso_t(70,69.97937_wp), &
                                 Iso_t(71,70.98368_wp)]

  ElementList(26)%Symbol = adjustl(PTab(26)) ! Fe
  ElementList(26)%Natural = 4
  call mma_Allocate(ElementList(26)%Isotopes,30)
  ElementList(26)%Isotopes(:) = [ &
                                 Iso_t(56,55.93493633_wp), &
                                 Iso_t(54,53.93960899_wp), &
                                 Iso_t(57,56.93539284_wp), &
                                 Iso_t(58,57.93327443_wp), &
                                 Iso_t(45,45.01442_wp), &
                                 Iso_t(46,46.00063_wp), &
                                 Iso_t(47,46.99185_wp), &
                                 Iso_t(48,47.98023_wp), &
                                 Iso_t(49,48.973429_wp), &
                                 Iso_t(50,49.962975_wp), &
                                 Iso_t(51,50.956841_wp), &
                                 Iso_t(52,51.9481131_wp), &
                                 Iso_t(53,52.9453064_wp), &
                                 Iso_t(55,54.93829199_wp), &
                                 Iso_t(59,58.93487434_wp), &
                                 Iso_t(60,59.9340711_wp), &
                                 Iso_t(61,60.9367462_wp), &
                                 Iso_t(62,61.9367918_wp), &
                                 Iso_t(63,62.9402727_wp), &
                                 Iso_t(64,63.9409878_wp), &
                                 Iso_t(65,64.9450115_wp), &
                                 Iso_t(66,65.94625_wp), &
                                 Iso_t(67,66.95054_wp), &
                                 Iso_t(68,67.95295_wp), &
                                 Iso_t(69,68.95807_wp), &
                                 Iso_t(70,69.96102_wp), &
                                 Iso_t(71,70.96672_wp), &
                                 Iso_t(72,71.96983_wp), &
                                 Iso_t(73,72.97572_wp), &
                                 Iso_t(74,73.97935_wp)]

  ElementList(27)%Symbol = adjustl(PTab(27)) ! Co
  ElementList(27)%Natural = 1
  call mma_Allocate(ElementList(27)%Isotopes,30)
  ElementList(27)%Isotopes(:) = [ &
                                 Iso_t(59,58.93319429_wp), &
                                 Iso_t(47,47.01057_wp), &
                                 Iso_t(48,48.00093_wp), &
                                 Iso_t(49,48.98891_wp), &
                                 Iso_t(50,49.98091_wp), &
                                 Iso_t(51,50.970647_wp), &
                                 Iso_t(52,51.96351_wp), &
                                 Iso_t(53,52.9542041_wp), &
                                 Iso_t(54,53.94845987_wp), &
                                 Iso_t(55,54.9419972_wp), &
                                 Iso_t(56,55.9398388_wp), &
                                 Iso_t(57,56.93629057_wp), &
                                 Iso_t(58,57.9357521_wp), &
                                 Iso_t(60,59.9338163_wp), &
                                 Iso_t(61,60.93247662_wp), &
                                 Iso_t(62,61.934059_wp), &
                                 Iso_t(63,62.9336_wp), &
                                 Iso_t(64,63.935811_wp), &
                                 Iso_t(65,64.9364621_wp), &
                                 Iso_t(66,65.939443_wp), &
                                 Iso_t(67,66.9406096_wp), &
                                 Iso_t(68,67.94426_wp), &
                                 Iso_t(69,68.94614_wp), &
                                 Iso_t(70,69.94963_wp), &
                                 Iso_t(71,70.95237_wp), &
                                 Iso_t(72,71.95729_wp), &
                                 Iso_t(73,72.96039_wp), &
                                 Iso_t(74,73.96515_wp), &
                                 Iso_t(75,74.96876_wp), &
                                 Iso_t(76,75.97413_wp)]

  ElementList(28)%Symbol = adjustl(PTab(28)) ! Ni
  ElementList(28)%Natural = 5
  call mma_Allocate(ElementList(28)%Isotopes,32)
  ElementList(28)%Isotopes(:) = [ &
                                 Iso_t(58,57.93534241_wp), &
                                 Iso_t(60,59.93078588_wp), &
                                 Iso_t(62,61.92834537_wp), &
                                 Iso_t(61,60.93105557_wp), &
                                 Iso_t(64,63.92796682_wp), &
                                 Iso_t(48,48.01769_wp), &
                                 Iso_t(49,49.0077_wp), &
                                 Iso_t(50,49.99474_wp), &
                                 Iso_t(51,50.98611_wp), &
                                 Iso_t(52,51.9748_wp), &
                                 Iso_t(53,52.96819_wp), &
                                 Iso_t(54,53.957892_wp), &
                                 Iso_t(55,54.95133063_wp), &
                                 Iso_t(56,55.94212855_wp), &
                                 Iso_t(57,56.93979218_wp), &
                                 Iso_t(59,58.9343462_wp), &
                                 Iso_t(63,62.92966963_wp), &
                                 Iso_t(65,64.93008517_wp), &
                                 Iso_t(66,65.9291393_wp), &
                                 Iso_t(67,66.9315694_wp), &
                                 Iso_t(68,67.9318688_wp), &
                                 Iso_t(69,68.9356103_wp), &
                                 Iso_t(70,69.9364313_wp), &
                                 Iso_t(71,70.940519_wp), &
                                 Iso_t(72,71.9417859_wp), &
                                 Iso_t(73,72.9462067_wp), &
                                 Iso_t(74,73.94798_wp), &
                                 Iso_t(75,74.9525_wp), &
                                 Iso_t(76,75.95533_wp), &
                                 Iso_t(77,76.96055_wp), &
                                 Iso_t(78,77.96336_wp), &
                                 Iso_t(79,78.97025_wp)]

  ElementList(29)%Symbol = adjustl(PTab(29)) ! Cu
  ElementList(29)%Natural = 2
  call mma_Allocate(ElementList(29)%Isotopes,31)
  ElementList(29)%Isotopes(:) = [ &
                                 Iso_t(63,62.92959772_wp), &
                                 Iso_t(65,64.9277897_wp), &
                                 Iso_t(52,51.99671_wp), &
                                 Iso_t(53,52.98459_wp), &
                                 Iso_t(54,53.97666_wp), &
                                 Iso_t(55,54.96604_wp), &
                                 Iso_t(56,55.95895_wp), &
                                 Iso_t(57,56.9492125_wp), &
                                 Iso_t(58,57.94453305_wp), &
                                 Iso_t(59,58.93949748_wp), &
                                 Iso_t(60,59.9373645_wp), &
                                 Iso_t(61,60.9334576_wp), &
                                 Iso_t(62,61.93259541_wp), &
                                 Iso_t(64,63.92976434_wp), &
                                 Iso_t(66,65.92886903_wp), &
                                 Iso_t(67,66.9277303_wp), &
                                 Iso_t(68,67.9296109_wp), &
                                 Iso_t(69,68.9294293_wp), &
                                 Iso_t(70,69.9323921_wp), &
                                 Iso_t(71,70.9326768_wp), &
                                 Iso_t(72,71.9358203_wp), &
                                 Iso_t(73,72.9366744_wp), &
                                 Iso_t(74,73.9398749_wp), &
                                 Iso_t(75,74.9415226_wp), &
                                 Iso_t(76,75.945275_wp), &
                                 Iso_t(77,76.94792_wp), &
                                 Iso_t(78,77.95223_wp), &
                                 Iso_t(79,78.95502_wp), &
                                 Iso_t(80,79.96089_wp), &
                                 Iso_t(81,80.96587_wp), &
                                 Iso_t(82,81.97244_wp)]

  ElementList(30)%Symbol = adjustl(PTab(30)) ! Zn
  ElementList(30)%Natural = 5
  call mma_Allocate(ElementList(30)%Isotopes,32)
  ElementList(30)%Isotopes(:) = [ &
                                 Iso_t(64,63.92914201_wp), &
                                 Iso_t(66,65.92603381_wp), &
                                 Iso_t(68,67.92484455_wp), &
                                 Iso_t(67,66.92712775_wp), &
                                 Iso_t(70,69.9253192_wp), &
                                 Iso_t(54,53.99204_wp), &
                                 Iso_t(55,54.98398_wp), &
                                 Iso_t(56,55.97254_wp), &
                                 Iso_t(57,56.96506_wp), &
                                 Iso_t(58,57.954591_wp), &
                                 Iso_t(59,58.94931266_wp), &
                                 Iso_t(60,59.9418421_wp), &
                                 Iso_t(61,60.939507_wp), &
                                 Iso_t(62,61.93433397_wp), &
                                 Iso_t(63,62.9332115_wp), &
                                 Iso_t(65,64.92924077_wp), &
                                 Iso_t(69,68.9265507_wp), &
                                 Iso_t(71,70.9277196_wp), &
                                 Iso_t(72,71.9268428_wp), &
                                 Iso_t(73,72.9295826_wp), &
                                 Iso_t(74,73.9294073_wp), &
                                 Iso_t(75,74.9328402_wp), &
                                 Iso_t(76,75.933115_wp), &
                                 Iso_t(77,76.9368872_wp), &
                                 Iso_t(78,77.9382892_wp), &
                                 Iso_t(79,78.9426381_wp), &
                                 Iso_t(80,79.9445529_wp), &
                                 Iso_t(81,80.9504026_wp), &
                                 Iso_t(82,81.95426_wp), &
                                 Iso_t(83,82.96056_wp), &
                                 Iso_t(84,83.96521_wp), &
                                 Iso_t(85,84.97226_wp)]

  ElementList(31)%Symbol = adjustl(PTab(31)) ! Ga
  ElementList(31)%Natural = 2
  call mma_Allocate(ElementList(31)%Isotopes,32)
  ElementList(31)%Isotopes(:) = [ &
                                 Iso_t(69,68.9255735_wp), &
                                 Iso_t(71,70.92470258_wp), &
                                 Iso_t(56,55.99536_wp), &
                                 Iso_t(57,56.9832_wp), &
                                 Iso_t(58,57.97478_wp), &
                                 Iso_t(59,58.96353_wp), &
                                 Iso_t(60,59.95729_wp), &
                                 Iso_t(61,60.949399_wp), &
                                 Iso_t(62,61.94419025_wp), &
                                 Iso_t(63,62.9392942_wp), &
                                 Iso_t(64,63.9368404_wp), &
                                 Iso_t(65,64.93273459_wp), &
                                 Iso_t(66,65.9315894_wp), &
                                 Iso_t(67,66.9282025_wp), &
                                 Iso_t(68,67.9279805_wp), &
                                 Iso_t(70,69.9260219_wp), &
                                 Iso_t(72,71.92636747_wp), &
                                 Iso_t(73,72.9251747_wp), &
                                 Iso_t(74,73.9269457_wp), &
                                 Iso_t(75,74.9265002_wp), &
                                 Iso_t(76,75.9288276_wp), &
                                 Iso_t(77,76.9291543_wp), &
                                 Iso_t(78,77.9316088_wp), &
                                 Iso_t(79,78.9328523_wp), &
                                 Iso_t(80,79.9364208_wp), &
                                 Iso_t(81,80.9381338_wp), &
                                 Iso_t(82,81.9431765_wp), &
                                 Iso_t(83,82.9471203_wp), &
                                 Iso_t(84,83.95246_wp), &
                                 Iso_t(85,84.95699_wp), &
                                 Iso_t(86,85.96301_wp), &
                                 Iso_t(87,86.96824_wp)]

  ElementList(32)%Symbol = adjustl(PTab(32)) ! Ge
  ElementList(32)%Natural = 5
  call mma_Allocate(ElementList(32)%Isotopes,33)
  ElementList(32)%Isotopes(:) = [ &
                                 Iso_t(74,73.921177761_wp), &
                                 Iso_t(72,71.922075826_wp), &
                                 Iso_t(70,69.92424875_wp), &
                                 Iso_t(73,72.923458956_wp), &
                                 Iso_t(76,75.921402726_wp), &
                                 Iso_t(58,57.99172_wp), &
                                 Iso_t(59,58.98249_wp), &
                                 Iso_t(60,59.97036_wp), &
                                 Iso_t(61,60.96379_wp), &
                                 Iso_t(62,61.95502_wp), &
                                 Iso_t(63,62.949628_wp), &
                                 Iso_t(64,63.9416899_wp), &
                                 Iso_t(65,64.9393681_wp), &
                                 Iso_t(66,65.9338621_wp), &
                                 Iso_t(67,66.9327339_wp), &
                                 Iso_t(68,67.9280953_wp), &
                                 Iso_t(69,68.9279645_wp), &
                                 Iso_t(71,70.92495233_wp), &
                                 Iso_t(75,74.92285837_wp), &
                                 Iso_t(77,76.923549843_wp), &
                                 Iso_t(78,77.9228529_wp), &
                                 Iso_t(79,78.92536_wp), &
                                 Iso_t(80,79.9253508_wp), &
                                 Iso_t(81,80.9288329_wp), &
                                 Iso_t(82,81.929774_wp), &
                                 Iso_t(83,82.9345391_wp), &
                                 Iso_t(84,83.9375751_wp), &
                                 Iso_t(85,84.9429697_wp), &
                                 Iso_t(86,85.94658_wp), &
                                 Iso_t(87,86.95268_wp), &
                                 Iso_t(88,87.95691_wp), &
                                 Iso_t(89,88.96379_wp), &
                                 Iso_t(90,89.96863_wp)]

  ElementList(33)%Symbol = adjustl(PTab(33)) ! As
  ElementList(33)%Natural = 1
  call mma_Allocate(ElementList(33)%Isotopes,33)
  ElementList(33)%Isotopes(:) = [ &
                                 Iso_t(75,74.92159457_wp), &
                                 Iso_t(60,59.99388_wp), &
                                 Iso_t(61,60.98112_wp), &
                                 Iso_t(62,61.97361_wp), &
                                 Iso_t(63,62.9639_wp), &
                                 Iso_t(64,63.95743_wp), &
                                 Iso_t(65,64.949611_wp), &
                                 Iso_t(66,65.9441488_wp), &
                                 Iso_t(67,66.93925111_wp), &
                                 Iso_t(68,67.9367741_wp), &
                                 Iso_t(69,68.932246_wp), &
                                 Iso_t(70,69.930926_wp), &
                                 Iso_t(71,70.9271138_wp), &
                                 Iso_t(72,71.9267523_wp), &
                                 Iso_t(73,72.9238291_wp), &
                                 Iso_t(74,73.9239286_wp), &
                                 Iso_t(76,75.92239202_wp), &
                                 Iso_t(77,76.9206476_wp), &
                                 Iso_t(78,77.921828_wp), &
                                 Iso_t(79,78.9209484_wp), &
                                 Iso_t(80,79.9224746_wp), &
                                 Iso_t(81,80.9221323_wp), &
                                 Iso_t(82,81.9247412_wp), &
                                 Iso_t(83,82.9252069_wp), &
                                 Iso_t(84,83.9293033_wp), &
                                 Iso_t(85,84.9321637_wp), &
                                 Iso_t(86,85.9367015_wp), &
                                 Iso_t(87,86.9402917_wp), &
                                 Iso_t(88,87.94555_wp), &
                                 Iso_t(89,88.94976_wp), &
                                 Iso_t(90,89.95563_wp), &
                                 Iso_t(91,90.96039_wp), &
                                 Iso_t(92,91.96674_wp)]

  ElementList(34)%Symbol = adjustl(PTab(34)) ! Se
  ElementList(34)%Natural = 6
  call mma_Allocate(ElementList(34)%Isotopes,32)
  ElementList(34)%Isotopes(:) = [ &
                                 Iso_t(80,79.9165218_wp), &
                                 Iso_t(78,77.91730928_wp), &
                                 Iso_t(76,75.919213704_wp), &
                                 Iso_t(82,81.9166995_wp), &
                                 Iso_t(77,76.919914154_wp), &
                                 Iso_t(74,73.922475934_wp), &
                                 Iso_t(64,63.97109_wp), &
                                 Iso_t(65,64.9644_wp), &
                                 Iso_t(66,65.95559_wp), &
                                 Iso_t(67,66.949994_wp), &
                                 Iso_t(68,67.94182524_wp), &
                                 Iso_t(69,68.9394148_wp), &
                                 Iso_t(70,69.9335155_wp), &
                                 Iso_t(71,70.9322094_wp), &
                                 Iso_t(72,71.9271405_wp), &
                                 Iso_t(73,72.9267549_wp), &
                                 Iso_t(75,74.92252287_wp), &
                                 Iso_t(79,78.91849929_wp), &
                                 Iso_t(81,80.917993_wp), &
                                 Iso_t(83,82.9191186_wp), &
                                 Iso_t(84,83.9184668_wp), &
                                 Iso_t(85,84.9222608_wp), &
                                 Iso_t(86,85.9243117_wp), &
                                 Iso_t(87,86.9286886_wp), &
                                 Iso_t(88,87.9314175_wp), &
                                 Iso_t(89,88.9366691_wp), &
                                 Iso_t(90,89.9401_wp), &
                                 Iso_t(91,90.94596_wp), &
                                 Iso_t(92,91.94984_wp), &
                                 Iso_t(93,92.95629_wp), &
                                 Iso_t(94,93.96049_wp), &
                                 Iso_t(95,94.9673_wp)]

  ElementList(35)%Symbol = adjustl(PTab(35)) ! Br
  ElementList(35)%Natural = 2
  call mma_Allocate(ElementList(35)%Isotopes,32)
  ElementList(35)%Isotopes(:) = [ &
                                 Iso_t(79,78.9183376_wp), &
                                 Iso_t(81,80.9162897_wp), &
                                 Iso_t(67,66.96465_wp), &
                                 Iso_t(68,67.95873_wp), &
                                 Iso_t(69,68.950497_wp), &
                                 Iso_t(70,69.944792_wp), &
                                 Iso_t(71,70.9393422_wp), &
                                 Iso_t(72,71.9365886_wp), &
                                 Iso_t(73,72.9316715_wp), &
                                 Iso_t(74,73.9299102_wp), &
                                 Iso_t(75,74.9258105_wp), &
                                 Iso_t(76,75.924542_wp), &
                                 Iso_t(77,76.9213792_wp), &
                                 Iso_t(78,77.9211459_wp), &
                                 Iso_t(80,79.9185298_wp), &
                                 Iso_t(82,81.9168032_wp), &
                                 Iso_t(83,82.9151756_wp), &
                                 Iso_t(84,83.916496_wp), &
                                 Iso_t(85,84.9156458_wp), &
                                 Iso_t(86,85.9188054_wp), &
                                 Iso_t(87,86.920674_wp), &
                                 Iso_t(88,87.9240833_wp), &
                                 Iso_t(89,88.9267046_wp), &
                                 Iso_t(90,89.9312928_wp), &
                                 Iso_t(91,90.9343986_wp), &
                                 Iso_t(92,91.9396316_wp), &
                                 Iso_t(93,92.94313_wp), &
                                 Iso_t(94,93.9489_wp), &
                                 Iso_t(95,94.95301_wp), &
                                 Iso_t(96,95.95903_wp), &
                                 Iso_t(97,96.96344_wp), &
                                 Iso_t(98,97.96946_wp)]

  ElementList(36)%Symbol = adjustl(PTab(36)) ! Kr
  ElementList(36)%Natural = 6
  call mma_Allocate(ElementList(36)%Isotopes,33)
  ElementList(36)%Isotopes(:) = [ &
                                 Iso_t(84,83.9114977282_wp), &
                                 Iso_t(86,85.9106106269_wp), &
                                 Iso_t(82,81.91348273_wp), &
                                 Iso_t(83,82.91412716_wp), &
                                 Iso_t(80,79.91637808_wp), &
                                 Iso_t(78,77.92036494_wp), &
                                 Iso_t(69,68.96518_wp), &
                                 Iso_t(70,69.95604_wp), &
                                 Iso_t(71,70.95027_wp), &
                                 Iso_t(72,71.9420924_wp), &
                                 Iso_t(73,72.9392892_wp), &
                                 Iso_t(74,73.933084_wp), &
                                 Iso_t(75,74.9309457_wp), &
                                 Iso_t(76,75.9259103_wp), &
                                 Iso_t(77,76.92467_wp), &
                                 Iso_t(79,78.9200829_wp), &
                                 Iso_t(81,80.9165912_wp), &
                                 Iso_t(85,84.9125273_wp), &
                                 Iso_t(87,86.91335476_wp), &
                                 Iso_t(88,87.9144479_wp), &
                                 Iso_t(89,88.9178355_wp), &
                                 Iso_t(90,89.9195279_wp), &
                                 Iso_t(91,90.9238063_wp), &
                                 Iso_t(92,91.9261731_wp), &
                                 Iso_t(93,92.9311472_wp), &
                                 Iso_t(94,93.93414_wp), &
                                 Iso_t(95,94.939711_wp), &
                                 Iso_t(96,95.943017_wp), &
                                 Iso_t(97,96.94909_wp), &
                                 Iso_t(98,97.95243_wp), &
                                 Iso_t(99,98.95839_wp), &
                                 Iso_t(100,99.96237_wp), &
                                 Iso_t(101,100.96873_wp)]

  ElementList(37)%Symbol = adjustl(PTab(37)) ! Rb
  ElementList(37)%Natural = 2
  call mma_Allocate(ElementList(37)%Isotopes,33)
  ElementList(37)%Isotopes(:) = [ &
                                 Iso_t(85,84.9117897379_wp), &
                                 Iso_t(87,86.909180531_wp), &
                                 Iso_t(71,70.96532_wp), &
                                 Iso_t(72,71.95908_wp), &
                                 Iso_t(73,72.95053_wp), &
                                 Iso_t(74,73.9442659_wp), &
                                 Iso_t(75,74.9385732_wp), &
                                 Iso_t(76,75.935073_wp), &
                                 Iso_t(77,76.9304016_wp), &
                                 Iso_t(78,77.9281419_wp), &
                                 Iso_t(79,78.9239899_wp), &
                                 Iso_t(80,79.9225164_wp), &
                                 Iso_t(81,80.9189939_wp), &
                                 Iso_t(82,81.918209_wp), &
                                 Iso_t(83,82.9151142_wp), &
                                 Iso_t(84,83.9143752_wp), &
                                 Iso_t(86,85.91116743_wp), &
                                 Iso_t(88,87.91131559_wp), &
                                 Iso_t(89,88.9122783_wp), &
                                 Iso_t(90,89.9147985_wp), &
                                 Iso_t(91,90.9165372_wp), &
                                 Iso_t(92,91.9197284_wp), &
                                 Iso_t(93,92.9220393_wp), &
                                 Iso_t(94,93.9263948_wp), &
                                 Iso_t(95,94.92926_wp), &
                                 Iso_t(96,95.9341334_wp), &
                                 Iso_t(97,96.9371771_wp), &
                                 Iso_t(98,97.9416869_wp), &
                                 Iso_t(99,98.94503_wp), &
                                 Iso_t(100,99.95003_wp), &
                                 Iso_t(101,100.95404_wp), &
                                 Iso_t(102,101.95952_wp), &
                                 Iso_t(103,102.96392_wp)]

  ElementList(38)%Symbol = adjustl(PTab(38)) ! Sr
  ElementList(38)%Natural = 4
  call mma_Allocate(ElementList(38)%Isotopes,35)
  ElementList(38)%Isotopes(:) = [ &
                                 Iso_t(88,87.9056125_wp), &
                                 Iso_t(86,85.9092606_wp), &
                                 Iso_t(87,86.9088775_wp), &
                                 Iso_t(84,83.9134191_wp), &
                                 Iso_t(73,72.9657_wp), &
                                 Iso_t(74,73.95617_wp), &
                                 Iso_t(75,74.94995_wp), &
                                 Iso_t(76,75.941763_wp), &
                                 Iso_t(77,76.9379455_wp), &
                                 Iso_t(78,77.93218_wp), &
                                 Iso_t(79,78.9297077_wp), &
                                 Iso_t(80,79.9245175_wp), &
                                 Iso_t(81,80.9232114_wp), &
                                 Iso_t(82,81.9183999_wp), &
                                 Iso_t(83,82.9175544_wp), &
                                 Iso_t(85,84.912932_wp), &
                                 Iso_t(89,88.9074511_wp), &
                                 Iso_t(90,89.90773_wp), &
                                 Iso_t(91,90.9101954_wp), &
                                 Iso_t(92,91.9110382_wp), &
                                 Iso_t(93,92.9140242_wp), &
                                 Iso_t(94,93.9153556_wp), &
                                 Iso_t(95,94.9193529_wp), &
                                 Iso_t(96,95.9217066_wp), &
                                 Iso_t(97,96.926374_wp), &
                                 Iso_t(98,97.9286888_wp), &
                                 Iso_t(99,98.9328907_wp), &
                                 Iso_t(100,99.93577_wp), &
                                 Iso_t(101,100.940352_wp), &
                                 Iso_t(102,101.943791_wp), &
                                 Iso_t(103,102.94909_wp), &
                                 Iso_t(104,103.95265_wp), &
                                 Iso_t(105,104.95855_wp), &
                                 Iso_t(106,105.96265_wp), &
                                 Iso_t(107,106.96897_wp)]

  ElementList(39)%Symbol = adjustl(PTab(39)) ! Y
  ElementList(39)%Natural = 1
  call mma_Allocate(ElementList(39)%Isotopes,34)
  ElementList(39)%Isotopes(:) = [ &
                                 Iso_t(89,88.9058403_wp), &
                                 Iso_t(76,75.95856_wp), &
                                 Iso_t(77,76.949781_wp), &
                                 Iso_t(78,77.94361_wp), &
                                 Iso_t(79,78.93735_wp), &
                                 Iso_t(80,79.9343561_wp), &
                                 Iso_t(81,80.9294556_wp), &
                                 Iso_t(82,81.9269314_wp), &
                                 Iso_t(83,82.922485_wp), &
                                 Iso_t(84,83.9206721_wp), &
                                 Iso_t(85,84.916433_wp), &
                                 Iso_t(86,85.914886_wp), &
                                 Iso_t(87,86.9108761_wp), &
                                 Iso_t(88,87.9095016_wp), &
                                 Iso_t(90,89.9071439_wp), &
                                 Iso_t(91,90.9072974_wp), &
                                 Iso_t(92,91.9089451_wp), &
                                 Iso_t(93,92.909578_wp), &
                                 Iso_t(94,93.9115906_wp), &
                                 Iso_t(95,94.9128161_wp), &
                                 Iso_t(96,95.9158968_wp), &
                                 Iso_t(97,96.9182741_wp), &
                                 Iso_t(98,97.9223821_wp), &
                                 Iso_t(99,98.924148_wp), &
                                 Iso_t(100,99.927715_wp), &
                                 Iso_t(101,100.9301477_wp), &
                                 Iso_t(102,101.9343277_wp), &
                                 Iso_t(103,102.937243_wp), &
                                 Iso_t(104,103.94196_wp), &
                                 Iso_t(105,104.94544_wp), &
                                 Iso_t(106,105.95056_wp), &
                                 Iso_t(107,106.95452_wp), &
                                 Iso_t(108,107.95996_wp), &
                                 Iso_t(109,108.96436_wp)]

  ElementList(40)%Symbol = adjustl(PTab(40)) ! Zr
  ElementList(40)%Natural = 5
  call mma_Allocate(ElementList(40)%Isotopes,35)
  ElementList(40)%Isotopes(:) = [ &
                                 Iso_t(90,89.9046977_wp), &
                                 Iso_t(94,93.9063108_wp), &
                                 Iso_t(92,91.9050347_wp), &
                                 Iso_t(91,90.9056396_wp), &
                                 Iso_t(96,95.9082714_wp), &
                                 Iso_t(78,77.95566_wp), &
                                 Iso_t(79,78.94948_wp), &
                                 Iso_t(80,79.9404_wp), &
                                 Iso_t(81,80.93731_wp), &
                                 Iso_t(82,81.93135_wp), &
                                 Iso_t(83,82.9292421_wp), &
                                 Iso_t(84,83.9233269_wp), &
                                 Iso_t(85,84.9214444_wp), &
                                 Iso_t(86,85.9162972_wp), &
                                 Iso_t(87,86.914818_wp), &
                                 Iso_t(88,87.9102213_wp), &
                                 Iso_t(89,88.9088814_wp), &
                                 Iso_t(93,92.9064699_wp), &
                                 Iso_t(95,94.9080385_wp), &
                                 Iso_t(97,96.9109512_wp), &
                                 Iso_t(98,97.9127289_wp), &
                                 Iso_t(99,98.916667_wp), &
                                 Iso_t(100,99.9180006_wp), &
                                 Iso_t(101,100.921448_wp), &
                                 Iso_t(102,101.9231409_wp), &
                                 Iso_t(103,102.927191_wp), &
                                 Iso_t(104,103.929436_wp), &
                                 Iso_t(105,104.934008_wp), &
                                 Iso_t(106,105.93676_wp), &
                                 Iso_t(107,106.94174_wp), &
                                 Iso_t(108,107.94487_wp), &
                                 Iso_t(109,108.95041_wp), &
                                 Iso_t(110,109.95396_wp), &
                                 Iso_t(111,110.95968_wp), &
                                 Iso_t(112,111.9637_wp)]

  ElementList(41)%Symbol = adjustl(PTab(41)) ! Nb
  ElementList(41)%Natural = 1
  call mma_Allocate(ElementList(41)%Isotopes,35)
  ElementList(41)%Isotopes(:) = [ &
                                 Iso_t(93,92.906373_wp), &
                                 Iso_t(81,80.9496_wp), &
                                 Iso_t(82,81.94396_wp), &
                                 Iso_t(83,82.93729_wp), &
                                 Iso_t(84,83.93449_wp), &
                                 Iso_t(85,84.9288458_wp), &
                                 Iso_t(86,85.9257828_wp), &
                                 Iso_t(87,86.9206937_wp), &
                                 Iso_t(88,87.918222_wp), &
                                 Iso_t(89,88.913445_wp), &
                                 Iso_t(90,89.9112584_wp), &
                                 Iso_t(91,90.9069897_wp), &
                                 Iso_t(92,91.9071881_wp), &
                                 Iso_t(94,93.9072788_wp), &
                                 Iso_t(95,94.9068324_wp), &
                                 Iso_t(96,95.9080973_wp), &
                                 Iso_t(97,96.9080959_wp), &
                                 Iso_t(98,97.9103265_wp), &
                                 Iso_t(99,98.911613_wp), &
                                 Iso_t(100,99.9143276_wp), &
                                 Iso_t(101,100.9153103_wp), &
                                 Iso_t(102,101.9180772_wp), &
                                 Iso_t(103,102.9194572_wp), &
                                 Iso_t(104,103.9228925_wp), &
                                 Iso_t(105,104.9249465_wp), &
                                 Iso_t(106,105.9289317_wp), &
                                 Iso_t(107,106.9315937_wp), &
                                 Iso_t(108,107.9360748_wp), &
                                 Iso_t(109,108.93922_wp), &
                                 Iso_t(110,109.94403_wp), &
                                 Iso_t(111,110.94753_wp), &
                                 Iso_t(112,111.95247_wp), &
                                 Iso_t(113,112.95651_wp), &
                                 Iso_t(114,113.96201_wp), &
                                 Iso_t(115,114.96634_wp)]

  ElementList(42)%Symbol = adjustl(PTab(42)) ! Mo
  ElementList(42)%Natural = 7
  call mma_Allocate(ElementList(42)%Isotopes,35)
  ElementList(42)%Isotopes(:) = [ &
                                 Iso_t(98,97.90540482_wp), &
                                 Iso_t(96,95.90467612_wp), &
                                 Iso_t(95,94.90583877_wp), &
                                 Iso_t(92,91.90680796_wp), &
                                 Iso_t(100,99.9074718_wp), &
                                 Iso_t(97,96.90601812_wp), &
                                 Iso_t(94,93.9050849_wp), &
                                 Iso_t(83,82.94988_wp), &
                                 Iso_t(84,83.94149_wp), &
                                 Iso_t(85,84.938261_wp), &
                                 Iso_t(86,85.9311748_wp), &
                                 Iso_t(87,86.9281962_wp), &
                                 Iso_t(88,87.9219678_wp), &
                                 Iso_t(89,88.9194682_wp), &
                                 Iso_t(90,89.9139309_wp), &
                                 Iso_t(91,90.9117453_wp), &
                                 Iso_t(93,92.90680958_wp), &
                                 Iso_t(99,98.90770851_wp), &
                                 Iso_t(101,100.9103414_wp), &
                                 Iso_t(102,101.9102834_wp), &
                                 Iso_t(103,102.913079_wp), &
                                 Iso_t(104,103.9137344_wp), &
                                 Iso_t(105,104.916969_wp), &
                                 Iso_t(106,105.918259_wp), &
                                 Iso_t(107,106.922106_wp), &
                                 Iso_t(108,107.924033_wp), &
                                 Iso_t(109,108.928424_wp), &
                                 Iso_t(110,109.930704_wp), &
                                 Iso_t(111,110.935654_wp), &
                                 Iso_t(112,111.93831_wp), &
                                 Iso_t(113,112.94335_wp), &
                                 Iso_t(114,113.94653_wp), &
                                 Iso_t(115,114.95196_wp), &
                                 Iso_t(116,115.95545_wp), &
                                 Iso_t(117,116.96117_wp)]

  ElementList(43)%Symbol = adjustl(PTab(43)) ! Tc
  ElementList(43)%Natural = 0
  call mma_Allocate(ElementList(43)%Isotopes,36)
  ElementList(43)%Isotopes(:) = [ &
                                 Iso_t(98,97.9072124_wp), &
                                 Iso_t(85,84.95058_wp), &
                                 Iso_t(86,85.94493_wp), &
                                 Iso_t(87,86.9380672_wp), &
                                 Iso_t(88,87.93378_wp), &
                                 Iso_t(89,88.9276487_wp), &
                                 Iso_t(90,89.9240739_wp), &
                                 Iso_t(91,90.9184254_wp), &
                                 Iso_t(92,91.9152698_wp), &
                                 Iso_t(93,92.910246_wp), &
                                 Iso_t(94,93.9096536_wp), &
                                 Iso_t(95,94.9076536_wp), &
                                 Iso_t(96,95.907868_wp), &
                                 Iso_t(97,96.9063667_wp), &
                                 Iso_t(99,98.9062508_wp), &
                                 Iso_t(100,99.9076539_wp), &
                                 Iso_t(101,100.907309_wp), &
                                 Iso_t(102,101.9092097_wp), &
                                 Iso_t(103,102.909176_wp), &
                                 Iso_t(104,103.911425_wp), &
                                 Iso_t(105,104.911655_wp), &
                                 Iso_t(106,105.914358_wp), &
                                 Iso_t(107,106.9154606_wp), &
                                 Iso_t(108,107.9184957_wp), &
                                 Iso_t(109,108.920256_wp), &
                                 Iso_t(110,109.923744_wp), &
                                 Iso_t(111,110.925901_wp), &
                                 Iso_t(112,111.9299458_wp), &
                                 Iso_t(113,112.932569_wp), &
                                 Iso_t(114,113.93691_wp), &
                                 Iso_t(115,114.93998_wp), &
                                 Iso_t(116,115.94476_wp), &
                                 Iso_t(117,116.94806_wp), &
                                 Iso_t(118,117.95299_wp), &
                                 Iso_t(119,118.95666_wp), &
                                 Iso_t(120,119.96187_wp)]

  ElementList(44)%Symbol = adjustl(PTab(44)) ! Ru
  ElementList(44)%Natural = 7
  call mma_Allocate(ElementList(44)%Isotopes,38)
  ElementList(44)%Isotopes(:) = [ &
                                 Iso_t(102,101.9043441_wp), &
                                 Iso_t(104,103.9054275_wp), &
                                 Iso_t(101,100.9055769_wp), &
                                 Iso_t(99,98.9059341_wp), &
                                 Iso_t(100,99.9042143_wp), &
                                 Iso_t(96,95.90759025_wp), &
                                 Iso_t(98,97.9052868_wp), &
                                 Iso_t(87,86.95069_wp), &
                                 Iso_t(88,87.9416_wp), &
                                 Iso_t(89,88.93762_wp), &
                                 Iso_t(90,89.9303444_wp), &
                                 Iso_t(91,90.9267419_wp), &
                                 Iso_t(92,91.9202344_wp), &
                                 Iso_t(93,92.9171044_wp), &
                                 Iso_t(94,93.9113429_wp), &
                                 Iso_t(95,94.910406_wp), &
                                 Iso_t(97,96.9075471_wp), &
                                 Iso_t(103,102.9063186_wp), &
                                 Iso_t(105,104.9077476_wp), &
                                 Iso_t(106,105.9073291_wp), &
                                 Iso_t(107,106.909972_wp), &
                                 Iso_t(108,107.910188_wp), &
                                 Iso_t(109,108.913326_wp), &
                                 Iso_t(110,109.9140407_wp), &
                                 Iso_t(111,110.91757_wp), &
                                 Iso_t(112,111.918809_wp), &
                                 Iso_t(113,112.922844_wp), &
                                 Iso_t(114,113.9246136_wp), &
                                 Iso_t(115,114.92882_wp), &
                                 Iso_t(116,115.9312192_wp), &
                                 Iso_t(117,116.9361_wp), &
                                 Iso_t(118,117.93853_wp), &
                                 Iso_t(119,118.94357_wp), &
                                 Iso_t(120,119.94631_wp), &
                                 Iso_t(121,120.95164_wp), &
                                 Iso_t(122,121.95447_wp), &
                                 Iso_t(123,122.95989_wp), &
                                 Iso_t(124,123.96305_wp)]

  ElementList(45)%Symbol = adjustl(PTab(45)) ! Rh
  ElementList(45)%Natural = 1
  call mma_Allocate(ElementList(45)%Isotopes,38)
  ElementList(45)%Isotopes(:) = [ &
                                 Iso_t(103,102.905498_wp), &
                                 Iso_t(89,88.95058_wp), &
                                 Iso_t(90,89.94422_wp), &
                                 Iso_t(91,90.93688_wp), &
                                 Iso_t(92,91.9323677_wp), &
                                 Iso_t(93,92.9259128_wp), &
                                 Iso_t(94,93.9217305_wp), &
                                 Iso_t(95,94.9158979_wp), &
                                 Iso_t(96,95.914453_wp), &
                                 Iso_t(97,96.911329_wp), &
                                 Iso_t(98,97.910708_wp), &
                                 Iso_t(99,98.9081282_wp), &
                                 Iso_t(100,99.908117_wp), &
                                 Iso_t(101,100.9061606_wp), &
                                 Iso_t(102,101.9068374_wp), &
                                 Iso_t(104,103.9066492_wp), &
                                 Iso_t(105,104.9056885_wp), &
                                 Iso_t(106,105.9072868_wp), &
                                 Iso_t(107,106.906748_wp), &
                                 Iso_t(108,107.908714_wp), &
                                 Iso_t(109,108.9087488_wp), &
                                 Iso_t(110,109.911079_wp), &
                                 Iso_t(111,110.9116423_wp), &
                                 Iso_t(112,111.914403_wp), &
                                 Iso_t(113,112.9154393_wp), &
                                 Iso_t(114,113.918718_wp), &
                                 Iso_t(115,114.9203116_wp), &
                                 Iso_t(116,115.924059_wp), &
                                 Iso_t(117,116.9260354_wp), &
                                 Iso_t(118,117.93034_wp), &
                                 Iso_t(119,118.932557_wp), &
                                 Iso_t(120,119.93686_wp), &
                                 Iso_t(121,120.93942_wp), &
                                 Iso_t(122,121.94399_wp), &
                                 Iso_t(123,122.94685_wp), &
                                 Iso_t(124,123.95151_wp), &
                                 Iso_t(125,124.95469_wp), &
                                 Iso_t(126,125.95946_wp)]

  ElementList(46)%Symbol = adjustl(PTab(46)) ! Pd
  ElementList(46)%Natural = 6
  call mma_Allocate(ElementList(46)%Isotopes,38)
  ElementList(46)%Isotopes(:) = [ &
                                 Iso_t(106,105.9034804_wp), &
                                 Iso_t(108,107.9038916_wp), &
                                 Iso_t(105,104.9050796_wp), &
                                 Iso_t(110,109.9051722_wp), &
                                 Iso_t(104,103.9040305_wp), &
                                 Iso_t(102,101.9056022_wp), &
                                 Iso_t(91,90.95032_wp), &
                                 Iso_t(92,91.94088_wp), &
                                 Iso_t(93,92.93651_wp), &
                                 Iso_t(94,93.9290376_wp), &
                                 Iso_t(95,94.9248898_wp), &
                                 Iso_t(96,95.9182151_wp), &
                                 Iso_t(97,96.916472_wp), &
                                 Iso_t(98,97.9126983_wp), &
                                 Iso_t(99,98.9117748_wp), &
                                 Iso_t(100,99.908505_wp), &
                                 Iso_t(101,100.9082864_wp), &
                                 Iso_t(103,102.9060809_wp), &
                                 Iso_t(107,106.9051282_wp), &
                                 Iso_t(109,108.9059504_wp), &
                                 Iso_t(111,110.90768968_wp), &
                                 Iso_t(112,111.9073297_wp), &
                                 Iso_t(113,112.910261_wp), &
                                 Iso_t(114,113.9103686_wp), &
                                 Iso_t(115,114.913659_wp), &
                                 Iso_t(116,115.914297_wp), &
                                 Iso_t(117,116.9179547_wp), &
                                 Iso_t(118,117.9190667_wp), &
                                 Iso_t(119,118.9233402_wp), &
                                 Iso_t(120,119.9245511_wp), &
                                 Iso_t(121,120.9289503_wp), &
                                 Iso_t(122,121.930632_wp), &
                                 Iso_t(123,122.93514_wp), &
                                 Iso_t(124,123.93714_wp), &
                                 Iso_t(125,124.94179_wp), &
                                 Iso_t(126,125.94416_wp), &
                                 Iso_t(127,126.94907_wp), &
                                 Iso_t(128,127.95183_wp)]

  ElementList(47)%Symbol = adjustl(PTab(47)) ! Ag
  ElementList(47)%Natural = 2
  call mma_Allocate(ElementList(47)%Isotopes,38)
  ElementList(47)%Isotopes(:) = [ &
                                 Iso_t(107,106.9050916_wp), &
                                 Iso_t(109,108.9047553_wp), &
                                 Iso_t(93,92.95033_wp), &
                                 Iso_t(94,93.94373_wp), &
                                 Iso_t(95,94.93602_wp), &
                                 Iso_t(96,95.930744_wp), &
                                 Iso_t(97,96.92397_wp), &
                                 Iso_t(98,97.92156_wp), &
                                 Iso_t(99,98.9176458_wp), &
                                 Iso_t(100,99.9161154_wp), &
                                 Iso_t(101,100.912684_wp), &
                                 Iso_t(102,101.9117047_wp), &
                                 Iso_t(103,102.9089631_wp), &
                                 Iso_t(104,103.9086239_wp), &
                                 Iso_t(105,104.9065256_wp), &
                                 Iso_t(106,105.9066636_wp), &
                                 Iso_t(108,107.9059503_wp), &
                                 Iso_t(110,109.9061102_wp), &
                                 Iso_t(111,110.9052959_wp), &
                                 Iso_t(112,111.9070486_wp), &
                                 Iso_t(113,112.906573_wp), &
                                 Iso_t(114,113.908823_wp), &
                                 Iso_t(115,114.908767_wp), &
                                 Iso_t(116,115.9113868_wp), &
                                 Iso_t(117,116.911774_wp), &
                                 Iso_t(118,117.9145955_wp), &
                                 Iso_t(119,118.91557_wp), &
                                 Iso_t(120,119.9187848_wp), &
                                 Iso_t(121,120.920125_wp), &
                                 Iso_t(122,121.923664_wp), &
                                 Iso_t(123,122.925337_wp), &
                                 Iso_t(124,123.92893_wp), &
                                 Iso_t(125,124.93105_wp), &
                                 Iso_t(126,125.93475_wp), &
                                 Iso_t(127,126.93711_wp), &
                                 Iso_t(128,127.94106_wp), &
                                 Iso_t(129,128.94395_wp), &
                                 Iso_t(130,129.9507_wp)]

  ElementList(48)%Symbol = adjustl(PTab(48)) ! Cd
  ElementList(48)%Natural = 8
  call mma_Allocate(ElementList(48)%Isotopes,39)
  ElementList(48)%Isotopes(:) = [ &
                                 Iso_t(114,113.90336509_wp), &
                                 Iso_t(112,111.90276287_wp), &
                                 Iso_t(111,110.90418287_wp), &
                                 Iso_t(110,109.90300661_wp), &
                                 Iso_t(113,112.90440813_wp), &
                                 Iso_t(116,115.90476315_wp), &
                                 Iso_t(106,105.9064599_wp), &
                                 Iso_t(108,107.9041834_wp), &
                                 Iso_t(95,94.94994_wp), &
                                 Iso_t(96,95.94034_wp), &
                                 Iso_t(97,96.9351_wp), &
                                 Iso_t(98,97.927389_wp), &
                                 Iso_t(99,98.9249258_wp), &
                                 Iso_t(100,99.9203488_wp), &
                                 Iso_t(101,100.9185862_wp), &
                                 Iso_t(102,101.914482_wp), &
                                 Iso_t(103,102.9134165_wp), &
                                 Iso_t(104,103.9098564_wp), &
                                 Iso_t(105,104.9094639_wp), &
                                 Iso_t(107,106.9066121_wp), &
                                 Iso_t(109,108.9049867_wp), &
                                 Iso_t(115,114.90543751_wp), &
                                 Iso_t(117,116.907226_wp), &
                                 Iso_t(118,117.906922_wp), &
                                 Iso_t(119,118.909847_wp), &
                                 Iso_t(120,119.9098681_wp), &
                                 Iso_t(121,120.9129637_wp), &
                                 Iso_t(122,121.9134591_wp), &
                                 Iso_t(123,122.9168925_wp), &
                                 Iso_t(124,123.9176574_wp), &
                                 Iso_t(125,124.9212576_wp), &
                                 Iso_t(126,125.9224291_wp), &
                                 Iso_t(127,126.926472_wp), &
                                 Iso_t(128,127.9278129_wp), &
                                 Iso_t(129,128.93182_wp), &
                                 Iso_t(130,129.93394_wp), &
                                 Iso_t(131,130.9406_wp), &
                                 Iso_t(132,131.94604_wp), &
                                 Iso_t(133,132.95285_wp)]

  ElementList(49)%Symbol = adjustl(PTab(49)) ! In
  ElementList(49)%Natural = 2
  call mma_Allocate(ElementList(49)%Isotopes,39)
  ElementList(49)%Isotopes(:) = [ &
                                 Iso_t(115,114.903878776_wp), &
                                 Iso_t(113,112.90406184_wp), &
                                 Iso_t(97,96.94934_wp), &
                                 Iso_t(98,97.94214_wp), &
                                 Iso_t(99,98.93411_wp), &
                                 Iso_t(100,99.93096_wp), &
                                 Iso_t(101,100.92634_wp), &
                                 Iso_t(102,101.9241071_wp), &
                                 Iso_t(103,102.9198819_wp), &
                                 Iso_t(104,103.9182145_wp), &
                                 Iso_t(105,104.914502_wp), &
                                 Iso_t(106,105.913464_wp), &
                                 Iso_t(107,106.91029_wp), &
                                 Iso_t(108,107.9096935_wp), &
                                 Iso_t(109,108.9071514_wp), &
                                 Iso_t(110,109.90717_wp), &
                                 Iso_t(111,110.9051085_wp), &
                                 Iso_t(112,111.9055377_wp), &
                                 Iso_t(114,113.90491791_wp), &
                                 Iso_t(116,115.90525999_wp), &
                                 Iso_t(117,116.9045157_wp), &
                                 Iso_t(118,117.9063566_wp), &
                                 Iso_t(119,118.9058507_wp), &
                                 Iso_t(120,119.907967_wp), &
                                 Iso_t(121,120.907851_wp), &
                                 Iso_t(122,121.910281_wp), &
                                 Iso_t(123,122.910434_wp), &
                                 Iso_t(124,123.913182_wp), &
                                 Iso_t(125,124.913605_wp), &
                                 Iso_t(126,125.916507_wp), &
                                 Iso_t(127,126.917446_wp), &
                                 Iso_t(128,127.9204_wp), &
                                 Iso_t(129,128.9218053_wp), &
                                 Iso_t(130,129.924977_wp), &
                                 Iso_t(131,130.9269715_wp), &
                                 Iso_t(132,131.933001_wp), &
                                 Iso_t(133,132.93831_wp), &
                                 Iso_t(134,133.94454_wp), &
                                 Iso_t(135,134.95005_wp)]

  ElementList(50)%Symbol = adjustl(PTab(50)) ! Sn
  ElementList(50)%Natural = 10
  call mma_Allocate(ElementList(50)%Isotopes,40)
  ElementList(50)%Isotopes(:) = [ &
                                 Iso_t(120,119.90220163_wp), &
                                 Iso_t(118,117.90160657_wp), &
                                 Iso_t(116,115.9017428_wp), &
                                 Iso_t(119,118.90331117_wp), &
                                 Iso_t(117,116.90295398_wp), &
                                 Iso_t(124,123.9052766_wp), &
                                 Iso_t(122,121.9034438_wp), &
                                 Iso_t(112,111.90482387_wp), &
                                 Iso_t(114,113.9027827_wp), &
                                 Iso_t(115,114.903344699_wp), &
                                 Iso_t(99,98.94853_wp), &
                                 Iso_t(100,99.9385_wp), &
                                 Iso_t(101,100.93526_wp), &
                                 Iso_t(102,101.93029_wp), &
                                 Iso_t(103,102.928105_wp), &
                                 Iso_t(104,103.9231052_wp), &
                                 Iso_t(105,104.9212684_wp), &
                                 Iso_t(106,105.9169574_wp), &
                                 Iso_t(107,106.9157137_wp), &
                                 Iso_t(108,107.9118943_wp), &
                                 Iso_t(109,108.9112921_wp), &
                                 Iso_t(110,109.907845_wp), &
                                 Iso_t(111,110.9077401_wp), &
                                 Iso_t(113,112.9051757_wp), &
                                 Iso_t(121,120.9042426_wp), &
                                 Iso_t(123,122.9057252_wp), &
                                 Iso_t(125,124.9077864_wp), &
                                 Iso_t(126,125.907659_wp), &
                                 Iso_t(127,126.91039_wp), &
                                 Iso_t(128,127.910507_wp), &
                                 Iso_t(129,128.913465_wp), &
                                 Iso_t(130,129.9139738_wp), &
                                 Iso_t(131,130.917045_wp), &
                                 Iso_t(132,131.9178267_wp), &
                                 Iso_t(133,132.9239134_wp), &
                                 Iso_t(134,133.9286821_wp), &
                                 Iso_t(135,134.9349086_wp), &
                                 Iso_t(136,135.93999_wp), &
                                 Iso_t(137,136.94655_wp), &
                                 Iso_t(138,137.95184_wp)]

  ElementList(51)%Symbol = adjustl(PTab(51)) ! Sb
  ElementList(51)%Natural = 2
  call mma_Allocate(ElementList(51)%Isotopes,38)
  ElementList(51)%Isotopes(:) = [ &
                                 Iso_t(121,120.903812_wp), &
                                 Iso_t(123,122.9042132_wp), &
                                 Iso_t(103,102.93969_wp), &
                                 Iso_t(104,103.93648_wp), &
                                 Iso_t(105,104.931276_wp), &
                                 Iso_t(106,105.928638_wp), &
                                 Iso_t(107,106.9241506_wp), &
                                 Iso_t(108,107.9222267_wp), &
                                 Iso_t(109,108.9181411_wp), &
                                 Iso_t(110,109.9168543_wp), &
                                 Iso_t(111,110.9132182_wp), &
                                 Iso_t(112,111.9124_wp), &
                                 Iso_t(113,112.909375_wp), &
                                 Iso_t(114,113.90929_wp), &
                                 Iso_t(115,114.906598_wp), &
                                 Iso_t(116,115.9067931_wp), &
                                 Iso_t(117,116.9048415_wp), &
                                 Iso_t(118,117.9055321_wp), &
                                 Iso_t(119,118.9039455_wp), &
                                 Iso_t(120,119.9050794_wp), &
                                 Iso_t(122,121.9051699_wp), &
                                 Iso_t(124,123.905935_wp), &
                                 Iso_t(125,124.905253_wp), &
                                 Iso_t(126,125.907253_wp), &
                                 Iso_t(127,126.9069243_wp), &
                                 Iso_t(128,127.909146_wp), &
                                 Iso_t(129,128.909147_wp), &
                                 Iso_t(130,129.911662_wp), &
                                 Iso_t(131,130.9119888_wp), &
                                 Iso_t(132,131.9145077_wp), &
                                 Iso_t(133,132.9152732_wp), &
                                 Iso_t(134,133.9205357_wp), &
                                 Iso_t(135,134.9251851_wp), &
                                 Iso_t(136,135.9307459_wp), &
                                 Iso_t(137,136.93555_wp), &
                                 Iso_t(138,137.94145_wp), &
                                 Iso_t(139,138.94655_wp), &
                                 Iso_t(140,139.95283_wp)]

  ElementList(52)%Symbol = adjustl(PTab(52)) ! Te
  ElementList(52)%Natural = 8
  call mma_Allocate(ElementList(52)%Isotopes,39)
  ElementList(52)%Isotopes(:) = [ &
                                 Iso_t(130,129.906222748_wp), &
                                 Iso_t(128,127.90446128_wp), &
                                 Iso_t(126,125.9033109_wp), &
                                 Iso_t(125,124.9044299_wp), &
                                 Iso_t(124,123.9028171_wp), &
                                 Iso_t(122,121.9030435_wp), &
                                 Iso_t(123,122.9042698_wp), &
                                 Iso_t(120,119.9040593_wp), &
                                 Iso_t(105,104.9433_wp), &
                                 Iso_t(106,105.9375_wp), &
                                 Iso_t(107,106.935012_wp), &
                                 Iso_t(108,107.9293805_wp), &
                                 Iso_t(109,108.9273045_wp), &
                                 Iso_t(110,109.9224581_wp), &
                                 Iso_t(111,110.9210006_wp), &
                                 Iso_t(112,111.9167279_wp), &
                                 Iso_t(113,112.915891_wp), &
                                 Iso_t(114,113.912089_wp), &
                                 Iso_t(115,114.911902_wp), &
                                 Iso_t(116,115.90846_wp), &
                                 Iso_t(117,116.908646_wp), &
                                 Iso_t(118,117.905854_wp), &
                                 Iso_t(119,118.9064071_wp), &
                                 Iso_t(121,120.904944_wp), &
                                 Iso_t(127,126.9052257_wp), &
                                 Iso_t(129,128.90659646_wp), &
                                 Iso_t(131,130.908522213_wp), &
                                 Iso_t(132,131.9085467_wp), &
                                 Iso_t(133,132.9109688_wp), &
                                 Iso_t(134,133.911394_wp), &
                                 Iso_t(135,134.9165557_wp), &
                                 Iso_t(136,135.9201006_wp), &
                                 Iso_t(137,136.9255989_wp), &
                                 Iso_t(138,137.9294722_wp), &
                                 Iso_t(139,138.9353672_wp), &
                                 Iso_t(140,139.939499_wp), &
                                 Iso_t(141,140.9458_wp), &
                                 Iso_t(142,141.95022_wp), &
                                 Iso_t(143,142.95676_wp)]

  ElementList(53)%Symbol = adjustl(PTab(53)) ! I
  ElementList(53)%Natural = 1
  call mma_Allocate(ElementList(53)%Isotopes,39)
  ElementList(53)%Isotopes(:) = [ &
                                 Iso_t(127,126.9044719_wp), &
                                 Iso_t(107,106.94678_wp), &
                                 Iso_t(108,107.94348_wp), &
                                 Iso_t(109,108.9380853_wp), &
                                 Iso_t(110,109.935089_wp), &
                                 Iso_t(111,110.9302692_wp), &
                                 Iso_t(112,111.928005_wp), &
                                 Iso_t(113,112.9236501_wp), &
                                 Iso_t(114,113.92185_wp), &
                                 Iso_t(115,114.918048_wp), &
                                 Iso_t(116,115.91681_wp), &
                                 Iso_t(117,116.913648_wp), &
                                 Iso_t(118,117.913074_wp), &
                                 Iso_t(119,118.910074_wp), &
                                 Iso_t(120,119.910087_wp), &
                                 Iso_t(121,120.9074051_wp), &
                                 Iso_t(122,121.9075888_wp), &
                                 Iso_t(123,122.9055885_wp), &
                                 Iso_t(124,123.906209_wp), &
                                 Iso_t(125,124.9046294_wp), &
                                 Iso_t(126,125.9056233_wp), &
                                 Iso_t(128,127.9058086_wp), &
                                 Iso_t(129,128.9049837_wp), &
                                 Iso_t(130,129.9066702_wp), &
                                 Iso_t(131,130.9061263_wp), &
                                 Iso_t(132,131.9079935_wp), &
                                 Iso_t(133,132.907797_wp), &
                                 Iso_t(134,133.9097588_wp), &
                                 Iso_t(135,134.9100488_wp), &
                                 Iso_t(136,135.914604_wp), &
                                 Iso_t(137,136.9180282_wp), &
                                 Iso_t(138,137.9227264_wp), &
                                 Iso_t(139,138.926506_wp), &
                                 Iso_t(140,139.93173_wp), &
                                 Iso_t(141,140.93569_wp), &
                                 Iso_t(142,141.9412_wp), &
                                 Iso_t(143,142.94565_wp), &
                                 Iso_t(144,143.95139_wp), &
                                 Iso_t(145,144.95605_wp)]

  ElementList(54)%Symbol = adjustl(PTab(54)) ! Xe
  ElementList(54)%Natural = 9
  call mma_Allocate(ElementList(54)%Isotopes,40)
  ElementList(54)%Isotopes(:) = [ &
                                 Iso_t(132,131.9041550856_wp), &
                                 Iso_t(129,128.9047808611_wp), &
                                 Iso_t(131,130.90508406_wp), &
                                 Iso_t(134,133.90539466_wp), &
                                 Iso_t(136,135.907214484_wp), &
                                 Iso_t(130,129.903509349_wp), &
                                 Iso_t(128,127.903531_wp), &
                                 Iso_t(124,123.905892_wp), &
                                 Iso_t(126,125.9042983_wp), &
                                 Iso_t(109,108.95043_wp), &
                                 Iso_t(110,109.94426_wp), &
                                 Iso_t(111,110.941607_wp), &
                                 Iso_t(112,111.935559_wp), &
                                 Iso_t(113,112.9332217_wp), &
                                 Iso_t(114,113.92798_wp), &
                                 Iso_t(115,114.926294_wp), &
                                 Iso_t(116,115.921581_wp), &
                                 Iso_t(117,116.920359_wp), &
                                 Iso_t(118,117.916179_wp), &
                                 Iso_t(119,118.915411_wp), &
                                 Iso_t(120,119.911784_wp), &
                                 Iso_t(121,120.911453_wp), &
                                 Iso_t(122,121.908368_wp), &
                                 Iso_t(123,122.908482_wp), &
                                 Iso_t(125,124.9063944_wp), &
                                 Iso_t(127,126.9051829_wp), &
                                 Iso_t(133,132.9059108_wp), &
                                 Iso_t(135,134.9072278_wp), &
                                 Iso_t(137,136.91155778_wp), &
                                 Iso_t(138,137.9141463_wp), &
                                 Iso_t(139,138.9187922_wp), &
                                 Iso_t(140,139.9216458_wp), &
                                 Iso_t(141,140.9267872_wp), &
                                 Iso_t(142,141.9299731_wp), &
                                 Iso_t(143,142.9353696_wp), &
                                 Iso_t(144,143.9389451_wp), &
                                 Iso_t(145,144.94472_wp), &
                                 Iso_t(146,145.948518_wp), &
                                 Iso_t(147,146.95426_wp), &
                                 Iso_t(148,147.95813_wp)]

  ElementList(55)%Symbol = adjustl(PTab(55)) ! Cs
  ElementList(55)%Natural = 1
  call mma_Allocate(ElementList(55)%Isotopes,40)
  ElementList(55)%Isotopes(:) = [ &
                                 Iso_t(133,132.905451961_wp), &
                                 Iso_t(112,111.950309_wp), &
                                 Iso_t(113,112.9444291_wp), &
                                 Iso_t(114,113.941296_wp), &
                                 Iso_t(115,114.93591_wp), &
                                 Iso_t(116,115.93337_wp), &
                                 Iso_t(117,116.928617_wp), &
                                 Iso_t(118,117.92656_wp), &
                                 Iso_t(119,118.922377_wp), &
                                 Iso_t(120,119.920677_wp), &
                                 Iso_t(121,120.917227_wp), &
                                 Iso_t(122,121.916108_wp), &
                                 Iso_t(123,122.912996_wp), &
                                 Iso_t(124,123.9122578_wp), &
                                 Iso_t(125,124.909728_wp), &
                                 Iso_t(126,125.909446_wp), &
                                 Iso_t(127,126.9074174_wp), &
                                 Iso_t(128,127.9077487_wp), &
                                 Iso_t(129,128.9060657_wp), &
                                 Iso_t(130,129.9067093_wp), &
                                 Iso_t(131,130.9054649_wp), &
                                 Iso_t(132,131.9064339_wp), &
                                 Iso_t(134,133.906718503_wp), &
                                 Iso_t(135,134.905977_wp), &
                                 Iso_t(136,135.9073114_wp), &
                                 Iso_t(137,136.90708923_wp), &
                                 Iso_t(138,137.9110171_wp), &
                                 Iso_t(139,138.9133638_wp), &
                                 Iso_t(140,139.9172831_wp), &
                                 Iso_t(141,140.9200455_wp), &
                                 Iso_t(142,141.924296_wp), &
                                 Iso_t(143,142.927349_wp), &
                                 Iso_t(144,143.932076_wp), &
                                 Iso_t(145,144.935527_wp), &
                                 Iso_t(146,145.940344_wp), &
                                 Iso_t(147,146.944156_wp), &
                                 Iso_t(148,147.94923_wp), &
                                 Iso_t(149,148.95302_wp), &
                                 Iso_t(150,149.95833_wp), &
                                 Iso_t(151,150.96258_wp)]

  ElementList(56)%Symbol = adjustl(PTab(56)) ! Ba
  ElementList(56)%Natural = 7
  call mma_Allocate(ElementList(56)%Isotopes,40)
  ElementList(56)%Isotopes(:) = [ &
                                 Iso_t(138,137.905247_wp), &
                                 Iso_t(137,136.90582714_wp), &
                                 Iso_t(136,135.90457573_wp), &
                                 Iso_t(135,134.90568838_wp), &
                                 Iso_t(134,133.90450818_wp), &
                                 Iso_t(130,129.9063207_wp), &
                                 Iso_t(132,131.9050611_wp), &
                                 Iso_t(114,113.95066_wp), &
                                 Iso_t(115,114.94737_wp), &
                                 Iso_t(116,115.94128_wp), &
                                 Iso_t(117,116.93814_wp), &
                                 Iso_t(118,117.93306_wp), &
                                 Iso_t(119,118.93066_wp), &
                                 Iso_t(120,119.92605_wp), &
                                 Iso_t(121,120.92405_wp), &
                                 Iso_t(122,121.919904_wp), &
                                 Iso_t(123,122.918781_wp), &
                                 Iso_t(124,123.915094_wp), &
                                 Iso_t(125,124.914472_wp), &
                                 Iso_t(126,125.91125_wp), &
                                 Iso_t(127,126.911091_wp), &
                                 Iso_t(128,127.908342_wp), &
                                 Iso_t(129,128.908681_wp), &
                                 Iso_t(131,130.906941_wp), &
                                 Iso_t(133,132.9060074_wp), &
                                 Iso_t(139,138.9088411_wp), &
                                 Iso_t(140,139.9106057_wp), &
                                 Iso_t(141,140.9144033_wp), &
                                 Iso_t(142,141.9164324_wp), &
                                 Iso_t(143,142.9206253_wp), &
                                 Iso_t(144,143.9229549_wp), &
                                 Iso_t(145,144.9275184_wp), &
                                 Iso_t(146,145.930284_wp), &
                                 Iso_t(147,146.935304_wp), &
                                 Iso_t(148,147.938171_wp), &
                                 Iso_t(149,148.94308_wp), &
                                 Iso_t(150,149.94605_wp), &
                                 Iso_t(151,150.95127_wp), &
                                 Iso_t(152,151.95481_wp), &
                                 Iso_t(153,152.96036_wp)]

  ElementList(57)%Symbol = adjustl(PTab(57)) ! La
  ElementList(57)%Natural = 2
  call mma_Allocate(ElementList(57)%Isotopes,40)
  ElementList(57)%Isotopes(:) = [ &
                                 Iso_t(139,138.9063563_wp), &
                                 Iso_t(138,137.9071149_wp), &
                                 Iso_t(116,115.9563_wp), &
                                 Iso_t(117,116.94999_wp), &
                                 Iso_t(118,117.94673_wp), &
                                 Iso_t(119,118.94099_wp), &
                                 Iso_t(120,119.93807_wp), &
                                 Iso_t(121,120.93315_wp), &
                                 Iso_t(122,121.93071_wp), &
                                 Iso_t(123,122.9263_wp), &
                                 Iso_t(124,123.924574_wp), &
                                 Iso_t(125,124.920816_wp), &
                                 Iso_t(126,125.919513_wp), &
                                 Iso_t(127,126.916375_wp), &
                                 Iso_t(128,127.915592_wp), &
                                 Iso_t(129,128.912694_wp), &
                                 Iso_t(130,129.912369_wp), &
                                 Iso_t(131,130.91007_wp), &
                                 Iso_t(132,131.910119_wp), &
                                 Iso_t(133,132.908218_wp), &
                                 Iso_t(134,133.908514_wp), &
                                 Iso_t(135,134.906984_wp), &
                                 Iso_t(136,135.907635_wp), &
                                 Iso_t(137,136.9064504_wp), &
                                 Iso_t(140,139.9094806_wp), &
                                 Iso_t(141,140.910966_wp), &
                                 Iso_t(142,141.9140909_wp), &
                                 Iso_t(143,142.9160795_wp), &
                                 Iso_t(144,143.919646_wp), &
                                 Iso_t(145,144.921808_wp), &
                                 Iso_t(146,145.925875_wp), &
                                 Iso_t(147,146.928418_wp), &
                                 Iso_t(148,147.932679_wp), &
                                 Iso_t(149,148.93535_wp), &
                                 Iso_t(150,149.93947_wp), &
                                 Iso_t(151,150.94232_wp), &
                                 Iso_t(152,151.94682_wp), &
                                 Iso_t(153,152.95036_wp), &
                                 Iso_t(154,153.95517_wp), &
                                 Iso_t(155,154.95901_wp)]

  ElementList(58)%Symbol = adjustl(PTab(58)) ! Ce
  ElementList(58)%Natural = 4
  call mma_Allocate(ElementList(58)%Isotopes,39)
  ElementList(58)%Isotopes(:) = [ &
                                 Iso_t(140,139.9054431_wp), &
                                 Iso_t(142,141.9092504_wp), &
                                 Iso_t(138,137.905991_wp), &
                                 Iso_t(136,135.90712921_wp), &
                                 Iso_t(119,118.95271_wp), &
                                 Iso_t(120,119.94654_wp), &
                                 Iso_t(121,120.94335_wp), &
                                 Iso_t(122,121.93787_wp), &
                                 Iso_t(123,122.93528_wp), &
                                 Iso_t(124,123.93031_wp), &
                                 Iso_t(125,124.92844_wp), &
                                 Iso_t(126,125.923971_wp), &
                                 Iso_t(127,126.922727_wp), &
                                 Iso_t(128,127.918911_wp), &
                                 Iso_t(129,128.918102_wp), &
                                 Iso_t(130,129.914736_wp), &
                                 Iso_t(131,130.914429_wp), &
                                 Iso_t(132,131.911464_wp), &
                                 Iso_t(133,132.91152_wp), &
                                 Iso_t(134,133.908928_wp), &
                                 Iso_t(135,134.909161_wp), &
                                 Iso_t(137,136.90776236_wp), &
                                 Iso_t(139,138.9066551_wp), &
                                 Iso_t(141,140.9082807_wp), &
                                 Iso_t(143,142.9123921_wp), &
                                 Iso_t(144,143.9136529_wp), &
                                 Iso_t(145,144.917265_wp), &
                                 Iso_t(146,145.918802_wp), &
                                 Iso_t(147,146.9226899_wp), &
                                 Iso_t(148,147.924424_wp), &
                                 Iso_t(149,148.928427_wp), &
                                 Iso_t(150,149.930384_wp), &
                                 Iso_t(151,150.934272_wp), &
                                 Iso_t(152,151.9366_wp), &
                                 Iso_t(153,152.94093_wp), &
                                 Iso_t(154,153.9438_wp), &
                                 Iso_t(155,154.94855_wp), &
                                 Iso_t(156,155.95183_wp), &
                                 Iso_t(157,156.95705_wp)]

  ElementList(59)%Symbol = adjustl(PTab(59)) ! Pr
  ElementList(59)%Natural = 1
  call mma_Allocate(ElementList(59)%Isotopes,39)
  ElementList(59)%Isotopes(:) = [ &
                                 Iso_t(141,140.9076576_wp), &
                                 Iso_t(121,120.95532_wp), &
                                 Iso_t(122,121.95175_wp), &
                                 Iso_t(123,122.94596_wp), &
                                 Iso_t(124,123.94294_wp), &
                                 Iso_t(125,124.9377_wp), &
                                 Iso_t(126,125.93524_wp), &
                                 Iso_t(127,126.93071_wp), &
                                 Iso_t(128,127.928791_wp), &
                                 Iso_t(129,128.925095_wp), &
                                 Iso_t(130,129.92359_wp), &
                                 Iso_t(131,130.920235_wp), &
                                 Iso_t(132,131.919255_wp), &
                                 Iso_t(133,132.916331_wp), &
                                 Iso_t(134,133.915697_wp), &
                                 Iso_t(135,134.913112_wp), &
                                 Iso_t(136,135.912677_wp), &
                                 Iso_t(137,136.9106792_wp), &
                                 Iso_t(138,137.910754_wp), &
                                 Iso_t(139,138.9089408_wp), &
                                 Iso_t(140,139.9090803_wp), &
                                 Iso_t(142,141.9100496_wp), &
                                 Iso_t(143,142.9108228_wp), &
                                 Iso_t(144,143.9133109_wp), &
                                 Iso_t(145,144.9145182_wp), &
                                 Iso_t(146,145.91768_wp), &
                                 Iso_t(147,146.919008_wp), &
                                 Iso_t(148,147.92213_wp), &
                                 Iso_t(149,148.923736_wp), &
                                 Iso_t(150,149.9266765_wp), &
                                 Iso_t(151,150.928309_wp), &
                                 Iso_t(152,151.931553_wp), &
                                 Iso_t(153,152.933904_wp), &
                                 Iso_t(154,153.93753_wp), &
                                 Iso_t(155,154.940509_wp), &
                                 Iso_t(156,155.94464_wp), &
                                 Iso_t(157,156.94789_wp), &
                                 Iso_t(158,157.95241_wp), &
                                 Iso_t(159,158.95589_wp)]

  ElementList(60)%Symbol = adjustl(PTab(60)) ! Nd
  ElementList(60)%Natural = 7
  call mma_Allocate(ElementList(60)%Isotopes,38)
  ElementList(60)%Isotopes(:) = [ &
                                 Iso_t(142,141.907729_wp), &
                                 Iso_t(144,143.910093_wp), &
                                 Iso_t(146,145.9131226_wp), &
                                 Iso_t(143,142.90982_wp), &
                                 Iso_t(145,144.9125793_wp), &
                                 Iso_t(148,147.9168993_wp), &
                                 Iso_t(150,149.9209022_wp), &
                                 Iso_t(124,123.9522_wp), &
                                 Iso_t(125,124.9489_wp), &
                                 Iso_t(126,125.94311_wp), &
                                 Iso_t(127,126.94038_wp), &
                                 Iso_t(128,127.93525_wp), &
                                 Iso_t(129,128.9331_wp), &
                                 Iso_t(130,129.928506_wp), &
                                 Iso_t(131,130.927248_wp), &
                                 Iso_t(132,131.923321_wp), &
                                 Iso_t(133,132.922348_wp), &
                                 Iso_t(134,133.91879_wp), &
                                 Iso_t(135,134.918181_wp), &
                                 Iso_t(136,135.914976_wp), &
                                 Iso_t(137,136.914562_wp), &
                                 Iso_t(138,137.91195_wp), &
                                 Iso_t(139,138.911954_wp), &
                                 Iso_t(140,139.90955_wp), &
                                 Iso_t(141,140.9096147_wp), &
                                 Iso_t(147,146.9161061_wp), &
                                 Iso_t(149,148.9201548_wp), &
                                 Iso_t(151,150.9238403_wp), &
                                 Iso_t(152,151.924692_wp), &
                                 Iso_t(153,152.927718_wp), &
                                 Iso_t(154,153.92948_wp), &
                                 Iso_t(155,154.9331357_wp), &
                                 Iso_t(156,155.93508_wp), &
                                 Iso_t(157,156.939386_wp), &
                                 Iso_t(158,157.94197_wp), &
                                 Iso_t(159,158.94653_wp), &
                                 Iso_t(160,159.9494_wp), &
                                 Iso_t(161,160.95428_wp)]

  ElementList(61)%Symbol = adjustl(PTab(61)) ! Pm
  ElementList(61)%Natural = 0
  call mma_Allocate(ElementList(61)%Isotopes,38)
  ElementList(61)%Isotopes(:) = [ &
                                 Iso_t(145,144.9127559_wp), &
                                 Iso_t(126,125.95792_wp), &
                                 Iso_t(127,126.95192_wp), &
                                 Iso_t(128,127.9487_wp), &
                                 Iso_t(129,128.94323_wp), &
                                 Iso_t(130,129.94053_wp), &
                                 Iso_t(131,130.93567_wp), &
                                 Iso_t(132,131.93384_wp), &
                                 Iso_t(133,132.929782_wp), &
                                 Iso_t(134,133.928353_wp), &
                                 Iso_t(135,134.924823_wp), &
                                 Iso_t(136,135.923585_wp), &
                                 Iso_t(137,136.92048_wp), &
                                 Iso_t(138,137.919548_wp), &
                                 Iso_t(139,138.9168_wp), &
                                 Iso_t(140,139.91604_wp), &
                                 Iso_t(141,140.913555_wp), &
                                 Iso_t(142,141.91289_wp), &
                                 Iso_t(143,142.9109383_wp), &
                                 Iso_t(144,143.9125964_wp), &
                                 Iso_t(146,145.9147024_wp), &
                                 Iso_t(147,146.915145_wp), &
                                 Iso_t(148,147.9174819_wp), &
                                 Iso_t(149,148.9183423_wp), &
                                 Iso_t(150,149.920991_wp), &
                                 Iso_t(151,150.9212175_wp), &
                                 Iso_t(152,151.923506_wp), &
                                 Iso_t(153,152.9241567_wp), &
                                 Iso_t(154,153.926472_wp), &
                                 Iso_t(155,154.928137_wp), &
                                 Iso_t(156,155.9311175_wp), &
                                 Iso_t(157,156.9331214_wp), &
                                 Iso_t(158,157.936565_wp), &
                                 Iso_t(159,158.939287_wp), &
                                 Iso_t(160,159.9431_wp), &
                                 Iso_t(161,160.94607_wp), &
                                 Iso_t(162,161.95022_wp), &
                                 Iso_t(163,162.95357_wp)]

  ElementList(62)%Symbol = adjustl(PTab(62)) ! Sm
  ElementList(62)%Natural = 7
  call mma_Allocate(ElementList(62)%Isotopes,38)
  ElementList(62)%Isotopes(:) = [ &
                                 Iso_t(152,151.9197397_wp), &
                                 Iso_t(154,153.9222169_wp), &
                                 Iso_t(147,146.9149044_wp), &
                                 Iso_t(149,148.9171921_wp), &
                                 Iso_t(148,147.9148292_wp), &
                                 Iso_t(150,149.9172829_wp), &
                                 Iso_t(144,143.9120065_wp), &
                                 Iso_t(128,127.95842_wp), &
                                 Iso_t(129,128.95476_wp), &
                                 Iso_t(130,129.949_wp), &
                                 Iso_t(131,130.94618_wp), &
                                 Iso_t(132,131.94087_wp), &
                                 Iso_t(133,132.93856_wp), &
                                 Iso_t(134,133.93411_wp), &
                                 Iso_t(135,134.93252_wp), &
                                 Iso_t(136,135.928276_wp), &
                                 Iso_t(137,136.926971_wp), &
                                 Iso_t(138,137.923244_wp), &
                                 Iso_t(139,138.922297_wp), &
                                 Iso_t(140,139.918995_wp), &
                                 Iso_t(141,140.9184816_wp), &
                                 Iso_t(142,141.9152044_wp), &
                                 Iso_t(143,142.9146353_wp), &
                                 Iso_t(145,144.9134173_wp), &
                                 Iso_t(146,145.913047_wp), &
                                 Iso_t(151,150.9199398_wp), &
                                 Iso_t(153,152.9221047_wp), &
                                 Iso_t(155,154.9246477_wp), &
                                 Iso_t(156,155.925536_wp), &
                                 Iso_t(157,156.9284187_wp), &
                                 Iso_t(158,157.929951_wp), &
                                 Iso_t(159,158.9332172_wp), &
                                 Iso_t(160,159.9353353_wp), &
                                 Iso_t(161,160.9391602_wp), &
                                 Iso_t(162,161.94146_wp), &
                                 Iso_t(163,162.94555_wp), &
                                 Iso_t(164,163.94836_wp), &
                                 Iso_t(165,164.95297_wp)]

  ElementList(63)%Symbol = adjustl(PTab(63)) ! Eu
  ElementList(63)%Natural = 2
  call mma_Allocate(ElementList(63)%Isotopes,38)
  ElementList(63)%Isotopes(:) = [ &
                                 Iso_t(153,152.921238_wp), &
                                 Iso_t(151,150.9198578_wp), &
                                 Iso_t(130,129.96369_wp), &
                                 Iso_t(131,130.95784_wp), &
                                 Iso_t(132,131.95467_wp), &
                                 Iso_t(133,132.94929_wp), &
                                 Iso_t(134,133.9464_wp), &
                                 Iso_t(135,134.94187_wp), &
                                 Iso_t(136,135.93962_wp), &
                                 Iso_t(137,136.93546_wp), &
                                 Iso_t(138,137.933709_wp), &
                                 Iso_t(139,138.929792_wp), &
                                 Iso_t(140,139.928088_wp), &
                                 Iso_t(141,140.924932_wp), &
                                 Iso_t(142,141.923442_wp), &
                                 Iso_t(143,142.920299_wp), &
                                 Iso_t(144,143.91882_wp), &
                                 Iso_t(145,144.9162726_wp), &
                                 Iso_t(146,145.917211_wp), &
                                 Iso_t(147,146.9167527_wp), &
                                 Iso_t(148,147.918089_wp), &
                                 Iso_t(149,148.9179378_wp), &
                                 Iso_t(150,149.9197077_wp), &
                                 Iso_t(152,151.9217522_wp), &
                                 Iso_t(154,153.922987_wp), &
                                 Iso_t(155,154.9229011_wp), &
                                 Iso_t(156,155.9247605_wp), &
                                 Iso_t(157,156.9254334_wp), &
                                 Iso_t(158,157.927799_wp), &
                                 Iso_t(159,158.9291001_wp), &
                                 Iso_t(160,159.931851_wp), &
                                 Iso_t(161,160.933664_wp), &
                                 Iso_t(162,161.936989_wp), &
                                 Iso_t(163,162.939196_wp), &
                                 Iso_t(164,163.94274_wp), &
                                 Iso_t(165,164.94559_wp), &
                                 Iso_t(166,165.94962_wp), &
                                 Iso_t(167,166.95289_wp)]

  ElementList(64)%Symbol = adjustl(PTab(64)) ! Gd
  ElementList(64)%Natural = 7
  call mma_Allocate(ElementList(64)%Isotopes,37)
  ElementList(64)%Isotopes(:) = [ &
                                 Iso_t(158,157.9241123_wp), &
                                 Iso_t(160,159.9270624_wp), &
                                 Iso_t(156,155.9221312_wp), &
                                 Iso_t(157,156.9239686_wp), &
                                 Iso_t(155,154.9226305_wp), &
                                 Iso_t(154,153.9208741_wp), &
                                 Iso_t(152,151.9197995_wp), &
                                 Iso_t(133,132.96133_wp), &
                                 Iso_t(134,133.95566_wp), &
                                 Iso_t(135,134.95245_wp), &
                                 Iso_t(136,135.9473_wp), &
                                 Iso_t(137,136.94502_wp), &
                                 Iso_t(138,137.94025_wp), &
                                 Iso_t(139,138.93813_wp), &
                                 Iso_t(140,139.933674_wp), &
                                 Iso_t(141,140.932126_wp), &
                                 Iso_t(142,141.928116_wp), &
                                 Iso_t(143,142.92675_wp), &
                                 Iso_t(144,143.922963_wp), &
                                 Iso_t(145,144.921713_wp), &
                                 Iso_t(146,145.9183188_wp), &
                                 Iso_t(147,146.9191014_wp), &
                                 Iso_t(148,147.9181215_wp), &
                                 Iso_t(149,148.9193481_wp), &
                                 Iso_t(150,149.9186644_wp), &
                                 Iso_t(151,150.920356_wp), &
                                 Iso_t(153,152.921758_wp), &
                                 Iso_t(159,158.926397_wp), &
                                 Iso_t(161,160.9296775_wp), &
                                 Iso_t(162,161.930993_wp), &
                                 Iso_t(163,162.9341769_wp), &
                                 Iso_t(164,163.93583_wp), &
                                 Iso_t(165,164.93936_wp), &
                                 Iso_t(166,165.94146_wp), &
                                 Iso_t(167,166.94545_wp), &
                                 Iso_t(168,167.94808_wp), &
                                 Iso_t(169,168.9526_wp)]

  ElementList(65)%Symbol = adjustl(PTab(65)) ! Tb
  ElementList(65)%Natural = 1
  call mma_Allocate(ElementList(65)%Isotopes,37)
  ElementList(65)%Isotopes(:) = [ &
                                 Iso_t(159,158.9253547_wp), &
                                 Iso_t(135,134.96476_wp), &
                                 Iso_t(136,135.96129_wp), &
                                 Iso_t(137,136.95602_wp), &
                                 Iso_t(138,137.95312_wp), &
                                 Iso_t(139,138.94833_wp), &
                                 Iso_t(140,139.94581_wp), &
                                 Iso_t(141,140.94145_wp), &
                                 Iso_t(142,141.93928_wp), &
                                 Iso_t(143,142.935137_wp), &
                                 Iso_t(144,143.933045_wp), &
                                 Iso_t(145,144.92882_wp), &
                                 Iso_t(146,145.927253_wp), &
                                 Iso_t(147,146.9240548_wp), &
                                 Iso_t(148,147.924282_wp), &
                                 Iso_t(149,148.9232535_wp), &
                                 Iso_t(150,149.9236649_wp), &
                                 Iso_t(151,150.9231096_wp), &
                                 Iso_t(152,151.924083_wp), &
                                 Iso_t(153,152.9234424_wp), &
                                 Iso_t(154,153.924685_wp), &
                                 Iso_t(155,154.923511_wp), &
                                 Iso_t(156,155.9247552_wp), &
                                 Iso_t(157,156.924033_wp), &
                                 Iso_t(158,157.9254209_wp), &
                                 Iso_t(160,159.9271756_wp), &
                                 Iso_t(161,160.9275778_wp), &
                                 Iso_t(162,161.929495_wp), &
                                 Iso_t(163,162.9306547_wp), &
                                 Iso_t(164,163.93336_wp), &
                                 Iso_t(165,164.93498_wp), &
                                 Iso_t(166,165.93786_wp), &
                                 Iso_t(167,166.93996_wp), &
                                 Iso_t(168,167.9434_wp), &
                                 Iso_t(169,168.94597_wp), &
                                 Iso_t(170,169.94984_wp), &
                                 Iso_t(171,170.95273_wp)]

  ElementList(66)%Symbol = adjustl(PTab(66)) ! Dy
  ElementList(66)%Natural = 7
  call mma_Allocate(ElementList(66)%Isotopes,36)
  ElementList(66)%Isotopes(:) = [ &
                                 Iso_t(164,163.9291819_wp), &
                                 Iso_t(162,161.9268056_wp), &
                                 Iso_t(163,162.9287383_wp), &
                                 Iso_t(161,160.9269405_wp), &
                                 Iso_t(160,159.9252046_wp), &
                                 Iso_t(158,157.9244159_wp), &
                                 Iso_t(156,155.9242847_wp), &
                                 Iso_t(138,137.9625_wp), &
                                 Iso_t(139,138.95959_wp), &
                                 Iso_t(140,139.95402_wp), &
                                 Iso_t(141,140.95128_wp), &
                                 Iso_t(142,141.94619_wp), &
                                 Iso_t(143,142.943994_wp), &
                                 Iso_t(144,143.9392695_wp), &
                                 Iso_t(145,144.937474_wp), &
                                 Iso_t(146,145.9328445_wp), &
                                 Iso_t(147,146.9310827_wp), &
                                 Iso_t(148,147.927157_wp), &
                                 Iso_t(149,148.927322_wp), &
                                 Iso_t(150,149.9255933_wp), &
                                 Iso_t(151,150.9261916_wp), &
                                 Iso_t(152,151.9247253_wp), &
                                 Iso_t(153,152.9257724_wp), &
                                 Iso_t(154,153.9244293_wp), &
                                 Iso_t(155,154.925759_wp), &
                                 Iso_t(157,156.9254707_wp), &
                                 Iso_t(159,158.925747_wp), &
                                 Iso_t(165,164.9317105_wp), &
                                 Iso_t(166,165.9328139_wp), &
                                 Iso_t(167,166.935661_wp), &
                                 Iso_t(168,167.93713_wp), &
                                 Iso_t(169,168.94031_wp), &
                                 Iso_t(170,169.94239_wp), &
                                 Iso_t(171,170.94612_wp), &
                                 Iso_t(172,171.94846_wp), &
                                 Iso_t(173,172.95283_wp)]

  ElementList(67)%Symbol = adjustl(PTab(67)) ! Ho
  ElementList(67)%Natural = 1
  call mma_Allocate(ElementList(67)%Isotopes,36)
  ElementList(67)%Isotopes(:) = [ &
                                 Iso_t(165,164.9303288_wp), &
                                 Iso_t(140,139.96859_wp), &
                                 Iso_t(141,140.96311_wp), &
                                 Iso_t(142,141.96001_wp), &
                                 Iso_t(143,142.95486_wp), &
                                 Iso_t(144,143.9521097_wp), &
                                 Iso_t(145,144.9472674_wp), &
                                 Iso_t(146,145.9449935_wp), &
                                 Iso_t(147,146.9401423_wp), &
                                 Iso_t(148,147.937744_wp), &
                                 Iso_t(149,148.933803_wp), &
                                 Iso_t(150,149.933498_wp), &
                                 Iso_t(151,150.9316983_wp), &
                                 Iso_t(152,151.931724_wp), &
                                 Iso_t(153,152.9302064_wp), &
                                 Iso_t(154,153.9306068_wp), &
                                 Iso_t(155,154.929104_wp), &
                                 Iso_t(156,155.929706_wp), &
                                 Iso_t(157,156.928254_wp), &
                                 Iso_t(158,157.928946_wp), &
                                 Iso_t(159,158.9277197_wp), &
                                 Iso_t(160,159.928737_wp), &
                                 Iso_t(161,160.9278615_wp), &
                                 Iso_t(162,161.9291023_wp), &
                                 Iso_t(163,162.928741_wp), &
                                 Iso_t(164,163.9302403_wp), &
                                 Iso_t(166,165.9322909_wp), &
                                 Iso_t(167,166.9331385_wp), &
                                 Iso_t(168,167.935522_wp), &
                                 Iso_t(169,168.936878_wp), &
                                 Iso_t(170,169.939625_wp), &
                                 Iso_t(171,170.94147_wp), &
                                 Iso_t(172,171.94473_wp), &
                                 Iso_t(173,172.94702_wp), &
                                 Iso_t(174,173.95095_wp), &
                                 Iso_t(175,174.95362_wp)]

  ElementList(68)%Symbol = adjustl(PTab(68)) ! Er
  ElementList(68)%Natural = 6
  call mma_Allocate(ElementList(68)%Isotopes,36)
  ElementList(68)%Isotopes(:) = [ &
                                 Iso_t(166,165.9302995_wp), &
                                 Iso_t(168,167.9323767_wp), &
                                 Iso_t(167,166.9320546_wp), &
                                 Iso_t(170,169.9354702_wp), &
                                 Iso_t(164,163.9292088_wp), &
                                 Iso_t(162,161.9287884_wp), &
                                 Iso_t(142,141.9701_wp), &
                                 Iso_t(143,142.96662_wp), &
                                 Iso_t(144,143.9607_wp), &
                                 Iso_t(145,144.95805_wp), &
                                 Iso_t(146,145.9524184_wp), &
                                 Iso_t(147,146.949964_wp), &
                                 Iso_t(148,147.944735_wp), &
                                 Iso_t(149,148.942306_wp), &
                                 Iso_t(150,149.937916_wp), &
                                 Iso_t(151,150.937449_wp), &
                                 Iso_t(152,151.935057_wp), &
                                 Iso_t(153,152.93508_wp), &
                                 Iso_t(154,153.9327908_wp), &
                                 Iso_t(155,154.9332159_wp), &
                                 Iso_t(156,155.931067_wp), &
                                 Iso_t(157,156.931949_wp), &
                                 Iso_t(158,157.929893_wp), &
                                 Iso_t(159,158.9306918_wp), &
                                 Iso_t(160,159.929077_wp), &
                                 Iso_t(161,160.9300046_wp), &
                                 Iso_t(163,162.9300408_wp), &
                                 Iso_t(165,164.9307345_wp), &
                                 Iso_t(169,168.9345968_wp), &
                                 Iso_t(171,170.9380357_wp), &
                                 Iso_t(172,171.9393619_wp), &
                                 Iso_t(173,172.9424_wp), &
                                 Iso_t(174,173.94423_wp), &
                                 Iso_t(175,174.94777_wp), &
                                 Iso_t(176,175.94994_wp), &
                                 Iso_t(177,176.95399_wp)]

  ElementList(69)%Symbol = adjustl(PTab(69)) ! Tm
  ElementList(69)%Natural = 1
  call mma_Allocate(ElementList(69)%Isotopes,36)
  ElementList(69)%Isotopes(:) = [ &
                                 Iso_t(169,168.9342179_wp), &
                                 Iso_t(144,143.97628_wp), &
                                 Iso_t(145,144.97039_wp), &
                                 Iso_t(146,145.96684_wp), &
                                 Iso_t(147,146.9613799_wp), &
                                 Iso_t(148,147.958384_wp), &
                                 Iso_t(149,148.95289_wp), &
                                 Iso_t(150,149.95009_wp), &
                                 Iso_t(151,150.945488_wp), &
                                 Iso_t(152,151.944422_wp), &
                                 Iso_t(153,152.94204_wp), &
                                 Iso_t(154,153.94157_wp), &
                                 Iso_t(155,154.93921_wp), &
                                 Iso_t(156,155.938992_wp), &
                                 Iso_t(157,156.936944_wp), &
                                 Iso_t(158,157.93698_wp), &
                                 Iso_t(159,158.934975_wp), &
                                 Iso_t(160,159.935263_wp), &
                                 Iso_t(161,160.933549_wp), &
                                 Iso_t(162,161.934002_wp), &
                                 Iso_t(163,162.9326592_wp), &
                                 Iso_t(164,163.933544_wp), &
                                 Iso_t(165,164.9324431_wp), &
                                 Iso_t(166,165.933561_wp), &
                                 Iso_t(167,166.9328562_wp), &
                                 Iso_t(168,167.9341774_wp), &
                                 Iso_t(170,169.935806_wp), &
                                 Iso_t(171,170.9364339_wp), &
                                 Iso_t(172,171.9384055_wp), &
                                 Iso_t(173,172.9396084_wp), &
                                 Iso_t(174,173.942173_wp), &
                                 Iso_t(175,174.943841_wp), &
                                 Iso_t(176,175.947_wp), &
                                 Iso_t(177,176.94904_wp), &
                                 Iso_t(178,177.95264_wp), &
                                 Iso_t(179,178.95534_wp)]

  ElementList(70)%Symbol = adjustl(PTab(70)) ! Yb
  ElementList(70)%Natural = 7
  call mma_Allocate(ElementList(70)%Isotopes,34)
  ElementList(70)%Isotopes(:) = [ &
                                 Iso_t(174,173.9388664_wp), &
                                 Iso_t(172,171.9363859_wp), &
                                 Iso_t(173,172.9382151_wp), &
                                 Iso_t(171,170.9363302_wp), &
                                 Iso_t(176,175.9425764_wp), &
                                 Iso_t(170,169.9347664_wp), &
                                 Iso_t(168,167.9338896_wp), &
                                 Iso_t(148,147.96758_wp), &
                                 Iso_t(149,148.96436_wp), &
                                 Iso_t(150,149.95852_wp), &
                                 Iso_t(151,150.9554_wp), &
                                 Iso_t(152,151.95027_wp), &
                                 Iso_t(153,152.94932_wp), &
                                 Iso_t(154,153.946396_wp), &
                                 Iso_t(155,154.945783_wp), &
                                 Iso_t(156,155.942825_wp), &
                                 Iso_t(157,156.942645_wp), &
                                 Iso_t(158,157.9398705_wp), &
                                 Iso_t(159,158.940055_wp), &
                                 Iso_t(160,159.937557_wp), &
                                 Iso_t(161,160.937907_wp), &
                                 Iso_t(162,161.935774_wp), &
                                 Iso_t(163,162.93634_wp), &
                                 Iso_t(164,163.934495_wp), &
                                 Iso_t(165,164.93527_wp), &
                                 Iso_t(166,165.9338747_wp), &
                                 Iso_t(167,166.934953_wp), &
                                 Iso_t(169,168.9351825_wp), &
                                 Iso_t(175,174.9412808_wp), &
                                 Iso_t(177,176.9452656_wp), &
                                 Iso_t(178,177.946651_wp), &
                                 Iso_t(179,178.95004_wp), &
                                 Iso_t(180,179.95212_wp), &
                                 Iso_t(181,180.95589_wp)]

  ElementList(71)%Symbol = adjustl(PTab(71)) ! Lu
  ElementList(71)%Natural = 2
  call mma_Allocate(ElementList(71)%Isotopes,36)
  ElementList(71)%Isotopes(:) = [ &
                                 Iso_t(175,174.9407752_wp), &
                                 Iso_t(176,175.9426897_wp), &
                                 Iso_t(150,149.97355_wp), &
                                 Iso_t(151,150.96768_wp), &
                                 Iso_t(152,151.96412_wp), &
                                 Iso_t(153,152.95875_wp), &
                                 Iso_t(154,153.95736_wp), &
                                 Iso_t(155,154.954321_wp), &
                                 Iso_t(156,155.953033_wp), &
                                 Iso_t(157,156.950127_wp), &
                                 Iso_t(158,157.949316_wp), &
                                 Iso_t(159,158.946636_wp), &
                                 Iso_t(160,159.946033_wp), &
                                 Iso_t(161,160.943572_wp), &
                                 Iso_t(162,161.943283_wp), &
                                 Iso_t(163,162.941179_wp), &
                                 Iso_t(164,163.941339_wp), &
                                 Iso_t(165,164.939407_wp), &
                                 Iso_t(166,165.939859_wp), &
                                 Iso_t(167,166.93827_wp), &
                                 Iso_t(168,167.938736_wp), &
                                 Iso_t(169,168.9376441_wp), &
                                 Iso_t(170,169.938478_wp), &
                                 Iso_t(171,170.937917_wp), &
                                 Iso_t(172,171.9390891_wp), &
                                 Iso_t(173,172.938934_wp), &
                                 Iso_t(174,173.9403409_wp), &
                                 Iso_t(177,176.9437615_wp), &
                                 Iso_t(178,177.945958_wp), &
                                 Iso_t(179,178.9473309_wp), &
                                 Iso_t(180,179.949888_wp), &
                                 Iso_t(181,180.95191_wp), &
                                 Iso_t(182,181.95504_wp), &
                                 Iso_t(183,182.957363_wp), &
                                 Iso_t(184,183.96091_wp), &
                                 Iso_t(185,184.96362_wp)]

  ElementList(72)%Symbol = adjustl(PTab(72)) ! Hf
  ElementList(72)%Natural = 6
  call mma_Allocate(ElementList(72)%Isotopes,37)
  ElementList(72)%Isotopes(:) = [ &
                                 Iso_t(180,179.946557_wp), &
                                 Iso_t(178,177.9437058_wp), &
                                 Iso_t(177,176.9432277_wp), &
                                 Iso_t(179,178.9458232_wp), &
                                 Iso_t(176,175.9414076_wp), &
                                 Iso_t(174,173.9400461_wp), &
                                 Iso_t(153,152.97069_wp), &
                                 Iso_t(154,153.96486_wp), &
                                 Iso_t(155,154.96311_wp), &
                                 Iso_t(156,155.95935_wp), &
                                 Iso_t(157,156.95824_wp), &
                                 Iso_t(158,157.954801_wp), &
                                 Iso_t(159,158.953996_wp), &
                                 Iso_t(160,159.950691_wp), &
                                 Iso_t(161,160.950278_wp), &
                                 Iso_t(162,161.9472148_wp), &
                                 Iso_t(163,162.947113_wp), &
                                 Iso_t(164,163.944371_wp), &
                                 Iso_t(165,164.944567_wp), &
                                 Iso_t(166,165.94218_wp), &
                                 Iso_t(167,166.9426_wp), &
                                 Iso_t(168,167.940568_wp), &
                                 Iso_t(169,168.941259_wp), &
                                 Iso_t(170,169.939609_wp), &
                                 Iso_t(171,170.940492_wp), &
                                 Iso_t(172,171.93945_wp), &
                                 Iso_t(173,172.940513_wp), &
                                 Iso_t(175,174.9415092_wp), &
                                 Iso_t(181,180.9491083_wp), &
                                 Iso_t(182,181.9505612_wp), &
                                 Iso_t(183,182.95353_wp), &
                                 Iso_t(184,183.955446_wp), &
                                 Iso_t(185,184.958862_wp), &
                                 Iso_t(186,185.960897_wp), &
                                 Iso_t(187,186.96477_wp), &
                                 Iso_t(188,187.96685_wp), &
                                 Iso_t(189,188.97084_wp)]

  ElementList(73)%Symbol = adjustl(PTab(73)) ! Ta
  ElementList(73)%Natural = 2
  call mma_Allocate(ElementList(73)%Isotopes,38)
  ElementList(73)%Isotopes(:) = [ &
                                 Iso_t(181,180.9479958_wp), &
                                 Iso_t(180,179.9474648_wp), &
                                 Iso_t(155,154.97424_wp), &
                                 Iso_t(156,155.97203_wp), &
                                 Iso_t(157,156.96818_wp), &
                                 Iso_t(158,157.96654_wp), &
                                 Iso_t(159,158.963023_wp), &
                                 Iso_t(160,159.961488_wp), &
                                 Iso_t(161,160.958452_wp), &
                                 Iso_t(162,161.957294_wp), &
                                 Iso_t(163,162.954337_wp), &
                                 Iso_t(164,163.953534_wp), &
                                 Iso_t(165,164.950781_wp), &
                                 Iso_t(166,165.950512_wp), &
                                 Iso_t(167,166.948093_wp), &
                                 Iso_t(168,167.948047_wp), &
                                 Iso_t(169,168.946011_wp), &
                                 Iso_t(170,169.946175_wp), &
                                 Iso_t(171,170.944476_wp), &
                                 Iso_t(172,171.944895_wp), &
                                 Iso_t(173,172.94375_wp), &
                                 Iso_t(174,173.944454_wp), &
                                 Iso_t(175,174.943737_wp), &
                                 Iso_t(176,175.944857_wp), &
                                 Iso_t(177,176.9444795_wp), &
                                 Iso_t(178,177.945678_wp), &
                                 Iso_t(179,178.9459366_wp), &
                                 Iso_t(182,181.9501519_wp), &
                                 Iso_t(183,182.9513726_wp), &
                                 Iso_t(184,183.954008_wp), &
                                 Iso_t(185,184.955559_wp), &
                                 Iso_t(186,185.958551_wp), &
                                 Iso_t(187,186.960386_wp), &
                                 Iso_t(188,187.963916_wp), &
                                 Iso_t(189,188.96583_wp), &
                                 Iso_t(190,189.96939_wp), &
                                 Iso_t(191,190.97156_wp), &
                                 Iso_t(192,191.97514_wp)]

  ElementList(74)%Symbol = adjustl(PTab(74)) ! W
  ElementList(74)%Natural = 5
  call mma_Allocate(ElementList(74)%Isotopes,38)
  ElementList(74)%Isotopes(:) = [ &
                                 Iso_t(184,183.95093092_wp), &
                                 Iso_t(186,185.9543628_wp), &
                                 Iso_t(182,181.94820394_wp), &
                                 Iso_t(183,182.95022275_wp), &
                                 Iso_t(180,179.9467108_wp), &
                                 Iso_t(157,156.97884_wp), &
                                 Iso_t(158,157.97456_wp), &
                                 Iso_t(159,158.97264_wp), &
                                 Iso_t(160,159.96846_wp), &
                                 Iso_t(161,160.9672_wp), &
                                 Iso_t(162,161.963499_wp), &
                                 Iso_t(163,162.962524_wp), &
                                 Iso_t(164,163.958961_wp), &
                                 Iso_t(165,164.958281_wp), &
                                 Iso_t(166,165.955031_wp), &
                                 Iso_t(167,166.954805_wp), &
                                 Iso_t(168,167.951806_wp), &
                                 Iso_t(169,168.951779_wp), &
                                 Iso_t(170,169.949232_wp), &
                                 Iso_t(171,170.949451_wp), &
                                 Iso_t(172,171.947292_wp), &
                                 Iso_t(173,172.947689_wp), &
                                 Iso_t(174,173.946079_wp), &
                                 Iso_t(175,174.946717_wp), &
                                 Iso_t(176,175.945634_wp), &
                                 Iso_t(177,176.946643_wp), &
                                 Iso_t(178,177.945883_wp), &
                                 Iso_t(179,178.947077_wp), &
                                 Iso_t(181,180.9481978_wp), &
                                 Iso_t(185,184.95341897_wp), &
                                 Iso_t(187,186.9571588_wp), &
                                 Iso_t(188,187.9584862_wp), &
                                 Iso_t(189,188.961763_wp), &
                                 Iso_t(190,189.963091_wp), &
                                 Iso_t(191,190.966531_wp), &
                                 Iso_t(192,191.96817_wp), &
                                 Iso_t(193,192.97178_wp), &
                                 Iso_t(194,193.97367_wp)]

  ElementList(75)%Symbol = adjustl(PTab(75)) ! Re
  ElementList(75)%Natural = 2
  call mma_Allocate(ElementList(75)%Isotopes,40)
  ElementList(75)%Isotopes(:) = [ &
                                 Iso_t(187,186.9557501_wp), &
                                 Iso_t(185,184.9529545_wp), &
                                 Iso_t(159,158.98418_wp), &
                                 Iso_t(160,159.98182_wp), &
                                 Iso_t(161,160.97757_wp), &
                                 Iso_t(162,161.97584_wp), &
                                 Iso_t(163,162.97208_wp), &
                                 Iso_t(164,163.970453_wp), &
                                 Iso_t(165,164.967103_wp), &
                                 Iso_t(166,165.965761_wp), &
                                 Iso_t(167,166.962595_wp), &
                                 Iso_t(168,167.961573_wp), &
                                 Iso_t(169,168.958766_wp), &
                                 Iso_t(170,169.95822_wp), &
                                 Iso_t(171,170.955716_wp), &
                                 Iso_t(172,171.95542_wp), &
                                 Iso_t(173,172.953243_wp), &
                                 Iso_t(174,173.953115_wp), &
                                 Iso_t(175,174.951381_wp), &
                                 Iso_t(176,175.951623_wp), &
                                 Iso_t(177,176.950328_wp), &
                                 Iso_t(178,177.950989_wp), &
                                 Iso_t(179,178.949989_wp), &
                                 Iso_t(180,179.950792_wp), &
                                 Iso_t(181,180.950058_wp), &
                                 Iso_t(182,181.95121_wp), &
                                 Iso_t(183,182.9508196_wp), &
                                 Iso_t(184,183.9525228_wp), &
                                 Iso_t(186,185.9549856_wp), &
                                 Iso_t(188,187.9581115_wp), &
                                 Iso_t(189,188.959226_wp), &
                                 Iso_t(190,189.961744_wp), &
                                 Iso_t(191,190.963122_wp), &
                                 Iso_t(192,191.966088_wp), &
                                 Iso_t(193,192.967541_wp), &
                                 Iso_t(194,193.97076_wp), &
                                 Iso_t(195,194.97254_wp), &
                                 Iso_t(196,195.9758_wp), &
                                 Iso_t(197,196.97799_wp), &
                                 Iso_t(198,197.9816_wp)]

  ElementList(76)%Symbol = adjustl(PTab(76)) ! Os
  ElementList(76)%Natural = 7
  call mma_Allocate(ElementList(76)%Isotopes,42)
  ElementList(76)%Isotopes(:) = [ &
                                 Iso_t(192,191.961477_wp), &
                                 Iso_t(190,189.9584437_wp), &
                                 Iso_t(189,188.9581442_wp), &
                                 Iso_t(188,187.9558352_wp), &
                                 Iso_t(187,186.9557474_wp), &
                                 Iso_t(186,185.953835_wp), &
                                 Iso_t(184,183.9524885_wp), &
                                 Iso_t(161,160.98903_wp), &
                                 Iso_t(162,161.98443_wp), &
                                 Iso_t(163,162.98241_wp), &
                                 Iso_t(164,163.97802_wp), &
                                 Iso_t(165,164.9766_wp), &
                                 Iso_t(166,165.972692_wp), &
                                 Iso_t(167,166.971549_wp), &
                                 Iso_t(168,167.967808_wp), &
                                 Iso_t(169,168.967018_wp), &
                                 Iso_t(170,169.963578_wp), &
                                 Iso_t(171,170.963174_wp), &
                                 Iso_t(172,171.960017_wp), &
                                 Iso_t(173,172.959808_wp), &
                                 Iso_t(174,173.957064_wp), &
                                 Iso_t(175,174.956945_wp), &
                                 Iso_t(176,175.954806_wp), &
                                 Iso_t(177,176.954966_wp), &
                                 Iso_t(178,177.953254_wp), &
                                 Iso_t(179,178.953817_wp), &
                                 Iso_t(180,179.952375_wp), &
                                 Iso_t(181,180.953247_wp), &
                                 Iso_t(182,181.95211_wp), &
                                 Iso_t(183,182.953125_wp), &
                                 Iso_t(185,184.9540417_wp), &
                                 Iso_t(191,190.9609264_wp), &
                                 Iso_t(193,192.9641479_wp), &
                                 Iso_t(194,193.9651772_wp), &
                                 Iso_t(195,194.968318_wp), &
                                 Iso_t(196,195.969641_wp), &
                                 Iso_t(197,196.97283_wp), &
                                 Iso_t(198,197.97441_wp), &
                                 Iso_t(199,198.97801_wp), &
                                 Iso_t(200,199.97984_wp), &
                                 Iso_t(201,200.98364_wp), &
                                 Iso_t(202,201.98595_wp)]

  ElementList(77)%Symbol = adjustl(PTab(77)) ! Ir
  ElementList(77)%Natural = 2
  call mma_Allocate(ElementList(77)%Isotopes,41)
  ElementList(77)%Isotopes(:) = [ &
                                 Iso_t(193,192.9629216_wp), &
                                 Iso_t(191,190.9605893_wp), &
                                 Iso_t(164,163.99191_wp), &
                                 Iso_t(165,164.9875_wp), &
                                 Iso_t(166,165.98566_wp), &
                                 Iso_t(167,166.981666_wp), &
                                 Iso_t(168,167.979907_wp), &
                                 Iso_t(169,168.976298_wp), &
                                 Iso_t(170,169.974922_wp), &
                                 Iso_t(171,170.97164_wp), &
                                 Iso_t(172,171.970607_wp), &
                                 Iso_t(173,172.967506_wp), &
                                 Iso_t(174,173.966861_wp), &
                                 Iso_t(175,174.96415_wp), &
                                 Iso_t(176,175.96365_wp), &
                                 Iso_t(177,176.961301_wp), &
                                 Iso_t(178,177.961082_wp), &
                                 Iso_t(179,178.95912_wp), &
                                 Iso_t(180,179.959229_wp), &
                                 Iso_t(181,180.957625_wp), &
                                 Iso_t(182,181.958076_wp), &
                                 Iso_t(183,182.95684_wp), &
                                 Iso_t(184,183.957476_wp), &
                                 Iso_t(185,184.956698_wp), &
                                 Iso_t(186,185.957944_wp), &
                                 Iso_t(187,186.957542_wp), &
                                 Iso_t(188,187.958828_wp), &
                                 Iso_t(189,188.958715_wp), &
                                 Iso_t(190,189.9605412_wp), &
                                 Iso_t(192,191.9626002_wp), &
                                 Iso_t(194,193.9650735_wp), &
                                 Iso_t(195,194.9659747_wp), &
                                 Iso_t(196,195.968397_wp), &
                                 Iso_t(197,196.969655_wp), &
                                 Iso_t(198,197.97228_wp), &
                                 Iso_t(199,198.973805_wp), &
                                 Iso_t(200,199.9768_wp), &
                                 Iso_t(201,200.97864_wp), &
                                 Iso_t(202,201.98199_wp), &
                                 Iso_t(203,202.98423_wp), &
                                 Iso_t(204,203.9896_wp)]

  ElementList(78)%Symbol = adjustl(PTab(78)) ! Pt
  ElementList(78)%Natural = 6
  call mma_Allocate(ElementList(78)%Isotopes,41)
  ElementList(78)%Isotopes(:) = [ &
                                 Iso_t(195,194.9647917_wp), &
                                 Iso_t(194,193.9626809_wp), &
                                 Iso_t(196,195.96495209_wp), &
                                 Iso_t(198,197.9678949_wp), &
                                 Iso_t(192,191.9610387_wp), &
                                 Iso_t(190,189.9599297_wp), &
                                 Iso_t(166,165.99486_wp), &
                                 Iso_t(167,166.99269_wp), &
                                 Iso_t(168,167.98813_wp), &
                                 Iso_t(169,168.98657_wp), &
                                 Iso_t(170,169.982496_wp), &
                                 Iso_t(171,170.981245_wp), &
                                 Iso_t(172,171.977351_wp), &
                                 Iso_t(173,172.976443_wp), &
                                 Iso_t(174,173.97282_wp), &
                                 Iso_t(175,174.97241_wp), &
                                 Iso_t(176,175.968938_wp), &
                                 Iso_t(177,176.96847_wp), &
                                 Iso_t(178,177.96565_wp), &
                                 Iso_t(179,178.965359_wp), &
                                 Iso_t(180,179.963032_wp), &
                                 Iso_t(181,180.963098_wp), &
                                 Iso_t(182,181.961172_wp), &
                                 Iso_t(183,182.961597_wp), &
                                 Iso_t(184,183.959915_wp), &
                                 Iso_t(185,184.960614_wp), &
                                 Iso_t(186,185.959351_wp), &
                                 Iso_t(187,186.960617_wp), &
                                 Iso_t(188,187.9593889_wp), &
                                 Iso_t(189,188.960831_wp), &
                                 Iso_t(191,190.9616729_wp), &
                                 Iso_t(193,192.9629824_wp), &
                                 Iso_t(197,196.96734069_wp), &
                                 Iso_t(199,198.9705952_wp), &
                                 Iso_t(200,199.971443_wp), &
                                 Iso_t(201,200.974513_wp), &
                                 Iso_t(202,201.975639_wp), &
                                 Iso_t(203,202.97893_wp), &
                                 Iso_t(204,203.98076_wp), &
                                 Iso_t(205,204.98608_wp), &
                                 Iso_t(206,205.98966_wp)]

  ElementList(79)%Symbol = adjustl(PTab(79)) ! Au
  ElementList(79)%Natural = 1
  call mma_Allocate(ElementList(79)%Isotopes,42)
  ElementList(79)%Isotopes(:) = [ &
                                 Iso_t(197,196.96656879_wp), &
                                 Iso_t(169,168.99808_wp), &
                                 Iso_t(170,169.99597_wp), &
                                 Iso_t(171,170.991876_wp), &
                                 Iso_t(172,171.989942_wp), &
                                 Iso_t(173,172.986241_wp), &
                                 Iso_t(174,173.984717_wp), &
                                 Iso_t(175,174.981304_wp), &
                                 Iso_t(176,175.98025_wp), &
                                 Iso_t(177,176.97687_wp), &
                                 Iso_t(178,177.976032_wp), &
                                 Iso_t(179,178.973174_wp), &
                                 Iso_t(180,179.972523_wp), &
                                 Iso_t(181,180.970079_wp), &
                                 Iso_t(182,181.969618_wp), &
                                 Iso_t(183,182.967591_wp), &
                                 Iso_t(184,183.967452_wp), &
                                 Iso_t(185,184.96579_wp), &
                                 Iso_t(186,185.965953_wp), &
                                 Iso_t(187,186.964543_wp), &
                                 Iso_t(188,187.965349_wp), &
                                 Iso_t(189,188.963948_wp), &
                                 Iso_t(190,189.964698_wp), &
                                 Iso_t(191,190.963702_wp), &
                                 Iso_t(192,191.964814_wp), &
                                 Iso_t(193,192.9641373_wp), &
                                 Iso_t(194,193.9654178_wp), &
                                 Iso_t(195,194.9650352_wp), &
                                 Iso_t(196,195.9665699_wp), &
                                 Iso_t(198,197.96824242_wp), &
                                 Iso_t(199,198.96876528_wp), &
                                 Iso_t(200,199.970756_wp), &
                                 Iso_t(201,200.9716575_wp), &
                                 Iso_t(202,201.973856_wp), &
                                 Iso_t(203,202.9751544_wp), &
                                 Iso_t(204,203.97783_wp), &
                                 Iso_t(205,204.97985_wp), &
                                 Iso_t(206,205.98474_wp), &
                                 Iso_t(207,206.9884_wp), &
                                 Iso_t(208,207.99345_wp), &
                                 Iso_t(209,208.99735_wp), &
                                 Iso_t(210,210.0025_wp)]

  ElementList(80)%Symbol = adjustl(PTab(80)) ! Hg
  ElementList(80)%Natural = 7
  call mma_Allocate(ElementList(80)%Isotopes,46)
  ElementList(80)%Isotopes(:) = [ &
                                 Iso_t(202,201.9706434_wp), &
                                 Iso_t(200,199.96832659_wp), &
                                 Iso_t(199,198.96828064_wp), &
                                 Iso_t(201,200.97030284_wp), &
                                 Iso_t(198,197.9667686_wp), &
                                 Iso_t(204,203.97349398_wp), &
                                 Iso_t(196,195.9658326_wp), &
                                 Iso_t(171,171.00353_wp), &
                                 Iso_t(172,171.99881_wp), &
                                 Iso_t(173,172.99709_wp), &
                                 Iso_t(174,173.992865_wp), &
                                 Iso_t(175,174.991441_wp), &
                                 Iso_t(176,175.987361_wp), &
                                 Iso_t(177,176.986277_wp), &
                                 Iso_t(178,177.982484_wp), &
                                 Iso_t(179,178.981831_wp), &
                                 Iso_t(180,179.97826_wp), &
                                 Iso_t(181,180.977819_wp), &
                                 Iso_t(182,181.974689_wp), &
                                 Iso_t(183,182.9744448_wp), &
                                 Iso_t(184,183.971714_wp), &
                                 Iso_t(185,184.971899_wp), &
                                 Iso_t(186,185.969362_wp), &
                                 Iso_t(187,186.969814_wp), &
                                 Iso_t(188,187.967567_wp), &
                                 Iso_t(189,188.968195_wp), &
                                 Iso_t(190,189.966323_wp), &
                                 Iso_t(191,190.967157_wp), &
                                 Iso_t(192,191.965635_wp), &
                                 Iso_t(193,192.966653_wp), &
                                 Iso_t(194,193.9654491_wp), &
                                 Iso_t(195,194.966721_wp), &
                                 Iso_t(197,196.9672128_wp), &
                                 Iso_t(203,202.9728728_wp), &
                                 Iso_t(205,204.9760734_wp), &
                                 Iso_t(206,205.977514_wp), &
                                 Iso_t(207,206.9823_wp), &
                                 Iso_t(208,207.985759_wp), &
                                 Iso_t(209,208.99072_wp), &
                                 Iso_t(210,209.99424_wp), &
                                 Iso_t(211,210.99933_wp), &
                                 Iso_t(212,212.00296_wp), &
                                 Iso_t(213,213.00823_wp), &
                                 Iso_t(214,214.012_wp), &
                                 Iso_t(215,215.0174_wp), &
                                 Iso_t(216,216.02132_wp)]

  ElementList(81)%Symbol = adjustl(PTab(81)) ! Tl
  ElementList(81)%Natural = 2
  call mma_Allocate(ElementList(81)%Isotopes,43)
  ElementList(81)%Isotopes(:) = [ &
                                 Iso_t(205,204.9744278_wp), &
                                 Iso_t(203,202.9723446_wp), &
                                 Iso_t(176,176.000624_wp), &
                                 Iso_t(177,176.996431_wp), &
                                 Iso_t(178,177.99485_wp), &
                                 Iso_t(179,178.991111_wp), &
                                 Iso_t(180,179.990057_wp), &
                                 Iso_t(181,180.98626_wp), &
                                 Iso_t(182,181.985713_wp), &
                                 Iso_t(183,182.982193_wp), &
                                 Iso_t(184,183.981886_wp), &
                                 Iso_t(185,184.978789_wp), &
                                 Iso_t(186,185.978651_wp), &
                                 Iso_t(187,186.9759063_wp), &
                                 Iso_t(188,187.976021_wp), &
                                 Iso_t(189,188.973588_wp), &
                                 Iso_t(190,189.973828_wp), &
                                 Iso_t(191,190.9717842_wp), &
                                 Iso_t(192,191.972225_wp), &
                                 Iso_t(193,192.970502_wp), &
                                 Iso_t(194,193.971081_wp), &
                                 Iso_t(195,194.969774_wp), &
                                 Iso_t(196,195.970481_wp), &
                                 Iso_t(197,196.969576_wp), &
                                 Iso_t(198,197.970483_wp), &
                                 Iso_t(199,198.969877_wp), &
                                 Iso_t(200,199.9709633_wp), &
                                 Iso_t(201,200.970822_wp), &
                                 Iso_t(202,201.972102_wp), &
                                 Iso_t(204,203.9738639_wp), &
                                 Iso_t(206,205.9761106_wp), &
                                 Iso_t(207,206.9774197_wp), &
                                 Iso_t(208,207.982019_wp), &
                                 Iso_t(209,208.9853594_wp), &
                                 Iso_t(210,209.990074_wp), &
                                 Iso_t(211,210.993475_wp), &
                                 Iso_t(212,211.99834_wp), &
                                 Iso_t(213,213.001915_wp), &
                                 Iso_t(214,214.00694_wp), &
                                 Iso_t(215,215.01064_wp), &
                                 Iso_t(216,216.0158_wp), &
                                 Iso_t(217,217.01966_wp), &
                                 Iso_t(218,218.02479_wp)]

  ElementList(82)%Symbol = adjustl(PTab(82)) ! Pb
  ElementList(82)%Natural = 4
  call mma_Allocate(ElementList(82)%Isotopes,43)
  ElementList(82)%Isotopes(:) = [ &
                                 Iso_t(208,207.9766525_wp), &
                                 Iso_t(206,205.9744657_wp), &
                                 Iso_t(207,206.9758973_wp), &
                                 Iso_t(204,203.973044_wp), &
                                 Iso_t(178,178.003831_wp), &
                                 Iso_t(179,179.002201_wp), &
                                 Iso_t(180,179.997928_wp), &
                                 Iso_t(181,180.996653_wp), &
                                 Iso_t(182,181.992672_wp), &
                                 Iso_t(183,182.991872_wp), &
                                 Iso_t(184,183.988136_wp), &
                                 Iso_t(185,184.98761_wp), &
                                 Iso_t(186,185.984238_wp), &
                                 Iso_t(187,186.9839109_wp), &
                                 Iso_t(188,187.980875_wp), &
                                 Iso_t(189,188.980807_wp), &
                                 Iso_t(190,189.978082_wp), &
                                 Iso_t(191,190.978276_wp), &
                                 Iso_t(192,191.975775_wp), &
                                 Iso_t(193,192.976173_wp), &
                                 Iso_t(194,193.974012_wp), &
                                 Iso_t(195,194.974543_wp), &
                                 Iso_t(196,195.972774_wp), &
                                 Iso_t(197,196.9734312_wp), &
                                 Iso_t(198,197.972034_wp), &
                                 Iso_t(199,198.972913_wp), &
                                 Iso_t(200,199.971819_wp), &
                                 Iso_t(201,200.972883_wp), &
                                 Iso_t(202,201.972152_wp), &
                                 Iso_t(203,202.9733911_wp), &
                                 Iso_t(205,204.9744822_wp), &
                                 Iso_t(209,208.9810905_wp), &
                                 Iso_t(210,209.9841889_wp), &
                                 Iso_t(211,210.9887371_wp), &
                                 Iso_t(212,211.9918977_wp), &
                                 Iso_t(213,212.9965629_wp), &
                                 Iso_t(214,213.9998059_wp), &
                                 Iso_t(215,215.00474_wp), &
                                 Iso_t(216,216.00803_wp), &
                                 Iso_t(217,217.01314_wp), &
                                 Iso_t(218,218.01659_wp), &
                                 Iso_t(219,219.02177_wp), &
                                 Iso_t(220,220.02541_wp)]

  ElementList(83)%Symbol = adjustl(PTab(83)) ! Bi
  ElementList(83)%Natural = 1
  call mma_Allocate(ElementList(83)%Isotopes,41)
  ElementList(83)%Isotopes(:) = [ &
                                 Iso_t(209,208.9803991_wp), &
                                 Iso_t(184,184.001275_wp), &
                                 Iso_t(185,184.9976_wp), &
                                 Iso_t(186,185.996644_wp), &
                                 Iso_t(187,186.993147_wp), &
                                 Iso_t(188,187.992287_wp), &
                                 Iso_t(189,188.989195_wp), &
                                 Iso_t(190,189.988622_wp), &
                                 Iso_t(191,190.9857866_wp), &
                                 Iso_t(192,191.985469_wp), &
                                 Iso_t(193,192.98296_wp), &
                                 Iso_t(194,193.982785_wp), &
                                 Iso_t(195,194.9806488_wp), &
                                 Iso_t(196,195.980667_wp), &
                                 Iso_t(197,196.9788651_wp), &
                                 Iso_t(198,197.979206_wp), &
                                 Iso_t(199,198.977673_wp), &
                                 Iso_t(200,199.978131_wp), &
                                 Iso_t(201,200.97701_wp), &
                                 Iso_t(202,201.977734_wp), &
                                 Iso_t(203,202.976893_wp), &
                                 Iso_t(204,203.9778361_wp), &
                                 Iso_t(205,204.9773867_wp), &
                                 Iso_t(206,205.9784993_wp), &
                                 Iso_t(207,206.978471_wp), &
                                 Iso_t(208,207.9797425_wp), &
                                 Iso_t(210,209.9841207_wp), &
                                 Iso_t(211,210.9872697_wp), &
                                 Iso_t(212,211.991286_wp), &
                                 Iso_t(213,212.9943851_wp), &
                                 Iso_t(214,213.998712_wp), &
                                 Iso_t(215,215.00177_wp), &
                                 Iso_t(216,216.006306_wp), &
                                 Iso_t(217,217.009372_wp), &
                                 Iso_t(218,218.014188_wp), &
                                 Iso_t(219,219.01748_wp), &
                                 Iso_t(220,220.02235_wp), &
                                 Iso_t(221,221.02587_wp), &
                                 Iso_t(222,222.03078_wp), &
                                 Iso_t(223,223.0345_wp), &
                                 Iso_t(224,224.03947_wp)]

  ElementList(84)%Symbol = adjustl(PTab(84)) ! Po
  ElementList(84)%Natural = 0
  call mma_Allocate(ElementList(84)%Isotopes,42)
  ElementList(84)%Isotopes(:) = [ &
                                 Iso_t(209,208.9824308_wp), &
                                 Iso_t(186,186.004393_wp), &
                                 Iso_t(187,187.003041_wp), &
                                 Iso_t(188,187.999416_wp), &
                                 Iso_t(189,188.998473_wp), &
                                 Iso_t(190,189.995101_wp), &
                                 Iso_t(191,190.9945585_wp), &
                                 Iso_t(192,191.991336_wp), &
                                 Iso_t(193,192.991026_wp), &
                                 Iso_t(194,193.988186_wp), &
                                 Iso_t(195,194.988126_wp), &
                                 Iso_t(196,195.985526_wp), &
                                 Iso_t(197,196.98566_wp), &
                                 Iso_t(198,197.983389_wp), &
                                 Iso_t(199,198.983667_wp), &
                                 Iso_t(200,199.981799_wp), &
                                 Iso_t(201,200.9822598_wp), &
                                 Iso_t(202,201.980758_wp), &
                                 Iso_t(203,202.9814161_wp), &
                                 Iso_t(204,203.98031_wp), &
                                 Iso_t(205,204.981203_wp), &
                                 Iso_t(206,205.980474_wp), &
                                 Iso_t(207,206.9815938_wp), &
                                 Iso_t(208,207.9812461_wp), &
                                 Iso_t(210,209.9828741_wp), &
                                 Iso_t(211,210.9866536_wp), &
                                 Iso_t(212,211.9888684_wp), &
                                 Iso_t(213,212.9928576_wp), &
                                 Iso_t(214,213.9952017_wp), &
                                 Iso_t(215,214.9994201_wp), &
                                 Iso_t(216,216.0019152_wp), &
                                 Iso_t(217,217.0063182_wp), &
                                 Iso_t(218,218.0089735_wp), &
                                 Iso_t(219,219.013614_wp), &
                                 Iso_t(220,220.016386_wp), &
                                 Iso_t(221,221.021228_wp), &
                                 Iso_t(222,222.02414_wp), &
                                 Iso_t(223,223.02907_wp), &
                                 Iso_t(224,224.03211_wp), &
                                 Iso_t(225,225.03707_wp), &
                                 Iso_t(226,226.04031_wp), &
                                 Iso_t(227,227.04539_wp)]

  ElementList(85)%Symbol = adjustl(PTab(85)) ! At
  ElementList(85)%Natural = 0
  call mma_Allocate(ElementList(85)%Isotopes,39)
  ElementList(85)%Isotopes(:) = [ &
                                 Iso_t(210,209.9871479_wp), &
                                 Iso_t(191,191.004148_wp), &
                                 Iso_t(192,192.003152_wp), &
                                 Iso_t(193,192.999927_wp), &
                                 Iso_t(194,193.999236_wp), &
                                 Iso_t(195,194.9962685_wp), &
                                 Iso_t(196,195.9958_wp), &
                                 Iso_t(197,196.993189_wp), &
                                 Iso_t(198,197.992784_wp), &
                                 Iso_t(199,198.9905277_wp), &
                                 Iso_t(200,199.990351_wp), &
                                 Iso_t(201,200.9884171_wp), &
                                 Iso_t(202,201.98863_wp), &
                                 Iso_t(203,202.986943_wp), &
                                 Iso_t(204,203.987251_wp), &
                                 Iso_t(205,204.986076_wp), &
                                 Iso_t(206,205.986657_wp), &
                                 Iso_t(207,206.9858_wp), &
                                 Iso_t(208,207.9866133_wp), &
                                 Iso_t(209,208.9861702_wp), &
                                 Iso_t(211,210.9874966_wp), &
                                 Iso_t(212,211.9907377_wp), &
                                 Iso_t(213,212.992937_wp), &
                                 Iso_t(214,213.9963721_wp), &
                                 Iso_t(215,214.9986528_wp), &
                                 Iso_t(216,216.0024236_wp), &
                                 Iso_t(217,217.0047192_wp), &
                                 Iso_t(218,218.008695_wp), &
                                 Iso_t(219,219.0111618_wp), &
                                 Iso_t(220,220.015433_wp), &
                                 Iso_t(221,221.018017_wp), &
                                 Iso_t(222,222.022494_wp), &
                                 Iso_t(223,223.025151_wp), &
                                 Iso_t(224,224.029749_wp), &
                                 Iso_t(225,225.03263_wp), &
                                 Iso_t(226,226.03716_wp), &
                                 Iso_t(227,227.04024_wp), &
                                 Iso_t(228,228.04475_wp), &
                                 Iso_t(229,229.04812_wp)]

  ElementList(86)%Symbol = adjustl(PTab(86)) ! Rn
  ElementList(86)%Natural = 0
  call mma_Allocate(ElementList(86)%Isotopes,39)
  ElementList(86)%Isotopes(:) = [ &
                                 Iso_t(222,222.0175782_wp), &
                                 Iso_t(193,193.009708_wp), &
                                 Iso_t(194,194.006144_wp), &
                                 Iso_t(195,195.005422_wp), &
                                 Iso_t(196,196.002116_wp), &
                                 Iso_t(197,197.001585_wp), &
                                 Iso_t(198,197.998679_wp), &
                                 Iso_t(199,198.99839_wp), &
                                 Iso_t(200,199.99569_wp), &
                                 Iso_t(201,200.995628_wp), &
                                 Iso_t(202,201.993264_wp), &
                                 Iso_t(203,202.993388_wp), &
                                 Iso_t(204,203.99143_wp), &
                                 Iso_t(205,204.991719_wp), &
                                 Iso_t(206,205.990214_wp), &
                                 Iso_t(207,206.9907303_wp), &
                                 Iso_t(208,207.989635_wp), &
                                 Iso_t(209,208.990415_wp), &
                                 Iso_t(210,209.9896891_wp), &
                                 Iso_t(211,210.9906011_wp), &
                                 Iso_t(212,211.9907039_wp), &
                                 Iso_t(213,212.9938831_wp), &
                                 Iso_t(214,213.995363_wp), &
                                 Iso_t(215,214.9987459_wp), &
                                 Iso_t(216,216.0002719_wp), &
                                 Iso_t(217,217.003928_wp), &
                                 Iso_t(218,218.0056016_wp), &
                                 Iso_t(219,219.0094804_wp), &
                                 Iso_t(220,220.0113941_wp), &
                                 Iso_t(221,221.0155371_wp), &
                                 Iso_t(223,223.0218893_wp), &
                                 Iso_t(224,224.024096_wp), &
                                 Iso_t(225,225.028486_wp), &
                                 Iso_t(226,226.030861_wp), &
                                 Iso_t(227,227.035304_wp), &
                                 Iso_t(228,228.037835_wp), &
                                 Iso_t(229,229.042257_wp), &
                                 Iso_t(230,230.04514_wp), &
                                 Iso_t(231,231.04987_wp)]

  ElementList(87)%Symbol = adjustl(PTab(87)) ! Fr
  ElementList(87)%Natural = 0
  call mma_Allocate(ElementList(87)%Isotopes,35)
  ElementList(87)%Isotopes(:) = [ &
                                 Iso_t(223,223.019736_wp), &
                                 Iso_t(199,199.007259_wp), &
                                 Iso_t(200,200.006586_wp), &
                                 Iso_t(201,201.003867_wp), &
                                 Iso_t(202,202.00332_wp), &
                                 Iso_t(203,203.0009407_wp), &
                                 Iso_t(204,204.000652_wp), &
                                 Iso_t(205,204.9985939_wp), &
                                 Iso_t(206,205.998666_wp), &
                                 Iso_t(207,206.996946_wp), &
                                 Iso_t(208,207.997138_wp), &
                                 Iso_t(209,208.995955_wp), &
                                 Iso_t(210,209.996422_wp), &
                                 Iso_t(211,210.995556_wp), &
                                 Iso_t(212,211.9962257_wp), &
                                 Iso_t(213,212.996186_wp), &
                                 Iso_t(214,213.9989713_wp), &
                                 Iso_t(215,215.0003418_wp), &
                                 Iso_t(216,216.0031899_wp), &
                                 Iso_t(217,217.0046323_wp), &
                                 Iso_t(218,218.0075787_wp), &
                                 Iso_t(219,219.0092524_wp), &
                                 Iso_t(220,220.0123277_wp), &
                                 Iso_t(221,221.0142552_wp), &
                                 Iso_t(222,222.017552_wp), &
                                 Iso_t(224,224.023398_wp), &
                                 Iso_t(225,225.025573_wp), &
                                 Iso_t(226,226.029566_wp), &
                                 Iso_t(227,227.031869_wp), &
                                 Iso_t(228,228.035823_wp), &
                                 Iso_t(229,229.038298_wp), &
                                 Iso_t(230,230.042416_wp), &
                                 Iso_t(231,231.045158_wp), &
                                 Iso_t(232,232.04937_wp), &
                                 Iso_t(233,233.05264_wp)]

  ElementList(88)%Symbol = adjustl(PTab(88)) ! Ra
  ElementList(88)%Natural = 0
  call mma_Allocate(ElementList(88)%Isotopes,35)
  ElementList(88)%Isotopes(:) = [ &
                                 Iso_t(226,226.0254103_wp), &
                                 Iso_t(201,201.01271_wp), &
                                 Iso_t(202,202.00976_wp), &
                                 Iso_t(203,203.009304_wp), &
                                 Iso_t(204,204.006492_wp), &
                                 Iso_t(205,205.006268_wp), &
                                 Iso_t(206,206.003828_wp), &
                                 Iso_t(207,207.003799_wp), &
                                 Iso_t(208,208.001841_wp), &
                                 Iso_t(209,209.00199_wp), &
                                 Iso_t(210,210.000494_wp), &
                                 Iso_t(211,211.0008932_wp), &
                                 Iso_t(212,211.999787_wp), &
                                 Iso_t(213,213.000384_wp), &
                                 Iso_t(214,214.0000997_wp), &
                                 Iso_t(215,215.0027204_wp), &
                                 Iso_t(216,216.0035334_wp), &
                                 Iso_t(217,217.0063207_wp), &
                                 Iso_t(218,218.007141_wp), &
                                 Iso_t(219,219.0100855_wp), &
                                 Iso_t(220,220.0110259_wp), &
                                 Iso_t(221,221.0139177_wp), &
                                 Iso_t(222,222.0153748_wp), &
                                 Iso_t(223,223.0185023_wp), &
                                 Iso_t(224,224.020212_wp), &
                                 Iso_t(225,225.0236119_wp), &
                                 Iso_t(227,227.0291783_wp), &
                                 Iso_t(228,228.0310707_wp), &
                                 Iso_t(229,229.034942_wp), &
                                 Iso_t(230,230.037055_wp), &
                                 Iso_t(231,231.041027_wp), &
                                 Iso_t(232,232.0434753_wp), &
                                 Iso_t(233,233.047582_wp), &
                                 Iso_t(234,234.050342_wp), &
                                 Iso_t(235,235.05497_wp)]

  ElementList(89)%Symbol = adjustl(PTab(89)) ! Ac
  ElementList(89)%Natural = 0
  call mma_Allocate(ElementList(89)%Isotopes,32)
  ElementList(89)%Isotopes(:) = [ &
                                 Iso_t(227,227.0277523_wp), &
                                 Iso_t(206,206.014452_wp), &
                                 Iso_t(207,207.011966_wp), &
                                 Iso_t(208,208.01155_wp), &
                                 Iso_t(209,209.009495_wp), &
                                 Iso_t(210,210.009436_wp), &
                                 Iso_t(211,211.007732_wp), &
                                 Iso_t(212,212.007813_wp), &
                                 Iso_t(213,213.006609_wp), &
                                 Iso_t(214,214.006918_wp), &
                                 Iso_t(215,215.006475_wp), &
                                 Iso_t(216,216.008743_wp), &
                                 Iso_t(217,217.009344_wp), &
                                 Iso_t(218,218.011642_wp), &
                                 Iso_t(219,219.012421_wp), &
                                 Iso_t(220,220.0147549_wp), &
                                 Iso_t(221,221.015592_wp), &
                                 Iso_t(222,222.0178442_wp), &
                                 Iso_t(223,223.0191377_wp), &
                                 Iso_t(224,224.0217232_wp), &
                                 Iso_t(225,225.02323_wp), &
                                 Iso_t(226,226.0260984_wp), &
                                 Iso_t(228,228.0310215_wp), &
                                 Iso_t(229,229.032956_wp), &
                                 Iso_t(230,230.036327_wp), &
                                 Iso_t(231,231.038393_wp), &
                                 Iso_t(232,232.042034_wp), &
                                 Iso_t(233,233.044346_wp), &
                                 Iso_t(234,234.048139_wp), &
                                 Iso_t(235,235.05084_wp), &
                                 Iso_t(236,236.054988_wp), &
                                 Iso_t(237,237.05827_wp)]

  ElementList(90)%Symbol = adjustl(PTab(90)) ! Th
  ElementList(90)%Natural = 1
  call mma_Allocate(ElementList(90)%Isotopes,32)
  ElementList(90)%Isotopes(:) = [ &
                                 Iso_t(232,232.0380558_wp), &
                                 Iso_t(208,208.0179_wp), &
                                 Iso_t(209,209.017753_wp), &
                                 Iso_t(210,210.015094_wp), &
                                 Iso_t(211,211.014929_wp), &
                                 Iso_t(212,212.012988_wp), &
                                 Iso_t(213,213.013009_wp), &
                                 Iso_t(214,214.0115_wp), &
                                 Iso_t(215,215.0117248_wp), &
                                 Iso_t(216,216.011056_wp), &
                                 Iso_t(217,217.013117_wp), &
                                 Iso_t(218,218.013276_wp), &
                                 Iso_t(219,219.015537_wp), &
                                 Iso_t(220,220.015748_wp), &
                                 Iso_t(221,221.018184_wp), &
                                 Iso_t(222,222.018469_wp), &
                                 Iso_t(223,223.0208119_wp), &
                                 Iso_t(224,224.021464_wp), &
                                 Iso_t(225,225.0239514_wp), &
                                 Iso_t(226,226.0249034_wp), &
                                 Iso_t(227,227.0277042_wp), &
                                 Iso_t(228,228.0287413_wp), &
                                 Iso_t(229,229.0317627_wp), &
                                 Iso_t(230,230.0331341_wp), &
                                 Iso_t(231,231.0363046_wp), &
                                 Iso_t(233,233.0415823_wp), &
                                 Iso_t(234,234.0436014_wp), &
                                 Iso_t(235,235.047255_wp), &
                                 Iso_t(236,236.049657_wp), &
                                 Iso_t(237,237.053629_wp), &
                                 Iso_t(238,238.0565_wp), &
                                 Iso_t(239,239.06077_wp)]

  ElementList(91)%Symbol = adjustl(PTab(91)) ! Pa
  ElementList(91)%Natural = 1
  call mma_Allocate(ElementList(91)%Isotopes,30)
  ElementList(91)%Isotopes(:) = [ &
                                 Iso_t(231,231.0358842_wp), &
                                 Iso_t(212,212.023203_wp), &
                                 Iso_t(213,213.021109_wp), &
                                 Iso_t(214,214.020918_wp), &
                                 Iso_t(215,215.019183_wp), &
                                 Iso_t(216,216.019109_wp), &
                                 Iso_t(217,217.018325_wp), &
                                 Iso_t(218,218.020059_wp), &
                                 Iso_t(219,219.019904_wp), &
                                 Iso_t(220,220.021705_wp), &
                                 Iso_t(221,221.021875_wp), &
                                 Iso_t(222,222.023784_wp), &
                                 Iso_t(223,223.023963_wp), &
                                 Iso_t(224,224.0256176_wp), &
                                 Iso_t(225,225.026131_wp), &
                                 Iso_t(226,226.027948_wp), &
                                 Iso_t(227,227.0288054_wp), &
                                 Iso_t(228,228.0310517_wp), &
                                 Iso_t(229,229.0320972_wp), &
                                 Iso_t(230,230.034541_wp), &
                                 Iso_t(232,232.0385917_wp), &
                                 Iso_t(233,233.0402472_wp), &
                                 Iso_t(234,234.0433072_wp), &
                                 Iso_t(235,235.045399_wp), &
                                 Iso_t(236,236.048668_wp), &
                                 Iso_t(237,237.051023_wp), &
                                 Iso_t(238,238.054637_wp), &
                                 Iso_t(239,239.05726_wp), &
                                 Iso_t(240,240.06098_wp), &
                                 Iso_t(241,241.06408_wp)]

  ElementList(92)%Symbol = adjustl(PTab(92)) ! U
  ElementList(92)%Natural = 3
  call mma_Allocate(ElementList(92)%Isotopes,27)
  ElementList(92)%Isotopes(:) = [ &
                                 Iso_t(238,238.0507884_wp), &
                                 Iso_t(235,235.0439301_wp), &
                                 Iso_t(234,234.0409523_wp), &
                                 Iso_t(217,217.02466_wp), &
                                 Iso_t(218,218.023523_wp), &
                                 Iso_t(219,219.024999_wp), &
                                 Iso_t(220,220.02462_wp), &
                                 Iso_t(221,221.02628_wp), &
                                 Iso_t(222,222.026_wp), &
                                 Iso_t(223,223.027739_wp), &
                                 Iso_t(224,224.027605_wp), &
                                 Iso_t(225,225.029391_wp), &
                                 Iso_t(226,226.029339_wp), &
                                 Iso_t(227,227.031157_wp), &
                                 Iso_t(228,228.031371_wp), &
                                 Iso_t(229,229.0335063_wp), &
                                 Iso_t(230,230.0339401_wp), &
                                 Iso_t(231,231.0362939_wp), &
                                 Iso_t(232,232.0371563_wp), &
                                 Iso_t(233,233.0396355_wp), &
                                 Iso_t(236,236.0455682_wp), &
                                 Iso_t(237,237.0487304_wp), &
                                 Iso_t(239,239.0542935_wp), &
                                 Iso_t(240,240.0565934_wp), &
                                 Iso_t(241,241.06033_wp), &
                                 Iso_t(242,242.06293_wp), &
                                 Iso_t(243,243.06699_wp)]

  ElementList(93)%Symbol = adjustl(PTab(93)) ! Np
  ElementList(93)%Natural = 0
  call mma_Allocate(ElementList(93)%Isotopes,27)
  ElementList(93)%Isotopes(:) = [ &
                                 Iso_t(237,237.0481736_wp), &
                                 Iso_t(219,219.03143_wp), &
                                 Iso_t(220,220.03254_wp), &
                                 Iso_t(221,221.03204_wp), &
                                 Iso_t(222,222.0333_wp), &
                                 Iso_t(223,223.03285_wp), &
                                 Iso_t(224,224.03422_wp), &
                                 Iso_t(225,225.033911_wp), &
                                 Iso_t(226,226.035188_wp), &
                                 Iso_t(227,227.034957_wp), &
                                 Iso_t(228,228.036067_wp), &
                                 Iso_t(229,229.036264_wp), &
                                 Iso_t(230,230.037828_wp), &
                                 Iso_t(231,231.038245_wp), &
                                 Iso_t(232,232.04011_wp), &
                                 Iso_t(233,233.040741_wp), &
                                 Iso_t(234,234.0428953_wp), &
                                 Iso_t(235,235.0440635_wp), &
                                 Iso_t(236,236.04657_wp), &
                                 Iso_t(238,238.0509466_wp), &
                                 Iso_t(239,239.0529392_wp), &
                                 Iso_t(240,240.056165_wp), &
                                 Iso_t(241,241.058253_wp), &
                                 Iso_t(242,242.06164_wp), &
                                 Iso_t(243,243.06428_wp), &
                                 Iso_t(244,244.06785_wp), &
                                 Iso_t(245,245.0708_wp)]

  ElementList(94)%Symbol = adjustl(PTab(94)) ! Pu
  ElementList(94)%Natural = 0
  call mma_Allocate(ElementList(94)%Isotopes,20)
  ElementList(94)%Isotopes(:) = [ &
                                 Iso_t(244,244.0642053_wp), &
                                 Iso_t(228,228.038732_wp), &
                                 Iso_t(229,229.040144_wp), &
                                 Iso_t(230,230.03965_wp), &
                                 Iso_t(231,231.041102_wp), &
                                 Iso_t(232,232.041185_wp), &
                                 Iso_t(233,233.042998_wp), &
                                 Iso_t(234,234.0433174_wp), &
                                 Iso_t(235,235.045286_wp), &
                                 Iso_t(236,236.0460581_wp), &
                                 Iso_t(237,237.0484098_wp), &
                                 Iso_t(238,238.0495601_wp), &
                                 Iso_t(239,239.0521636_wp), &
                                 Iso_t(240,240.0538138_wp), &
                                 Iso_t(241,241.0568517_wp), &
                                 Iso_t(242,242.0587428_wp), &
                                 Iso_t(243,243.0620036_wp), &
                                 Iso_t(245,245.067826_wp), &
                                 Iso_t(246,246.070205_wp), &
                                 Iso_t(247,247.07419_wp)]

  ElementList(95)%Symbol = adjustl(PTab(95)) ! Am
  ElementList(95)%Natural = 0
  call mma_Allocate(ElementList(95)%Isotopes,20)
  ElementList(95)%Isotopes(:) = [ &
                                 Iso_t(243,243.0613813_wp), &
                                 Iso_t(230,230.04609_wp), &
                                 Iso_t(231,231.04556_wp), &
                                 Iso_t(232,232.04645_wp), &
                                 Iso_t(233,233.04644_wp), &
                                 Iso_t(234,234.04773_wp), &
                                 Iso_t(235,235.047908_wp), &
                                 Iso_t(236,236.04943_wp), &
                                 Iso_t(237,237.049996_wp), &
                                 Iso_t(238,238.051985_wp), &
                                 Iso_t(239,239.0530247_wp), &
                                 Iso_t(240,240.0553_wp), &
                                 Iso_t(241,241.0568293_wp), &
                                 Iso_t(242,242.0595494_wp), &
                                 Iso_t(244,244.0642851_wp), &
                                 Iso_t(245,245.0664548_wp), &
                                 Iso_t(246,246.069775_wp), &
                                 Iso_t(247,247.07209_wp), &
                                 Iso_t(248,248.07575_wp), &
                                 Iso_t(249,249.07848_wp)]

  ElementList(96)%Symbol = adjustl(PTab(96)) ! Cm
  ElementList(96)%Natural = 0
  call mma_Allocate(ElementList(96)%Isotopes,21)
  ElementList(96)%Isotopes(:) = [ &
                                 Iso_t(247,247.0703541_wp), &
                                 Iso_t(232,232.04982_wp), &
                                 Iso_t(233,233.05077_wp), &
                                 Iso_t(234,234.05016_wp), &
                                 Iso_t(235,235.05154_wp), &
                                 Iso_t(236,236.051374_wp), &
                                 Iso_t(237,237.052869_wp), &
                                 Iso_t(238,238.053081_wp), &
                                 Iso_t(239,239.05491_wp), &
                                 Iso_t(240,240.0555297_wp), &
                                 Iso_t(241,241.0576532_wp), &
                                 Iso_t(242,242.058836_wp), &
                                 Iso_t(243,243.0613893_wp), &
                                 Iso_t(244,244.0627528_wp), &
                                 Iso_t(245,245.0654915_wp), &
                                 Iso_t(246,246.0672238_wp), &
                                 Iso_t(248,248.0723499_wp), &
                                 Iso_t(249,249.0759548_wp), &
                                 Iso_t(250,250.078358_wp), &
                                 Iso_t(251,251.082286_wp), &
                                 Iso_t(252,252.08487_wp)]

  ElementList(97)%Symbol = adjustl(PTab(97)) ! Bk
  ElementList(97)%Natural = 0
  call mma_Allocate(ElementList(97)%Isotopes,21)
  ElementList(97)%Isotopes(:) = [ &
                                 Iso_t(247,247.0703073_wp), &
                                 Iso_t(234,234.05727_wp), &
                                 Iso_t(235,235.05658_wp), &
                                 Iso_t(236,236.05748_wp), &
                                 Iso_t(237,237.0571_wp), &
                                 Iso_t(238,238.0582_wp), &
                                 Iso_t(239,239.05824_wp), &
                                 Iso_t(240,240.05976_wp), &
                                 Iso_t(241,241.06016_wp), &
                                 Iso_t(242,242.06198_wp), &
                                 Iso_t(243,243.0630078_wp), &
                                 Iso_t(244,244.065181_wp), &
                                 Iso_t(245,245.0663618_wp), &
                                 Iso_t(246,246.068673_wp), &
                                 Iso_t(248,248.073088_wp), &
                                 Iso_t(249,249.0749877_wp), &
                                 Iso_t(250,250.0783167_wp), &
                                 Iso_t(251,251.080762_wp), &
                                 Iso_t(252,252.08431_wp), &
                                 Iso_t(253,253.08688_wp), &
                                 Iso_t(254,254.0906_wp)]

  ElementList(98)%Symbol = adjustl(PTab(98)) ! Cf
  ElementList(98)%Natural = 0
  call mma_Allocate(ElementList(98)%Isotopes,20)
  ElementList(98)%Isotopes(:) = [ &
                                 Iso_t(251,251.0795886_wp), &
                                 Iso_t(237,237.062198_wp), &
                                 Iso_t(238,238.06149_wp), &
                                 Iso_t(239,239.06253_wp), &
                                 Iso_t(240,240.062256_wp), &
                                 Iso_t(241,241.06369_wp), &
                                 Iso_t(242,242.063754_wp), &
                                 Iso_t(243,243.06548_wp), &
                                 Iso_t(244,244.0660008_wp), &
                                 Iso_t(245,245.0680487_wp), &
                                 Iso_t(246,246.0688055_wp), &
                                 Iso_t(247,247.070965_wp), &
                                 Iso_t(248,248.0721851_wp), &
                                 Iso_t(249,249.0748539_wp), &
                                 Iso_t(250,250.0764062_wp), &
                                 Iso_t(252,252.0816272_wp), &
                                 Iso_t(253,253.0851345_wp), &
                                 Iso_t(254,254.087324_wp), &
                                 Iso_t(255,255.09105_wp), &
                                 Iso_t(256,256.09344_wp)]

  ElementList(99)%Symbol = adjustl(PTab(99)) ! Es
  ElementList(99)%Natural = 0
  call mma_Allocate(ElementList(99)%Isotopes,20)
  ElementList(99)%Isotopes(:) = [ &
                                 Iso_t(252,252.08298_wp), &
                                 Iso_t(239,239.06823_wp), &
                                 Iso_t(240,240.06892_wp), &
                                 Iso_t(241,241.06856_wp), &
                                 Iso_t(242,242.06957_wp), &
                                 Iso_t(243,243.06951_wp), &
                                 Iso_t(244,244.07088_wp), &
                                 Iso_t(245,245.07125_wp), &
                                 Iso_t(246,246.0729_wp), &
                                 Iso_t(247,247.073622_wp), &
                                 Iso_t(248,248.075471_wp), &
                                 Iso_t(249,249.076411_wp), &
                                 Iso_t(250,250.07861_wp), &
                                 Iso_t(251,251.0799936_wp), &
                                 Iso_t(253,253.0848257_wp), &
                                 Iso_t(254,254.0880222_wp), &
                                 Iso_t(255,255.090275_wp), &
                                 Iso_t(256,256.0936_wp), &
                                 Iso_t(257,257.09598_wp), &
                                 Iso_t(258,258.09952_wp)]

  ElementList(100)%Symbol = adjustl(PTab(100)) ! Fm
  ElementList(100)%Natural = 0
  call mma_Allocate(ElementList(100)%Isotopes,20)
  ElementList(100)%Isotopes(:) = [ &
                                  Iso_t(257,257.0951061_wp), &
                                  Iso_t(241,241.07421_wp), &
                                  Iso_t(242,242.07343_wp), &
                                  Iso_t(243,243.07446_wp), &
                                  Iso_t(244,244.07404_wp), &
                                  Iso_t(245,245.07535_wp), &
                                  Iso_t(246,246.07535_wp), &
                                  Iso_t(247,247.07694_wp), &
                                  Iso_t(248,248.0771865_wp), &
                                  Iso_t(249,249.0789275_wp), &
                                  Iso_t(250,250.079521_wp), &
                                  Iso_t(251,251.08154_wp), &
                                  Iso_t(252,252.0824671_wp), &
                                  Iso_t(253,253.0851846_wp), &
                                  Iso_t(254,254.0868544_wp), &
                                  Iso_t(255,255.089964_wp), &
                                  Iso_t(256,256.0917745_wp), &
                                  Iso_t(258,258.09708_wp), &
                                  Iso_t(259,259.1006_wp), &
                                  Iso_t(260,260.10281_wp)]

  ElementList(101)%Symbol = adjustl(PTab(101)) ! Md
  ElementList(101)%Natural = 0
  call mma_Allocate(ElementList(101)%Isotopes,18)
  ElementList(101)%Isotopes(:) = [ &
                                  Iso_t(258,258.0984315_wp), &
                                  Iso_t(245,245.08081_wp), &
                                  Iso_t(246,246.08171_wp), &
                                  Iso_t(247,247.08152_wp), &
                                  Iso_t(248,248.08282_wp), &
                                  Iso_t(249,249.08291_wp), &
                                  Iso_t(250,250.08441_wp), &
                                  Iso_t(251,251.084774_wp), &
                                  Iso_t(252,252.08643_wp), &
                                  Iso_t(253,253.087144_wp), &
                                  Iso_t(254,254.08959_wp), &
                                  Iso_t(255,255.0910841_wp), &
                                  Iso_t(256,256.09389_wp), &
                                  Iso_t(257,257.0955424_wp), &
                                  Iso_t(259,259.10051_wp), &
                                  Iso_t(260,260.10365_wp), &
                                  Iso_t(261,261.10583_wp), &
                                  Iso_t(262,262.1091_wp)]

  ElementList(102)%Symbol = adjustl(PTab(102)) ! No
  ElementList(102)%Natural = 0
  call mma_Allocate(ElementList(102)%Isotopes,17)
  ElementList(102)%Isotopes(:) = [ &
                                  Iso_t(259,259.10103_wp), &
                                  Iso_t(248,248.08655_wp), &
                                  Iso_t(249,249.0878_wp), &
                                  Iso_t(250,250.08756_wp), &
                                  Iso_t(251,251.08894_wp), &
                                  Iso_t(252,252.088967_wp), &
                                  Iso_t(253,253.0905641_wp), &
                                  Iso_t(254,254.090956_wp), &
                                  Iso_t(255,255.093191_wp), &
                                  Iso_t(256,256.0942829_wp), &
                                  Iso_t(257,257.0968878_wp), &
                                  Iso_t(258,258.09821_wp), &
                                  Iso_t(260,260.10264_wp), &
                                  Iso_t(261,261.1057_wp), &
                                  Iso_t(262,262.10746_wp), &
                                  Iso_t(263,263.11071_wp), &
                                  Iso_t(264,264.11273_wp)]

  ElementList(103)%Symbol = adjustl(PTab(103)) ! Lr
  ElementList(103)%Natural = 0
  call mma_Allocate(ElementList(103)%Isotopes,16)
  ElementList(103)%Isotopes(:) = [ &
                                  Iso_t(262,262.10961_wp), &
                                  Iso_t(251,251.09418_wp), &
                                  Iso_t(252,252.09526_wp), &
                                  Iso_t(253,253.09509_wp), &
                                  Iso_t(254,254.09648_wp), &
                                  Iso_t(255,255.096562_wp), &
                                  Iso_t(256,256.098494_wp), &
                                  Iso_t(257,257.099418_wp), &
                                  Iso_t(258,258.10176_wp), &
                                  Iso_t(259,259.102902_wp), &
                                  Iso_t(260,260.1055_wp), &
                                  Iso_t(261,261.10688_wp), &
                                  Iso_t(263,263.11136_wp), &
                                  Iso_t(264,264.1142_wp), &
                                  Iso_t(265,265.11619_wp), &
                                  Iso_t(266,266.11983_wp)]

  ElementList(104)%Symbol = adjustl(PTab(104)) ! Rf
  ElementList(104)%Natural = 0
  call mma_Allocate(ElementList(104)%Isotopes,16)
  ElementList(104)%Isotopes(:) = [ &
                                  Iso_t(267,267.12179_wp), &
                                  Iso_t(253,253.10044_wp), &
                                  Iso_t(254,254.10005_wp), &
                                  Iso_t(255,255.10127_wp), &
                                  Iso_t(256,256.101152_wp), &
                                  Iso_t(257,257.102918_wp), &
                                  Iso_t(258,258.103428_wp), &
                                  Iso_t(259,259.105596_wp), &
                                  Iso_t(260,260.10644_wp), &
                                  Iso_t(261,261.108773_wp), &
                                  Iso_t(262,262.10992_wp), &
                                  Iso_t(263,263.11249_wp), &
                                  Iso_t(264,264.11388_wp), &
                                  Iso_t(265,265.11668_wp), &
                                  Iso_t(266,266.11817_wp), &
                                  Iso_t(268,268.12397_wp)]

  ElementList(105)%Symbol = adjustl(PTab(105)) ! Db
  ElementList(105)%Natural = 0
  call mma_Allocate(ElementList(105)%Isotopes,16)
  ElementList(105)%Isotopes(:) = [ &
                                  Iso_t(268,268.12567_wp), &
                                  Iso_t(255,255.10707_wp), &
                                  Iso_t(256,256.10789_wp), &
                                  Iso_t(257,257.10758_wp), &
                                  Iso_t(258,258.10928_wp), &
                                  Iso_t(259,259.109492_wp), &
                                  Iso_t(260,260.1113_wp), &
                                  Iso_t(261,261.11192_wp), &
                                  Iso_t(262,262.11407_wp), &
                                  Iso_t(263,263.11499_wp), &
                                  Iso_t(264,264.11741_wp), &
                                  Iso_t(265,265.11861_wp), &
                                  Iso_t(266,266.12103_wp), &
                                  Iso_t(267,267.12247_wp), &
                                  Iso_t(269,269.12791_wp), &
                                  Iso_t(270,270.13136_wp)]

  ElementList(106)%Symbol = adjustl(PTab(106)) ! Sg
  ElementList(106)%Natural = 0
  call mma_Allocate(ElementList(106)%Isotopes,16)
  ElementList(106)%Isotopes(:) = [ &
                                  Iso_t(269,269.12863_wp), &
                                  Iso_t(258,258.11298_wp), &
                                  Iso_t(259,259.1144_wp), &
                                  Iso_t(260,260.114384_wp), &
                                  Iso_t(261,261.115949_wp), &
                                  Iso_t(262,262.116337_wp), &
                                  Iso_t(263,263.11829_wp), &
                                  Iso_t(264,264.11893_wp), &
                                  Iso_t(265,265.12109_wp), &
                                  Iso_t(266,266.12198_wp), &
                                  Iso_t(267,267.12436_wp), &
                                  Iso_t(268,268.12539_wp), &
                                  Iso_t(270,270.13043_wp), &
                                  Iso_t(271,271.13393_wp), &
                                  Iso_t(272,272.13589_wp), &
                                  Iso_t(273,273.13958_wp)]

  ElementList(107)%Symbol = adjustl(PTab(107)) ! Bh
  ElementList(107)%Natural = 0
  call mma_Allocate(ElementList(107)%Isotopes,16)
  ElementList(107)%Isotopes(:) = [ &
                                  Iso_t(270,270.13336_wp), &
                                  Iso_t(260,260.12166_wp), &
                                  Iso_t(261,261.12145_wp), &
                                  Iso_t(262,262.12297_wp), &
                                  Iso_t(263,263.12292_wp), &
                                  Iso_t(264,264.12459_wp), &
                                  Iso_t(265,265.12491_wp), &
                                  Iso_t(266,266.12679_wp), &
                                  Iso_t(267,267.1275_wp), &
                                  Iso_t(268,268.12969_wp), &
                                  Iso_t(269,269.13042_wp), &
                                  Iso_t(271,271.13526_wp), &
                                  Iso_t(272,272.13826_wp), &
                                  Iso_t(273,273.14024_wp), &
                                  Iso_t(274,274.14355_wp), &
                                  Iso_t(275,275.14567_wp)]

  ElementList(108)%Symbol = adjustl(PTab(108)) ! Hs
  ElementList(108)%Natural = 0
  call mma_Allocate(ElementList(108)%Isotopes,15)
  ElementList(108)%Isotopes(:) = [ &
                                  Iso_t(269,269.13375_wp), &
                                  Iso_t(263,263.12852_wp), &
                                  Iso_t(264,264.128357_wp), &
                                  Iso_t(265,265.129793_wp), &
                                  Iso_t(266,266.130046_wp), &
                                  Iso_t(267,267.13167_wp), &
                                  Iso_t(268,268.13186_wp), &
                                  Iso_t(270,270.13429_wp), &
                                  Iso_t(271,271.13717_wp), &
                                  Iso_t(272,272.1385_wp), &
                                  Iso_t(273,273.14168_wp), &
                                  Iso_t(274,274.1433_wp), &
                                  Iso_t(275,275.14667_wp), &
                                  Iso_t(276,276.14846_wp), &
                                  Iso_t(277,277.1519_wp)]

  ElementList(109)%Symbol = adjustl(PTab(109)) ! Mt
  ElementList(109)%Natural = 0
  call mma_Allocate(ElementList(109)%Isotopes,15)
  ElementList(109)%Isotopes(:) = [ &
                                  Iso_t(278,278.15631_wp), &
                                  Iso_t(265,265.136_wp), &
                                  Iso_t(266,266.13737_wp), &
                                  Iso_t(267,267.13719_wp), &
                                  Iso_t(268,268.13865_wp), &
                                  Iso_t(269,269.13882_wp), &
                                  Iso_t(270,270.14033_wp), &
                                  Iso_t(271,271.14074_wp), &
                                  Iso_t(272,272.14341_wp), &
                                  Iso_t(273,273.1444_wp), &
                                  Iso_t(274,274.14724_wp), &
                                  Iso_t(275,275.14882_wp), &
                                  Iso_t(276,276.15159_wp), &
                                  Iso_t(277,277.15327_wp), &
                                  Iso_t(279,279.15808_wp)]

  ElementList(110)%Symbol = adjustl(PTab(110)) ! Ds
  ElementList(110)%Natural = 0
  call mma_Allocate(ElementList(110)%Isotopes,15)
  ElementList(110)%Isotopes(:) = [ &
                                  Iso_t(281,281.16451_wp), &
                                  Iso_t(267,267.14377_wp), &
                                  Iso_t(268,268.14348_wp), &
                                  Iso_t(269,269.144752_wp), &
                                  Iso_t(270,270.144584_wp), &
                                  Iso_t(271,271.14595_wp), &
                                  Iso_t(272,272.14602_wp), &
                                  Iso_t(273,273.14856_wp), &
                                  Iso_t(274,274.14941_wp), &
                                  Iso_t(275,275.15203_wp), &
                                  Iso_t(276,276.15303_wp), &
                                  Iso_t(277,277.15591_wp), &
                                  Iso_t(278,278.15704_wp), &
                                  Iso_t(279,279.1601_wp), &
                                  Iso_t(280,280.16131_wp)]

  ElementList(111)%Symbol = adjustl(PTab(111)) ! Rg
  ElementList(111)%Natural = 0
  call mma_Allocate(ElementList(111)%Isotopes,12)
  ElementList(111)%Isotopes(:) = [ &
                                  Iso_t(281,281.16636_wp), &
                                  Iso_t(272,272.15327_wp), &
                                  Iso_t(273,273.15313_wp), &
                                  Iso_t(274,274.15525_wp), &
                                  Iso_t(275,275.15594_wp), &
                                  Iso_t(276,276.15833_wp), &
                                  Iso_t(277,277.15907_wp), &
                                  Iso_t(278,278.16149_wp), &
                                  Iso_t(279,279.16272_wp), &
                                  Iso_t(280,280.16514_wp), &
                                  Iso_t(282,282.16912_wp), &
                                  Iso_t(283,283.17054_wp)]

  ElementList(112)%Symbol = adjustl(PTab(112)) ! Cn
  ElementList(112)%Natural = 0
  call mma_Allocate(ElementList(112)%Isotopes,10)
  ElementList(112)%Isotopes(:) = [ &
                                  Iso_t(283,283.17327_wp), &
                                  Iso_t(276,276.16141_wp), &
                                  Iso_t(277,277.16364_wp), &
                                  Iso_t(278,278.16416_wp), &
                                  Iso_t(279,279.16654_wp), &
                                  Iso_t(280,280.16715_wp), &
                                  Iso_t(281,281.16975_wp), &
                                  Iso_t(282,282.1705_wp), &
                                  Iso_t(284,284.17416_wp), &
                                  Iso_t(285,285.17712_wp)]

  ElementList(113)%Symbol = adjustl(PTab(113)) ! Nh
  ElementList(113)%Natural = 0
  call mma_Allocate(ElementList(113)%Isotopes,10)
  ElementList(113)%Isotopes(:) = [ &
                                  Iso_t(287,287.18339_wp), &
                                  Iso_t(278,278.17058_wp), &
                                  Iso_t(279,279.17095_wp), &
                                  Iso_t(280,280.17293_wp), &
                                  Iso_t(281,281.17348_wp), &
                                  Iso_t(282,282.17567_wp), &
                                  Iso_t(283,283.17657_wp), &
                                  Iso_t(284,284.17873_wp), &
                                  Iso_t(285,285.17973_wp), &
                                  Iso_t(286,286.18221_wp)]

  ElementList(114)%Symbol = adjustl(PTab(114)) ! Fl
  ElementList(114)%Natural = 0
  call mma_Allocate(ElementList(114)%Isotopes,5)
  ElementList(114)%Isotopes(:) = [ &
                                  Iso_t(289,289.19042_wp), &
                                  Iso_t(285,285.18364_wp), &
                                  Iso_t(286,286.18423_wp), &
                                  Iso_t(287,287.18678_wp), &
                                  Iso_t(288,288.18757_wp)]

  ElementList(115)%Symbol = adjustl(PTab(115)) ! Mc
  ElementList(115)%Natural = 0
  call mma_Allocate(ElementList(115)%Isotopes,5)
  ElementList(115)%Isotopes(:) = [ &
                                  Iso_t(288,288.19274_wp), &
                                  Iso_t(287,287.1907_wp), &
                                  Iso_t(289,289.19363_wp), &
                                  Iso_t(290,290.19598_wp), &
                                  Iso_t(291,291.19707_wp)]

  ElementList(116)%Symbol = adjustl(PTab(116)) ! Lv
  ElementList(116)%Natural = 0
  call mma_Allocate(ElementList(116)%Isotopes,5)
  ElementList(116)%Isotopes(:) = [ &
                                  Iso_t(293,293.20449_wp), &
                                  Iso_t(289,289.19816_wp), &
                                  Iso_t(290,290.19864_wp), &
                                  Iso_t(291,291.20108_wp), &
                                  Iso_t(292,292.20174_wp)]

  ElementList(117)%Symbol = adjustl(PTab(117)) ! Ts
  ElementList(117)%Natural = 0
  call mma_Allocate(ElementList(117)%Isotopes,4)
  ElementList(117)%Isotopes(:) = [ &
                                  Iso_t(294,294.21046_wp), &
                                  Iso_t(291,291.20553_wp), &
                                  Iso_t(292,292.20746_wp), &
                                  Iso_t(293,293.20824_wp)]

  ElementList(118)%Symbol = adjustl(PTab(118)) ! Og
  ElementList(118)%Natural = 0
  call mma_Allocate(ElementList(118)%Isotopes,3)
  ElementList(118)%Isotopes(:) = [ &
                                  Iso_t(294,294.21392_wp), &
                                  Iso_t(293,293.21356_wp), &
                                  Iso_t(295,295.21624_wp)]

end subroutine Initialize_Isotopes

! This subroutine frees up the memory

subroutine Free_Isotopes()
  integer(kind=iwp) :: i
  if (.not. allocated(ElementList)) return
  do i=1,size(ElementList,1)
    call mma_Deallocate(ElementList(i)%Isotopes)
  end do
  call mma_Deallocate(ElementList)
end subroutine Free_Isotopes

! Subroutine(s) to get the Mass of the isotope IsNr belonging to the
! element Atom. If IsNr=0, the most abundant isotope (or the most
! stable if all are radioactive) is selected. The mass is returned
! in atomic units (m_e).
! Atom can be an atomic symbol or an atomic number.

subroutine Isotope_sym(IsNr,Atom,Mass)
  integer(kind=iwp), intent(inout) :: IsNr
  character(len=2), intent(in) :: Atom
  real(kind=wp), intent(out) :: Mass
  integer(kind=iwp) :: i, This
  character(len=2) :: Sym, Sym2

  call Initialize_Isotopes()

  Sym2 = adjustl(Atom)
  call UpCase(Sym2)
  if ((Sym2 == 'D') .or. (Sym2 == 'T')) Sym2 = 'H'
  This = 0
  do i=1,MaxAtomNum
    Sym = adjustl(ElementList(i)%Symbol)
    call UpCase(Sym)
    if (Sym == Sym2) then
      This = i
      exit
    end if
  end do

  if (This == 0) then
    write(u6,*) 'Isotope: Did not find atom!'
    write(u6,*) 'Atom=',Atom
    call Abend()
  end if

  if (IsNr == 0) IsNr = ElementList(This)%Isotopes(1)%A
  if (Sym2 == 'D') IsNr = 2
  if (Sym2 == 'T') IsNr = 3
  do i=1,size(ElementList(This)%Isotopes,1)
    if (ElementList(This)%Isotopes(i)%A == IsNr) then
      Mass = uToau*ElementList(This)%Isotopes(i)%m
      return
    end if
  end do

  write(u6,*) 'Isotope: Did not find isotope!'
  write(u6,*) 'IsNr=',IsNr
  write(u6,*) 'Atom=',Atom
  call Abend()

end subroutine Isotope_sym

subroutine Isotope_num(IsNr,Atom,Mass)
  integer(kind=iwp), intent(inout) :: IsNr
  integer(kind=iwp), intent(in) :: Atom
  real(kind=wp), intent(out) :: Mass
  integer(kind=iwp) :: i

  call Initialize_Isotopes()

  if ((Atom < 0) .or. (Atom > MaxAtomNum)) then
    write(u6,*) 'Isotope: Did not find atom!'
    write(u6,*) 'Atom=',Atom
    call Abend()
  end if

  if (IsNr == 0) IsNr = ElementList(Atom)%Isotopes(1)%A
  do i=1,size(ElementList(Atom)%Isotopes,1)
    if (ElementList(Atom)%Isotopes(i)%A == IsNr) then
      Mass = uToau*ElementList(Atom)%Isotopes(i)%m
      return
    end if
  end do

  write(u6,*) 'Isotope: Did not find isotope!'
  write(u6,*) 'IsNr=',IsNr
  write(u6,*) 'Atom=',Atom
  call Abend()

end subroutine Isotope_num

! Function that returns the mass in atomic units (m_e) of a particular
! nuclide with Z protons and A-Z neutrons. Returns -1.0 if the nuclide
! is unknown.

function NuclideMass(Z,A)
  use Constants, only: One
  real(kind=wp) :: NuclideMass
  integer(kind=iwp), intent(in) :: Z, A
  integer(kind=iwp) :: i

  call Initialize_Isotopes()

  NuclideMass = -One
  if ((Z < 1) .or. (Z > MaxAtomNum)) return
  do i=1,size(ElementList(Z)%Isotopes,1)
    if (ElementList(Z)%Isotopes(i)%A == A) then
      NuclideMass = uToau*ElementList(Z)%Isotopes(i)%m
      exit
    end if
  end do

end function NuclideMass

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define elm_cptr2loff, element_mma_allo_1D, element_mma_allo_1D_lim, element_mma_free_1D
#define _TYPE_ type(element_t)
#  define _FUNC_NAME_ elm_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ element_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'elm_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

! Define iso_cptr2loff, isotope_mma_allo_1D, isotope_mma_allo_1D_lim, isotope_mma_free_1D
#define _TYPE_ type(iso_t)
#  define _FUNC_NAME_ iso_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ isotope_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'iso_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module Isotopes
