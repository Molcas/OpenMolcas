************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************
*
* Isotope numbers and masses
*
* Taken from NIST: "Atomic Weights and Isotopic Compositions with
* Relative Atomic Masses" (https://www.nist.gov/pml/data/comp.cfm)
* Last updated: March 2017 (Version 4.1, July 2015)
*
* Each item in the the array ElementList contains data for an element:
*  - %Symbol: symbol
*  - %Nat: number of natural ocurring isotopes
*  - %Isotopes: array with all the isotopes for the element, sorted
*               by order of abundance (most to least, artificial
*               sorted by increasing mass number); elements with no
*               natural isotopes have the most stable first. Each item
*               in this array contains:
*    - %A: mass number (protons + neutrons)
*    - %m: isotopic mass in Da
*
* The "default" isotope for each element is simply the first item in
* the %Isotopes member.
*
* Manual changes from NIST data:
*  - Definitive symbols for all elements
*  - Most stable isotope (from Wikipedia) selected for Z > 94
*  - 3H and 14C included as natural
*
      Module Isotopes
      Use stdalloc, Only: mma_Allocate, mma_Deallocate
      Implicit None
      Private
      Type Iso
        Integer :: A
        Real*8 :: m
      End Type Iso
      Type Element
        Character(Len=2) :: Symbol
        Integer :: Z, Nat
        Type(Iso), Allocatable :: Isotopes(:)
      End Type Element
      Integer, Parameter :: MaxAtomNum=118
      Type(Element), Allocatable :: ElementList(:)
#include "constants2.fh"

      Interface Isotope
        Module Procedure Isotope_sym, Isotope_num
      End Interface Isotope

      Public :: MaxAtomNum, ElementList, Isotope, Initialize_Isotopes,
     &          Free_Isotopes, NuclideMass

*
* Private extensions to mma interfaces
*
      Interface cptr2loff
        Module Procedure elm_cptr2loff
        Module Procedure iso_cptr2loff
      End Interface
      Interface mma_Allocate
        Module Procedure element_mma_allo_1D, element_mma_allo_1D_lim
        Module Procedure isotope_mma_allo_1D, isotope_mma_allo_1D_lim
      End Interface
      Interface mma_Deallocate
        Module Procedure element_mma_free_1D
        Module Procedure isotope_mma_free_1D
      End Interface

      Contains

*
* This subroutine allocates and fills the data in ElementList
* Since each array has a different size, it has to be done dynamically
*
      Subroutine Initialize_Isotopes()
#ifdef _GARBLE_
      Interface
        Subroutine c_null_alloc(A)
          Import Iso
          Type(Iso), Allocatable :: A(:)
        End Subroutine c_null_alloc
      End Interface
      Integer :: i
#endif
      If (Allocated(ElementList)) Return
      Call mma_Allocate(ElementList,MaxAtomNum,'ElmList')
#ifdef _GARBLE_
*     Garbling corrupts the allocation status of allocatable
*     components, use a hack to reset it
      Do i=1,Size(ElementList,1)
        Call c_null_alloc(ElementList(i)%Isotopes)
      End Do
#endif

      ElementList(1)%Symbol = 'H'
      ElementList(1)%nat = 3
      Call mma_Allocate(ElementList(1)%Isotopes,7)
      ElementList(1)%Isotopes(:) = [
     &  Iso(1, 1.00782503223d0),
     &  Iso(2, 2.01410177812d0),
     &  Iso(3, 3.0160492779d0),
     &  Iso(4, 4.02643d0),
     &  Iso(5, 5.035311d0),
     &  Iso(6, 6.04496d0),
     &  Iso(7, 7.0527d0) ]

      ElementList(2)%Symbol = 'He'
      ElementList(2)%nat = 2
      Call mma_Allocate(ElementList(2)%Isotopes,8)
      ElementList(2)%Isotopes(:) = [
     &  Iso(4, 4.00260325413d0),
     &  Iso(3, 3.0160293201d0),
     &  Iso(5, 5.012057d0),
     &  Iso(6, 6.018885891d0),
     &  Iso(7, 7.0279907d0),
     &  Iso(8, 8.03393439d0),
     &  Iso(9, 9.043946d0),
     &  Iso(10, 10.05279d0) ]

      ElementList(3)%Symbol = 'Li'
      ElementList(3)%nat = 2
      Call mma_Allocate(ElementList(3)%Isotopes,11)
      ElementList(3)%Isotopes(:) = [
     &  Iso(7, 7.0160034366d0),
     &  Iso(6, 6.0151228874d0),
     &  Iso(3, 3.0308d0),
     &  Iso(4, 4.02719d0),
     &  Iso(5, 5.012538d0),
     &  Iso(8, 8.022486246d0),
     &  Iso(9, 9.02679019d0),
     &  Iso(10, 10.035483d0),
     &  Iso(11, 11.04372358d0),
     &  Iso(12, 12.052517d0),
     &  Iso(13, 13.06263d0) ]

      ElementList(4)%Symbol = 'Be'
      ElementList(4)%nat = 1
      Call mma_Allocate(ElementList(4)%Isotopes,12)
      ElementList(4)%Isotopes(:) = [
     &  Iso(9, 9.012183065d0),
     &  Iso(5, 5.0399d0),
     &  Iso(6, 6.0197264d0),
     &  Iso(7, 7.016928717d0),
     &  Iso(8, 8.005305102d0),
     &  Iso(10, 10.013534695d0),
     &  Iso(11, 11.02166108d0),
     &  Iso(12, 12.0269221d0),
     &  Iso(13, 13.036135d0),
     &  Iso(14, 14.04289d0),
     &  Iso(15, 15.05342d0),
     &  Iso(16, 16.06167d0) ]

      ElementList(5)%Symbol = 'B'
      ElementList(5)%nat = 2
      Call mma_Allocate(ElementList(5)%Isotopes,16)
      ElementList(5)%Isotopes(:) = [
     &  Iso(11, 11.00930536d0),
     &  Iso(10, 10.01293695d0),
     &  Iso(6, 6.0508d0),
     &  Iso(7, 7.029712d0),
     &  Iso(8, 8.0246073d0),
     &  Iso(9, 9.01332965d0),
     &  Iso(12, 12.0143527d0),
     &  Iso(13, 13.0177802d0),
     &  Iso(14, 14.025404d0),
     &  Iso(15, 15.031088d0),
     &  Iso(16, 16.039842d0),
     &  Iso(17, 17.04699d0),
     &  Iso(18, 18.05566d0),
     &  Iso(19, 19.0631d0),
     &  Iso(20, 20.07207d0),
     &  Iso(21, 21.08129d0) ]

      ElementList(6)%Symbol = 'C'
      ElementList(6)%nat = 3
      Call mma_Allocate(ElementList(6)%Isotopes,16)
      ElementList(6)%Isotopes(:) = [
     &  Iso(12, 12.0d0),
     &  Iso(13, 13.00335483507d0),
     &  Iso(14, 14.0032419884d0),
     &  Iso(8, 8.037643d0),
     &  Iso(9, 9.0310372d0),
     &  Iso(10, 10.01685331d0),
     &  Iso(11, 11.0114336d0),
     &  Iso(15, 15.01059926d0),
     &  Iso(16, 16.0147013d0),
     &  Iso(17, 17.022577d0),
     &  Iso(18, 18.026751d0),
     &  Iso(19, 19.0348d0),
     &  Iso(20, 20.04032d0),
     &  Iso(21, 21.049d0),
     &  Iso(22, 22.05753d0),
     &  Iso(23, 23.0689d0) ]

      ElementList(7)%Symbol = 'N'
      ElementList(7)%nat = 2
      Call mma_Allocate(ElementList(7)%Isotopes,16)
      ElementList(7)%Isotopes(:) = [
     &  Iso(14, 14.00307400443d0),
     &  Iso(15, 15.00010889888d0),
     &  Iso(10, 10.04165d0),
     &  Iso(11, 11.026091d0),
     &  Iso(12, 12.0186132d0),
     &  Iso(13, 13.00573861d0),
     &  Iso(16, 16.0061019d0),
     &  Iso(17, 17.008449d0),
     &  Iso(18, 18.014078d0),
     &  Iso(19, 19.017022d0),
     &  Iso(20, 20.023366d0),
     &  Iso(21, 21.02711d0),
     &  Iso(22, 22.03439d0),
     &  Iso(23, 23.04114d0),
     &  Iso(24, 24.05039d0),
     &  Iso(25, 25.0601d0) ]

      ElementList(8)%Symbol = 'O'
      ElementList(8)%nat = 3
      Call mma_Allocate(ElementList(8)%Isotopes,17)
      ElementList(8)%Isotopes(:) = [
     &  Iso(16, 15.99491461957d0),
     &  Iso(18, 17.99915961286d0),
     &  Iso(17, 16.9991317565d0),
     &  Iso(12, 12.034262d0),
     &  Iso(13, 13.024815d0),
     &  Iso(14, 14.00859636d0),
     &  Iso(15, 15.00306562d0),
     &  Iso(19, 19.003578d0),
     &  Iso(20, 20.00407535d0),
     &  Iso(21, 21.008655d0),
     &  Iso(22, 22.009966d0),
     &  Iso(23, 23.015696d0),
     &  Iso(24, 24.01986d0),
     &  Iso(25, 25.02936d0),
     &  Iso(26, 26.03729d0),
     &  Iso(27, 27.04772d0),
     &  Iso(28, 28.05591d0) ]

      ElementList(9)%Symbol = 'F'
      ElementList(9)%nat = 1
      Call mma_Allocate(ElementList(9)%Isotopes,18)
      ElementList(9)%Isotopes(:) = [
     &  Iso(19, 18.99840316273d0),
     &  Iso(14, 14.034315d0),
     &  Iso(15, 15.018043d0),
     &  Iso(16, 16.0114657d0),
     &  Iso(17, 17.00209524d0),
     &  Iso(18, 18.00093733d0),
     &  Iso(20, 19.999981252d0),
     &  Iso(21, 20.9999489d0),
     &  Iso(22, 22.002999d0),
     &  Iso(23, 23.003557d0),
     &  Iso(24, 24.008115d0),
     &  Iso(25, 25.012199d0),
     &  Iso(26, 26.020038d0),
     &  Iso(27, 27.02644d0),
     &  Iso(28, 28.03534d0),
     &  Iso(29, 29.04254d0),
     &  Iso(30, 30.05165d0),
     &  Iso(31, 31.05971d0) ]

      ElementList(10)%Symbol = 'Ne'
      ElementList(10)%nat = 3
      Call mma_Allocate(ElementList(10)%Isotopes,19)
      ElementList(10)%Isotopes(:) = [
     &  Iso(20, 19.9924401762d0),
     &  Iso(22, 21.991385114d0),
     &  Iso(21, 20.993846685d0),
     &  Iso(16, 16.02575d0),
     &  Iso(17, 17.01771396d0),
     &  Iso(18, 18.0057087d0),
     &  Iso(19, 19.00188091d0),
     &  Iso(23, 22.99446691d0),
     &  Iso(24, 23.99361065d0),
     &  Iso(25, 24.997789d0),
     &  Iso(26, 26.000515d0),
     &  Iso(27, 27.007553d0),
     &  Iso(28, 28.01212d0),
     &  Iso(29, 29.01975d0),
     &  Iso(30, 30.02473d0),
     &  Iso(31, 31.0331d0),
     &  Iso(32, 32.03972d0),
     &  Iso(33, 33.04938d0),
     &  Iso(34, 34.05673d0) ]

      ElementList(11)%Symbol = 'Na'
      ElementList(11)%nat = 1
      Call mma_Allocate(ElementList(11)%Isotopes,20)
      ElementList(11)%Isotopes(:) = [
     &  Iso(23, 22.989769282d0),
     &  Iso(18, 18.02688d0),
     &  Iso(19, 19.01388d0),
     &  Iso(20, 20.0073544d0),
     &  Iso(21, 20.99765469d0),
     &  Iso(22, 21.99443741d0),
     &  Iso(24, 23.99096295d0),
     &  Iso(25, 24.989954d0),
     &  Iso(26, 25.9926346d0),
     &  Iso(27, 26.9940765d0),
     &  Iso(28, 27.998939d0),
     &  Iso(29, 29.0028771d0),
     &  Iso(30, 30.0090979d0),
     &  Iso(31, 31.013163d0),
     &  Iso(32, 32.02019d0),
     &  Iso(33, 33.02573d0),
     &  Iso(34, 34.03359d0),
     &  Iso(35, 35.04062d0),
     &  Iso(36, 36.04929d0),
     &  Iso(37, 37.05705d0) ]

      ElementList(12)%Symbol = 'Mg'
      ElementList(12)%nat = 3
      Call mma_Allocate(ElementList(12)%Isotopes,22)
      ElementList(12)%Isotopes(:) = [
     &  Iso(24, 23.985041697d0),
     &  Iso(26, 25.982592968d0),
     &  Iso(25, 24.985836976d0),
     &  Iso(19, 19.034169d0),
     &  Iso(20, 20.01885d0),
     &  Iso(21, 21.011716d0),
     &  Iso(22, 21.99957065d0),
     &  Iso(23, 22.99412421d0),
     &  Iso(27, 26.984340624d0),
     &  Iso(28, 27.9838767d0),
     &  Iso(29, 28.988617d0),
     &  Iso(30, 29.9904629d0),
     &  Iso(31, 30.996648d0),
     &  Iso(32, 31.9991102d0),
     &  Iso(33, 33.0053271d0),
     &  Iso(34, 34.008935d0),
     &  Iso(35, 35.01679d0),
     &  Iso(36, 36.02188d0),
     &  Iso(37, 37.03037d0),
     &  Iso(38, 38.03658d0),
     &  Iso(39, 39.04538d0),
     &  Iso(40, 40.05218d0) ]

      ElementList(13)%Symbol = 'Al'
      ElementList(13)%nat = 1
      Call mma_Allocate(ElementList(13)%Isotopes,23)
      ElementList(13)%Isotopes(:) = [
     &  Iso(27, 26.98153853d0),
     &  Iso(21, 21.02897d0),
     &  Iso(22, 22.01954d0),
     &  Iso(23, 23.00724435d0),
     &  Iso(24, 23.9999489d0),
     &  Iso(25, 24.9904281d0),
     &  Iso(26, 25.986891904d0),
     &  Iso(28, 27.98191021d0),
     &  Iso(29, 28.9804565d0),
     &  Iso(30, 29.98296d0),
     &  Iso(31, 30.983945d0),
     &  Iso(32, 31.988085d0),
     &  Iso(33, 32.990909d0),
     &  Iso(34, 33.996705d0),
     &  Iso(35, 34.999764d0),
     &  Iso(36, 36.00639d0),
     &  Iso(37, 37.01053d0),
     &  Iso(38, 38.0174d0),
     &  Iso(39, 39.02254d0),
     &  Iso(40, 40.03003d0),
     &  Iso(41, 41.03638d0),
     &  Iso(42, 42.04384d0),
     &  Iso(43, 43.05147d0) ]

      ElementList(14)%Symbol = 'Si'
      ElementList(14)%nat = 3
      Call mma_Allocate(ElementList(14)%Isotopes,24)
      ElementList(14)%Isotopes(:) = [
     &  Iso(28, 27.97692653465d0),
     &  Iso(29, 28.9764946649d0),
     &  Iso(30, 29.973770136d0),
     &  Iso(22, 22.03579d0),
     &  Iso(23, 23.02544d0),
     &  Iso(24, 24.011535d0),
     &  Iso(25, 25.004109d0),
     &  Iso(26, 25.99233384d0),
     &  Iso(27, 26.98670481d0),
     &  Iso(31, 30.975363194d0),
     &  Iso(32, 31.97415154d0),
     &  Iso(33, 32.97797696d0),
     &  Iso(34, 33.978576d0),
     &  Iso(35, 34.984583d0),
     &  Iso(36, 35.986695d0),
     &  Iso(37, 36.992921d0),
     &  Iso(38, 37.995523d0),
     &  Iso(39, 39.002491d0),
     &  Iso(40, 40.00583d0),
     &  Iso(41, 41.01301d0),
     &  Iso(42, 42.01778d0),
     &  Iso(43, 43.0248d0),
     &  Iso(44, 44.03061d0),
     &  Iso(45, 45.03995d0) ]

      ElementList(15)%Symbol = 'P'
      ElementList(15)%nat = 1
      Call mma_Allocate(ElementList(15)%Isotopes,24)
      ElementList(15)%Isotopes(:) = [
     &  Iso(31, 30.97376199842d0),
     &  Iso(24, 24.03577d0),
     &  Iso(25, 25.02119d0),
     &  Iso(26, 26.01178d0),
     &  Iso(27, 26.999224d0),
     &  Iso(28, 27.9923266d0),
     &  Iso(29, 28.98180079d0),
     &  Iso(30, 29.97831375d0),
     &  Iso(32, 31.973907643d0),
     &  Iso(33, 32.9717257d0),
     &  Iso(34, 33.97364589d0),
     &  Iso(35, 34.9733141d0),
     &  Iso(36, 35.97826d0),
     &  Iso(37, 36.979607d0),
     &  Iso(38, 37.984252d0),
     &  Iso(39, 38.986227d0),
     &  Iso(40, 39.99133d0),
     &  Iso(41, 40.994654d0),
     &  Iso(42, 42.00108d0),
     &  Iso(43, 43.00502d0),
     &  Iso(44, 44.01121d0),
     &  Iso(45, 45.01645d0),
     &  Iso(46, 46.02446d0),
     &  Iso(47, 47.03139d0) ]

      ElementList(16)%Symbol = 'S'
      ElementList(16)%nat = 4
      Call mma_Allocate(ElementList(16)%Isotopes,24)
      ElementList(16)%Isotopes(:) = [
     &  Iso(32, 31.9720711744d0),
     &  Iso(34, 33.967867004d0),
     &  Iso(33, 32.9714589098d0),
     &  Iso(36, 35.96708071d0),
     &  Iso(26, 26.02907d0),
     &  Iso(27, 27.01828d0),
     &  Iso(28, 28.00437d0),
     &  Iso(29, 28.996611d0),
     &  Iso(30, 29.98490703d0),
     &  Iso(31, 30.97955701d0),
     &  Iso(35, 34.96903231d0),
     &  Iso(37, 36.97112551d0),
     &  Iso(38, 37.9711633d0),
     &  Iso(39, 38.975134d0),
     &  Iso(40, 39.9754826d0),
     &  Iso(41, 40.9795935d0),
     &  Iso(42, 41.9810651d0),
     &  Iso(43, 42.9869076d0),
     &  Iso(44, 43.9901188d0),
     &  Iso(45, 44.99572d0),
     &  Iso(46, 46.00004d0),
     &  Iso(47, 47.00795d0),
     &  Iso(48, 48.0137d0),
     &  Iso(49, 49.02276d0) ]

      ElementList(17)%Symbol = 'Cl'
      ElementList(17)%nat = 2
      Call mma_Allocate(ElementList(17)%Isotopes,24)
      ElementList(17)%Isotopes(:) = [
     &  Iso(35, 34.968852682d0),
     &  Iso(37, 36.965902602d0),
     &  Iso(28, 28.02954d0),
     &  Iso(29, 29.01478d0),
     &  Iso(30, 30.00477d0),
     &  Iso(31, 30.992414d0),
     &  Iso(32, 31.98568464d0),
     &  Iso(33, 32.97745199d0),
     &  Iso(34, 33.973762485d0),
     &  Iso(36, 35.968306809d0),
     &  Iso(38, 37.96801044d0),
     &  Iso(39, 38.9680082d0),
     &  Iso(40, 39.970415d0),
     &  Iso(41, 40.970685d0),
     &  Iso(42, 41.97325d0),
     &  Iso(43, 42.97389d0),
     &  Iso(44, 43.97787d0),
     &  Iso(45, 44.98029d0),
     &  Iso(46, 45.98517d0),
     &  Iso(47, 46.98916d0),
     &  Iso(48, 47.99564d0),
     &  Iso(49, 49.00123d0),
     &  Iso(50, 50.00905d0),
     &  Iso(51, 51.01554d0) ]

      ElementList(18)%Symbol = 'Ar'
      ElementList(18)%nat = 3
      Call mma_Allocate(ElementList(18)%Isotopes,24)
      ElementList(18)%Isotopes(:) = [
     &  Iso(40, 39.9623831237d0),
     &  Iso(36, 35.967545105d0),
     &  Iso(38, 37.96273211d0),
     &  Iso(30, 30.02307d0),
     &  Iso(31, 31.01212d0),
     &  Iso(32, 31.9976378d0),
     &  Iso(33, 32.98992555d0),
     &  Iso(34, 33.98027009d0),
     &  Iso(35, 34.97525759d0),
     &  Iso(37, 36.96677633d0),
     &  Iso(39, 38.964313d0),
     &  Iso(41, 40.96450057d0),
     &  Iso(42, 41.9630457d0),
     &  Iso(43, 42.9656361d0),
     &  Iso(44, 43.9649238d0),
     &  Iso(45, 44.96803973d0),
     &  Iso(46, 45.968083d0),
     &  Iso(47, 46.972935d0),
     &  Iso(48, 47.97591d0),
     &  Iso(49, 48.9819d0),
     &  Iso(50, 49.98613d0),
     &  Iso(51, 50.9937d0),
     &  Iso(52, 51.99896d0),
     &  Iso(53, 53.00729d0) ]

      ElementList(19)%Symbol = 'K'
      ElementList(19)%nat = 3
      Call mma_Allocate(ElementList(19)%Isotopes,25)
      ElementList(19)%Isotopes(:) = [
     &  Iso(39, 38.9637064864d0),
     &  Iso(41, 40.9618252579d0),
     &  Iso(40, 39.963998166d0),
     &  Iso(32, 32.02265d0),
     &  Iso(33, 33.00756d0),
     &  Iso(34, 33.99869d0),
     &  Iso(35, 34.98800541d0),
     &  Iso(36, 35.98130201d0),
     &  Iso(37, 36.97337589d0),
     &  Iso(38, 37.96908112d0),
     &  Iso(42, 41.96240231d0),
     &  Iso(43, 42.9607347d0),
     &  Iso(44, 43.96158699d0),
     &  Iso(45, 44.96069149d0),
     &  Iso(46, 45.96198159d0),
     &  Iso(47, 46.9616616d0),
     &  Iso(48, 47.96534119d0),
     &  Iso(49, 48.96821075d0),
     &  Iso(50, 49.97238d0),
     &  Iso(51, 50.975828d0),
     &  Iso(52, 51.98224d0),
     &  Iso(53, 52.98746d0),
     &  Iso(54, 53.99463d0),
     &  Iso(55, 55.00076d0),
     &  Iso(56, 56.00851d0) ]

      ElementList(20)%Symbol = 'Ca'
      ElementList(20)%nat = 6
      Call mma_Allocate(ElementList(20)%Isotopes,25)
      ElementList(20)%Isotopes(:) = [
     &  Iso(40, 39.962590863d0),
     &  Iso(44, 43.95548156d0),
     &  Iso(42, 41.95861783d0),
     &  Iso(48, 47.95252276d0),
     &  Iso(43, 42.95876644d0),
     &  Iso(46, 45.953689d0),
     &  Iso(34, 34.01487d0),
     &  Iso(35, 35.00514d0),
     &  Iso(36, 35.993074d0),
     &  Iso(37, 36.98589785d0),
     &  Iso(38, 37.97631922d0),
     &  Iso(39, 38.97071081d0),
     &  Iso(41, 40.96227792d0),
     &  Iso(45, 44.95618635d0),
     &  Iso(47, 46.9545424d0),
     &  Iso(49, 48.95566274d0),
     &  Iso(50, 49.9574992d0),
     &  Iso(51, 50.960989d0),
     &  Iso(52, 51.963217d0),
     &  Iso(53, 52.96945d0),
     &  Iso(54, 53.9734d0),
     &  Iso(55, 54.9803d0),
     &  Iso(56, 55.98508d0),
     &  Iso(57, 56.99262d0),
     &  Iso(58, 57.99794d0) ]

      ElementList(21)%Symbol = 'Sc'
      ElementList(21)%nat = 1
      Call mma_Allocate(ElementList(21)%Isotopes,26)
      ElementList(21)%Isotopes(:) = [
     &  Iso(45, 44.95590828d0),
     &  Iso(36, 36.01648d0),
     &  Iso(37, 37.00374d0),
     &  Iso(38, 37.99512d0),
     &  Iso(39, 38.984785d0),
     &  Iso(40, 39.9779673d0),
     &  Iso(41, 40.969251105d0),
     &  Iso(42, 41.96551653d0),
     &  Iso(43, 42.9611505d0),
     &  Iso(44, 43.9594029d0),
     &  Iso(46, 45.95516826d0),
     &  Iso(47, 46.9524037d0),
     &  Iso(48, 47.9522236d0),
     &  Iso(49, 48.9500146d0),
     &  Iso(50, 49.952176d0),
     &  Iso(51, 50.953592d0),
     &  Iso(52, 51.95688d0),
     &  Iso(53, 52.95909d0),
     &  Iso(54, 53.96393d0),
     &  Iso(55, 54.96782d0),
     &  Iso(56, 55.97345d0),
     &  Iso(57, 56.97777d0),
     &  Iso(58, 57.98403d0),
     &  Iso(59, 58.98894d0),
     &  Iso(60, 59.99565d0),
     &  Iso(61, 61.001d0) ]

      ElementList(22)%Symbol = 'Ti'
      ElementList(22)%nat = 5
      Call mma_Allocate(ElementList(22)%Isotopes,26)
      ElementList(22)%Isotopes(:) = [
     &  Iso(48, 47.94794198d0),
     &  Iso(46, 45.95262772d0),
     &  Iso(47, 46.95175879d0),
     &  Iso(49, 48.94786568d0),
     &  Iso(50, 49.94478689d0),
     &  Iso(38, 38.01145d0),
     &  Iso(39, 39.00236d0),
     &  Iso(40, 39.9905d0),
     &  Iso(41, 40.983148d0),
     &  Iso(42, 41.97304903d0),
     &  Iso(43, 42.9685225d0),
     &  Iso(44, 43.95968995d0),
     &  Iso(45, 44.95812198d0),
     &  Iso(51, 50.94661065d0),
     &  Iso(52, 51.946893d0),
     &  Iso(53, 52.94973d0),
     &  Iso(54, 53.95105d0),
     &  Iso(55, 54.95527d0),
     &  Iso(56, 55.95791d0),
     &  Iso(57, 56.96364d0),
     &  Iso(58, 57.9666d0),
     &  Iso(59, 58.97247d0),
     &  Iso(60, 59.97603d0),
     &  Iso(61, 60.98245d0),
     &  Iso(62, 61.98651d0),
     &  Iso(63, 62.99375d0) ]

      ElementList(23)%Symbol = 'V'
      ElementList(23)%nat = 2
      Call mma_Allocate(ElementList(23)%Isotopes,27)
      ElementList(23)%Isotopes(:) = [
     &  Iso(51, 50.94395704d0),
     &  Iso(50, 49.94715601d0),
     &  Iso(40, 40.01276d0),
     &  Iso(41, 41.00021d0),
     &  Iso(42, 41.99182d0),
     &  Iso(43, 42.980766d0),
     &  Iso(44, 43.97411d0),
     &  Iso(45, 44.9657748d0),
     &  Iso(46, 45.96019878d0),
     &  Iso(47, 46.95490491d0),
     &  Iso(48, 47.9522522d0),
     &  Iso(49, 48.9485118d0),
     &  Iso(52, 51.94477301d0),
     &  Iso(53, 52.9443367d0),
     &  Iso(54, 53.946439d0),
     &  Iso(55, 54.94724d0),
     &  Iso(56, 55.95048d0),
     &  Iso(57, 56.95252d0),
     &  Iso(58, 57.95672d0),
     &  Iso(59, 58.95939d0),
     &  Iso(60, 59.96431d0),
     &  Iso(61, 60.96725d0),
     &  Iso(62, 61.97265d0),
     &  Iso(63, 62.97639d0),
     &  Iso(64, 63.98264d0),
     &  Iso(65, 64.9875d0),
     &  Iso(66, 65.99398d0) ]

      ElementList(24)%Symbol = 'Cr'
      ElementList(24)%nat = 4
      Call mma_Allocate(ElementList(24)%Isotopes,27)
      ElementList(24)%Isotopes(:) = [
     &  Iso(52, 51.94050623d0),
     &  Iso(53, 52.94064815d0),
     &  Iso(50, 49.94604183d0),
     &  Iso(54, 53.93887916d0),
     &  Iso(42, 42.0067d0),
     &  Iso(43, 42.99753d0),
     &  Iso(44, 43.98536d0),
     &  Iso(45, 44.97905d0),
     &  Iso(46, 45.968359d0),
     &  Iso(47, 46.9628974d0),
     &  Iso(48, 47.9540291d0),
     &  Iso(49, 48.9513333d0),
     &  Iso(51, 50.94476502d0),
     &  Iso(55, 54.94083843d0),
     &  Iso(56, 55.9406531d0),
     &  Iso(57, 56.943613d0),
     &  Iso(58, 57.94435d0),
     &  Iso(59, 58.94859d0),
     &  Iso(60, 59.95008d0),
     &  Iso(61, 60.95442d0),
     &  Iso(62, 61.9561d0),
     &  Iso(63, 62.96165d0),
     &  Iso(64, 63.96408d0),
     &  Iso(65, 64.96996d0),
     &  Iso(66, 65.97366d0),
     &  Iso(67, 66.98016d0),
     &  Iso(68, 67.98403d0) ]

      ElementList(25)%Symbol = 'Mn'
      ElementList(25)%nat = 1
      Call mma_Allocate(ElementList(25)%Isotopes,28)
      ElementList(25)%Isotopes(:) = [
     &  Iso(55, 54.93804391d0),
     &  Iso(44, 44.00715d0),
     &  Iso(45, 44.99449d0),
     &  Iso(46, 45.98609d0),
     &  Iso(47, 46.975775d0),
     &  Iso(48, 47.96852d0),
     &  Iso(49, 48.959595d0),
     &  Iso(50, 49.95423778d0),
     &  Iso(51, 50.94820847d0),
     &  Iso(52, 51.9455639d0),
     &  Iso(53, 52.94128889d0),
     &  Iso(54, 53.9403576d0),
     &  Iso(56, 55.93890369d0),
     &  Iso(57, 56.9382861d0),
     &  Iso(58, 57.9400666d0),
     &  Iso(59, 58.9403911d0),
     &  Iso(60, 59.9431366d0),
     &  Iso(61, 60.9444525d0),
     &  Iso(62, 61.94795d0),
     &  Iso(63, 62.9496647d0),
     &  Iso(64, 63.9538494d0),
     &  Iso(65, 64.9560198d0),
     &  Iso(66, 65.960547d0),
     &  Iso(67, 66.96424d0),
     &  Iso(68, 67.96962d0),
     &  Iso(69, 68.97366d0),
     &  Iso(70, 69.97937d0),
     &  Iso(71, 70.98368d0) ]

      ElementList(26)%Symbol = 'Fe'
      ElementList(26)%nat = 4
      Call mma_Allocate(ElementList(26)%Isotopes,30)
      ElementList(26)%Isotopes(:) = [
     &  Iso(56, 55.93493633d0),
     &  Iso(54, 53.93960899d0),
     &  Iso(57, 56.93539284d0),
     &  Iso(58, 57.93327443d0),
     &  Iso(45, 45.01442d0),
     &  Iso(46, 46.00063d0),
     &  Iso(47, 46.99185d0),
     &  Iso(48, 47.98023d0),
     &  Iso(49, 48.973429d0),
     &  Iso(50, 49.962975d0),
     &  Iso(51, 50.956841d0),
     &  Iso(52, 51.9481131d0),
     &  Iso(53, 52.9453064d0),
     &  Iso(55, 54.93829199d0),
     &  Iso(59, 58.93487434d0),
     &  Iso(60, 59.9340711d0),
     &  Iso(61, 60.9367462d0),
     &  Iso(62, 61.9367918d0),
     &  Iso(63, 62.9402727d0),
     &  Iso(64, 63.9409878d0),
     &  Iso(65, 64.9450115d0),
     &  Iso(66, 65.94625d0),
     &  Iso(67, 66.95054d0),
     &  Iso(68, 67.95295d0),
     &  Iso(69, 68.95807d0),
     &  Iso(70, 69.96102d0),
     &  Iso(71, 70.96672d0),
     &  Iso(72, 71.96983d0),
     &  Iso(73, 72.97572d0),
     &  Iso(74, 73.97935d0) ]

      ElementList(27)%Symbol = 'Co'
      ElementList(27)%nat = 1
      Call mma_Allocate(ElementList(27)%Isotopes,30)
      ElementList(27)%Isotopes(:) = [
     &  Iso(59, 58.93319429d0),
     &  Iso(47, 47.01057d0),
     &  Iso(48, 48.00093d0),
     &  Iso(49, 48.98891d0),
     &  Iso(50, 49.98091d0),
     &  Iso(51, 50.970647d0),
     &  Iso(52, 51.96351d0),
     &  Iso(53, 52.9542041d0),
     &  Iso(54, 53.94845987d0),
     &  Iso(55, 54.9419972d0),
     &  Iso(56, 55.9398388d0),
     &  Iso(57, 56.93629057d0),
     &  Iso(58, 57.9357521d0),
     &  Iso(60, 59.9338163d0),
     &  Iso(61, 60.93247662d0),
     &  Iso(62, 61.934059d0),
     &  Iso(63, 62.9336d0),
     &  Iso(64, 63.935811d0),
     &  Iso(65, 64.9364621d0),
     &  Iso(66, 65.939443d0),
     &  Iso(67, 66.9406096d0),
     &  Iso(68, 67.94426d0),
     &  Iso(69, 68.94614d0),
     &  Iso(70, 69.94963d0),
     &  Iso(71, 70.95237d0),
     &  Iso(72, 71.95729d0),
     &  Iso(73, 72.96039d0),
     &  Iso(74, 73.96515d0),
     &  Iso(75, 74.96876d0),
     &  Iso(76, 75.97413d0) ]

      ElementList(28)%Symbol = 'Ni'
      ElementList(28)%nat = 5
      Call mma_Allocate(ElementList(28)%Isotopes,32)
      ElementList(28)%Isotopes(:) = [
     &  Iso(58, 57.93534241d0),
     &  Iso(60, 59.93078588d0),
     &  Iso(62, 61.92834537d0),
     &  Iso(61, 60.93105557d0),
     &  Iso(64, 63.92796682d0),
     &  Iso(48, 48.01769d0),
     &  Iso(49, 49.0077d0),
     &  Iso(50, 49.99474d0),
     &  Iso(51, 50.98611d0),
     &  Iso(52, 51.9748d0),
     &  Iso(53, 52.96819d0),
     &  Iso(54, 53.957892d0),
     &  Iso(55, 54.95133063d0),
     &  Iso(56, 55.94212855d0),
     &  Iso(57, 56.93979218d0),
     &  Iso(59, 58.9343462d0),
     &  Iso(63, 62.92966963d0),
     &  Iso(65, 64.93008517d0),
     &  Iso(66, 65.9291393d0),
     &  Iso(67, 66.9315694d0),
     &  Iso(68, 67.9318688d0),
     &  Iso(69, 68.9356103d0),
     &  Iso(70, 69.9364313d0),
     &  Iso(71, 70.940519d0),
     &  Iso(72, 71.9417859d0),
     &  Iso(73, 72.9462067d0),
     &  Iso(74, 73.94798d0),
     &  Iso(75, 74.9525d0),
     &  Iso(76, 75.95533d0),
     &  Iso(77, 76.96055d0),
     &  Iso(78, 77.96336d0),
     &  Iso(79, 78.97025d0) ]

      ElementList(29)%Symbol = 'Cu'
      ElementList(29)%nat = 2
      Call mma_Allocate(ElementList(29)%Isotopes,31)
      ElementList(29)%Isotopes(:) = [
     &  Iso(63, 62.92959772d0),
     &  Iso(65, 64.9277897d0),
     &  Iso(52, 51.99671d0),
     &  Iso(53, 52.98459d0),
     &  Iso(54, 53.97666d0),
     &  Iso(55, 54.96604d0),
     &  Iso(56, 55.95895d0),
     &  Iso(57, 56.9492125d0),
     &  Iso(58, 57.94453305d0),
     &  Iso(59, 58.93949748d0),
     &  Iso(60, 59.9373645d0),
     &  Iso(61, 60.9334576d0),
     &  Iso(62, 61.93259541d0),
     &  Iso(64, 63.92976434d0),
     &  Iso(66, 65.92886903d0),
     &  Iso(67, 66.9277303d0),
     &  Iso(68, 67.9296109d0),
     &  Iso(69, 68.9294293d0),
     &  Iso(70, 69.9323921d0),
     &  Iso(71, 70.9326768d0),
     &  Iso(72, 71.9358203d0),
     &  Iso(73, 72.9366744d0),
     &  Iso(74, 73.9398749d0),
     &  Iso(75, 74.9415226d0),
     &  Iso(76, 75.945275d0),
     &  Iso(77, 76.94792d0),
     &  Iso(78, 77.95223d0),
     &  Iso(79, 78.95502d0),
     &  Iso(80, 79.96089d0),
     &  Iso(81, 80.96587d0),
     &  Iso(82, 81.97244d0) ]

      ElementList(30)%Symbol = 'Zn'
      ElementList(30)%nat = 5
      Call mma_Allocate(ElementList(30)%Isotopes,32)
      ElementList(30)%Isotopes(:) = [
     &  Iso(64, 63.92914201d0),
     &  Iso(66, 65.92603381d0),
     &  Iso(68, 67.92484455d0),
     &  Iso(67, 66.92712775d0),
     &  Iso(70, 69.9253192d0),
     &  Iso(54, 53.99204d0),
     &  Iso(55, 54.98398d0),
     &  Iso(56, 55.97254d0),
     &  Iso(57, 56.96506d0),
     &  Iso(58, 57.954591d0),
     &  Iso(59, 58.94931266d0),
     &  Iso(60, 59.9418421d0),
     &  Iso(61, 60.939507d0),
     &  Iso(62, 61.93433397d0),
     &  Iso(63, 62.9332115d0),
     &  Iso(65, 64.92924077d0),
     &  Iso(69, 68.9265507d0),
     &  Iso(71, 70.9277196d0),
     &  Iso(72, 71.9268428d0),
     &  Iso(73, 72.9295826d0),
     &  Iso(74, 73.9294073d0),
     &  Iso(75, 74.9328402d0),
     &  Iso(76, 75.933115d0),
     &  Iso(77, 76.9368872d0),
     &  Iso(78, 77.9382892d0),
     &  Iso(79, 78.9426381d0),
     &  Iso(80, 79.9445529d0),
     &  Iso(81, 80.9504026d0),
     &  Iso(82, 81.95426d0),
     &  Iso(83, 82.96056d0),
     &  Iso(84, 83.96521d0),
     &  Iso(85, 84.97226d0) ]

      ElementList(31)%Symbol = 'Ga'
      ElementList(31)%nat = 2
      Call mma_Allocate(ElementList(31)%Isotopes,32)
      ElementList(31)%Isotopes(:) = [
     &  Iso(69, 68.9255735d0),
     &  Iso(71, 70.92470258d0),
     &  Iso(56, 55.99536d0),
     &  Iso(57, 56.9832d0),
     &  Iso(58, 57.97478d0),
     &  Iso(59, 58.96353d0),
     &  Iso(60, 59.95729d0),
     &  Iso(61, 60.949399d0),
     &  Iso(62, 61.94419025d0),
     &  Iso(63, 62.9392942d0),
     &  Iso(64, 63.9368404d0),
     &  Iso(65, 64.93273459d0),
     &  Iso(66, 65.9315894d0),
     &  Iso(67, 66.9282025d0),
     &  Iso(68, 67.9279805d0),
     &  Iso(70, 69.9260219d0),
     &  Iso(72, 71.92636747d0),
     &  Iso(73, 72.9251747d0),
     &  Iso(74, 73.9269457d0),
     &  Iso(75, 74.9265002d0),
     &  Iso(76, 75.9288276d0),
     &  Iso(77, 76.9291543d0),
     &  Iso(78, 77.9316088d0),
     &  Iso(79, 78.9328523d0),
     &  Iso(80, 79.9364208d0),
     &  Iso(81, 80.9381338d0),
     &  Iso(82, 81.9431765d0),
     &  Iso(83, 82.9471203d0),
     &  Iso(84, 83.95246d0),
     &  Iso(85, 84.95699d0),
     &  Iso(86, 85.96301d0),
     &  Iso(87, 86.96824d0) ]

      ElementList(32)%Symbol = 'Ge'
      ElementList(32)%nat = 5
      Call mma_Allocate(ElementList(32)%Isotopes,33)
      ElementList(32)%Isotopes(:) = [
     &  Iso(74, 73.921177761d0),
     &  Iso(72, 71.922075826d0),
     &  Iso(70, 69.92424875d0),
     &  Iso(73, 72.923458956d0),
     &  Iso(76, 75.921402726d0),
     &  Iso(58, 57.99172d0),
     &  Iso(59, 58.98249d0),
     &  Iso(60, 59.97036d0),
     &  Iso(61, 60.96379d0),
     &  Iso(62, 61.95502d0),
     &  Iso(63, 62.949628d0),
     &  Iso(64, 63.9416899d0),
     &  Iso(65, 64.9393681d0),
     &  Iso(66, 65.9338621d0),
     &  Iso(67, 66.9327339d0),
     &  Iso(68, 67.9280953d0),
     &  Iso(69, 68.9279645d0),
     &  Iso(71, 70.92495233d0),
     &  Iso(75, 74.92285837d0),
     &  Iso(77, 76.923549843d0),
     &  Iso(78, 77.9228529d0),
     &  Iso(79, 78.92536d0),
     &  Iso(80, 79.9253508d0),
     &  Iso(81, 80.9288329d0),
     &  Iso(82, 81.929774d0),
     &  Iso(83, 82.9345391d0),
     &  Iso(84, 83.9375751d0),
     &  Iso(85, 84.9429697d0),
     &  Iso(86, 85.94658d0),
     &  Iso(87, 86.95268d0),
     &  Iso(88, 87.95691d0),
     &  Iso(89, 88.96379d0),
     &  Iso(90, 89.96863d0) ]

      ElementList(33)%Symbol = 'As'
      ElementList(33)%nat = 1
      Call mma_Allocate(ElementList(33)%Isotopes,33)
      ElementList(33)%Isotopes(:) = [
     &  Iso(75, 74.92159457d0),
     &  Iso(60, 59.99388d0),
     &  Iso(61, 60.98112d0),
     &  Iso(62, 61.97361d0),
     &  Iso(63, 62.9639d0),
     &  Iso(64, 63.95743d0),
     &  Iso(65, 64.949611d0),
     &  Iso(66, 65.9441488d0),
     &  Iso(67, 66.93925111d0),
     &  Iso(68, 67.9367741d0),
     &  Iso(69, 68.932246d0),
     &  Iso(70, 69.930926d0),
     &  Iso(71, 70.9271138d0),
     &  Iso(72, 71.9267523d0),
     &  Iso(73, 72.9238291d0),
     &  Iso(74, 73.9239286d0),
     &  Iso(76, 75.92239202d0),
     &  Iso(77, 76.9206476d0),
     &  Iso(78, 77.921828d0),
     &  Iso(79, 78.9209484d0),
     &  Iso(80, 79.9224746d0),
     &  Iso(81, 80.9221323d0),
     &  Iso(82, 81.9247412d0),
     &  Iso(83, 82.9252069d0),
     &  Iso(84, 83.9293033d0),
     &  Iso(85, 84.9321637d0),
     &  Iso(86, 85.9367015d0),
     &  Iso(87, 86.9402917d0),
     &  Iso(88, 87.94555d0),
     &  Iso(89, 88.94976d0),
     &  Iso(90, 89.95563d0),
     &  Iso(91, 90.96039d0),
     &  Iso(92, 91.96674d0) ]

      ElementList(34)%Symbol = 'Se'
      ElementList(34)%nat = 6
      Call mma_Allocate(ElementList(34)%Isotopes,32)
      ElementList(34)%Isotopes(:) = [
     &  Iso(80, 79.9165218d0),
     &  Iso(78, 77.91730928d0),
     &  Iso(76, 75.919213704d0),
     &  Iso(82, 81.9166995d0),
     &  Iso(77, 76.919914154d0),
     &  Iso(74, 73.922475934d0),
     &  Iso(64, 63.97109d0),
     &  Iso(65, 64.9644d0),
     &  Iso(66, 65.95559d0),
     &  Iso(67, 66.949994d0),
     &  Iso(68, 67.94182524d0),
     &  Iso(69, 68.9394148d0),
     &  Iso(70, 69.9335155d0),
     &  Iso(71, 70.9322094d0),
     &  Iso(72, 71.9271405d0),
     &  Iso(73, 72.9267549d0),
     &  Iso(75, 74.92252287d0),
     &  Iso(79, 78.91849929d0),
     &  Iso(81, 80.917993d0),
     &  Iso(83, 82.9191186d0),
     &  Iso(84, 83.9184668d0),
     &  Iso(85, 84.9222608d0),
     &  Iso(86, 85.9243117d0),
     &  Iso(87, 86.9286886d0),
     &  Iso(88, 87.9314175d0),
     &  Iso(89, 88.9366691d0),
     &  Iso(90, 89.9401d0),
     &  Iso(91, 90.94596d0),
     &  Iso(92, 91.94984d0),
     &  Iso(93, 92.95629d0),
     &  Iso(94, 93.96049d0),
     &  Iso(95, 94.9673d0) ]

      ElementList(35)%Symbol = 'Br'
      ElementList(35)%nat = 2
      Call mma_Allocate(ElementList(35)%Isotopes,32)
      ElementList(35)%Isotopes(:) = [
     &  Iso(79, 78.9183376d0),
     &  Iso(81, 80.9162897d0),
     &  Iso(67, 66.96465d0),
     &  Iso(68, 67.95873d0),
     &  Iso(69, 68.950497d0),
     &  Iso(70, 69.944792d0),
     &  Iso(71, 70.9393422d0),
     &  Iso(72, 71.9365886d0),
     &  Iso(73, 72.9316715d0),
     &  Iso(74, 73.9299102d0),
     &  Iso(75, 74.9258105d0),
     &  Iso(76, 75.924542d0),
     &  Iso(77, 76.9213792d0),
     &  Iso(78, 77.9211459d0),
     &  Iso(80, 79.9185298d0),
     &  Iso(82, 81.9168032d0),
     &  Iso(83, 82.9151756d0),
     &  Iso(84, 83.916496d0),
     &  Iso(85, 84.9156458d0),
     &  Iso(86, 85.9188054d0),
     &  Iso(87, 86.920674d0),
     &  Iso(88, 87.9240833d0),
     &  Iso(89, 88.9267046d0),
     &  Iso(90, 89.9312928d0),
     &  Iso(91, 90.9343986d0),
     &  Iso(92, 91.9396316d0),
     &  Iso(93, 92.94313d0),
     &  Iso(94, 93.9489d0),
     &  Iso(95, 94.95301d0),
     &  Iso(96, 95.95903d0),
     &  Iso(97, 96.96344d0),
     &  Iso(98, 97.96946d0) ]

      ElementList(36)%Symbol = 'Kr'
      ElementList(36)%nat = 6
      Call mma_Allocate(ElementList(36)%Isotopes,33)
      ElementList(36)%Isotopes(:) = [
     &  Iso(84, 83.9114977282d0),
     &  Iso(86, 85.9106106269d0),
     &  Iso(82, 81.91348273d0),
     &  Iso(83, 82.91412716d0),
     &  Iso(80, 79.91637808d0),
     &  Iso(78, 77.92036494d0),
     &  Iso(69, 68.96518d0),
     &  Iso(70, 69.95604d0),
     &  Iso(71, 70.95027d0),
     &  Iso(72, 71.9420924d0),
     &  Iso(73, 72.9392892d0),
     &  Iso(74, 73.933084d0),
     &  Iso(75, 74.9309457d0),
     &  Iso(76, 75.9259103d0),
     &  Iso(77, 76.92467d0),
     &  Iso(79, 78.9200829d0),
     &  Iso(81, 80.9165912d0),
     &  Iso(85, 84.9125273d0),
     &  Iso(87, 86.91335476d0),
     &  Iso(88, 87.9144479d0),
     &  Iso(89, 88.9178355d0),
     &  Iso(90, 89.9195279d0),
     &  Iso(91, 90.9238063d0),
     &  Iso(92, 91.9261731d0),
     &  Iso(93, 92.9311472d0),
     &  Iso(94, 93.93414d0),
     &  Iso(95, 94.939711d0),
     &  Iso(96, 95.943017d0),
     &  Iso(97, 96.94909d0),
     &  Iso(98, 97.95243d0),
     &  Iso(99, 98.95839d0),
     &  Iso(100, 99.96237d0),
     &  Iso(101, 100.96873d0) ]

      ElementList(37)%Symbol = 'Rb'
      ElementList(37)%nat = 2
      Call mma_Allocate(ElementList(37)%Isotopes,33)
      ElementList(37)%Isotopes(:) = [
     &  Iso(85, 84.9117897379d0),
     &  Iso(87, 86.909180531d0),
     &  Iso(71, 70.96532d0),
     &  Iso(72, 71.95908d0),
     &  Iso(73, 72.95053d0),
     &  Iso(74, 73.9442659d0),
     &  Iso(75, 74.9385732d0),
     &  Iso(76, 75.935073d0),
     &  Iso(77, 76.9304016d0),
     &  Iso(78, 77.9281419d0),
     &  Iso(79, 78.9239899d0),
     &  Iso(80, 79.9225164d0),
     &  Iso(81, 80.9189939d0),
     &  Iso(82, 81.918209d0),
     &  Iso(83, 82.9151142d0),
     &  Iso(84, 83.9143752d0),
     &  Iso(86, 85.91116743d0),
     &  Iso(88, 87.91131559d0),
     &  Iso(89, 88.9122783d0),
     &  Iso(90, 89.9147985d0),
     &  Iso(91, 90.9165372d0),
     &  Iso(92, 91.9197284d0),
     &  Iso(93, 92.9220393d0),
     &  Iso(94, 93.9263948d0),
     &  Iso(95, 94.92926d0),
     &  Iso(96, 95.9341334d0),
     &  Iso(97, 96.9371771d0),
     &  Iso(98, 97.9416869d0),
     &  Iso(99, 98.94503d0),
     &  Iso(100, 99.95003d0),
     &  Iso(101, 100.95404d0),
     &  Iso(102, 101.95952d0),
     &  Iso(103, 102.96392d0) ]

      ElementList(38)%Symbol = 'Sr'
      ElementList(38)%nat = 4
      Call mma_Allocate(ElementList(38)%Isotopes,35)
      ElementList(38)%Isotopes(:) = [
     &  Iso(88, 87.9056125d0),
     &  Iso(86, 85.9092606d0),
     &  Iso(87, 86.9088775d0),
     &  Iso(84, 83.9134191d0),
     &  Iso(73, 72.9657d0),
     &  Iso(74, 73.95617d0),
     &  Iso(75, 74.94995d0),
     &  Iso(76, 75.941763d0),
     &  Iso(77, 76.9379455d0),
     &  Iso(78, 77.93218d0),
     &  Iso(79, 78.9297077d0),
     &  Iso(80, 79.9245175d0),
     &  Iso(81, 80.9232114d0),
     &  Iso(82, 81.9183999d0),
     &  Iso(83, 82.9175544d0),
     &  Iso(85, 84.912932d0),
     &  Iso(89, 88.9074511d0),
     &  Iso(90, 89.90773d0),
     &  Iso(91, 90.9101954d0),
     &  Iso(92, 91.9110382d0),
     &  Iso(93, 92.9140242d0),
     &  Iso(94, 93.9153556d0),
     &  Iso(95, 94.9193529d0),
     &  Iso(96, 95.9217066d0),
     &  Iso(97, 96.926374d0),
     &  Iso(98, 97.9286888d0),
     &  Iso(99, 98.9328907d0),
     &  Iso(100, 99.93577d0),
     &  Iso(101, 100.940352d0),
     &  Iso(102, 101.943791d0),
     &  Iso(103, 102.94909d0),
     &  Iso(104, 103.95265d0),
     &  Iso(105, 104.95855d0),
     &  Iso(106, 105.96265d0),
     &  Iso(107, 106.96897d0) ]

      ElementList(39)%Symbol = 'Y'
      ElementList(39)%nat = 1
      Call mma_Allocate(ElementList(39)%Isotopes,34)
      ElementList(39)%Isotopes(:) = [
     &  Iso(89, 88.9058403d0),
     &  Iso(76, 75.95856d0),
     &  Iso(77, 76.949781d0),
     &  Iso(78, 77.94361d0),
     &  Iso(79, 78.93735d0),
     &  Iso(80, 79.9343561d0),
     &  Iso(81, 80.9294556d0),
     &  Iso(82, 81.9269314d0),
     &  Iso(83, 82.922485d0),
     &  Iso(84, 83.9206721d0),
     &  Iso(85, 84.916433d0),
     &  Iso(86, 85.914886d0),
     &  Iso(87, 86.9108761d0),
     &  Iso(88, 87.9095016d0),
     &  Iso(90, 89.9071439d0),
     &  Iso(91, 90.9072974d0),
     &  Iso(92, 91.9089451d0),
     &  Iso(93, 92.909578d0),
     &  Iso(94, 93.9115906d0),
     &  Iso(95, 94.9128161d0),
     &  Iso(96, 95.9158968d0),
     &  Iso(97, 96.9182741d0),
     &  Iso(98, 97.9223821d0),
     &  Iso(99, 98.924148d0),
     &  Iso(100, 99.927715d0),
     &  Iso(101, 100.9301477d0),
     &  Iso(102, 101.9343277d0),
     &  Iso(103, 102.937243d0),
     &  Iso(104, 103.94196d0),
     &  Iso(105, 104.94544d0),
     &  Iso(106, 105.95056d0),
     &  Iso(107, 106.95452d0),
     &  Iso(108, 107.95996d0),
     &  Iso(109, 108.96436d0) ]

      ElementList(40)%Symbol = 'Zr'
      ElementList(40)%nat = 5
      Call mma_Allocate(ElementList(40)%Isotopes,35)
      ElementList(40)%Isotopes(:) = [
     &  Iso(90, 89.9046977d0),
     &  Iso(94, 93.9063108d0),
     &  Iso(92, 91.9050347d0),
     &  Iso(91, 90.9056396d0),
     &  Iso(96, 95.9082714d0),
     &  Iso(78, 77.95566d0),
     &  Iso(79, 78.94948d0),
     &  Iso(80, 79.9404d0),
     &  Iso(81, 80.93731d0),
     &  Iso(82, 81.93135d0),
     &  Iso(83, 82.9292421d0),
     &  Iso(84, 83.9233269d0),
     &  Iso(85, 84.9214444d0),
     &  Iso(86, 85.9162972d0),
     &  Iso(87, 86.914818d0),
     &  Iso(88, 87.9102213d0),
     &  Iso(89, 88.9088814d0),
     &  Iso(93, 92.9064699d0),
     &  Iso(95, 94.9080385d0),
     &  Iso(97, 96.9109512d0),
     &  Iso(98, 97.9127289d0),
     &  Iso(99, 98.916667d0),
     &  Iso(100, 99.9180006d0),
     &  Iso(101, 100.921448d0),
     &  Iso(102, 101.9231409d0),
     &  Iso(103, 102.927191d0),
     &  Iso(104, 103.929436d0),
     &  Iso(105, 104.934008d0),
     &  Iso(106, 105.93676d0),
     &  Iso(107, 106.94174d0),
     &  Iso(108, 107.94487d0),
     &  Iso(109, 108.95041d0),
     &  Iso(110, 109.95396d0),
     &  Iso(111, 110.95968d0),
     &  Iso(112, 111.9637d0) ]

      ElementList(41)%Symbol = 'Nb'
      ElementList(41)%nat = 1
      Call mma_Allocate(ElementList(41)%Isotopes,35)
      ElementList(41)%Isotopes(:) = [
     &  Iso(93, 92.906373d0),
     &  Iso(81, 80.9496d0),
     &  Iso(82, 81.94396d0),
     &  Iso(83, 82.93729d0),
     &  Iso(84, 83.93449d0),
     &  Iso(85, 84.9288458d0),
     &  Iso(86, 85.9257828d0),
     &  Iso(87, 86.9206937d0),
     &  Iso(88, 87.918222d0),
     &  Iso(89, 88.913445d0),
     &  Iso(90, 89.9112584d0),
     &  Iso(91, 90.9069897d0),
     &  Iso(92, 91.9071881d0),
     &  Iso(94, 93.9072788d0),
     &  Iso(95, 94.9068324d0),
     &  Iso(96, 95.9080973d0),
     &  Iso(97, 96.9080959d0),
     &  Iso(98, 97.9103265d0),
     &  Iso(99, 98.911613d0),
     &  Iso(100, 99.9143276d0),
     &  Iso(101, 100.9153103d0),
     &  Iso(102, 101.9180772d0),
     &  Iso(103, 102.9194572d0),
     &  Iso(104, 103.9228925d0),
     &  Iso(105, 104.9249465d0),
     &  Iso(106, 105.9289317d0),
     &  Iso(107, 106.9315937d0),
     &  Iso(108, 107.9360748d0),
     &  Iso(109, 108.93922d0),
     &  Iso(110, 109.94403d0),
     &  Iso(111, 110.94753d0),
     &  Iso(112, 111.95247d0),
     &  Iso(113, 112.95651d0),
     &  Iso(114, 113.96201d0),
     &  Iso(115, 114.96634d0) ]

      ElementList(42)%Symbol = 'Mo'
      ElementList(42)%nat = 7
      Call mma_Allocate(ElementList(42)%Isotopes,35)
      ElementList(42)%Isotopes(:) = [
     &  Iso(98, 97.90540482d0),
     &  Iso(96, 95.90467612d0),
     &  Iso(95, 94.90583877d0),
     &  Iso(92, 91.90680796d0),
     &  Iso(100, 99.9074718d0),
     &  Iso(97, 96.90601812d0),
     &  Iso(94, 93.9050849d0),
     &  Iso(83, 82.94988d0),
     &  Iso(84, 83.94149d0),
     &  Iso(85, 84.938261d0),
     &  Iso(86, 85.9311748d0),
     &  Iso(87, 86.9281962d0),
     &  Iso(88, 87.9219678d0),
     &  Iso(89, 88.9194682d0),
     &  Iso(90, 89.9139309d0),
     &  Iso(91, 90.9117453d0),
     &  Iso(93, 92.90680958d0),
     &  Iso(99, 98.90770851d0),
     &  Iso(101, 100.9103414d0),
     &  Iso(102, 101.9102834d0),
     &  Iso(103, 102.913079d0),
     &  Iso(104, 103.9137344d0),
     &  Iso(105, 104.916969d0),
     &  Iso(106, 105.918259d0),
     &  Iso(107, 106.922106d0),
     &  Iso(108, 107.924033d0),
     &  Iso(109, 108.928424d0),
     &  Iso(110, 109.930704d0),
     &  Iso(111, 110.935654d0),
     &  Iso(112, 111.93831d0),
     &  Iso(113, 112.94335d0),
     &  Iso(114, 113.94653d0),
     &  Iso(115, 114.95196d0),
     &  Iso(116, 115.95545d0),
     &  Iso(117, 116.96117d0) ]

      ElementList(43)%Symbol = 'Tc'
      ElementList(43)%nat = 0
      Call mma_Allocate(ElementList(43)%Isotopes,36)
      ElementList(43)%Isotopes(:) = [
     &  Iso(98, 97.9072124d0),
     &  Iso(85, 84.95058d0),
     &  Iso(86, 85.94493d0),
     &  Iso(87, 86.9380672d0),
     &  Iso(88, 87.93378d0),
     &  Iso(89, 88.9276487d0),
     &  Iso(90, 89.9240739d0),
     &  Iso(91, 90.9184254d0),
     &  Iso(92, 91.9152698d0),
     &  Iso(93, 92.910246d0),
     &  Iso(94, 93.9096536d0),
     &  Iso(95, 94.9076536d0),
     &  Iso(96, 95.907868d0),
     &  Iso(97, 96.9063667d0),
     &  Iso(99, 98.9062508d0),
     &  Iso(100, 99.9076539d0),
     &  Iso(101, 100.907309d0),
     &  Iso(102, 101.9092097d0),
     &  Iso(103, 102.909176d0),
     &  Iso(104, 103.911425d0),
     &  Iso(105, 104.911655d0),
     &  Iso(106, 105.914358d0),
     &  Iso(107, 106.9154606d0),
     &  Iso(108, 107.9184957d0),
     &  Iso(109, 108.920256d0),
     &  Iso(110, 109.923744d0),
     &  Iso(111, 110.925901d0),
     &  Iso(112, 111.9299458d0),
     &  Iso(113, 112.932569d0),
     &  Iso(114, 113.93691d0),
     &  Iso(115, 114.93998d0),
     &  Iso(116, 115.94476d0),
     &  Iso(117, 116.94806d0),
     &  Iso(118, 117.95299d0),
     &  Iso(119, 118.95666d0),
     &  Iso(120, 119.96187d0) ]

      ElementList(44)%Symbol = 'Ru'
      ElementList(44)%nat = 7
      Call mma_Allocate(ElementList(44)%Isotopes,38)
      ElementList(44)%Isotopes(:) = [
     &  Iso(102, 101.9043441d0),
     &  Iso(104, 103.9054275d0),
     &  Iso(101, 100.9055769d0),
     &  Iso(99, 98.9059341d0),
     &  Iso(100, 99.9042143d0),
     &  Iso(96, 95.90759025d0),
     &  Iso(98, 97.9052868d0),
     &  Iso(87, 86.95069d0),
     &  Iso(88, 87.9416d0),
     &  Iso(89, 88.93762d0),
     &  Iso(90, 89.9303444d0),
     &  Iso(91, 90.9267419d0),
     &  Iso(92, 91.9202344d0),
     &  Iso(93, 92.9171044d0),
     &  Iso(94, 93.9113429d0),
     &  Iso(95, 94.910406d0),
     &  Iso(97, 96.9075471d0),
     &  Iso(103, 102.9063186d0),
     &  Iso(105, 104.9077476d0),
     &  Iso(106, 105.9073291d0),
     &  Iso(107, 106.909972d0),
     &  Iso(108, 107.910188d0),
     &  Iso(109, 108.913326d0),
     &  Iso(110, 109.9140407d0),
     &  Iso(111, 110.91757d0),
     &  Iso(112, 111.918809d0),
     &  Iso(113, 112.922844d0),
     &  Iso(114, 113.9246136d0),
     &  Iso(115, 114.92882d0),
     &  Iso(116, 115.9312192d0),
     &  Iso(117, 116.9361d0),
     &  Iso(118, 117.93853d0),
     &  Iso(119, 118.94357d0),
     &  Iso(120, 119.94631d0),
     &  Iso(121, 120.95164d0),
     &  Iso(122, 121.95447d0),
     &  Iso(123, 122.95989d0),
     &  Iso(124, 123.96305d0) ]

      ElementList(45)%Symbol = 'Rh'
      ElementList(45)%nat = 1
      Call mma_Allocate(ElementList(45)%Isotopes,38)
      ElementList(45)%Isotopes(:) = [
     &  Iso(103, 102.905498d0),
     &  Iso(89, 88.95058d0),
     &  Iso(90, 89.94422d0),
     &  Iso(91, 90.93688d0),
     &  Iso(92, 91.9323677d0),
     &  Iso(93, 92.9259128d0),
     &  Iso(94, 93.9217305d0),
     &  Iso(95, 94.9158979d0),
     &  Iso(96, 95.914453d0),
     &  Iso(97, 96.911329d0),
     &  Iso(98, 97.910708d0),
     &  Iso(99, 98.9081282d0),
     &  Iso(100, 99.908117d0),
     &  Iso(101, 100.9061606d0),
     &  Iso(102, 101.9068374d0),
     &  Iso(104, 103.9066492d0),
     &  Iso(105, 104.9056885d0),
     &  Iso(106, 105.9072868d0),
     &  Iso(107, 106.906748d0),
     &  Iso(108, 107.908714d0),
     &  Iso(109, 108.9087488d0),
     &  Iso(110, 109.911079d0),
     &  Iso(111, 110.9116423d0),
     &  Iso(112, 111.914403d0),
     &  Iso(113, 112.9154393d0),
     &  Iso(114, 113.918718d0),
     &  Iso(115, 114.9203116d0),
     &  Iso(116, 115.924059d0),
     &  Iso(117, 116.9260354d0),
     &  Iso(118, 117.93034d0),
     &  Iso(119, 118.932557d0),
     &  Iso(120, 119.93686d0),
     &  Iso(121, 120.93942d0),
     &  Iso(122, 121.94399d0),
     &  Iso(123, 122.94685d0),
     &  Iso(124, 123.95151d0),
     &  Iso(125, 124.95469d0),
     &  Iso(126, 125.95946d0) ]

      ElementList(46)%Symbol = 'Pd'
      ElementList(46)%nat = 6
      Call mma_Allocate(ElementList(46)%Isotopes,38)
      ElementList(46)%Isotopes(:) = [
     &  Iso(106, 105.9034804d0),
     &  Iso(108, 107.9038916d0),
     &  Iso(105, 104.9050796d0),
     &  Iso(110, 109.9051722d0),
     &  Iso(104, 103.9040305d0),
     &  Iso(102, 101.9056022d0),
     &  Iso(91, 90.95032d0),
     &  Iso(92, 91.94088d0),
     &  Iso(93, 92.93651d0),
     &  Iso(94, 93.9290376d0),
     &  Iso(95, 94.9248898d0),
     &  Iso(96, 95.9182151d0),
     &  Iso(97, 96.916472d0),
     &  Iso(98, 97.9126983d0),
     &  Iso(99, 98.9117748d0),
     &  Iso(100, 99.908505d0),
     &  Iso(101, 100.9082864d0),
     &  Iso(103, 102.9060809d0),
     &  Iso(107, 106.9051282d0),
     &  Iso(109, 108.9059504d0),
     &  Iso(111, 110.90768968d0),
     &  Iso(112, 111.9073297d0),
     &  Iso(113, 112.910261d0),
     &  Iso(114, 113.9103686d0),
     &  Iso(115, 114.913659d0),
     &  Iso(116, 115.914297d0),
     &  Iso(117, 116.9179547d0),
     &  Iso(118, 117.9190667d0),
     &  Iso(119, 118.9233402d0),
     &  Iso(120, 119.9245511d0),
     &  Iso(121, 120.9289503d0),
     &  Iso(122, 121.930632d0),
     &  Iso(123, 122.93514d0),
     &  Iso(124, 123.93714d0),
     &  Iso(125, 124.94179d0),
     &  Iso(126, 125.94416d0),
     &  Iso(127, 126.94907d0),
     &  Iso(128, 127.95183d0) ]

      ElementList(47)%Symbol = 'Ag'
      ElementList(47)%nat = 2
      Call mma_Allocate(ElementList(47)%Isotopes,38)
      ElementList(47)%Isotopes(:) = [
     &  Iso(107, 106.9050916d0),
     &  Iso(109, 108.9047553d0),
     &  Iso(93, 92.95033d0),
     &  Iso(94, 93.94373d0),
     &  Iso(95, 94.93602d0),
     &  Iso(96, 95.930744d0),
     &  Iso(97, 96.92397d0),
     &  Iso(98, 97.92156d0),
     &  Iso(99, 98.9176458d0),
     &  Iso(100, 99.9161154d0),
     &  Iso(101, 100.912684d0),
     &  Iso(102, 101.9117047d0),
     &  Iso(103, 102.9089631d0),
     &  Iso(104, 103.9086239d0),
     &  Iso(105, 104.9065256d0),
     &  Iso(106, 105.9066636d0),
     &  Iso(108, 107.9059503d0),
     &  Iso(110, 109.9061102d0),
     &  Iso(111, 110.9052959d0),
     &  Iso(112, 111.9070486d0),
     &  Iso(113, 112.906573d0),
     &  Iso(114, 113.908823d0),
     &  Iso(115, 114.908767d0),
     &  Iso(116, 115.9113868d0),
     &  Iso(117, 116.911774d0),
     &  Iso(118, 117.9145955d0),
     &  Iso(119, 118.91557d0),
     &  Iso(120, 119.9187848d0),
     &  Iso(121, 120.920125d0),
     &  Iso(122, 121.923664d0),
     &  Iso(123, 122.925337d0),
     &  Iso(124, 123.92893d0),
     &  Iso(125, 124.93105d0),
     &  Iso(126, 125.93475d0),
     &  Iso(127, 126.93711d0),
     &  Iso(128, 127.94106d0),
     &  Iso(129, 128.94395d0),
     &  Iso(130, 129.9507d0) ]

      ElementList(48)%Symbol = 'Cd'
      ElementList(48)%nat = 8
      Call mma_Allocate(ElementList(48)%Isotopes,39)
      ElementList(48)%Isotopes(:) = [
     &  Iso(114, 113.90336509d0),
     &  Iso(112, 111.90276287d0),
     &  Iso(111, 110.90418287d0),
     &  Iso(110, 109.90300661d0),
     &  Iso(113, 112.90440813d0),
     &  Iso(116, 115.90476315d0),
     &  Iso(106, 105.9064599d0),
     &  Iso(108, 107.9041834d0),
     &  Iso(95, 94.94994d0),
     &  Iso(96, 95.94034d0),
     &  Iso(97, 96.9351d0),
     &  Iso(98, 97.927389d0),
     &  Iso(99, 98.9249258d0),
     &  Iso(100, 99.9203488d0),
     &  Iso(101, 100.9185862d0),
     &  Iso(102, 101.914482d0),
     &  Iso(103, 102.9134165d0),
     &  Iso(104, 103.9098564d0),
     &  Iso(105, 104.9094639d0),
     &  Iso(107, 106.9066121d0),
     &  Iso(109, 108.9049867d0),
     &  Iso(115, 114.90543751d0),
     &  Iso(117, 116.907226d0),
     &  Iso(118, 117.906922d0),
     &  Iso(119, 118.909847d0),
     &  Iso(120, 119.9098681d0),
     &  Iso(121, 120.9129637d0),
     &  Iso(122, 121.9134591d0),
     &  Iso(123, 122.9168925d0),
     &  Iso(124, 123.9176574d0),
     &  Iso(125, 124.9212576d0),
     &  Iso(126, 125.9224291d0),
     &  Iso(127, 126.926472d0),
     &  Iso(128, 127.9278129d0),
     &  Iso(129, 128.93182d0),
     &  Iso(130, 129.93394d0),
     &  Iso(131, 130.9406d0),
     &  Iso(132, 131.94604d0),
     &  Iso(133, 132.95285d0) ]

      ElementList(49)%Symbol = 'In'
      ElementList(49)%nat = 2
      Call mma_Allocate(ElementList(49)%Isotopes,39)
      ElementList(49)%Isotopes(:) = [
     &  Iso(115, 114.903878776d0),
     &  Iso(113, 112.90406184d0),
     &  Iso(97, 96.94934d0),
     &  Iso(98, 97.94214d0),
     &  Iso(99, 98.93411d0),
     &  Iso(100, 99.93096d0),
     &  Iso(101, 100.92634d0),
     &  Iso(102, 101.9241071d0),
     &  Iso(103, 102.9198819d0),
     &  Iso(104, 103.9182145d0),
     &  Iso(105, 104.914502d0),
     &  Iso(106, 105.913464d0),
     &  Iso(107, 106.91029d0),
     &  Iso(108, 107.9096935d0),
     &  Iso(109, 108.9071514d0),
     &  Iso(110, 109.90717d0),
     &  Iso(111, 110.9051085d0),
     &  Iso(112, 111.9055377d0),
     &  Iso(114, 113.90491791d0),
     &  Iso(116, 115.90525999d0),
     &  Iso(117, 116.9045157d0),
     &  Iso(118, 117.9063566d0),
     &  Iso(119, 118.9058507d0),
     &  Iso(120, 119.907967d0),
     &  Iso(121, 120.907851d0),
     &  Iso(122, 121.910281d0),
     &  Iso(123, 122.910434d0),
     &  Iso(124, 123.913182d0),
     &  Iso(125, 124.913605d0),
     &  Iso(126, 125.916507d0),
     &  Iso(127, 126.917446d0),
     &  Iso(128, 127.9204d0),
     &  Iso(129, 128.9218053d0),
     &  Iso(130, 129.924977d0),
     &  Iso(131, 130.9269715d0),
     &  Iso(132, 131.933001d0),
     &  Iso(133, 132.93831d0),
     &  Iso(134, 133.94454d0),
     &  Iso(135, 134.95005d0) ]

      ElementList(50)%Symbol = 'Sn'
      ElementList(50)%nat = 10
      Call mma_Allocate(ElementList(50)%Isotopes,40)
      ElementList(50)%Isotopes(:) = [
     &  Iso(120, 119.90220163d0),
     &  Iso(118, 117.90160657d0),
     &  Iso(116, 115.9017428d0),
     &  Iso(119, 118.90331117d0),
     &  Iso(117, 116.90295398d0),
     &  Iso(124, 123.9052766d0),
     &  Iso(122, 121.9034438d0),
     &  Iso(112, 111.90482387d0),
     &  Iso(114, 113.9027827d0),
     &  Iso(115, 114.903344699d0),
     &  Iso(99, 98.94853d0),
     &  Iso(100, 99.9385d0),
     &  Iso(101, 100.93526d0),
     &  Iso(102, 101.93029d0),
     &  Iso(103, 102.928105d0),
     &  Iso(104, 103.9231052d0),
     &  Iso(105, 104.9212684d0),
     &  Iso(106, 105.9169574d0),
     &  Iso(107, 106.9157137d0),
     &  Iso(108, 107.9118943d0),
     &  Iso(109, 108.9112921d0),
     &  Iso(110, 109.907845d0),
     &  Iso(111, 110.9077401d0),
     &  Iso(113, 112.9051757d0),
     &  Iso(121, 120.9042426d0),
     &  Iso(123, 122.9057252d0),
     &  Iso(125, 124.9077864d0),
     &  Iso(126, 125.907659d0),
     &  Iso(127, 126.91039d0),
     &  Iso(128, 127.910507d0),
     &  Iso(129, 128.913465d0),
     &  Iso(130, 129.9139738d0),
     &  Iso(131, 130.917045d0),
     &  Iso(132, 131.9178267d0),
     &  Iso(133, 132.9239134d0),
     &  Iso(134, 133.9286821d0),
     &  Iso(135, 134.9349086d0),
     &  Iso(136, 135.93999d0),
     &  Iso(137, 136.94655d0),
     &  Iso(138, 137.95184d0) ]

      ElementList(51)%Symbol = 'Sb'
      ElementList(51)%nat = 2
      Call mma_Allocate(ElementList(51)%Isotopes,38)
      ElementList(51)%Isotopes(:) = [
     &  Iso(121, 120.903812d0),
     &  Iso(123, 122.9042132d0),
     &  Iso(103, 102.93969d0),
     &  Iso(104, 103.93648d0),
     &  Iso(105, 104.931276d0),
     &  Iso(106, 105.928638d0),
     &  Iso(107, 106.9241506d0),
     &  Iso(108, 107.9222267d0),
     &  Iso(109, 108.9181411d0),
     &  Iso(110, 109.9168543d0),
     &  Iso(111, 110.9132182d0),
     &  Iso(112, 111.9124d0),
     &  Iso(113, 112.909375d0),
     &  Iso(114, 113.90929d0),
     &  Iso(115, 114.906598d0),
     &  Iso(116, 115.9067931d0),
     &  Iso(117, 116.9048415d0),
     &  Iso(118, 117.9055321d0),
     &  Iso(119, 118.9039455d0),
     &  Iso(120, 119.9050794d0),
     &  Iso(122, 121.9051699d0),
     &  Iso(124, 123.905935d0),
     &  Iso(125, 124.905253d0),
     &  Iso(126, 125.907253d0),
     &  Iso(127, 126.9069243d0),
     &  Iso(128, 127.909146d0),
     &  Iso(129, 128.909147d0),
     &  Iso(130, 129.911662d0),
     &  Iso(131, 130.9119888d0),
     &  Iso(132, 131.9145077d0),
     &  Iso(133, 132.9152732d0),
     &  Iso(134, 133.9205357d0),
     &  Iso(135, 134.9251851d0),
     &  Iso(136, 135.9307459d0),
     &  Iso(137, 136.93555d0),
     &  Iso(138, 137.94145d0),
     &  Iso(139, 138.94655d0),
     &  Iso(140, 139.95283d0) ]

      ElementList(52)%Symbol = 'Te'
      ElementList(52)%nat = 8
      Call mma_Allocate(ElementList(52)%Isotopes,39)
      ElementList(52)%Isotopes(:) = [
     &  Iso(130, 129.906222748d0),
     &  Iso(128, 127.90446128d0),
     &  Iso(126, 125.9033109d0),
     &  Iso(125, 124.9044299d0),
     &  Iso(124, 123.9028171d0),
     &  Iso(122, 121.9030435d0),
     &  Iso(123, 122.9042698d0),
     &  Iso(120, 119.9040593d0),
     &  Iso(105, 104.9433d0),
     &  Iso(106, 105.9375d0),
     &  Iso(107, 106.935012d0),
     &  Iso(108, 107.9293805d0),
     &  Iso(109, 108.9273045d0),
     &  Iso(110, 109.9224581d0),
     &  Iso(111, 110.9210006d0),
     &  Iso(112, 111.9167279d0),
     &  Iso(113, 112.915891d0),
     &  Iso(114, 113.912089d0),
     &  Iso(115, 114.911902d0),
     &  Iso(116, 115.90846d0),
     &  Iso(117, 116.908646d0),
     &  Iso(118, 117.905854d0),
     &  Iso(119, 118.9064071d0),
     &  Iso(121, 120.904944d0),
     &  Iso(127, 126.9052257d0),
     &  Iso(129, 128.90659646d0),
     &  Iso(131, 130.908522213d0),
     &  Iso(132, 131.9085467d0),
     &  Iso(133, 132.9109688d0),
     &  Iso(134, 133.911394d0),
     &  Iso(135, 134.9165557d0),
     &  Iso(136, 135.9201006d0),
     &  Iso(137, 136.9255989d0),
     &  Iso(138, 137.9294722d0),
     &  Iso(139, 138.9353672d0),
     &  Iso(140, 139.939499d0),
     &  Iso(141, 140.9458d0),
     &  Iso(142, 141.95022d0),
     &  Iso(143, 142.95676d0) ]

      ElementList(53)%Symbol = 'I'
      ElementList(53)%nat = 1
      Call mma_Allocate(ElementList(53)%Isotopes,39)
      ElementList(53)%Isotopes(:) = [
     &  Iso(127, 126.9044719d0),
     &  Iso(107, 106.94678d0),
     &  Iso(108, 107.94348d0),
     &  Iso(109, 108.9380853d0),
     &  Iso(110, 109.935089d0),
     &  Iso(111, 110.9302692d0),
     &  Iso(112, 111.928005d0),
     &  Iso(113, 112.9236501d0),
     &  Iso(114, 113.92185d0),
     &  Iso(115, 114.918048d0),
     &  Iso(116, 115.91681d0),
     &  Iso(117, 116.913648d0),
     &  Iso(118, 117.913074d0),
     &  Iso(119, 118.910074d0),
     &  Iso(120, 119.910087d0),
     &  Iso(121, 120.9074051d0),
     &  Iso(122, 121.9075888d0),
     &  Iso(123, 122.9055885d0),
     &  Iso(124, 123.906209d0),
     &  Iso(125, 124.9046294d0),
     &  Iso(126, 125.9056233d0),
     &  Iso(128, 127.9058086d0),
     &  Iso(129, 128.9049837d0),
     &  Iso(130, 129.9066702d0),
     &  Iso(131, 130.9061263d0),
     &  Iso(132, 131.9079935d0),
     &  Iso(133, 132.907797d0),
     &  Iso(134, 133.9097588d0),
     &  Iso(135, 134.9100488d0),
     &  Iso(136, 135.914604d0),
     &  Iso(137, 136.9180282d0),
     &  Iso(138, 137.9227264d0),
     &  Iso(139, 138.926506d0),
     &  Iso(140, 139.93173d0),
     &  Iso(141, 140.93569d0),
     &  Iso(142, 141.9412d0),
     &  Iso(143, 142.94565d0),
     &  Iso(144, 143.95139d0),
     &  Iso(145, 144.95605d0) ]

      ElementList(54)%Symbol = 'Xe'
      ElementList(54)%nat = 9
      Call mma_Allocate(ElementList(54)%Isotopes,40)
      ElementList(54)%Isotopes(:) = [
     &  Iso(132, 131.9041550856d0),
     &  Iso(129, 128.9047808611d0),
     &  Iso(131, 130.90508406d0),
     &  Iso(134, 133.90539466d0),
     &  Iso(136, 135.907214484d0),
     &  Iso(130, 129.903509349d0),
     &  Iso(128, 127.903531d0),
     &  Iso(124, 123.905892d0),
     &  Iso(126, 125.9042983d0),
     &  Iso(109, 108.95043d0),
     &  Iso(110, 109.94426d0),
     &  Iso(111, 110.941607d0),
     &  Iso(112, 111.935559d0),
     &  Iso(113, 112.9332217d0),
     &  Iso(114, 113.92798d0),
     &  Iso(115, 114.926294d0),
     &  Iso(116, 115.921581d0),
     &  Iso(117, 116.920359d0),
     &  Iso(118, 117.916179d0),
     &  Iso(119, 118.915411d0),
     &  Iso(120, 119.911784d0),
     &  Iso(121, 120.911453d0),
     &  Iso(122, 121.908368d0),
     &  Iso(123, 122.908482d0),
     &  Iso(125, 124.9063944d0),
     &  Iso(127, 126.9051829d0),
     &  Iso(133, 132.9059108d0),
     &  Iso(135, 134.9072278d0),
     &  Iso(137, 136.91155778d0),
     &  Iso(138, 137.9141463d0),
     &  Iso(139, 138.9187922d0),
     &  Iso(140, 139.9216458d0),
     &  Iso(141, 140.9267872d0),
     &  Iso(142, 141.9299731d0),
     &  Iso(143, 142.9353696d0),
     &  Iso(144, 143.9389451d0),
     &  Iso(145, 144.94472d0),
     &  Iso(146, 145.948518d0),
     &  Iso(147, 146.95426d0),
     &  Iso(148, 147.95813d0) ]

      ElementList(55)%Symbol = 'Cs'
      ElementList(55)%nat = 1
      Call mma_Allocate(ElementList(55)%Isotopes,40)
      ElementList(55)%Isotopes(:) = [
     &  Iso(133, 132.905451961d0),
     &  Iso(112, 111.950309d0),
     &  Iso(113, 112.9444291d0),
     &  Iso(114, 113.941296d0),
     &  Iso(115, 114.93591d0),
     &  Iso(116, 115.93337d0),
     &  Iso(117, 116.928617d0),
     &  Iso(118, 117.92656d0),
     &  Iso(119, 118.922377d0),
     &  Iso(120, 119.920677d0),
     &  Iso(121, 120.917227d0),
     &  Iso(122, 121.916108d0),
     &  Iso(123, 122.912996d0),
     &  Iso(124, 123.9122578d0),
     &  Iso(125, 124.909728d0),
     &  Iso(126, 125.909446d0),
     &  Iso(127, 126.9074174d0),
     &  Iso(128, 127.9077487d0),
     &  Iso(129, 128.9060657d0),
     &  Iso(130, 129.9067093d0),
     &  Iso(131, 130.9054649d0),
     &  Iso(132, 131.9064339d0),
     &  Iso(134, 133.906718503d0),
     &  Iso(135, 134.905977d0),
     &  Iso(136, 135.9073114d0),
     &  Iso(137, 136.90708923d0),
     &  Iso(138, 137.9110171d0),
     &  Iso(139, 138.9133638d0),
     &  Iso(140, 139.9172831d0),
     &  Iso(141, 140.9200455d0),
     &  Iso(142, 141.924296d0),
     &  Iso(143, 142.927349d0),
     &  Iso(144, 143.932076d0),
     &  Iso(145, 144.935527d0),
     &  Iso(146, 145.940344d0),
     &  Iso(147, 146.944156d0),
     &  Iso(148, 147.94923d0),
     &  Iso(149, 148.95302d0),
     &  Iso(150, 149.95833d0),
     &  Iso(151, 150.96258d0) ]

      ElementList(56)%Symbol = 'Ba'
      ElementList(56)%nat = 7
      Call mma_Allocate(ElementList(56)%Isotopes,40)
      ElementList(56)%Isotopes(:) = [
     &  Iso(138, 137.905247d0),
     &  Iso(137, 136.90582714d0),
     &  Iso(136, 135.90457573d0),
     &  Iso(135, 134.90568838d0),
     &  Iso(134, 133.90450818d0),
     &  Iso(130, 129.9063207d0),
     &  Iso(132, 131.9050611d0),
     &  Iso(114, 113.95066d0),
     &  Iso(115, 114.94737d0),
     &  Iso(116, 115.94128d0),
     &  Iso(117, 116.93814d0),
     &  Iso(118, 117.93306d0),
     &  Iso(119, 118.93066d0),
     &  Iso(120, 119.92605d0),
     &  Iso(121, 120.92405d0),
     &  Iso(122, 121.919904d0),
     &  Iso(123, 122.918781d0),
     &  Iso(124, 123.915094d0),
     &  Iso(125, 124.914472d0),
     &  Iso(126, 125.91125d0),
     &  Iso(127, 126.911091d0),
     &  Iso(128, 127.908342d0),
     &  Iso(129, 128.908681d0),
     &  Iso(131, 130.906941d0),
     &  Iso(133, 132.9060074d0),
     &  Iso(139, 138.9088411d0),
     &  Iso(140, 139.9106057d0),
     &  Iso(141, 140.9144033d0),
     &  Iso(142, 141.9164324d0),
     &  Iso(143, 142.9206253d0),
     &  Iso(144, 143.9229549d0),
     &  Iso(145, 144.9275184d0),
     &  Iso(146, 145.930284d0),
     &  Iso(147, 146.935304d0),
     &  Iso(148, 147.938171d0),
     &  Iso(149, 148.94308d0),
     &  Iso(150, 149.94605d0),
     &  Iso(151, 150.95127d0),
     &  Iso(152, 151.95481d0),
     &  Iso(153, 152.96036d0) ]

      ElementList(57)%Symbol = 'La'
      ElementList(57)%nat = 2
      Call mma_Allocate(ElementList(57)%Isotopes,40)
      ElementList(57)%Isotopes(:) = [
     &  Iso(139, 138.9063563d0),
     &  Iso(138, 137.9071149d0),
     &  Iso(116, 115.9563d0),
     &  Iso(117, 116.94999d0),
     &  Iso(118, 117.94673d0),
     &  Iso(119, 118.94099d0),
     &  Iso(120, 119.93807d0),
     &  Iso(121, 120.93315d0),
     &  Iso(122, 121.93071d0),
     &  Iso(123, 122.9263d0),
     &  Iso(124, 123.924574d0),
     &  Iso(125, 124.920816d0),
     &  Iso(126, 125.919513d0),
     &  Iso(127, 126.916375d0),
     &  Iso(128, 127.915592d0),
     &  Iso(129, 128.912694d0),
     &  Iso(130, 129.912369d0),
     &  Iso(131, 130.91007d0),
     &  Iso(132, 131.910119d0),
     &  Iso(133, 132.908218d0),
     &  Iso(134, 133.908514d0),
     &  Iso(135, 134.906984d0),
     &  Iso(136, 135.907635d0),
     &  Iso(137, 136.9064504d0),
     &  Iso(140, 139.9094806d0),
     &  Iso(141, 140.910966d0),
     &  Iso(142, 141.9140909d0),
     &  Iso(143, 142.9160795d0),
     &  Iso(144, 143.919646d0),
     &  Iso(145, 144.921808d0),
     &  Iso(146, 145.925875d0),
     &  Iso(147, 146.928418d0),
     &  Iso(148, 147.932679d0),
     &  Iso(149, 148.93535d0),
     &  Iso(150, 149.93947d0),
     &  Iso(151, 150.94232d0),
     &  Iso(152, 151.94682d0),
     &  Iso(153, 152.95036d0),
     &  Iso(154, 153.95517d0),
     &  Iso(155, 154.95901d0) ]

      ElementList(58)%Symbol = 'Ce'
      ElementList(58)%nat = 4
      Call mma_Allocate(ElementList(58)%Isotopes,39)
      ElementList(58)%Isotopes(:) = [
     &  Iso(140, 139.9054431d0),
     &  Iso(142, 141.9092504d0),
     &  Iso(138, 137.905991d0),
     &  Iso(136, 135.90712921d0),
     &  Iso(119, 118.95271d0),
     &  Iso(120, 119.94654d0),
     &  Iso(121, 120.94335d0),
     &  Iso(122, 121.93787d0),
     &  Iso(123, 122.93528d0),
     &  Iso(124, 123.93031d0),
     &  Iso(125, 124.92844d0),
     &  Iso(126, 125.923971d0),
     &  Iso(127, 126.922727d0),
     &  Iso(128, 127.918911d0),
     &  Iso(129, 128.918102d0),
     &  Iso(130, 129.914736d0),
     &  Iso(131, 130.914429d0),
     &  Iso(132, 131.911464d0),
     &  Iso(133, 132.91152d0),
     &  Iso(134, 133.908928d0),
     &  Iso(135, 134.909161d0),
     &  Iso(137, 136.90776236d0),
     &  Iso(139, 138.9066551d0),
     &  Iso(141, 140.9082807d0),
     &  Iso(143, 142.9123921d0),
     &  Iso(144, 143.9136529d0),
     &  Iso(145, 144.917265d0),
     &  Iso(146, 145.918802d0),
     &  Iso(147, 146.9226899d0),
     &  Iso(148, 147.924424d0),
     &  Iso(149, 148.928427d0),
     &  Iso(150, 149.930384d0),
     &  Iso(151, 150.934272d0),
     &  Iso(152, 151.9366d0),
     &  Iso(153, 152.94093d0),
     &  Iso(154, 153.9438d0),
     &  Iso(155, 154.94855d0),
     &  Iso(156, 155.95183d0),
     &  Iso(157, 156.95705d0) ]

      ElementList(59)%Symbol = 'Pr'
      ElementList(59)%nat = 1
      Call mma_Allocate(ElementList(59)%Isotopes,39)
      ElementList(59)%Isotopes(:) = [
     &  Iso(141, 140.9076576d0),
     &  Iso(121, 120.95532d0),
     &  Iso(122, 121.95175d0),
     &  Iso(123, 122.94596d0),
     &  Iso(124, 123.94294d0),
     &  Iso(125, 124.9377d0),
     &  Iso(126, 125.93524d0),
     &  Iso(127, 126.93071d0),
     &  Iso(128, 127.928791d0),
     &  Iso(129, 128.925095d0),
     &  Iso(130, 129.92359d0),
     &  Iso(131, 130.920235d0),
     &  Iso(132, 131.919255d0),
     &  Iso(133, 132.916331d0),
     &  Iso(134, 133.915697d0),
     &  Iso(135, 134.913112d0),
     &  Iso(136, 135.912677d0),
     &  Iso(137, 136.9106792d0),
     &  Iso(138, 137.910754d0),
     &  Iso(139, 138.9089408d0),
     &  Iso(140, 139.9090803d0),
     &  Iso(142, 141.9100496d0),
     &  Iso(143, 142.9108228d0),
     &  Iso(144, 143.9133109d0),
     &  Iso(145, 144.9145182d0),
     &  Iso(146, 145.91768d0),
     &  Iso(147, 146.919008d0),
     &  Iso(148, 147.92213d0),
     &  Iso(149, 148.923736d0),
     &  Iso(150, 149.9266765d0),
     &  Iso(151, 150.928309d0),
     &  Iso(152, 151.931553d0),
     &  Iso(153, 152.933904d0),
     &  Iso(154, 153.93753d0),
     &  Iso(155, 154.940509d0),
     &  Iso(156, 155.94464d0),
     &  Iso(157, 156.94789d0),
     &  Iso(158, 157.95241d0),
     &  Iso(159, 158.95589d0) ]

      ElementList(60)%Symbol = 'Nd'
      ElementList(60)%nat = 7
      Call mma_Allocate(ElementList(60)%Isotopes,38)
      ElementList(60)%Isotopes(:) = [
     &  Iso(142, 141.907729d0),
     &  Iso(144, 143.910093d0),
     &  Iso(146, 145.9131226d0),
     &  Iso(143, 142.90982d0),
     &  Iso(145, 144.9125793d0),
     &  Iso(148, 147.9168993d0),
     &  Iso(150, 149.9209022d0),
     &  Iso(124, 123.9522d0),
     &  Iso(125, 124.9489d0),
     &  Iso(126, 125.94311d0),
     &  Iso(127, 126.94038d0),
     &  Iso(128, 127.93525d0),
     &  Iso(129, 128.9331d0),
     &  Iso(130, 129.928506d0),
     &  Iso(131, 130.927248d0),
     &  Iso(132, 131.923321d0),
     &  Iso(133, 132.922348d0),
     &  Iso(134, 133.91879d0),
     &  Iso(135, 134.918181d0),
     &  Iso(136, 135.914976d0),
     &  Iso(137, 136.914562d0),
     &  Iso(138, 137.91195d0),
     &  Iso(139, 138.911954d0),
     &  Iso(140, 139.90955d0),
     &  Iso(141, 140.9096147d0),
     &  Iso(147, 146.9161061d0),
     &  Iso(149, 148.9201548d0),
     &  Iso(151, 150.9238403d0),
     &  Iso(152, 151.924692d0),
     &  Iso(153, 152.927718d0),
     &  Iso(154, 153.92948d0),
     &  Iso(155, 154.9331357d0),
     &  Iso(156, 155.93508d0),
     &  Iso(157, 156.939386d0),
     &  Iso(158, 157.94197d0),
     &  Iso(159, 158.94653d0),
     &  Iso(160, 159.9494d0),
     &  Iso(161, 160.95428d0) ]

      ElementList(61)%Symbol = 'Pm'
      ElementList(61)%nat = 0
      Call mma_Allocate(ElementList(61)%Isotopes,38)
      ElementList(61)%Isotopes(:) = [
     &  Iso(145, 144.9127559d0),
     &  Iso(126, 125.95792d0),
     &  Iso(127, 126.95192d0),
     &  Iso(128, 127.9487d0),
     &  Iso(129, 128.94323d0),
     &  Iso(130, 129.94053d0),
     &  Iso(131, 130.93567d0),
     &  Iso(132, 131.93384d0),
     &  Iso(133, 132.929782d0),
     &  Iso(134, 133.928353d0),
     &  Iso(135, 134.924823d0),
     &  Iso(136, 135.923585d0),
     &  Iso(137, 136.92048d0),
     &  Iso(138, 137.919548d0),
     &  Iso(139, 138.9168d0),
     &  Iso(140, 139.91604d0),
     &  Iso(141, 140.913555d0),
     &  Iso(142, 141.91289d0),
     &  Iso(143, 142.9109383d0),
     &  Iso(144, 143.9125964d0),
     &  Iso(146, 145.9147024d0),
     &  Iso(147, 146.915145d0),
     &  Iso(148, 147.9174819d0),
     &  Iso(149, 148.9183423d0),
     &  Iso(150, 149.920991d0),
     &  Iso(151, 150.9212175d0),
     &  Iso(152, 151.923506d0),
     &  Iso(153, 152.9241567d0),
     &  Iso(154, 153.926472d0),
     &  Iso(155, 154.928137d0),
     &  Iso(156, 155.9311175d0),
     &  Iso(157, 156.9331214d0),
     &  Iso(158, 157.936565d0),
     &  Iso(159, 158.939287d0),
     &  Iso(160, 159.9431d0),
     &  Iso(161, 160.94607d0),
     &  Iso(162, 161.95022d0),
     &  Iso(163, 162.95357d0) ]

      ElementList(62)%Symbol = 'Sm'
      ElementList(62)%nat = 7
      Call mma_Allocate(ElementList(62)%Isotopes,38)
      ElementList(62)%Isotopes(:) = [
     &  Iso(152, 151.9197397d0),
     &  Iso(154, 153.9222169d0),
     &  Iso(147, 146.9149044d0),
     &  Iso(149, 148.9171921d0),
     &  Iso(148, 147.9148292d0),
     &  Iso(150, 149.9172829d0),
     &  Iso(144, 143.9120065d0),
     &  Iso(128, 127.95842d0),
     &  Iso(129, 128.95476d0),
     &  Iso(130, 129.949d0),
     &  Iso(131, 130.94618d0),
     &  Iso(132, 131.94087d0),
     &  Iso(133, 132.93856d0),
     &  Iso(134, 133.93411d0),
     &  Iso(135, 134.93252d0),
     &  Iso(136, 135.928276d0),
     &  Iso(137, 136.926971d0),
     &  Iso(138, 137.923244d0),
     &  Iso(139, 138.922297d0),
     &  Iso(140, 139.918995d0),
     &  Iso(141, 140.9184816d0),
     &  Iso(142, 141.9152044d0),
     &  Iso(143, 142.9146353d0),
     &  Iso(145, 144.9134173d0),
     &  Iso(146, 145.913047d0),
     &  Iso(151, 150.9199398d0),
     &  Iso(153, 152.9221047d0),
     &  Iso(155, 154.9246477d0),
     &  Iso(156, 155.925536d0),
     &  Iso(157, 156.9284187d0),
     &  Iso(158, 157.929951d0),
     &  Iso(159, 158.9332172d0),
     &  Iso(160, 159.9353353d0),
     &  Iso(161, 160.9391602d0),
     &  Iso(162, 161.94146d0),
     &  Iso(163, 162.94555d0),
     &  Iso(164, 163.94836d0),
     &  Iso(165, 164.95297d0) ]

      ElementList(63)%Symbol = 'Eu'
      ElementList(63)%nat = 2
      Call mma_Allocate(ElementList(63)%Isotopes,38)
      ElementList(63)%Isotopes(:) = [
     &  Iso(153, 152.921238d0),
     &  Iso(151, 150.9198578d0),
     &  Iso(130, 129.96369d0),
     &  Iso(131, 130.95784d0),
     &  Iso(132, 131.95467d0),
     &  Iso(133, 132.94929d0),
     &  Iso(134, 133.9464d0),
     &  Iso(135, 134.94187d0),
     &  Iso(136, 135.93962d0),
     &  Iso(137, 136.93546d0),
     &  Iso(138, 137.933709d0),
     &  Iso(139, 138.929792d0),
     &  Iso(140, 139.928088d0),
     &  Iso(141, 140.924932d0),
     &  Iso(142, 141.923442d0),
     &  Iso(143, 142.920299d0),
     &  Iso(144, 143.91882d0),
     &  Iso(145, 144.9162726d0),
     &  Iso(146, 145.917211d0),
     &  Iso(147, 146.9167527d0),
     &  Iso(148, 147.918089d0),
     &  Iso(149, 148.9179378d0),
     &  Iso(150, 149.9197077d0),
     &  Iso(152, 151.9217522d0),
     &  Iso(154, 153.922987d0),
     &  Iso(155, 154.9229011d0),
     &  Iso(156, 155.9247605d0),
     &  Iso(157, 156.9254334d0),
     &  Iso(158, 157.927799d0),
     &  Iso(159, 158.9291001d0),
     &  Iso(160, 159.931851d0),
     &  Iso(161, 160.933664d0),
     &  Iso(162, 161.936989d0),
     &  Iso(163, 162.939196d0),
     &  Iso(164, 163.94274d0),
     &  Iso(165, 164.94559d0),
     &  Iso(166, 165.94962d0),
     &  Iso(167, 166.95289d0) ]

      ElementList(64)%Symbol = 'Gd'
      ElementList(64)%nat = 7
      Call mma_Allocate(ElementList(64)%Isotopes,37)
      ElementList(64)%Isotopes(:) = [
     &  Iso(158, 157.9241123d0),
     &  Iso(160, 159.9270624d0),
     &  Iso(156, 155.9221312d0),
     &  Iso(157, 156.9239686d0),
     &  Iso(155, 154.9226305d0),
     &  Iso(154, 153.9208741d0),
     &  Iso(152, 151.9197995d0),
     &  Iso(133, 132.96133d0),
     &  Iso(134, 133.95566d0),
     &  Iso(135, 134.95245d0),
     &  Iso(136, 135.9473d0),
     &  Iso(137, 136.94502d0),
     &  Iso(138, 137.94025d0),
     &  Iso(139, 138.93813d0),
     &  Iso(140, 139.933674d0),
     &  Iso(141, 140.932126d0),
     &  Iso(142, 141.928116d0),
     &  Iso(143, 142.92675d0),
     &  Iso(144, 143.922963d0),
     &  Iso(145, 144.921713d0),
     &  Iso(146, 145.9183188d0),
     &  Iso(147, 146.9191014d0),
     &  Iso(148, 147.9181215d0),
     &  Iso(149, 148.9193481d0),
     &  Iso(150, 149.9186644d0),
     &  Iso(151, 150.920356d0),
     &  Iso(153, 152.921758d0),
     &  Iso(159, 158.926397d0),
     &  Iso(161, 160.9296775d0),
     &  Iso(162, 161.930993d0),
     &  Iso(163, 162.9341769d0),
     &  Iso(164, 163.93583d0),
     &  Iso(165, 164.93936d0),
     &  Iso(166, 165.94146d0),
     &  Iso(167, 166.94545d0),
     &  Iso(168, 167.94808d0),
     &  Iso(169, 168.9526d0) ]

      ElementList(65)%Symbol = 'Tb'
      ElementList(65)%nat = 1
      Call mma_Allocate(ElementList(65)%Isotopes,37)
      ElementList(65)%Isotopes(:) = [
     &  Iso(159, 158.9253547d0),
     &  Iso(135, 134.96476d0),
     &  Iso(136, 135.96129d0),
     &  Iso(137, 136.95602d0),
     &  Iso(138, 137.95312d0),
     &  Iso(139, 138.94833d0),
     &  Iso(140, 139.94581d0),
     &  Iso(141, 140.94145d0),
     &  Iso(142, 141.93928d0),
     &  Iso(143, 142.935137d0),
     &  Iso(144, 143.933045d0),
     &  Iso(145, 144.92882d0),
     &  Iso(146, 145.927253d0),
     &  Iso(147, 146.9240548d0),
     &  Iso(148, 147.924282d0),
     &  Iso(149, 148.9232535d0),
     &  Iso(150, 149.9236649d0),
     &  Iso(151, 150.9231096d0),
     &  Iso(152, 151.924083d0),
     &  Iso(153, 152.9234424d0),
     &  Iso(154, 153.924685d0),
     &  Iso(155, 154.923511d0),
     &  Iso(156, 155.9247552d0),
     &  Iso(157, 156.924033d0),
     &  Iso(158, 157.9254209d0),
     &  Iso(160, 159.9271756d0),
     &  Iso(161, 160.9275778d0),
     &  Iso(162, 161.929495d0),
     &  Iso(163, 162.9306547d0),
     &  Iso(164, 163.93336d0),
     &  Iso(165, 164.93498d0),
     &  Iso(166, 165.93786d0),
     &  Iso(167, 166.93996d0),
     &  Iso(168, 167.9434d0),
     &  Iso(169, 168.94597d0),
     &  Iso(170, 169.94984d0),
     &  Iso(171, 170.95273d0) ]

      ElementList(66)%Symbol = 'Dy'
      ElementList(66)%nat = 7
      Call mma_Allocate(ElementList(66)%Isotopes,36)
      ElementList(66)%Isotopes(:) = [
     &  Iso(164, 163.9291819d0),
     &  Iso(162, 161.9268056d0),
     &  Iso(163, 162.9287383d0),
     &  Iso(161, 160.9269405d0),
     &  Iso(160, 159.9252046d0),
     &  Iso(158, 157.9244159d0),
     &  Iso(156, 155.9242847d0),
     &  Iso(138, 137.9625d0),
     &  Iso(139, 138.95959d0),
     &  Iso(140, 139.95402d0),
     &  Iso(141, 140.95128d0),
     &  Iso(142, 141.94619d0),
     &  Iso(143, 142.943994d0),
     &  Iso(144, 143.9392695d0),
     &  Iso(145, 144.937474d0),
     &  Iso(146, 145.9328445d0),
     &  Iso(147, 146.9310827d0),
     &  Iso(148, 147.927157d0),
     &  Iso(149, 148.927322d0),
     &  Iso(150, 149.9255933d0),
     &  Iso(151, 150.9261916d0),
     &  Iso(152, 151.9247253d0),
     &  Iso(153, 152.9257724d0),
     &  Iso(154, 153.9244293d0),
     &  Iso(155, 154.925759d0),
     &  Iso(157, 156.9254707d0),
     &  Iso(159, 158.925747d0),
     &  Iso(165, 164.9317105d0),
     &  Iso(166, 165.9328139d0),
     &  Iso(167, 166.935661d0),
     &  Iso(168, 167.93713d0),
     &  Iso(169, 168.94031d0),
     &  Iso(170, 169.94239d0),
     &  Iso(171, 170.94612d0),
     &  Iso(172, 171.94846d0),
     &  Iso(173, 172.95283d0) ]

      ElementList(67)%Symbol = 'Ho'
      ElementList(67)%nat = 1
      Call mma_Allocate(ElementList(67)%Isotopes,36)
      ElementList(67)%Isotopes(:) = [
     &  Iso(165, 164.9303288d0),
     &  Iso(140, 139.96859d0),
     &  Iso(141, 140.96311d0),
     &  Iso(142, 141.96001d0),
     &  Iso(143, 142.95486d0),
     &  Iso(144, 143.9521097d0),
     &  Iso(145, 144.9472674d0),
     &  Iso(146, 145.9449935d0),
     &  Iso(147, 146.9401423d0),
     &  Iso(148, 147.937744d0),
     &  Iso(149, 148.933803d0),
     &  Iso(150, 149.933498d0),
     &  Iso(151, 150.9316983d0),
     &  Iso(152, 151.931724d0),
     &  Iso(153, 152.9302064d0),
     &  Iso(154, 153.9306068d0),
     &  Iso(155, 154.929104d0),
     &  Iso(156, 155.929706d0),
     &  Iso(157, 156.928254d0),
     &  Iso(158, 157.928946d0),
     &  Iso(159, 158.9277197d0),
     &  Iso(160, 159.928737d0),
     &  Iso(161, 160.9278615d0),
     &  Iso(162, 161.9291023d0),
     &  Iso(163, 162.928741d0),
     &  Iso(164, 163.9302403d0),
     &  Iso(166, 165.9322909d0),
     &  Iso(167, 166.9331385d0),
     &  Iso(168, 167.935522d0),
     &  Iso(169, 168.936878d0),
     &  Iso(170, 169.939625d0),
     &  Iso(171, 170.94147d0),
     &  Iso(172, 171.94473d0),
     &  Iso(173, 172.94702d0),
     &  Iso(174, 173.95095d0),
     &  Iso(175, 174.95362d0) ]

      ElementList(68)%Symbol = 'Er'
      ElementList(68)%nat = 6
      Call mma_Allocate(ElementList(68)%Isotopes,36)
      ElementList(68)%Isotopes(:) = [
     &  Iso(166, 165.9302995d0),
     &  Iso(168, 167.9323767d0),
     &  Iso(167, 166.9320546d0),
     &  Iso(170, 169.9354702d0),
     &  Iso(164, 163.9292088d0),
     &  Iso(162, 161.9287884d0),
     &  Iso(142, 141.9701d0),
     &  Iso(143, 142.96662d0),
     &  Iso(144, 143.9607d0),
     &  Iso(145, 144.95805d0),
     &  Iso(146, 145.9524184d0),
     &  Iso(147, 146.949964d0),
     &  Iso(148, 147.944735d0),
     &  Iso(149, 148.942306d0),
     &  Iso(150, 149.937916d0),
     &  Iso(151, 150.937449d0),
     &  Iso(152, 151.935057d0),
     &  Iso(153, 152.93508d0),
     &  Iso(154, 153.9327908d0),
     &  Iso(155, 154.9332159d0),
     &  Iso(156, 155.931067d0),
     &  Iso(157, 156.931949d0),
     &  Iso(158, 157.929893d0),
     &  Iso(159, 158.9306918d0),
     &  Iso(160, 159.929077d0),
     &  Iso(161, 160.9300046d0),
     &  Iso(163, 162.9300408d0),
     &  Iso(165, 164.9307345d0),
     &  Iso(169, 168.9345968d0),
     &  Iso(171, 170.9380357d0),
     &  Iso(172, 171.9393619d0),
     &  Iso(173, 172.9424d0),
     &  Iso(174, 173.94423d0),
     &  Iso(175, 174.94777d0),
     &  Iso(176, 175.94994d0),
     &  Iso(177, 176.95399d0) ]

      ElementList(69)%Symbol = 'Tm'
      ElementList(69)%nat = 1
      Call mma_Allocate(ElementList(69)%Isotopes,36)
      ElementList(69)%Isotopes(:) = [
     &  Iso(169, 168.9342179d0),
     &  Iso(144, 143.97628d0),
     &  Iso(145, 144.97039d0),
     &  Iso(146, 145.96684d0),
     &  Iso(147, 146.9613799d0),
     &  Iso(148, 147.958384d0),
     &  Iso(149, 148.95289d0),
     &  Iso(150, 149.95009d0),
     &  Iso(151, 150.945488d0),
     &  Iso(152, 151.944422d0),
     &  Iso(153, 152.94204d0),
     &  Iso(154, 153.94157d0),
     &  Iso(155, 154.93921d0),
     &  Iso(156, 155.938992d0),
     &  Iso(157, 156.936944d0),
     &  Iso(158, 157.93698d0),
     &  Iso(159, 158.934975d0),
     &  Iso(160, 159.935263d0),
     &  Iso(161, 160.933549d0),
     &  Iso(162, 161.934002d0),
     &  Iso(163, 162.9326592d0),
     &  Iso(164, 163.933544d0),
     &  Iso(165, 164.9324431d0),
     &  Iso(166, 165.933561d0),
     &  Iso(167, 166.9328562d0),
     &  Iso(168, 167.9341774d0),
     &  Iso(170, 169.935806d0),
     &  Iso(171, 170.9364339d0),
     &  Iso(172, 171.9384055d0),
     &  Iso(173, 172.9396084d0),
     &  Iso(174, 173.942173d0),
     &  Iso(175, 174.943841d0),
     &  Iso(176, 175.947d0),
     &  Iso(177, 176.94904d0),
     &  Iso(178, 177.95264d0),
     &  Iso(179, 178.95534d0) ]

      ElementList(70)%Symbol = 'Yb'
      ElementList(70)%nat = 7
      Call mma_Allocate(ElementList(70)%Isotopes,34)
      ElementList(70)%Isotopes(:) = [
     &  Iso(174, 173.9388664d0),
     &  Iso(172, 171.9363859d0),
     &  Iso(173, 172.9382151d0),
     &  Iso(171, 170.9363302d0),
     &  Iso(176, 175.9425764d0),
     &  Iso(170, 169.9347664d0),
     &  Iso(168, 167.9338896d0),
     &  Iso(148, 147.96758d0),
     &  Iso(149, 148.96436d0),
     &  Iso(150, 149.95852d0),
     &  Iso(151, 150.9554d0),
     &  Iso(152, 151.95027d0),
     &  Iso(153, 152.94932d0),
     &  Iso(154, 153.946396d0),
     &  Iso(155, 154.945783d0),
     &  Iso(156, 155.942825d0),
     &  Iso(157, 156.942645d0),
     &  Iso(158, 157.9398705d0),
     &  Iso(159, 158.940055d0),
     &  Iso(160, 159.937557d0),
     &  Iso(161, 160.937907d0),
     &  Iso(162, 161.935774d0),
     &  Iso(163, 162.93634d0),
     &  Iso(164, 163.934495d0),
     &  Iso(165, 164.93527d0),
     &  Iso(166, 165.9338747d0),
     &  Iso(167, 166.934953d0),
     &  Iso(169, 168.9351825d0),
     &  Iso(175, 174.9412808d0),
     &  Iso(177, 176.9452656d0),
     &  Iso(178, 177.946651d0),
     &  Iso(179, 178.95004d0),
     &  Iso(180, 179.95212d0),
     &  Iso(181, 180.95589d0) ]

      ElementList(71)%Symbol = 'Lu'
      ElementList(71)%nat = 2
      Call mma_Allocate(ElementList(71)%Isotopes,36)
      ElementList(71)%Isotopes(:) = [
     &  Iso(175, 174.9407752d0),
     &  Iso(176, 175.9426897d0),
     &  Iso(150, 149.97355d0),
     &  Iso(151, 150.96768d0),
     &  Iso(152, 151.96412d0),
     &  Iso(153, 152.95875d0),
     &  Iso(154, 153.95736d0),
     &  Iso(155, 154.954321d0),
     &  Iso(156, 155.953033d0),
     &  Iso(157, 156.950127d0),
     &  Iso(158, 157.949316d0),
     &  Iso(159, 158.946636d0),
     &  Iso(160, 159.946033d0),
     &  Iso(161, 160.943572d0),
     &  Iso(162, 161.943283d0),
     &  Iso(163, 162.941179d0),
     &  Iso(164, 163.941339d0),
     &  Iso(165, 164.939407d0),
     &  Iso(166, 165.939859d0),
     &  Iso(167, 166.93827d0),
     &  Iso(168, 167.938736d0),
     &  Iso(169, 168.9376441d0),
     &  Iso(170, 169.938478d0),
     &  Iso(171, 170.937917d0),
     &  Iso(172, 171.9390891d0),
     &  Iso(173, 172.938934d0),
     &  Iso(174, 173.9403409d0),
     &  Iso(177, 176.9437615d0),
     &  Iso(178, 177.945958d0),
     &  Iso(179, 178.9473309d0),
     &  Iso(180, 179.949888d0),
     &  Iso(181, 180.95191d0),
     &  Iso(182, 181.95504d0),
     &  Iso(183, 182.957363d0),
     &  Iso(184, 183.96091d0),
     &  Iso(185, 184.96362d0) ]

      ElementList(72)%Symbol = 'Hf'
      ElementList(72)%nat = 6
      Call mma_Allocate(ElementList(72)%Isotopes,37)
      ElementList(72)%Isotopes(:) = [
     &  Iso(180, 179.946557d0),
     &  Iso(178, 177.9437058d0),
     &  Iso(177, 176.9432277d0),
     &  Iso(179, 178.9458232d0),
     &  Iso(176, 175.9414076d0),
     &  Iso(174, 173.9400461d0),
     &  Iso(153, 152.97069d0),
     &  Iso(154, 153.96486d0),
     &  Iso(155, 154.96311d0),
     &  Iso(156, 155.95935d0),
     &  Iso(157, 156.95824d0),
     &  Iso(158, 157.954801d0),
     &  Iso(159, 158.953996d0),
     &  Iso(160, 159.950691d0),
     &  Iso(161, 160.950278d0),
     &  Iso(162, 161.9472148d0),
     &  Iso(163, 162.947113d0),
     &  Iso(164, 163.944371d0),
     &  Iso(165, 164.944567d0),
     &  Iso(166, 165.94218d0),
     &  Iso(167, 166.9426d0),
     &  Iso(168, 167.940568d0),
     &  Iso(169, 168.941259d0),
     &  Iso(170, 169.939609d0),
     &  Iso(171, 170.940492d0),
     &  Iso(172, 171.93945d0),
     &  Iso(173, 172.940513d0),
     &  Iso(175, 174.9415092d0),
     &  Iso(181, 180.9491083d0),
     &  Iso(182, 181.9505612d0),
     &  Iso(183, 182.95353d0),
     &  Iso(184, 183.955446d0),
     &  Iso(185, 184.958862d0),
     &  Iso(186, 185.960897d0),
     &  Iso(187, 186.96477d0),
     &  Iso(188, 187.96685d0),
     &  Iso(189, 188.97084d0) ]

      ElementList(73)%Symbol = 'Ta'
      ElementList(73)%nat = 2
      Call mma_Allocate(ElementList(73)%Isotopes,38)
      ElementList(73)%Isotopes(:) = [
     &  Iso(181, 180.9479958d0),
     &  Iso(180, 179.9474648d0),
     &  Iso(155, 154.97424d0),
     &  Iso(156, 155.97203d0),
     &  Iso(157, 156.96818d0),
     &  Iso(158, 157.96654d0),
     &  Iso(159, 158.963023d0),
     &  Iso(160, 159.961488d0),
     &  Iso(161, 160.958452d0),
     &  Iso(162, 161.957294d0),
     &  Iso(163, 162.954337d0),
     &  Iso(164, 163.953534d0),
     &  Iso(165, 164.950781d0),
     &  Iso(166, 165.950512d0),
     &  Iso(167, 166.948093d0),
     &  Iso(168, 167.948047d0),
     &  Iso(169, 168.946011d0),
     &  Iso(170, 169.946175d0),
     &  Iso(171, 170.944476d0),
     &  Iso(172, 171.944895d0),
     &  Iso(173, 172.94375d0),
     &  Iso(174, 173.944454d0),
     &  Iso(175, 174.943737d0),
     &  Iso(176, 175.944857d0),
     &  Iso(177, 176.9444795d0),
     &  Iso(178, 177.945678d0),
     &  Iso(179, 178.9459366d0),
     &  Iso(182, 181.9501519d0),
     &  Iso(183, 182.9513726d0),
     &  Iso(184, 183.954008d0),
     &  Iso(185, 184.955559d0),
     &  Iso(186, 185.958551d0),
     &  Iso(187, 186.960386d0),
     &  Iso(188, 187.963916d0),
     &  Iso(189, 188.96583d0),
     &  Iso(190, 189.96939d0),
     &  Iso(191, 190.97156d0),
     &  Iso(192, 191.97514d0) ]

      ElementList(74)%Symbol = 'W'
      ElementList(74)%nat = 5
      Call mma_Allocate(ElementList(74)%Isotopes,38)
      ElementList(74)%Isotopes(:) = [
     &  Iso(184, 183.95093092d0),
     &  Iso(186, 185.9543628d0),
     &  Iso(182, 181.94820394d0),
     &  Iso(183, 182.95022275d0),
     &  Iso(180, 179.9467108d0),
     &  Iso(157, 156.97884d0),
     &  Iso(158, 157.97456d0),
     &  Iso(159, 158.97264d0),
     &  Iso(160, 159.96846d0),
     &  Iso(161, 160.9672d0),
     &  Iso(162, 161.963499d0),
     &  Iso(163, 162.962524d0),
     &  Iso(164, 163.958961d0),
     &  Iso(165, 164.958281d0),
     &  Iso(166, 165.955031d0),
     &  Iso(167, 166.954805d0),
     &  Iso(168, 167.951806d0),
     &  Iso(169, 168.951779d0),
     &  Iso(170, 169.949232d0),
     &  Iso(171, 170.949451d0),
     &  Iso(172, 171.947292d0),
     &  Iso(173, 172.947689d0),
     &  Iso(174, 173.946079d0),
     &  Iso(175, 174.946717d0),
     &  Iso(176, 175.945634d0),
     &  Iso(177, 176.946643d0),
     &  Iso(178, 177.945883d0),
     &  Iso(179, 178.947077d0),
     &  Iso(181, 180.9481978d0),
     &  Iso(185, 184.95341897d0),
     &  Iso(187, 186.9571588d0),
     &  Iso(188, 187.9584862d0),
     &  Iso(189, 188.961763d0),
     &  Iso(190, 189.963091d0),
     &  Iso(191, 190.966531d0),
     &  Iso(192, 191.96817d0),
     &  Iso(193, 192.97178d0),
     &  Iso(194, 193.97367d0) ]

      ElementList(75)%Symbol = 'Re'
      ElementList(75)%nat = 2
      Call mma_Allocate(ElementList(75)%Isotopes,40)
      ElementList(75)%Isotopes(:) = [
     &  Iso(187, 186.9557501d0),
     &  Iso(185, 184.9529545d0),
     &  Iso(159, 158.98418d0),
     &  Iso(160, 159.98182d0),
     &  Iso(161, 160.97757d0),
     &  Iso(162, 161.97584d0),
     &  Iso(163, 162.97208d0),
     &  Iso(164, 163.970453d0),
     &  Iso(165, 164.967103d0),
     &  Iso(166, 165.965761d0),
     &  Iso(167, 166.962595d0),
     &  Iso(168, 167.961573d0),
     &  Iso(169, 168.958766d0),
     &  Iso(170, 169.95822d0),
     &  Iso(171, 170.955716d0),
     &  Iso(172, 171.95542d0),
     &  Iso(173, 172.953243d0),
     &  Iso(174, 173.953115d0),
     &  Iso(175, 174.951381d0),
     &  Iso(176, 175.951623d0),
     &  Iso(177, 176.950328d0),
     &  Iso(178, 177.950989d0),
     &  Iso(179, 178.949989d0),
     &  Iso(180, 179.950792d0),
     &  Iso(181, 180.950058d0),
     &  Iso(182, 181.95121d0),
     &  Iso(183, 182.9508196d0),
     &  Iso(184, 183.9525228d0),
     &  Iso(186, 185.9549856d0),
     &  Iso(188, 187.9581115d0),
     &  Iso(189, 188.959226d0),
     &  Iso(190, 189.961744d0),
     &  Iso(191, 190.963122d0),
     &  Iso(192, 191.966088d0),
     &  Iso(193, 192.967541d0),
     &  Iso(194, 193.97076d0),
     &  Iso(195, 194.97254d0),
     &  Iso(196, 195.9758d0),
     &  Iso(197, 196.97799d0),
     &  Iso(198, 197.9816d0) ]

      ElementList(76)%Symbol = 'Os'
      ElementList(76)%nat = 7
      Call mma_Allocate(ElementList(76)%Isotopes,42)
      ElementList(76)%Isotopes(:) = [
     &  Iso(192, 191.961477d0),
     &  Iso(190, 189.9584437d0),
     &  Iso(189, 188.9581442d0),
     &  Iso(188, 187.9558352d0),
     &  Iso(187, 186.9557474d0),
     &  Iso(186, 185.953835d0),
     &  Iso(184, 183.9524885d0),
     &  Iso(161, 160.98903d0),
     &  Iso(162, 161.98443d0),
     &  Iso(163, 162.98241d0),
     &  Iso(164, 163.97802d0),
     &  Iso(165, 164.9766d0),
     &  Iso(166, 165.972692d0),
     &  Iso(167, 166.971549d0),
     &  Iso(168, 167.967808d0),
     &  Iso(169, 168.967018d0),
     &  Iso(170, 169.963578d0),
     &  Iso(171, 170.963174d0),
     &  Iso(172, 171.960017d0),
     &  Iso(173, 172.959808d0),
     &  Iso(174, 173.957064d0),
     &  Iso(175, 174.956945d0),
     &  Iso(176, 175.954806d0),
     &  Iso(177, 176.954966d0),
     &  Iso(178, 177.953254d0),
     &  Iso(179, 178.953817d0),
     &  Iso(180, 179.952375d0),
     &  Iso(181, 180.953247d0),
     &  Iso(182, 181.95211d0),
     &  Iso(183, 182.953125d0),
     &  Iso(185, 184.9540417d0),
     &  Iso(191, 190.9609264d0),
     &  Iso(193, 192.9641479d0),
     &  Iso(194, 193.9651772d0),
     &  Iso(195, 194.968318d0),
     &  Iso(196, 195.969641d0),
     &  Iso(197, 196.97283d0),
     &  Iso(198, 197.97441d0),
     &  Iso(199, 198.97801d0),
     &  Iso(200, 199.97984d0),
     &  Iso(201, 200.98364d0),
     &  Iso(202, 201.98595d0) ]

      ElementList(77)%Symbol = 'Ir'
      ElementList(77)%nat = 2
      Call mma_Allocate(ElementList(77)%Isotopes,41)
      ElementList(77)%Isotopes(:) = [
     &  Iso(193, 192.9629216d0),
     &  Iso(191, 190.9605893d0),
     &  Iso(164, 163.99191d0),
     &  Iso(165, 164.9875d0),
     &  Iso(166, 165.98566d0),
     &  Iso(167, 166.981666d0),
     &  Iso(168, 167.979907d0),
     &  Iso(169, 168.976298d0),
     &  Iso(170, 169.974922d0),
     &  Iso(171, 170.97164d0),
     &  Iso(172, 171.970607d0),
     &  Iso(173, 172.967506d0),
     &  Iso(174, 173.966861d0),
     &  Iso(175, 174.96415d0),
     &  Iso(176, 175.96365d0),
     &  Iso(177, 176.961301d0),
     &  Iso(178, 177.961082d0),
     &  Iso(179, 178.95912d0),
     &  Iso(180, 179.959229d0),
     &  Iso(181, 180.957625d0),
     &  Iso(182, 181.958076d0),
     &  Iso(183, 182.95684d0),
     &  Iso(184, 183.957476d0),
     &  Iso(185, 184.956698d0),
     &  Iso(186, 185.957944d0),
     &  Iso(187, 186.957542d0),
     &  Iso(188, 187.958828d0),
     &  Iso(189, 188.958715d0),
     &  Iso(190, 189.9605412d0),
     &  Iso(192, 191.9626002d0),
     &  Iso(194, 193.9650735d0),
     &  Iso(195, 194.9659747d0),
     &  Iso(196, 195.968397d0),
     &  Iso(197, 196.969655d0),
     &  Iso(198, 197.97228d0),
     &  Iso(199, 198.973805d0),
     &  Iso(200, 199.9768d0),
     &  Iso(201, 200.97864d0),
     &  Iso(202, 201.98199d0),
     &  Iso(203, 202.98423d0),
     &  Iso(204, 203.9896d0) ]

      ElementList(78)%Symbol = 'Pt'
      ElementList(78)%nat = 6
      Call mma_Allocate(ElementList(78)%Isotopes,41)
      ElementList(78)%Isotopes(:) = [
     &  Iso(195, 194.9647917d0),
     &  Iso(194, 193.9626809d0),
     &  Iso(196, 195.96495209d0),
     &  Iso(198, 197.9678949d0),
     &  Iso(192, 191.9610387d0),
     &  Iso(190, 189.9599297d0),
     &  Iso(166, 165.99486d0),
     &  Iso(167, 166.99269d0),
     &  Iso(168, 167.98813d0),
     &  Iso(169, 168.98657d0),
     &  Iso(170, 169.982496d0),
     &  Iso(171, 170.981245d0),
     &  Iso(172, 171.977351d0),
     &  Iso(173, 172.976443d0),
     &  Iso(174, 173.97282d0),
     &  Iso(175, 174.97241d0),
     &  Iso(176, 175.968938d0),
     &  Iso(177, 176.96847d0),
     &  Iso(178, 177.96565d0),
     &  Iso(179, 178.965359d0),
     &  Iso(180, 179.963032d0),
     &  Iso(181, 180.963098d0),
     &  Iso(182, 181.961172d0),
     &  Iso(183, 182.961597d0),
     &  Iso(184, 183.959915d0),
     &  Iso(185, 184.960614d0),
     &  Iso(186, 185.959351d0),
     &  Iso(187, 186.960617d0),
     &  Iso(188, 187.9593889d0),
     &  Iso(189, 188.960831d0),
     &  Iso(191, 190.9616729d0),
     &  Iso(193, 192.9629824d0),
     &  Iso(197, 196.96734069d0),
     &  Iso(199, 198.9705952d0),
     &  Iso(200, 199.971443d0),
     &  Iso(201, 200.974513d0),
     &  Iso(202, 201.975639d0),
     &  Iso(203, 202.97893d0),
     &  Iso(204, 203.98076d0),
     &  Iso(205, 204.98608d0),
     &  Iso(206, 205.98966d0) ]

      ElementList(79)%Symbol = 'Au'
      ElementList(79)%nat = 1
      Call mma_Allocate(ElementList(79)%Isotopes,42)
      ElementList(79)%Isotopes(:) = [
     &  Iso(197, 196.96656879d0),
     &  Iso(169, 168.99808d0),
     &  Iso(170, 169.99597d0),
     &  Iso(171, 170.991876d0),
     &  Iso(172, 171.989942d0),
     &  Iso(173, 172.986241d0),
     &  Iso(174, 173.984717d0),
     &  Iso(175, 174.981304d0),
     &  Iso(176, 175.98025d0),
     &  Iso(177, 176.97687d0),
     &  Iso(178, 177.976032d0),
     &  Iso(179, 178.973174d0),
     &  Iso(180, 179.972523d0),
     &  Iso(181, 180.970079d0),
     &  Iso(182, 181.969618d0),
     &  Iso(183, 182.967591d0),
     &  Iso(184, 183.967452d0),
     &  Iso(185, 184.96579d0),
     &  Iso(186, 185.965953d0),
     &  Iso(187, 186.964543d0),
     &  Iso(188, 187.965349d0),
     &  Iso(189, 188.963948d0),
     &  Iso(190, 189.964698d0),
     &  Iso(191, 190.963702d0),
     &  Iso(192, 191.964814d0),
     &  Iso(193, 192.9641373d0),
     &  Iso(194, 193.9654178d0),
     &  Iso(195, 194.9650352d0),
     &  Iso(196, 195.9665699d0),
     &  Iso(198, 197.96824242d0),
     &  Iso(199, 198.96876528d0),
     &  Iso(200, 199.970756d0),
     &  Iso(201, 200.9716575d0),
     &  Iso(202, 201.973856d0),
     &  Iso(203, 202.9751544d0),
     &  Iso(204, 203.97783d0),
     &  Iso(205, 204.97985d0),
     &  Iso(206, 205.98474d0),
     &  Iso(207, 206.9884d0),
     &  Iso(208, 207.99345d0),
     &  Iso(209, 208.99735d0),
     &  Iso(210, 210.0025d0) ]

      ElementList(80)%Symbol = 'Hg'
      ElementList(80)%nat = 7
      Call mma_Allocate(ElementList(80)%Isotopes,46)
      ElementList(80)%Isotopes(:) = [
     &  Iso(202, 201.9706434d0),
     &  Iso(200, 199.96832659d0),
     &  Iso(199, 198.96828064d0),
     &  Iso(201, 200.97030284d0),
     &  Iso(198, 197.9667686d0),
     &  Iso(204, 203.97349398d0),
     &  Iso(196, 195.9658326d0),
     &  Iso(171, 171.00353d0),
     &  Iso(172, 171.99881d0),
     &  Iso(173, 172.99709d0),
     &  Iso(174, 173.992865d0),
     &  Iso(175, 174.991441d0),
     &  Iso(176, 175.987361d0),
     &  Iso(177, 176.986277d0),
     &  Iso(178, 177.982484d0),
     &  Iso(179, 178.981831d0),
     &  Iso(180, 179.97826d0),
     &  Iso(181, 180.977819d0),
     &  Iso(182, 181.974689d0),
     &  Iso(183, 182.9744448d0),
     &  Iso(184, 183.971714d0),
     &  Iso(185, 184.971899d0),
     &  Iso(186, 185.969362d0),
     &  Iso(187, 186.969814d0),
     &  Iso(188, 187.967567d0),
     &  Iso(189, 188.968195d0),
     &  Iso(190, 189.966323d0),
     &  Iso(191, 190.967157d0),
     &  Iso(192, 191.965635d0),
     &  Iso(193, 192.966653d0),
     &  Iso(194, 193.9654491d0),
     &  Iso(195, 194.966721d0),
     &  Iso(197, 196.9672128d0),
     &  Iso(203, 202.9728728d0),
     &  Iso(205, 204.9760734d0),
     &  Iso(206, 205.977514d0),
     &  Iso(207, 206.9823d0),
     &  Iso(208, 207.985759d0),
     &  Iso(209, 208.99072d0),
     &  Iso(210, 209.99424d0),
     &  Iso(211, 210.99933d0),
     &  Iso(212, 212.00296d0),
     &  Iso(213, 213.00823d0),
     &  Iso(214, 214.012d0),
     &  Iso(215, 215.0174d0),
     &  Iso(216, 216.02132d0) ]

      ElementList(81)%Symbol = 'Tl'
      ElementList(81)%nat = 2
      Call mma_Allocate(ElementList(81)%Isotopes,43)
      ElementList(81)%Isotopes(:) = [
     &  Iso(205, 204.9744278d0),
     &  Iso(203, 202.9723446d0),
     &  Iso(176, 176.000624d0),
     &  Iso(177, 176.996431d0),
     &  Iso(178, 177.99485d0),
     &  Iso(179, 178.991111d0),
     &  Iso(180, 179.990057d0),
     &  Iso(181, 180.98626d0),
     &  Iso(182, 181.985713d0),
     &  Iso(183, 182.982193d0),
     &  Iso(184, 183.981886d0),
     &  Iso(185, 184.978789d0),
     &  Iso(186, 185.978651d0),
     &  Iso(187, 186.9759063d0),
     &  Iso(188, 187.976021d0),
     &  Iso(189, 188.973588d0),
     &  Iso(190, 189.973828d0),
     &  Iso(191, 190.9717842d0),
     &  Iso(192, 191.972225d0),
     &  Iso(193, 192.970502d0),
     &  Iso(194, 193.971081d0),
     &  Iso(195, 194.969774d0),
     &  Iso(196, 195.970481d0),
     &  Iso(197, 196.969576d0),
     &  Iso(198, 197.970483d0),
     &  Iso(199, 198.969877d0),
     &  Iso(200, 199.9709633d0),
     &  Iso(201, 200.970822d0),
     &  Iso(202, 201.972102d0),
     &  Iso(204, 203.9738639d0),
     &  Iso(206, 205.9761106d0),
     &  Iso(207, 206.9774197d0),
     &  Iso(208, 207.982019d0),
     &  Iso(209, 208.9853594d0),
     &  Iso(210, 209.990074d0),
     &  Iso(211, 210.993475d0),
     &  Iso(212, 211.99834d0),
     &  Iso(213, 213.001915d0),
     &  Iso(214, 214.00694d0),
     &  Iso(215, 215.01064d0),
     &  Iso(216, 216.0158d0),
     &  Iso(217, 217.01966d0),
     &  Iso(218, 218.02479d0) ]

      ElementList(82)%Symbol = 'Pb'
      ElementList(82)%nat = 4
      Call mma_Allocate(ElementList(82)%Isotopes,43)
      ElementList(82)%Isotopes(:) = [
     &  Iso(208, 207.9766525d0),
     &  Iso(206, 205.9744657d0),
     &  Iso(207, 206.9758973d0),
     &  Iso(204, 203.973044d0),
     &  Iso(178, 178.003831d0),
     &  Iso(179, 179.002201d0),
     &  Iso(180, 179.997928d0),
     &  Iso(181, 180.996653d0),
     &  Iso(182, 181.992672d0),
     &  Iso(183, 182.991872d0),
     &  Iso(184, 183.988136d0),
     &  Iso(185, 184.98761d0),
     &  Iso(186, 185.984238d0),
     &  Iso(187, 186.9839109d0),
     &  Iso(188, 187.980875d0),
     &  Iso(189, 188.980807d0),
     &  Iso(190, 189.978082d0),
     &  Iso(191, 190.978276d0),
     &  Iso(192, 191.975775d0),
     &  Iso(193, 192.976173d0),
     &  Iso(194, 193.974012d0),
     &  Iso(195, 194.974543d0),
     &  Iso(196, 195.972774d0),
     &  Iso(197, 196.9734312d0),
     &  Iso(198, 197.972034d0),
     &  Iso(199, 198.972913d0),
     &  Iso(200, 199.971819d0),
     &  Iso(201, 200.972883d0),
     &  Iso(202, 201.972152d0),
     &  Iso(203, 202.9733911d0),
     &  Iso(205, 204.9744822d0),
     &  Iso(209, 208.9810905d0),
     &  Iso(210, 209.9841889d0),
     &  Iso(211, 210.9887371d0),
     &  Iso(212, 211.9918977d0),
     &  Iso(213, 212.9965629d0),
     &  Iso(214, 213.9998059d0),
     &  Iso(215, 215.00474d0),
     &  Iso(216, 216.00803d0),
     &  Iso(217, 217.01314d0),
     &  Iso(218, 218.01659d0),
     &  Iso(219, 219.02177d0),
     &  Iso(220, 220.02541d0) ]

      ElementList(83)%Symbol = 'Bi'
      ElementList(83)%nat = 1
      Call mma_Allocate(ElementList(83)%Isotopes,41)
      ElementList(83)%Isotopes(:) = [
     &  Iso(209, 208.9803991d0),
     &  Iso(184, 184.001275d0),
     &  Iso(185, 184.9976d0),
     &  Iso(186, 185.996644d0),
     &  Iso(187, 186.993147d0),
     &  Iso(188, 187.992287d0),
     &  Iso(189, 188.989195d0),
     &  Iso(190, 189.988622d0),
     &  Iso(191, 190.9857866d0),
     &  Iso(192, 191.985469d0),
     &  Iso(193, 192.98296d0),
     &  Iso(194, 193.982785d0),
     &  Iso(195, 194.9806488d0),
     &  Iso(196, 195.980667d0),
     &  Iso(197, 196.9788651d0),
     &  Iso(198, 197.979206d0),
     &  Iso(199, 198.977673d0),
     &  Iso(200, 199.978131d0),
     &  Iso(201, 200.97701d0),
     &  Iso(202, 201.977734d0),
     &  Iso(203, 202.976893d0),
     &  Iso(204, 203.9778361d0),
     &  Iso(205, 204.9773867d0),
     &  Iso(206, 205.9784993d0),
     &  Iso(207, 206.978471d0),
     &  Iso(208, 207.9797425d0),
     &  Iso(210, 209.9841207d0),
     &  Iso(211, 210.9872697d0),
     &  Iso(212, 211.991286d0),
     &  Iso(213, 212.9943851d0),
     &  Iso(214, 213.998712d0),
     &  Iso(215, 215.00177d0),
     &  Iso(216, 216.006306d0),
     &  Iso(217, 217.009372d0),
     &  Iso(218, 218.014188d0),
     &  Iso(219, 219.01748d0),
     &  Iso(220, 220.02235d0),
     &  Iso(221, 221.02587d0),
     &  Iso(222, 222.03078d0),
     &  Iso(223, 223.0345d0),
     &  Iso(224, 224.03947d0) ]

      ElementList(84)%Symbol = 'Po'
      ElementList(84)%nat = 0
      Call mma_Allocate(ElementList(84)%Isotopes,42)
      ElementList(84)%Isotopes(:) = [
     &  Iso(209, 208.9824308d0),
     &  Iso(186, 186.004393d0),
     &  Iso(187, 187.003041d0),
     &  Iso(188, 187.999416d0),
     &  Iso(189, 188.998473d0),
     &  Iso(190, 189.995101d0),
     &  Iso(191, 190.9945585d0),
     &  Iso(192, 191.991336d0),
     &  Iso(193, 192.991026d0),
     &  Iso(194, 193.988186d0),
     &  Iso(195, 194.988126d0),
     &  Iso(196, 195.985526d0),
     &  Iso(197, 196.98566d0),
     &  Iso(198, 197.983389d0),
     &  Iso(199, 198.983667d0),
     &  Iso(200, 199.981799d0),
     &  Iso(201, 200.9822598d0),
     &  Iso(202, 201.980758d0),
     &  Iso(203, 202.9814161d0),
     &  Iso(204, 203.98031d0),
     &  Iso(205, 204.981203d0),
     &  Iso(206, 205.980474d0),
     &  Iso(207, 206.9815938d0),
     &  Iso(208, 207.9812461d0),
     &  Iso(210, 209.9828741d0),
     &  Iso(211, 210.9866536d0),
     &  Iso(212, 211.9888684d0),
     &  Iso(213, 212.9928576d0),
     &  Iso(214, 213.9952017d0),
     &  Iso(215, 214.9994201d0),
     &  Iso(216, 216.0019152d0),
     &  Iso(217, 217.0063182d0),
     &  Iso(218, 218.0089735d0),
     &  Iso(219, 219.013614d0),
     &  Iso(220, 220.016386d0),
     &  Iso(221, 221.021228d0),
     &  Iso(222, 222.02414d0),
     &  Iso(223, 223.02907d0),
     &  Iso(224, 224.03211d0),
     &  Iso(225, 225.03707d0),
     &  Iso(226, 226.04031d0),
     &  Iso(227, 227.04539d0) ]

      ElementList(85)%Symbol = 'At'
      ElementList(85)%nat = 0
      Call mma_Allocate(ElementList(85)%Isotopes,39)
      ElementList(85)%Isotopes(:) = [
     &  Iso(210, 209.9871479d0),
     &  Iso(191, 191.004148d0),
     &  Iso(192, 192.003152d0),
     &  Iso(193, 192.999927d0),
     &  Iso(194, 193.999236d0),
     &  Iso(195, 194.9962685d0),
     &  Iso(196, 195.9958d0),
     &  Iso(197, 196.993189d0),
     &  Iso(198, 197.992784d0),
     &  Iso(199, 198.9905277d0),
     &  Iso(200, 199.990351d0),
     &  Iso(201, 200.9884171d0),
     &  Iso(202, 201.98863d0),
     &  Iso(203, 202.986943d0),
     &  Iso(204, 203.987251d0),
     &  Iso(205, 204.986076d0),
     &  Iso(206, 205.986657d0),
     &  Iso(207, 206.9858d0),
     &  Iso(208, 207.9866133d0),
     &  Iso(209, 208.9861702d0),
     &  Iso(211, 210.9874966d0),
     &  Iso(212, 211.9907377d0),
     &  Iso(213, 212.992937d0),
     &  Iso(214, 213.9963721d0),
     &  Iso(215, 214.9986528d0),
     &  Iso(216, 216.0024236d0),
     &  Iso(217, 217.0047192d0),
     &  Iso(218, 218.008695d0),
     &  Iso(219, 219.0111618d0),
     &  Iso(220, 220.015433d0),
     &  Iso(221, 221.018017d0),
     &  Iso(222, 222.022494d0),
     &  Iso(223, 223.025151d0),
     &  Iso(224, 224.029749d0),
     &  Iso(225, 225.03263d0),
     &  Iso(226, 226.03716d0),
     &  Iso(227, 227.04024d0),
     &  Iso(228, 228.04475d0),
     &  Iso(229, 229.04812d0) ]

      ElementList(86)%Symbol = 'Rn'
      ElementList(86)%nat = 0
      Call mma_Allocate(ElementList(86)%Isotopes,39)
      ElementList(86)%Isotopes(:) = [
     &  Iso(222, 222.0175782d0),
     &  Iso(193, 193.009708d0),
     &  Iso(194, 194.006144d0),
     &  Iso(195, 195.005422d0),
     &  Iso(196, 196.002116d0),
     &  Iso(197, 197.001585d0),
     &  Iso(198, 197.998679d0),
     &  Iso(199, 198.99839d0),
     &  Iso(200, 199.99569d0),
     &  Iso(201, 200.995628d0),
     &  Iso(202, 201.993264d0),
     &  Iso(203, 202.993388d0),
     &  Iso(204, 203.99143d0),
     &  Iso(205, 204.991719d0),
     &  Iso(206, 205.990214d0),
     &  Iso(207, 206.9907303d0),
     &  Iso(208, 207.989635d0),
     &  Iso(209, 208.990415d0),
     &  Iso(210, 209.9896891d0),
     &  Iso(211, 210.9906011d0),
     &  Iso(212, 211.9907039d0),
     &  Iso(213, 212.9938831d0),
     &  Iso(214, 213.995363d0),
     &  Iso(215, 214.9987459d0),
     &  Iso(216, 216.0002719d0),
     &  Iso(217, 217.003928d0),
     &  Iso(218, 218.0056016d0),
     &  Iso(219, 219.0094804d0),
     &  Iso(220, 220.0113941d0),
     &  Iso(221, 221.0155371d0),
     &  Iso(223, 223.0218893d0),
     &  Iso(224, 224.024096d0),
     &  Iso(225, 225.028486d0),
     &  Iso(226, 226.030861d0),
     &  Iso(227, 227.035304d0),
     &  Iso(228, 228.037835d0),
     &  Iso(229, 229.042257d0),
     &  Iso(230, 230.04514d0),
     &  Iso(231, 231.04987d0) ]

      ElementList(87)%Symbol = 'Fr'
      ElementList(87)%nat = 0
      Call mma_Allocate(ElementList(87)%Isotopes,35)
      ElementList(87)%Isotopes(:) = [
     &  Iso(223, 223.019736d0),
     &  Iso(199, 199.007259d0),
     &  Iso(200, 200.006586d0),
     &  Iso(201, 201.003867d0),
     &  Iso(202, 202.00332d0),
     &  Iso(203, 203.0009407d0),
     &  Iso(204, 204.000652d0),
     &  Iso(205, 204.9985939d0),
     &  Iso(206, 205.998666d0),
     &  Iso(207, 206.996946d0),
     &  Iso(208, 207.997138d0),
     &  Iso(209, 208.995955d0),
     &  Iso(210, 209.996422d0),
     &  Iso(211, 210.995556d0),
     &  Iso(212, 211.9962257d0),
     &  Iso(213, 212.996186d0),
     &  Iso(214, 213.9989713d0),
     &  Iso(215, 215.0003418d0),
     &  Iso(216, 216.0031899d0),
     &  Iso(217, 217.0046323d0),
     &  Iso(218, 218.0075787d0),
     &  Iso(219, 219.0092524d0),
     &  Iso(220, 220.0123277d0),
     &  Iso(221, 221.0142552d0),
     &  Iso(222, 222.017552d0),
     &  Iso(224, 224.023398d0),
     &  Iso(225, 225.025573d0),
     &  Iso(226, 226.029566d0),
     &  Iso(227, 227.031869d0),
     &  Iso(228, 228.035823d0),
     &  Iso(229, 229.038298d0),
     &  Iso(230, 230.042416d0),
     &  Iso(231, 231.045158d0),
     &  Iso(232, 232.04937d0),
     &  Iso(233, 233.05264d0) ]

      ElementList(88)%Symbol = 'Ra'
      ElementList(88)%nat = 0
      Call mma_Allocate(ElementList(88)%Isotopes,35)
      ElementList(88)%Isotopes(:) = [
     &  Iso(226, 226.0254103d0),
     &  Iso(201, 201.01271d0),
     &  Iso(202, 202.00976d0),
     &  Iso(203, 203.009304d0),
     &  Iso(204, 204.006492d0),
     &  Iso(205, 205.006268d0),
     &  Iso(206, 206.003828d0),
     &  Iso(207, 207.003799d0),
     &  Iso(208, 208.001841d0),
     &  Iso(209, 209.00199d0),
     &  Iso(210, 210.000494d0),
     &  Iso(211, 211.0008932d0),
     &  Iso(212, 211.999787d0),
     &  Iso(213, 213.000384d0),
     &  Iso(214, 214.0000997d0),
     &  Iso(215, 215.0027204d0),
     &  Iso(216, 216.0035334d0),
     &  Iso(217, 217.0063207d0),
     &  Iso(218, 218.007141d0),
     &  Iso(219, 219.0100855d0),
     &  Iso(220, 220.0110259d0),
     &  Iso(221, 221.0139177d0),
     &  Iso(222, 222.0153748d0),
     &  Iso(223, 223.0185023d0),
     &  Iso(224, 224.020212d0),
     &  Iso(225, 225.0236119d0),
     &  Iso(227, 227.0291783d0),
     &  Iso(228, 228.0310707d0),
     &  Iso(229, 229.034942d0),
     &  Iso(230, 230.037055d0),
     &  Iso(231, 231.041027d0),
     &  Iso(232, 232.0434753d0),
     &  Iso(233, 233.047582d0),
     &  Iso(234, 234.050342d0),
     &  Iso(235, 235.05497d0) ]

      ElementList(89)%Symbol = 'Ac'
      ElementList(89)%nat = 0
      Call mma_Allocate(ElementList(89)%Isotopes,32)
      ElementList(89)%Isotopes(:) = [
     &  Iso(227, 227.0277523d0),
     &  Iso(206, 206.014452d0),
     &  Iso(207, 207.011966d0),
     &  Iso(208, 208.01155d0),
     &  Iso(209, 209.009495d0),
     &  Iso(210, 210.009436d0),
     &  Iso(211, 211.007732d0),
     &  Iso(212, 212.007813d0),
     &  Iso(213, 213.006609d0),
     &  Iso(214, 214.006918d0),
     &  Iso(215, 215.006475d0),
     &  Iso(216, 216.008743d0),
     &  Iso(217, 217.009344d0),
     &  Iso(218, 218.011642d0),
     &  Iso(219, 219.012421d0),
     &  Iso(220, 220.0147549d0),
     &  Iso(221, 221.015592d0),
     &  Iso(222, 222.0178442d0),
     &  Iso(223, 223.0191377d0),
     &  Iso(224, 224.0217232d0),
     &  Iso(225, 225.02323d0),
     &  Iso(226, 226.0260984d0),
     &  Iso(228, 228.0310215d0),
     &  Iso(229, 229.032956d0),
     &  Iso(230, 230.036327d0),
     &  Iso(231, 231.038393d0),
     &  Iso(232, 232.042034d0),
     &  Iso(233, 233.044346d0),
     &  Iso(234, 234.048139d0),
     &  Iso(235, 235.05084d0),
     &  Iso(236, 236.054988d0),
     &  Iso(237, 237.05827d0) ]

      ElementList(90)%Symbol = 'Th'
      ElementList(90)%nat = 1
      Call mma_Allocate(ElementList(90)%Isotopes,32)
      ElementList(90)%Isotopes(:) = [
     &  Iso(232, 232.0380558d0),
     &  Iso(208, 208.0179d0),
     &  Iso(209, 209.017753d0),
     &  Iso(210, 210.015094d0),
     &  Iso(211, 211.014929d0),
     &  Iso(212, 212.012988d0),
     &  Iso(213, 213.013009d0),
     &  Iso(214, 214.0115d0),
     &  Iso(215, 215.0117248d0),
     &  Iso(216, 216.011056d0),
     &  Iso(217, 217.013117d0),
     &  Iso(218, 218.013276d0),
     &  Iso(219, 219.015537d0),
     &  Iso(220, 220.015748d0),
     &  Iso(221, 221.018184d0),
     &  Iso(222, 222.018469d0),
     &  Iso(223, 223.0208119d0),
     &  Iso(224, 224.021464d0),
     &  Iso(225, 225.0239514d0),
     &  Iso(226, 226.0249034d0),
     &  Iso(227, 227.0277042d0),
     &  Iso(228, 228.0287413d0),
     &  Iso(229, 229.0317627d0),
     &  Iso(230, 230.0331341d0),
     &  Iso(231, 231.0363046d0),
     &  Iso(233, 233.0415823d0),
     &  Iso(234, 234.0436014d0),
     &  Iso(235, 235.047255d0),
     &  Iso(236, 236.049657d0),
     &  Iso(237, 237.053629d0),
     &  Iso(238, 238.0565d0),
     &  Iso(239, 239.06077d0) ]

      ElementList(91)%Symbol = 'Pa'
      ElementList(91)%nat = 1
      Call mma_Allocate(ElementList(91)%Isotopes,30)
      ElementList(91)%Isotopes(:) = [
     &  Iso(231, 231.0358842d0),
     &  Iso(212, 212.023203d0),
     &  Iso(213, 213.021109d0),
     &  Iso(214, 214.020918d0),
     &  Iso(215, 215.019183d0),
     &  Iso(216, 216.019109d0),
     &  Iso(217, 217.018325d0),
     &  Iso(218, 218.020059d0),
     &  Iso(219, 219.019904d0),
     &  Iso(220, 220.021705d0),
     &  Iso(221, 221.021875d0),
     &  Iso(222, 222.023784d0),
     &  Iso(223, 223.023963d0),
     &  Iso(224, 224.0256176d0),
     &  Iso(225, 225.026131d0),
     &  Iso(226, 226.027948d0),
     &  Iso(227, 227.0288054d0),
     &  Iso(228, 228.0310517d0),
     &  Iso(229, 229.0320972d0),
     &  Iso(230, 230.034541d0),
     &  Iso(232, 232.0385917d0),
     &  Iso(233, 233.0402472d0),
     &  Iso(234, 234.0433072d0),
     &  Iso(235, 235.045399d0),
     &  Iso(236, 236.048668d0),
     &  Iso(237, 237.051023d0),
     &  Iso(238, 238.054637d0),
     &  Iso(239, 239.05726d0),
     &  Iso(240, 240.06098d0),
     &  Iso(241, 241.06408d0) ]

      ElementList(92)%Symbol = 'U'
      ElementList(92)%nat = 3
      Call mma_Allocate(ElementList(92)%Isotopes,27)
      ElementList(92)%Isotopes(:) = [
     &  Iso(238, 238.0507884d0),
     &  Iso(235, 235.0439301d0),
     &  Iso(234, 234.0409523d0),
     &  Iso(217, 217.02466d0),
     &  Iso(218, 218.023523d0),
     &  Iso(219, 219.024999d0),
     &  Iso(220, 220.02462d0),
     &  Iso(221, 221.02628d0),
     &  Iso(222, 222.026d0),
     &  Iso(223, 223.027739d0),
     &  Iso(224, 224.027605d0),
     &  Iso(225, 225.029391d0),
     &  Iso(226, 226.029339d0),
     &  Iso(227, 227.031157d0),
     &  Iso(228, 228.031371d0),
     &  Iso(229, 229.0335063d0),
     &  Iso(230, 230.0339401d0),
     &  Iso(231, 231.0362939d0),
     &  Iso(232, 232.0371563d0),
     &  Iso(233, 233.0396355d0),
     &  Iso(236, 236.0455682d0),
     &  Iso(237, 237.0487304d0),
     &  Iso(239, 239.0542935d0),
     &  Iso(240, 240.0565934d0),
     &  Iso(241, 241.06033d0),
     &  Iso(242, 242.06293d0),
     &  Iso(243, 243.06699d0) ]

      ElementList(93)%Symbol = 'Np'
      ElementList(93)%nat = 0
      Call mma_Allocate(ElementList(93)%Isotopes,27)
      ElementList(93)%Isotopes(:) = [
     &  Iso(237, 237.0481736d0),
     &  Iso(219, 219.03143d0),
     &  Iso(220, 220.03254d0),
     &  Iso(221, 221.03204d0),
     &  Iso(222, 222.0333d0),
     &  Iso(223, 223.03285d0),
     &  Iso(224, 224.03422d0),
     &  Iso(225, 225.033911d0),
     &  Iso(226, 226.035188d0),
     &  Iso(227, 227.034957d0),
     &  Iso(228, 228.036067d0),
     &  Iso(229, 229.036264d0),
     &  Iso(230, 230.037828d0),
     &  Iso(231, 231.038245d0),
     &  Iso(232, 232.04011d0),
     &  Iso(233, 233.040741d0),
     &  Iso(234, 234.0428953d0),
     &  Iso(235, 235.0440635d0),
     &  Iso(236, 236.04657d0),
     &  Iso(238, 238.0509466d0),
     &  Iso(239, 239.0529392d0),
     &  Iso(240, 240.056165d0),
     &  Iso(241, 241.058253d0),
     &  Iso(242, 242.06164d0),
     &  Iso(243, 243.06428d0),
     &  Iso(244, 244.06785d0),
     &  Iso(245, 245.0708d0) ]

      ElementList(94)%Symbol = 'Pu'
      ElementList(94)%nat = 0
      Call mma_Allocate(ElementList(94)%Isotopes,20)
      ElementList(94)%Isotopes(:) = [
     &  Iso(244, 244.0642053d0),
     &  Iso(228, 228.038732d0),
     &  Iso(229, 229.040144d0),
     &  Iso(230, 230.03965d0),
     &  Iso(231, 231.041102d0),
     &  Iso(232, 232.041185d0),
     &  Iso(233, 233.042998d0),
     &  Iso(234, 234.0433174d0),
     &  Iso(235, 235.045286d0),
     &  Iso(236, 236.0460581d0),
     &  Iso(237, 237.0484098d0),
     &  Iso(238, 238.0495601d0),
     &  Iso(239, 239.0521636d0),
     &  Iso(240, 240.0538138d0),
     &  Iso(241, 241.0568517d0),
     &  Iso(242, 242.0587428d0),
     &  Iso(243, 243.0620036d0),
     &  Iso(245, 245.067826d0),
     &  Iso(246, 246.070205d0),
     &  Iso(247, 247.07419d0) ]

      ElementList(95)%Symbol = 'Am'
      ElementList(95)%nat = 0
      Call mma_Allocate(ElementList(95)%Isotopes,20)
      ElementList(95)%Isotopes(:) = [
     &  Iso(243, 243.0613813d0),
     &  Iso(230, 230.04609d0),
     &  Iso(231, 231.04556d0),
     &  Iso(232, 232.04645d0),
     &  Iso(233, 233.04644d0),
     &  Iso(234, 234.04773d0),
     &  Iso(235, 235.047908d0),
     &  Iso(236, 236.04943d0),
     &  Iso(237, 237.049996d0),
     &  Iso(238, 238.051985d0),
     &  Iso(239, 239.0530247d0),
     &  Iso(240, 240.0553d0),
     &  Iso(241, 241.0568293d0),
     &  Iso(242, 242.0595494d0),
     &  Iso(244, 244.0642851d0),
     &  Iso(245, 245.0664548d0),
     &  Iso(246, 246.069775d0),
     &  Iso(247, 247.07209d0),
     &  Iso(248, 248.07575d0),
     &  Iso(249, 249.07848d0) ]

      ElementList(96)%Symbol = 'Cm'
      ElementList(96)%nat = 0
      Call mma_Allocate(ElementList(96)%Isotopes,21)
      ElementList(96)%Isotopes(:) = [
     &  Iso(247, 247.0703541d0),
     &  Iso(232, 232.04982d0),
     &  Iso(233, 233.05077d0),
     &  Iso(234, 234.05016d0),
     &  Iso(235, 235.05154d0),
     &  Iso(236, 236.051374d0),
     &  Iso(237, 237.052869d0),
     &  Iso(238, 238.053081d0),
     &  Iso(239, 239.05491d0),
     &  Iso(240, 240.0555297d0),
     &  Iso(241, 241.0576532d0),
     &  Iso(242, 242.058836d0),
     &  Iso(243, 243.0613893d0),
     &  Iso(244, 244.0627528d0),
     &  Iso(245, 245.0654915d0),
     &  Iso(246, 246.0672238d0),
     &  Iso(248, 248.0723499d0),
     &  Iso(249, 249.0759548d0),
     &  Iso(250, 250.078358d0),
     &  Iso(251, 251.082286d0),
     &  Iso(252, 252.08487d0) ]

      ElementList(97)%Symbol = 'Bk'
      ElementList(97)%nat = 0
      Call mma_Allocate(ElementList(97)%Isotopes,21)
      ElementList(97)%Isotopes(:) = [
     &  Iso(247, 247.0703073d0),
     &  Iso(234, 234.05727d0),
     &  Iso(235, 235.05658d0),
     &  Iso(236, 236.05748d0),
     &  Iso(237, 237.0571d0),
     &  Iso(238, 238.0582d0),
     &  Iso(239, 239.05824d0),
     &  Iso(240, 240.05976d0),
     &  Iso(241, 241.06016d0),
     &  Iso(242, 242.06198d0),
     &  Iso(243, 243.0630078d0),
     &  Iso(244, 244.065181d0),
     &  Iso(245, 245.0663618d0),
     &  Iso(246, 246.068673d0),
     &  Iso(248, 248.073088d0),
     &  Iso(249, 249.0749877d0),
     &  Iso(250, 250.0783167d0),
     &  Iso(251, 251.080762d0),
     &  Iso(252, 252.08431d0),
     &  Iso(253, 253.08688d0),
     &  Iso(254, 254.0906d0) ]

      ElementList(98)%Symbol = 'Cf'
      ElementList(98)%nat = 0
      Call mma_Allocate(ElementList(98)%Isotopes,20)
      ElementList(98)%Isotopes(:) = [
     &  Iso(251, 251.0795886d0),
     &  Iso(237, 237.062198d0),
     &  Iso(238, 238.06149d0),
     &  Iso(239, 239.06253d0),
     &  Iso(240, 240.062256d0),
     &  Iso(241, 241.06369d0),
     &  Iso(242, 242.063754d0),
     &  Iso(243, 243.06548d0),
     &  Iso(244, 244.0660008d0),
     &  Iso(245, 245.0680487d0),
     &  Iso(246, 246.0688055d0),
     &  Iso(247, 247.070965d0),
     &  Iso(248, 248.0721851d0),
     &  Iso(249, 249.0748539d0),
     &  Iso(250, 250.0764062d0),
     &  Iso(252, 252.0816272d0),
     &  Iso(253, 253.0851345d0),
     &  Iso(254, 254.087324d0),
     &  Iso(255, 255.09105d0),
     &  Iso(256, 256.09344d0) ]

      ElementList(99)%Symbol = 'Es'
      ElementList(99)%nat = 0
      Call mma_Allocate(ElementList(99)%Isotopes,20)
      ElementList(99)%Isotopes(:) = [
     &  Iso(252, 252.08298d0),
     &  Iso(239, 239.06823d0),
     &  Iso(240, 240.06892d0),
     &  Iso(241, 241.06856d0),
     &  Iso(242, 242.06957d0),
     &  Iso(243, 243.06951d0),
     &  Iso(244, 244.07088d0),
     &  Iso(245, 245.07125d0),
     &  Iso(246, 246.0729d0),
     &  Iso(247, 247.073622d0),
     &  Iso(248, 248.075471d0),
     &  Iso(249, 249.076411d0),
     &  Iso(250, 250.07861d0),
     &  Iso(251, 251.0799936d0),
     &  Iso(253, 253.0848257d0),
     &  Iso(254, 254.0880222d0),
     &  Iso(255, 255.090275d0),
     &  Iso(256, 256.0936d0),
     &  Iso(257, 257.09598d0),
     &  Iso(258, 258.09952d0) ]

      ElementList(100)%Symbol = 'Fm'
      ElementList(100)%nat = 0
      Call mma_Allocate(ElementList(100)%Isotopes,20)
      ElementList(100)%Isotopes(:) = [
     &  Iso(257, 257.0951061d0),
     &  Iso(241, 241.07421d0),
     &  Iso(242, 242.07343d0),
     &  Iso(243, 243.07446d0),
     &  Iso(244, 244.07404d0),
     &  Iso(245, 245.07535d0),
     &  Iso(246, 246.07535d0),
     &  Iso(247, 247.07694d0),
     &  Iso(248, 248.0771865d0),
     &  Iso(249, 249.0789275d0),
     &  Iso(250, 250.079521d0),
     &  Iso(251, 251.08154d0),
     &  Iso(252, 252.0824671d0),
     &  Iso(253, 253.0851846d0),
     &  Iso(254, 254.0868544d0),
     &  Iso(255, 255.089964d0),
     &  Iso(256, 256.0917745d0),
     &  Iso(258, 258.09708d0),
     &  Iso(259, 259.1006d0),
     &  Iso(260, 260.10281d0) ]

      ElementList(101)%Symbol = 'Md'
      ElementList(101)%nat = 0
      Call mma_Allocate(ElementList(101)%Isotopes,18)
      ElementList(101)%Isotopes(:) = [
     &  Iso(258, 258.0984315d0),
     &  Iso(245, 245.08081d0),
     &  Iso(246, 246.08171d0),
     &  Iso(247, 247.08152d0),
     &  Iso(248, 248.08282d0),
     &  Iso(249, 249.08291d0),
     &  Iso(250, 250.08441d0),
     &  Iso(251, 251.084774d0),
     &  Iso(252, 252.08643d0),
     &  Iso(253, 253.087144d0),
     &  Iso(254, 254.08959d0),
     &  Iso(255, 255.0910841d0),
     &  Iso(256, 256.09389d0),
     &  Iso(257, 257.0955424d0),
     &  Iso(259, 259.10051d0),
     &  Iso(260, 260.10365d0),
     &  Iso(261, 261.10583d0),
     &  Iso(262, 262.1091d0) ]

      ElementList(102)%Symbol = 'No'
      ElementList(102)%nat = 0
      Call mma_Allocate(ElementList(102)%Isotopes,17)
      ElementList(102)%Isotopes(:) = [
     &  Iso(259, 259.10103d0),
     &  Iso(248, 248.08655d0),
     &  Iso(249, 249.0878d0),
     &  Iso(250, 250.08756d0),
     &  Iso(251, 251.08894d0),
     &  Iso(252, 252.088967d0),
     &  Iso(253, 253.0905641d0),
     &  Iso(254, 254.090956d0),
     &  Iso(255, 255.093191d0),
     &  Iso(256, 256.0942829d0),
     &  Iso(257, 257.0968878d0),
     &  Iso(258, 258.09821d0),
     &  Iso(260, 260.10264d0),
     &  Iso(261, 261.1057d0),
     &  Iso(262, 262.10746d0),
     &  Iso(263, 263.11071d0),
     &  Iso(264, 264.11273d0) ]

      ElementList(103)%Symbol = 'Lr'
      ElementList(103)%nat = 0
      Call mma_Allocate(ElementList(103)%Isotopes,16)
      ElementList(103)%Isotopes(:) = [
     &  Iso(262, 262.10961d0),
     &  Iso(251, 251.09418d0),
     &  Iso(252, 252.09526d0),
     &  Iso(253, 253.09509d0),
     &  Iso(254, 254.09648d0),
     &  Iso(255, 255.096562d0),
     &  Iso(256, 256.098494d0),
     &  Iso(257, 257.099418d0),
     &  Iso(258, 258.10176d0),
     &  Iso(259, 259.102902d0),
     &  Iso(260, 260.1055d0),
     &  Iso(261, 261.10688d0),
     &  Iso(263, 263.11136d0),
     &  Iso(264, 264.1142d0),
     &  Iso(265, 265.11619d0),
     &  Iso(266, 266.11983d0) ]

      ElementList(104)%Symbol = 'Rf'
      ElementList(104)%nat = 0
      Call mma_Allocate(ElementList(104)%Isotopes,16)
      ElementList(104)%Isotopes(:) = [
     &  Iso(267, 267.12179d0),
     &  Iso(253, 253.10044d0),
     &  Iso(254, 254.10005d0),
     &  Iso(255, 255.10127d0),
     &  Iso(256, 256.101152d0),
     &  Iso(257, 257.102918d0),
     &  Iso(258, 258.103428d0),
     &  Iso(259, 259.105596d0),
     &  Iso(260, 260.10644d0),
     &  Iso(261, 261.108773d0),
     &  Iso(262, 262.10992d0),
     &  Iso(263, 263.11249d0),
     &  Iso(264, 264.11388d0),
     &  Iso(265, 265.11668d0),
     &  Iso(266, 266.11817d0),
     &  Iso(268, 268.12397d0) ]

      ElementList(105)%Symbol = 'Db'
      ElementList(105)%nat = 0
      Call mma_Allocate(ElementList(105)%Isotopes,16)
      ElementList(105)%Isotopes(:) = [
     &  Iso(268, 268.12567d0),
     &  Iso(255, 255.10707d0),
     &  Iso(256, 256.10789d0),
     &  Iso(257, 257.10758d0),
     &  Iso(258, 258.10928d0),
     &  Iso(259, 259.109492d0),
     &  Iso(260, 260.1113d0),
     &  Iso(261, 261.11192d0),
     &  Iso(262, 262.11407d0),
     &  Iso(263, 263.11499d0),
     &  Iso(264, 264.11741d0),
     &  Iso(265, 265.11861d0),
     &  Iso(266, 266.12103d0),
     &  Iso(267, 267.12247d0),
     &  Iso(269, 269.12791d0),
     &  Iso(270, 270.13136d0) ]

      ElementList(106)%Symbol = 'Sg'
      ElementList(106)%nat = 0
      Call mma_Allocate(ElementList(106)%Isotopes,16)
      ElementList(106)%Isotopes(:) = [
     &  Iso(269, 269.12863d0),
     &  Iso(258, 258.11298d0),
     &  Iso(259, 259.1144d0),
     &  Iso(260, 260.114384d0),
     &  Iso(261, 261.115949d0),
     &  Iso(262, 262.116337d0),
     &  Iso(263, 263.11829d0),
     &  Iso(264, 264.11893d0),
     &  Iso(265, 265.12109d0),
     &  Iso(266, 266.12198d0),
     &  Iso(267, 267.12436d0),
     &  Iso(268, 268.12539d0),
     &  Iso(270, 270.13043d0),
     &  Iso(271, 271.13393d0),
     &  Iso(272, 272.13589d0),
     &  Iso(273, 273.13958d0) ]

      ElementList(107)%Symbol = 'Bh'
      ElementList(107)%nat = 0
      Call mma_Allocate(ElementList(107)%Isotopes,16)
      ElementList(107)%Isotopes(:) = [
     &  Iso(270, 270.13336d0),
     &  Iso(260, 260.12166d0),
     &  Iso(261, 261.12145d0),
     &  Iso(262, 262.12297d0),
     &  Iso(263, 263.12292d0),
     &  Iso(264, 264.12459d0),
     &  Iso(265, 265.12491d0),
     &  Iso(266, 266.12679d0),
     &  Iso(267, 267.1275d0),
     &  Iso(268, 268.12969d0),
     &  Iso(269, 269.13042d0),
     &  Iso(271, 271.13526d0),
     &  Iso(272, 272.13826d0),
     &  Iso(273, 273.14024d0),
     &  Iso(274, 274.14355d0),
     &  Iso(275, 275.14567d0) ]

      ElementList(108)%Symbol = 'Hs'
      ElementList(108)%nat = 0
      Call mma_Allocate(ElementList(108)%Isotopes,15)
      ElementList(108)%Isotopes(:) = [
     &  Iso(269, 269.13375d0),
     &  Iso(263, 263.12852d0),
     &  Iso(264, 264.128357d0),
     &  Iso(265, 265.129793d0),
     &  Iso(266, 266.130046d0),
     &  Iso(267, 267.13167d0),
     &  Iso(268, 268.13186d0),
     &  Iso(270, 270.13429d0),
     &  Iso(271, 271.13717d0),
     &  Iso(272, 272.1385d0),
     &  Iso(273, 273.14168d0),
     &  Iso(274, 274.1433d0),
     &  Iso(275, 275.14667d0),
     &  Iso(276, 276.14846d0),
     &  Iso(277, 277.1519d0) ]

      ElementList(109)%Symbol = 'Mt'
      ElementList(109)%nat = 0
      Call mma_Allocate(ElementList(109)%Isotopes,15)
      ElementList(109)%Isotopes(:) = [
     &  Iso(278, 278.15631d0),
     &  Iso(265, 265.136d0),
     &  Iso(266, 266.13737d0),
     &  Iso(267, 267.13719d0),
     &  Iso(268, 268.13865d0),
     &  Iso(269, 269.13882d0),
     &  Iso(270, 270.14033d0),
     &  Iso(271, 271.14074d0),
     &  Iso(272, 272.14341d0),
     &  Iso(273, 273.1444d0),
     &  Iso(274, 274.14724d0),
     &  Iso(275, 275.14882d0),
     &  Iso(276, 276.15159d0),
     &  Iso(277, 277.15327d0),
     &  Iso(279, 279.15808d0) ]

      ElementList(110)%Symbol = 'Ds'
      ElementList(110)%nat = 0
      Call mma_Allocate(ElementList(110)%Isotopes,15)
      ElementList(110)%Isotopes(:) = [
     &  Iso(281, 281.16451d0),
     &  Iso(267, 267.14377d0),
     &  Iso(268, 268.14348d0),
     &  Iso(269, 269.144752d0),
     &  Iso(270, 270.144584d0),
     &  Iso(271, 271.14595d0),
     &  Iso(272, 272.14602d0),
     &  Iso(273, 273.14856d0),
     &  Iso(274, 274.14941d0),
     &  Iso(275, 275.15203d0),
     &  Iso(276, 276.15303d0),
     &  Iso(277, 277.15591d0),
     &  Iso(278, 278.15704d0),
     &  Iso(279, 279.1601d0),
     &  Iso(280, 280.16131d0) ]

      ElementList(111)%Symbol = 'Rg'
      ElementList(111)%nat = 0
      Call mma_Allocate(ElementList(111)%Isotopes,12)
      ElementList(111)%Isotopes(:) = [
     &  Iso(281, 281.16636d0),
     &  Iso(272, 272.15327d0),
     &  Iso(273, 273.15313d0),
     &  Iso(274, 274.15525d0),
     &  Iso(275, 275.15594d0),
     &  Iso(276, 276.15833d0),
     &  Iso(277, 277.15907d0),
     &  Iso(278, 278.16149d0),
     &  Iso(279, 279.16272d0),
     &  Iso(280, 280.16514d0),
     &  Iso(282, 282.16912d0),
     &  Iso(283, 283.17054d0) ]

      ElementList(112)%Symbol = 'Cn'
      ElementList(112)%nat = 0
      Call mma_Allocate(ElementList(112)%Isotopes,10)
      ElementList(112)%Isotopes(:) = [
     &  Iso(283, 283.17327d0),
     &  Iso(276, 276.16141d0),
     &  Iso(277, 277.16364d0),
     &  Iso(278, 278.16416d0),
     &  Iso(279, 279.16654d0),
     &  Iso(280, 280.16715d0),
     &  Iso(281, 281.16975d0),
     &  Iso(282, 282.1705d0),
     &  Iso(284, 284.17416d0),
     &  Iso(285, 285.17712d0) ]

      ElementList(113)%Symbol = 'Nh'
      ElementList(113)%nat = 0
      Call mma_Allocate(ElementList(113)%Isotopes,10)
      ElementList(113)%Isotopes(:) = [
     &  Iso(287, 287.18339d0),
     &  Iso(278, 278.17058d0),
     &  Iso(279, 279.17095d0),
     &  Iso(280, 280.17293d0),
     &  Iso(281, 281.17348d0),
     &  Iso(282, 282.17567d0),
     &  Iso(283, 283.17657d0),
     &  Iso(284, 284.17873d0),
     &  Iso(285, 285.17973d0),
     &  Iso(286, 286.18221d0) ]

      ElementList(114)%Symbol = 'Fl'
      ElementList(114)%nat = 0
      Call mma_Allocate(ElementList(114)%Isotopes,5)
      ElementList(114)%Isotopes(:) = [
     &  Iso(289, 289.19042d0),
     &  Iso(285, 285.18364d0),
     &  Iso(286, 286.18423d0),
     &  Iso(287, 287.18678d0),
     &  Iso(288, 288.18757d0) ]

      ElementList(115)%Symbol = 'Mc'
      ElementList(115)%nat = 0
      Call mma_Allocate(ElementList(115)%Isotopes,5)
      ElementList(115)%Isotopes(:) = [
     &  Iso(288, 288.19274d0),
     &  Iso(287, 287.1907d0),
     &  Iso(289, 289.19363d0),
     &  Iso(290, 290.19598d0),
     &  Iso(291, 291.19707d0) ]

      ElementList(116)%Symbol = 'Lv'
      ElementList(116)%nat = 0
      Call mma_Allocate(ElementList(116)%Isotopes,5)
      ElementList(116)%Isotopes(:) = [
     &  Iso(293, 293.20449d0),
     &  Iso(289, 289.19816d0),
     &  Iso(290, 290.19864d0),
     &  Iso(291, 291.20108d0),
     &  Iso(292, 292.20174d0) ]

      ElementList(117)%Symbol = 'Ts'
      ElementList(117)%nat = 0
      Call mma_Allocate(ElementList(117)%Isotopes,4)
      ElementList(117)%Isotopes(:) = [
     &  Iso(294, 294.21046d0),
     &  Iso(291, 291.20553d0),
     &  Iso(292, 292.20746d0),
     &  Iso(293, 293.20824d0) ]

      ElementList(118)%Symbol = 'Og'
      ElementList(118)%nat = 0
      Call mma_Allocate(ElementList(118)%Isotopes,3)
      ElementList(118)%Isotopes(:) = [
     &  Iso(294, 294.21392d0),
     &  Iso(293, 293.21356d0),
     &  Iso(295, 295.21624d0) ]

      End Subroutine Initialize_Isotopes

*
* This subroutine frees up the memory
*
      Subroutine Free_Isotopes()
      Integer :: i
      If (.Not. Allocated(ElementList)) Return
      Do i=1,Size(ElementList,1)
        Call mma_Deallocate(ElementList(i)%Isotopes)
      End Do
      Call mma_Deallocate(ElementList)
#ifdef _WARNING_WORKAROUND_
      If (.False.) Then
*       Since this should never be executed, don't deallocate
        Call mma_Allocate(ElementList(1)%Isotopes,[0,0])
        Call mma_Allocate(ElementList,[0,0])
      End If
#endif
      End Subroutine Free_Isotopes

*
* Subroutine(s) to get the Mass of the isotope IsNr belonging to the
* element Atom. If IsNr=0, the most abundant isotope (or the most
* stable if all are radioactive) is selected. The mass is returned
* in atomic units (m_e).
* Atom can be an atomic symbol or an atomic number.
*
      Subroutine Isotope_sym(IsNr, Atom, Mass)
      Integer, Intent(InOut) :: IsNr
      Character(Len=2), Intent(In) :: Atom
      Real*8, Intent(Out) :: Mass
      Integer :: i, This
      Character(Len=2) :: Sym, Sym2

      Call Initialize_Isotopes()

      Sym2 = AdjustL(Atom)
      Call UpCase(Sym2)
      If ((Sym2 .eq. 'D') .or. (Sym2 .eq. 'T')) Sym2 = 'H'
      This = 0
      Do i=1,MaxAtomNum
        Sym = AdjustL(ElementList(i)%Symbol)
        Call UpCase(Sym)
        If (Sym .eq. Sym2) Then
          This = i
          Exit
        End If
      End Do

      If (This .eq. 0) Then
        Write (6,*) 'Isotope: Did not find atom!'
        Write (6,*) 'Atom=',Atom
        Call Abend()
      End If

      If (IsNr .eq. 0) IsNr = ElementList(This)%Isotopes(1)%A
      If (Sym2 .eq. 'D') IsNr = 2
      If (Sym2 .eq. 'T') IsNr = 3
      Do i=1,Size(ElementList(This)%Isotopes,1)
        If (ElementList(This)%Isotopes(i)%A .eq. IsNr) Then
          Mass = uToau*ElementList(This)%Isotopes(i)%m
          Return
        End If
      End Do

      Write (6,*) 'Isotope: Did not find isotope!'
      Write (6,*) 'IsNr=',IsNr
      Write (6,*) 'Atom=',Atom
      Call Abend()

      End Subroutine Isotope_sym

      Subroutine Isotope_num(IsNr, Atom, Mass)
      Integer, Intent(InOut) :: IsNr
      Integer, Intent(In) :: Atom
      Real*8, Intent(Out) :: Mass
      Integer :: i

      Call Initialize_Isotopes()

      If ((Atom .lt. 0) .or. (Atom .gt. MaxAtomNum)) Then
        Write (6,*) 'Isotope: Did not find atom!'
        Write (6,*) 'Atom=',Atom
        Call Abend()
      End If

      If (IsNr .eq. 0) IsNr = ElementList(Atom)%Isotopes(1)%A
      Do i=1,Size(ElementList(Atom)%Isotopes,1)
        If (ElementList(Atom)%Isotopes(i)%A .eq. IsNr) Then
          Mass = uToau*ElementList(Atom)%Isotopes(i)%m
          Return
        End If
      End Do

      Write (6,*) 'Isotope: Did not find isotope!'
      Write (6,*) 'IsNr=',IsNr
      Write (6,*) 'Atom=',Atom
      Call Abend()

      End Subroutine Isotope_num

*
* Function that returns the mass in atomic units (m_e) of a particular
* nuclide with Z protons and A-Z neutrons. Returns -1.0 if the nuclide
* is unknown.
*
      Function NuclideMass(Z, A)
      Integer, Intent(In) :: Z, A
      Real*8 :: NuclideMass
      Integer :: i

      Call Initialize_Isotopes()

      NuclideMass = -1.0d0
      If ((Z < 1) .or. (Z .gt. MaxAtomNum)) Return
      Do i=1,Size(ElementList(Z)%Isotopes,1)
        If (ElementList(Z)%Isotopes(i)%A .ne. A) Cycle
        NuclideMass = uToau*ElementList(Z)%Isotopes(i)%m
        Exit
      End Do

      End Function NuclideMass

*
* Private extensions to mma_interfaces, using preprocessor templates
* (see src/mma_util/stdalloc.f)
*

* Define elm_cptr2loff, element_mma_allo_1D, element_mma_allo_1D_lim, element_mma_free_1D
#define _TYPE_ type(element)
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

* Define iso_cptr2loff, isotope_mma_allo_1D, isotope_mma_allo_1D_lim, isotope_mma_free_1D
#define _TYPE_ type(iso)
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

      End Module Isotopes
