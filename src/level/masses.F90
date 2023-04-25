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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine MASSES(IAN,IMN,ANAME,GELGS,GNS,MASS,ABUND)
!***********************************************************************
!** For isotope with (input) atomic number IAN and mass number IMN,
!  return (output):  (i) as the right-adjusted 2-character variable ANAME
!  the alphabetic symbol for that element,  (ii) the ground state
!  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy GNS,
!  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
!  abundance ABUND [in percent].   GELGS values based on atomic states
!  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
!  from the 2003 mass table [Audi, Wapstra & Thibault, Nucl.Phys. A729,
!  337-676 (2003)] and other quantities from Tables 6.2 and 6.3 of
!  "Quantities, Units and Symbols in Physical Chemistry", by Mills et
!  al. (Blackwell, 2'nd Edition, Oxford, 1993).
!** If the input value of IMN does not equal one of the tabulated values
!  for atomic species IAN, return the abundance-averaged standard atomic
!  weight of that atom and set GNS=-1 and ABUND=-1.
!** Before Git tracking, contributions were from:
!** N Dattani, G T Kraemer, R J Le Roy, J Y Seto
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ian, imn, gelgs, gns
character(len=2) :: ANAME
real(kind=wp) :: mass, abund
integer(kind=iwp), parameter :: numel = 123, numiso = 10
integer(kind=iwp) :: gel(numel), i, mn(numel,numiso), nmn(numel), ns2(numel,numiso)
real(kind=wp) :: ab(numel,numiso), zm(numel,0:numiso)
character(len=2) :: AT(numel)

data at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
data(zm(1,i),i=0,3)/1.00794_wp,1.00782503207_wp,2.0141017778_wp,3.0160492777_wp/
data(ns2(1,i),i=1,3)/1,2,1/
data(ab(1,i),i=1,3)/99.985_wp,0.015_wp,Zero/

data at(2),gel(2),nmn(2),(mn(2,i),i=1,2)/'He',1,2,3,4/
data(zm(2,i),i=0,2)/4.002602_wp,3.0160293191_wp,4.00260325415_wp/
data(ns2(2,i),i=1,2)/1,0/
data(ab(2,i),i=1,2)/0.000137_wp,99.999863_wp/

data at(3),gel(3),nmn(3),(mn(3,i),i=1,2)/'Li',2,2,6,7/
data(zm(3,i),i=0,2)/6.941_wp,6.015122795_wp,7.01600455_wp/
data(ns2(3,i),i=1,2)/2,3/
data(ab(3,i),i=1,2)/7.5_wp,92.5_wp/

data at(4),gel(4),nmn(4),(mn(4,i),i=1,1)/'Be',1,1,9/
data(zm(4,i),i=0,1)/9.012182_wp,9.0121822_wp/
data(ns2(4,i),i=1,1)/3/
data(ab(4,i),i=1,1)/100.0_wp/

data at(5),gel(5),nmn(5),(mn(5,i),i=1,2)/' B',2,2,10,11/
data(zm(5,i),i=0,2)/10.811_wp,10.0129370_wp,11.0093054_wp/
data(ns2(5,i),i=1,2)/6,3/
data(ab(5,i),i=1,2)/19.9_wp,80.1_wp/

data at(6),gel(6),nmn(6),(mn(6,i),i=1,3)/' C',1,3,12,13,14/
data(zm(6,i),i=0,3)/12.011_wp,12.0_wp,13.0033548378_wp,14.003241989_wp/
data(ns2(6,i),i=1,3)/0,1,0/
data(ab(6,i),i=1,3)/98.90_wp,1.10_wp,Zero/

data at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
data(zm(7,i),i=0,2)/14.00674_wp,14.0030740048_wp,15.0001088982_wp/
data(ns2(7,i),i=1,2)/2,1/
data(ab(7,i),i=1,2)/99.634_wp,0.366_wp/

data at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
data(zm(8,i),i=0,3)/15.9994_wp,15.99491461956_wp,16.99913170_wp,17.9991610_wp/
data(ns2(8,i),i=1,3)/0,5,0/
data(ab(8,i),i=1,3)/99.762_wp,0.038_wp,0.200_wp/

data at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
data(zm(9,i),i=0,1)/18.9984032_wp,18.99840322_wp/
data(ns2(9,i),i=1,1)/1/
data(ab(9,i),i=1,1)/100.0_wp/

data at(10),gel(10),nmn(10),(mn(10,i),i=1,3)/'Ne',1,3,20,21,22/
data(zm(10,i),i=0,3)/20.1797_wp,19.9924401754_wp,20.99384668_wp,21.991385114_wp/
data(ns2(10,i),i=1,3)/0,3,0/
data(ab(10,i),i=1,3)/90.48_wp,0.27_wp,9.25_wp/

data at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
data(zm(11,i),i=0,1)/22.989768_wp,22.9897692809_wp/
data(ns2(11,i),i=1,1)/3/
data(ab(11,i),i=1,1)/100.0_wp/

data at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
data(zm(12,i),i=0,3)/24.3050_wp,23.985041700_wp,24.98583692_wp,25.982592929_wp/
data(ns2(12,i),i=1,3)/0,5,0/
data(ab(12,i),i=1,3)/78.99_wp,10.00_wp,11.01_wp/

data at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
data(zm(13,i),i=0,1)/26.981539_wp,26.98153863_wp/
data(ns2(13,i),i=1,1)/5/
data(ab(13,i),i=1,1)/100.0_wp/

data at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
data(zm(14,i),i=0,3)/28.0855_wp,27.9769265325_wp,28.976494700_wp,29.97377017_wp/
data(ns2(14,i),i=1,3)/0,1,0/
data(ab(14,i),i=1,3)/92.23_wp,4.67_wp,3.10_wp/

data at(15),gel(15),nmn(15),(mn(15,i),i=1,1)/' P',4,1,31/
data(zm(15,i),i=0,1)/30.973762_wp,30.97376163_wp/
data(ns2(15,i),i=1,1)/1/
data(ab(15,i),i=1,1)/100.0_wp/

data at(16),gel(16),nmn(16),(mn(16,i),i=1,4)/' S',5,4,32,33,34,36/
data(zm(16,i),i=0,4)/32.066_wp,31.97207100_wp,32.97145876_wp,33.96786690_wp,35.96708076_wp/
data(ns2(16,i),i=1,4)/0,3,0,0/
data(ab(16,i),i=1,4)/95.02_wp,0.75_wp,4.21_wp,0.02_wp/

data at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
data(zm(17,i),i=0,2)/35.4527_wp,34.96885268_wp,36.96590259_wp/
data(ns2(17,i),i=1,2)/3,3/
data(ab(17,i),i=1,2)/75.77_wp,24.23_wp/

data at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
data(zm(18,i),i=0,3)/39.948_wp,35.967545106_wp,37.9627324_wp,39.9623831225_wp/
data(ns2(18,i),i=1,3)/0,0,0/
data(ab(18,i),i=1,3)/0.337_wp,0.063_wp,99.600_wp/

data at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
data(zm(19,i),i=0,3)/39.0983_wp,38.96370668_wp,39.96399848_wp,40.96182576_wp/
data(ns2(19,i),i=1,3)/3,8,3/
data(ab(19,i),i=1,3)/93.2581_wp,0.0117_wp,6.7302_wp/

data at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,46,48/
data(zm(20,i),i=0,6)/40.078_wp,39.96259098_wp,41.95861801_wp,42.9587666_wp,43.9554818_wp,45.9536926_wp,47.952534_wp/
data(ns2(20,i),i=1,6)/0,0,7,0,0,0/
data(ab(20,i),i=1,6)/96.941_wp,0.647_wp,0.135_wp,2.086_wp,0.004_wp,0.187_wp/

data at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
data(zm(21,i),i=0,1)/44.955910_wp,44.9559119_wp/
data(ns2(21,i),i=1,1)/7/
data(ab(21,i),i=1,1)/100.0_wp/

data at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,50/
data(zm(22,i),i=0,5)/47.88_wp,45.9526316_wp,46.9517631_wp,47.9479463_wp,48.9478700_wp,49.9447912_wp/
data(ns2(22,i),i=1,5)/0,5,0,7,0/
data(ab(22,i),i=1,5)/8.0_wp,7.3_wp,73.8_wp,5.5_wp,5.4_wp/

data at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
data(zm(23,i),i=0,2)/50.9415_wp,49.9471585_wp,50.9439595_wp/
data(ns2(23,i),i=1,2)/12,7/
data(ab(23,i),i=1,2)/0.250_wp,99.750_wp/

data at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
data(zm(24,i),i=0,4)/51.9961_wp,49.9460442_wp,51.9405075_wp,52.9406494_wp,53.9388804_wp/
data(ns2(24,i),i=1,4)/0,0,3,0/
data(ab(24,i),i=1,4)/4.345_wp,83.789_wp,9.501_wp,2.365_wp/

data at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
data(zm(25,i),i=0,1)/54.93805_wp,54.9380451_wp/
data(ns2(25,i),i=1,1)/5/
data(ab(25,i),i=1,1)/100.0_wp/

data at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
data(zm(26,i),i=0,4)/55.847_wp,53.9396105_wp,55.9349375_wp,56.9353940_wp,57.9332756_wp/
data(ns2(26,i),i=1,4)/0,0,1,0/
data(ab(26,i),i=1,4)/5.8_wp,91.72_wp,2.2_wp,0.28_wp/

data at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
data(zm(27,i),i=0,1)/58.93320_wp,58.9331950_wp/
data(ns2(27,i),i=1,1)/7/
data(ab(27,i),i=1,1)/100.0_wp/

data at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,64/
data(zm(28,i),i=0,5)/58.69_wp,57.9353429_wp,59.9307864_wp,60.9310560_wp,61.9283451_wp,63.9279660_wp/
data(ns2(28,i),i=1,5)/0,0,3,0,0/
data(ab(28,i),i=1,5)/68.077_wp,26.223_wp,1.140_wp,3.634_wp,0.926_wp/

data at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
data(zm(29,i),i=0,2)/63.546_wp,62.9295975_wp,64.9277895_wp/
data(ns2(29,i),i=1,2)/3,3/
data(ab(29,i),i=1,2)/69.17_wp,30.83_wp/

data at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,70/
data(zm(30,i),i=0,5)/65.40_wp,63.9291422_wp,65.9260334_wp,66.9271273_wp,67.9248442_wp,69.9253193_wp/
data(ns2(30,i),i=1,5)/0,0,5,0,0/
data(ab(30,i),i=1,5)/48.6_wp,27.9_wp,4.1_wp,18.8_wp,0.6_wp/

data at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
data(zm(31,i),i=0,2)/69.723_wp,68.9255736_wp,70.9247013_wp/
data(ns2(31,i),i=1,2)/3,3/
data(ab(31,i),i=1,2)/60.108_wp,39.892_wp/

data at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,76/
data(zm(32,i),i=0,5)/72.61_wp,69.9242474_wp,71.9220758_wp,72.9234589_wp,73.9211778_wp,75.9214026_wp/
data(ns2(32,i),i=1,5)/0,0,9,0,0/
data(ab(32,i),i=1,5)/21.23_wp,27.66_wp,7.73_wp,35.94_wp,7.44_wp/

data at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
data(zm(33,i),i=0,1)/74.92159_wp,74.9215965_wp/
data(ns2(33,i),i=1,1)/3/
data(ab(33,i),i=1,1)/100.0_wp/

data at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,80,82/
data(zm(34,i),i=0,6)/78.96_wp,73.9224764_wp,75.9192136_wp,76.9199140_wp,77.9173091_wp,79.9165213_wp,81.9166994_wp/
data(ns2(34,i),i=1,6)/0,0,1,0,0,0/
data(ab(34,i),i=1,6)/0.89_wp,9.36_wp,7.63_wp,23.78_wp,49.61_wp,8.73_wp/

data at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
data(zm(35,i),i=0,2)/79.904_wp,78.9183371_wp,80.9162906_wp/
data(ns2(35,i),i=1,2)/3,3/
data(ab(35,i),i=1,2)/50.69_wp,49.31_wp/

data at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,84,86/
data(zm(36,i),i=0,6)/83.80_wp,77.9203648_wp,79.9163790_wp,81.9134836_wp,82.914136_wp,83.911507_wp,85.91061073_wp/
data(ns2(36,i),i=1,6)/0,0,0,9,0,0/
data(ab(36,i),i=1,6)/0.35_wp,2.25_wp,11.6_wp,11.5_wp,57.0_wp,17.3_wp/

data at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
data(zm(37,i),i=0,2)/85.4678_wp,84.911789738_wp,86.909180527_wp/
data(ns2(37,i),i=1,2)/5,3/
data(ab(37,i),i=1,2)/72.165_wp,27.835_wp/

data at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
data(zm(38,i),i=0,4)/87.62_wp,83.913425_wp,85.9092602_wp,86.9088771_wp,87.9056121_wp/
data(ns2(38,i),i=1,4)/0,0,9,0/
data(ab(38,i),i=1,4)/0.56_wp,9.86_wp,7.00_wp,82.58_wp/

data at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
data(zm(39,i),i=0,1)/88.90585_wp,88.9058483_wp/
data(ns2(39,i),i=1,1)/1/
data(ab(39,i),i=1,1)/100.0_wp/

data at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,96/
data(zm(40,i),i=0,5)/91.224_wp,89.9047044_wp,90.9056458_wp,91.9050408_wp,93.9063152_wp,95.9082734_wp/
data(ns2(40,i),i=1,5)/0,5,0,0,0/
data(ab(40,i),i=1,5)/51.45_wp,11.22_wp,17.15_wp,17.38_wp,2.80_wp/

data at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
data(zm(41,i),i=0,1)/92.90638_wp,92.9063781_wp/
data(ns2(41,i),i=1,1)/9/
data(ab(41,i),i=1,1)/100.0_wp/

data at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,97,98,100/
data(zm(42,i),i=0,7)/95.94_wp,91.906811_wp,93.9050883_wp,94.9058421_wp,95.9046795_wp,96.9060215_wp,97.9054082_wp,99.907477_wp/
data(ns2(42,i),i=1,7)/0,0,5,0,5,0,0/
data(ab(42,i),i=1,7)/14.84_wp,9.25_wp,15.92_wp,16.68_wp,9.55_wp,24.13_wp,9.63_wp/

data at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
data(zm(43,i),i=0,1)/97.907215_wp,97.907216_wp/
data(ns2(43,i),i=1,1)/12/
data(ab(43,i),i=1,1)/100.0_wp/

data at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,100,101,102,104/
data(zm(44,i),i=0,7)/101.07_wp,95.907598_wp,97.905287_wp,98.9059393_wp,99.9042195_wp,100.9055821_wp,101.9043493_wp,103.905433_wp/
data(ns2(44,i),i=1,7)/0,0,5,0,5,0,0/
data(ab(44,i),i=1,7)/5.52_wp,1.88_wp,12.7_wp,12.6_wp,17.0_wp,31.6_wp,18.7_wp/

data at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
data(zm(45,i),i=0,1)/102.90550_wp,102.905504_wp/
data(ns2(45,i),i=1,1)/1/
data(ab(45,i),i=1,1)/100.0_wp/

data at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,106,108,110/
data(zm(46,i),i=0,6)/106.42_wp,101.905609_wp,103.904036_wp,104.905085_wp,105.903486_wp,107.903892_wp,109.905153_wp/
data(ns2(46,i),i=1,6)/0,0,5,0,0,0/
data(ab(46,i),i=1,6)/1.02_wp,11.14_wp,22.33_wp,27.33_wp,26.46_wp,11.72_wp/

data at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
data(zm(47,i),i=0,2)/107.8682_wp,106.905097_wp,108.904752_wp/
data(ns2(47,i),i=1,2)/1,1/
data(ab(47,i),i=1,2)/51.839_wp,48.161_wp/

data at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,111,112,113,114,116/
data(zm(48,i),i=0,8)/112.411_wp,105.906459_wp,107.904184_wp,109.9030021_wp,110.9041781_wp,111.9027578_wp,112.9044017_wp, &
                     113.9033585_wp,115.904756_wp/
data(ns2(48,i),i=1,8)/0,0,0,1,0,1,0,0/
data(ab(48,i),i=1,8)/1.25_wp,0.89_wp,12.49_wp,12.80_wp,24.13_wp,12.22_wp,28.73_wp,7.49_wp/

data at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
data(zm(49,i),i=0,2)/114.818_wp,112.904058_wp,114.903878_wp/
data(ns2(49,i),i=1,2)/9,9/
data(ab(49,i),i=1,2)/4.3_wp,95.7_wp/

data at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,115,116,117,118,119,120,122,124/
data(zm(50,i),i=0,10)/118.710_wp,111.904818_wp,113.902779_wp,114.903342_wp,115.901741_wp,116.902952_wp,117.901603_wp, &
                      118.903308_wp,119.9021947_wp,121.9034390_wp,123.9052739_wp/
data(ns2(50,i),i=1,10)/0,0,1,0,1,0,1,0,0,0/
data(ab(50,i),i=1,10)/0.97_wp,0.65_wp,0.34_wp,14.53_wp,7.68_wp,24.23_wp,8.59_wp,32.59_wp,4.63_wp,5.79_wp/

data at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
data(zm(51,i),i=0,2)/121.757_wp,120.9038157_wp,122.9042140_wp/
data(ns2(51,i),i=1,2)/5,7/
data(ab(51,i),i=1,2)/57.36_wp,42.64_wp/

data at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,124,125,126,128,130/
data(zm(52,i),i=0,8)/127.60_wp,119.904020_wp,121.9030439_wp,122.9042700_wp,123.9028179_wp,124.9044307_wp,125.9033117_wp, &
                     127.9044631_wp,129.9062244_wp/
data(ns2(52,i),i=1,8)/0,0,1,0,1,0,0,0/
data(ab(52,i),i=1,8)/0.096_wp,2.603_wp,0.908_wp,4.816_wp,7.139_wp,18.95_wp,31.69_wp,33.80_wp/

data at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
data(zm(53,i),i=0,2)/126.90447_wp,126.904473_wp,128.904988_wp/
data(ns2(53,i),i=1,2)/5,7/
data(ab(53,i),i=1,2)/100.0_wp,Zero/

data at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,129,130,131,132,134,136/
data(zm(54,i),i=0,9)/131.29_wp,123.9058930_wp,125.904274_wp,127.9035313_wp,128.9047794_wp,129.9035080_wp,130.9050824_wp, &
                     131.9041535_wp,133.9053945_wp,135.907219_wp/
data(ns2(54,i),i=1,9)/0,0,0,1,0,3,0,0,0/
data(ab(54,i),i=1,9)/0.10_wp,0.09_wp,1.91_wp,26.4_wp,4.1_wp,21.2_wp,26.9_wp,10.4_wp,8.9_wp/

data at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
data(zm(55,i),i=0,1)/132.90543_wp,132.905451933_wp/
data(ns2(55,i),i=1,1)/7/
data(ab(55,i),i=1,1)/100.0_wp/

data at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,135,136,137,138/
data(zm(56,i),i=0,7)/137.327_wp,129.9063208_wp,131.9050613_wp,133.9045084_wp,134.9056886_wp,135.9045759_wp,136.9058274_wp, &
                     137.9052472_wp/
data(ns2(56,i),i=1,7)/0,0,0,3,0,3,0/
data(ab(56,i),i=1,7)/0.106_wp,0.101_wp,2.417_wp,6.592_wp,7.854_wp,11.23_wp,71.70_wp/

data at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
data(zm(57,i),i=0,2)/138.9055_wp,137.907112_wp,138.9063533_wp/
data(ns2(57,i),i=1,2)/10,7/
data(ab(57,i),i=1,2)/0.0902_wp,99.9098_wp/

data at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,142/
data(zm(58,i),i=0,4)/140.115_wp,135.907172_wp,137.905991_wp,139.9054387_wp,141.909244_wp/
data(ns2(58,i),i=1,4)/0,0,0,0/
data(ab(58,i),i=1,4)/0.19_wp,0.25_wp,88.48_wp,11.08_wp/

data at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
data(zm(59,i),i=0,1)/140.90765_wp,140.9076528_wp/
data(ns2(59,i),i=1,1)/5/
data(ab(59,i),i=1,1)/100.0_wp/

data at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,145,146,148,150/
data(zm(60,i),i=0,7)/144.24_wp,141.9077233_wp,142.9098143_wp,143.9100873_wp,144.9125736_wp,145.9131169_wp,147.916893_wp, &
                     149.920891_wp/
data(ns2(60,i),i=1,7)/0,7,0,7,0,0,0/
data(ab(60,i),i=1,7)/27.13_wp,12.18_wp,23.80_wp,8.30_wp,17.19_wp,5.76_wp,5.64_wp/

data at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
data(zm(61,i),i=0,1)/144.912743_wp,144.912749_wp/
data(ns2(61,i),i=1,1)/5/
data(ab(61,i),i=1,1)/100.0_wp/

data at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,149,150,152,154/
data(zm(62,i),i=0,7)/150.36_wp,143.911999_wp,146.9148979_wp,147.9148227_wp,148.9171847_wp,149.9172755_wp,151.9197324_wp, &
                     153.9222093_wp/
data(ns2(62,i),i=1,7)/0,7,0,7,0,0,0/
data(ab(62,i),i=1,7)/3.1_wp,15.0_wp,11.3_wp,13.8_wp,7.4_wp,26.7_wp,22.7_wp/

data at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
data(zm(63,i),i=0,2)/151.965_wp,150.9198502_wp,152.9212303_wp/
data(ns2(63,i),i=1,2)/5,5/
data(ab(63,i),i=1,2)/47.8_wp,52.2_wp/

data at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,156,157,158,160/
data(zm(64,i),i=0,7)/157.25_wp,151.9197910_wp,153.92086560_wp,154.9226220_wp,155.9221227_wp,156.9239601_wp,157.9241039_wp, &
                     159.9270541_wp/
data(ns2(64,i),i=1,7)/0,0,3,0,3,0,0/
data(ab(64,i),i=1,7)/0.20_wp,2.18_wp,14.80_wp,20.47_wp,15.65_wp,24.84_wp,21.86_wp/

data at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
data(zm(65,i),i=0,1)/158.92534_wp,158.9253468_wp/
data(ns2(65,i),i=1,1)/3/
data(ab(65,i),i=1,1)/100.0_wp/

data at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,160,161,162,163,164/
data(zm(66,i),i=0,7)/162.50_wp,155.924283_wp,157.924409_wp,159.9251975_wp,160.9269334_wp,161.9267984_wp,162.9287312_wp, &
                     163.9291748_wp/
data(ns2(66,i),i=1,7)/0,0,0,5,0,5,0/
data(ab(66,i),i=1,7)/0.06_wp,0.10_wp,2.34_wp,18.9_wp,25.5_wp,24.9_wp,28.2_wp/

data at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
data(zm(67,i),i=0,1)/164.93032_wp,164.9303221_wp/
data(ns2(67,i),i=1,1)/7/
data(ab(67,i),i=1,1)/100.0_wp/

data at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,166,167,168,170/
data(zm(68,i),i=0,6)/167.26_wp,161.928778_wp,163.929200_wp,165.9302931_wp,166.9320482_wp,167.9323702_wp,169.9354643_wp/
data(ns2(68,i),i=1,6)/0,0,0,7,0,0/
data(ab(68,i),i=1,6)/0.14_wp,1.61_wp,33.6_wp,22.95_wp,26.8_wp,14.9_wp/

data at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/
data(zm(69,i),i=0,1)/168.93421_wp,168.9342133_wp/
data(ns2(69,i),i=1,1)/1/
data(ab(69,i),i=1,1)/100.0_wp/

data at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,172,173,174,176/
data(zm(70,i),i=0,7)/173.04_wp,167.933897_wp,169.9347618_wp,170.936323580_wp,171.9363815_wp,172.9382108_wp,173.9388621_wp, &
                     175.9425717_wp/
data(ns2(70,i),i=1,7)/0,0,1,0,5,0,0/
data(ab(70,i),i=1,7)/0.13_wp,3.05_wp,14.3_wp,21.9_wp,16.12_wp,31.8_wp,12.7_wp/

data at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
data(zm(71,i),i=0,2)/174.967_wp,174.9407718_wp,175.9426863_wp/
data(ns2(71,i),i=1,2)/7,14/
data(ab(71,i),i=1,2)/97.41_wp,2.59_wp/

data at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,178,179,180/
data(zm(72,i),i=0,6)/178.49_wp,173.940046_wp,175.9414086_wp,176.9432207_wp,177.9436988_wp,178.9458161_wp,179.9465500_wp/
data(ns2(72,i),i=1,6)/0,0,7,0,9,0/
data(ab(72,i),i=1,6)/0.162_wp,5.206_wp,18.606_wp,27.297_wp,13.629_wp,35.100_wp/

data at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
data(zm(73,i),i=0,2)/180.9479_wp,179.9474648_wp,180.9479958_wp/
data(ns2(73,i),i=1,2)/16,7/
data(ab(73,i),i=1,2)/0.012_wp,99.988_wp/

data at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,184,186/
data(zm(74,i),i=0,5)/183.84_wp,179.946704_wp,181.9482042_wp,182.9502230_wp,183.9509312_wp,185.9543641_wp/
data(ns2(74,i),i=1,5)/0,0,1,0,0/
data(ab(74,i),i=1,5)/0.13_wp,26.3_wp,14.3_wp,30.67_wp,28.6_wp/

data at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
data(zm(75,i),i=0,2)/186.207_wp,184.9529550_wp,186.9557531_wp/
data(ns2(75,i),i=1,2)/5,5/
data(ab(75,i),i=1,2)/37.40_wp,62.60_wp/

data at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,188,189,190,192/
data(zm(76,i),i=0,7)/190.23_wp,183.9524891_wp,185.9538382_wp,186.9557505_wp,187.9558382_wp,188.9581475_wp,189.9584470_wp, &
                     191.9614807_wp/
data(ns2(76,i),i=1,7)/0,0,1,0,3,0,0/
data(ab(76,i),i=1,7)/0.02_wp,1.58_wp,1.6_wp,13.3_wp,16.1_wp,26.4_wp,41.0_wp/

data at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
data(zm(77,i),i=0,2)/192.22_wp,190.9605940_wp,192.9629264_wp/
data(ns2(77,i),i=1,2)/3,3/
data(ab(77,i),i=1,2)/37.3_wp,62.7_wp/

data at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,195,196,198/
data(zm(78,i),i=0,6)/195.08_wp,189.959932_wp,191.9610380_wp,193.9626803_wp,194.9647911_wp,195.9649515_wp,197.967893_wp/
data(ns2(78,i),i=1,6)/0,0,0,1,0,0/
data(ab(78,i),i=1,6)/0.01_wp,0.79_wp,32.9_wp,33.8_wp,25.3_wp,7.2_wp/

data at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
data(zm(79,i),i=0,1)/196.96654_wp,196.9665687_wp/
data(ns2(79,i),i=1,1)/3/
data(ab(79,i),i=1,1)/100.0_wp/

data at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,200,201,202,204/
data(zm(80,i),i=0,7)/200.59_wp,195.965833_wp,197.9667690_wp,198.9682799_wp,199.9683260_wp,200.9703023_wp,201.9706430_wp, &
                     203.9734939_wp/
data(ns2(80,i),i=1,7)/0,0,1,0,3,0,0/
data(ab(80,i),i=1,7)/0.15_wp,9.97_wp,16.87_wp,23.10_wp,13.18_wp,29.86_wp,6.87_wp/

data at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
data(zm(81,i),i=0,2)/204.3833_wp,202.9723442_wp,204.9744275_wp/
data(ns2(81,i),i=1,2)/1,1/
data(ab(81,i),i=1,2)/29.524_wp,70.476_wp/

data at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,208/
data(zm(82,i),i=0,4)/207.2_wp,203.9730436_wp,205.9744653_wp,206.9758969_wp,207.9766521_wp/
data(ns2(82,i),i=1,4)/0,0,1,0/
data(ab(82,i),i=1,4)/1.4_wp,24.1_wp,22.1_wp,52.4_wp/

data at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
data(zm(83,i),i=0,1)/208.98037_wp,208.9803987_wp/
data(ns2(83,i),i=1,1)/9/
data(ab(83,i),i=1,1)/100.0_wp/

data at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
data(zm(84,i),i=0,1)/208.982404_wp,208.9824304_wp/
data(ns2(84,i),i=1,1)/1/
data(ab(84,i),i=1,1)/100.0_wp/

data at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
data(zm(85,i),i=0,1)/209.987126_wp,209.987148_wp/
data(ns2(85,i),i=1,1)/10/
data(ab(85,i),i=1,1)/100.0_wp/

data at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
data(zm(86,i),i=0,1)/222.017571_wp,222.0175777_wp/
data(ns2(86,i),i=1,1)/0/
data(ab(86,i),i=1,1)/100.0_wp/

data at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
data(zm(87,i),i=0,1)/223.019733_wp,223.0197359_wp/
data(ns2(87,i),i=1,1)/3/
data(ab(87,i),i=1,1)/100.0_wp/

data at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
data(zm(88,i),i=0,1)/226.025403_wp,226.0254098_wp/
data(ns2(88,i),i=1,1)/0/
data(ab(88,i),i=1,1)/100.0_wp/

data at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
data(zm(89,i),i=0,1)/227.027750_wp,227.0277521_wp/
data(ns2(89,i),i=1,1)/3/
data(ab(89,i),i=1,1)/100.0_wp/

data at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
data(zm(90,i),i=0,1)/232.038_wp,232.0380553_wp/
data(ns2(90,i),i=1,1)/0/
data(ab(90,i),i=1,1)/100.0_wp/

data at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
data(zm(91,i),i=0,1)/231.03588_wp,231.0358840_wp/
data(ns2(91,i),i=1,1)/3/
data(ab(91,i),i=1,1)/100.0_wp/

data at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,235,238/
data(zm(92,i),i=0,4)/238.0289_wp,233.0396352_wp,234.0409521_wp,235.0439299_wp,238.0507882_wp/
data(ns2(92,i),i=1,4)/5,0,7,0/
data(ab(92,i),i=1,4)/Zero,0.0055_wp,0.7200_wp,99.2745_wp/

data at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
data(zm(93,i),i=0,1)/237.0481678_wp,237.0481734_wp/
data(ns2(93,i),i=1,1)/5/
data(ab(93,i),i=1,1)/100.0_wp/

data at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
data(zm(94,i),i=0,1)/244.064199_wp,244.064204_wp/
data(ns2(94,i),i=1,1)/0/
data(ab(94,i),i=1,1)/100.0_wp/

data at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
data(zm(95,i),i=0,1)/243.061375_wp,243.0613811_wp/
data(ns2(95,i),i=1,1)/5/
data(ab(95,i),i=1,1)/100.0_wp/

data at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
data(zm(96,i),i=0,1)/247.070347_wp,247.070354_wp/
data(ns2(96,i),i=1,1)/9/
data(ab(96,i),i=1,1)/100.0_wp/

data at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
data(zm(97,i),i=0,1)/247.070300_wp,247.070307_wp/
data(ns2(97,i),i=1,1)/3/
data(ab(97,i),i=1,1)/100.0_wp/

data at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
data(zm(98,i),i=0,1)/251.079580_wp,251.079587_wp/
data(ns2(98,i),i=1,1)/1/
data(ab(98,i),i=1,1)/100.0_wp/

data at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
data(zm(99,i),i=0,1)/252.082944_wp,252.082980_wp/
data(ns2(99,i),i=1,1)/10/
data(ab(99,i),i=1,1)/100.0_wp/

data at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
data(zm(100,i),i=0,1)/257.095099_wp,257.095105_wp/
data(ns2(100,i),i=1,1)/9/
data(ab(100,i),i=1,1)/100.0_wp/

data at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
data(zm(101,i),i=0,1)/258.09857_wp,258.098431_wp/
data(ns2(101,i),i=1,1)/16/
data(ab(101,i),i=1,1)/100.0_wp/

data at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
data(zm(102,i),i=0,1)/259.100931_wp,259.101030_wp/
data(ns2(102,i),i=1,1)/9/
data(ab(102,i),i=1,1)/100.0_wp/

data at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
data(zm(103,i),i=0,1)/260.105320_wp,260.105500_wp/
data(ns2(103,i),i=1,1)/-1/
data(ab(103,i),i=1,1)/100.0_wp/

data at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
data(zm(104,i),i=0,1)/261.10869_wp,261.108770_wp/
data(ns2(104,i),i=1,1)/-1/
data(ab(104,i),i=1,1)/100.0_wp/

data at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
data(zm(105,i),i=0,1)/262.11376_wp,262.114080_wp/
data(ns2(105,i),i=1,1)/-1/
data(ab(105,i),i=1,1)/100.0_wp/

data at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
data(zm(106,i),i=0,1)/263.11822_wp,263.118320_wp/
data(ns2(106,i),i=1,1)/-1/
data(ab(106,i),i=1,1)/100.0_wp/

data at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
data(zm(107,i),i=0,1)/262.12293_wp,262.122890_wp/
data(ns2(107,i),i=1,1)/-1/
data(ab(107,i),i=1,1)/100.0_wp/

data at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
data(zm(108,i),i=0,1)/265.13016_wp,265.130090_wp/
data(ns2(108,i),i=1,1)/-1/
data(ab(108,i),i=1,1)/100.0_wp/

data at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
data(zm(109,i),i=0,1)/266.13764_wp,266.137300_wp/
data(ns2(109,i),i=1,1)/-1/
data(ab(109,i),i=1,1)/100.0_wp/

if ((IAN <= 0) .or. (IAN > 109)) then
  MASS = Zero
  ANAME = 'XX'
  IMN = 0
  write(u6,601) IAN
  return
else
  ANAME = AT(IAN)
end if
if ((IAN == 1) .and. (IMN /= 1)) then
  ! Special case: insert common name for deuterium or tritium
  if (IMN == 2) ANAME = ' D'
  if (IMN == 3) ANAME = ' T'
end if
GELGS = GEL(IAN)
MASS = -One
GNS = -1
ABUND = -One
do I=1,NMN(IAN)
  if (i > 10) write(u6,606) ian,imn,nmn(ian)
  606 format(3i9)
  if (IMN == MN(IAN,I)) then
    MASS = ZM(IAN,I)
    GNS = NS2(IAN,I)+1
    ABUND = AB(IAN,I)
  end if
end do
if (MASS < Zero) then
  MASS = ZM(IAN,0)
  if (IMN /= 0) write(u6,602) AT(IAN),IMN
  IMN = 0
end if

return

601 format(' *** MASSES Data base does not include Atomic Number=',i4)
602 format(' *** MASSES Data base does not include ',A2,'(',i3,'), so use average atomic mass.')

end subroutine MASSES
