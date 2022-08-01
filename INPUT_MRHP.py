################################################################################
#|----------------------------------------------------------------------------|#
#|                         OpenMC Python API Input                            |#
#|            Micro Reactor Heat Pipe (MRHP) with Uranium Carbide             |#
#|----------------------------------------------------------------------------|#
################################################################################

import openmc
import matplotlib.pyplot as plt
import math

################################################################################
#                                  Materials                                   #
################################################################################

###################################################
Enrichment = 11.1 #enrichment of UC in %wt U     ##
###################################################

EC = 0.048049
EU = 1 - EC
E24 = 0.000254
E26 = 0.000131
E25 = ( Enrichment / 100 )*EU
E28 = 1-(E24+E25+E26+EC)

#TEMPERATURE VARIABLES ON VARIOUS COMPONENTS(in K, kelvin)
TMP0 = 350 # temperature of concrete, Nitrogen, and n aux (unchanged)
TMPA = 950 # temperature of fuel cladding and moderator (unchanged)
TMPB = 900 # temperature of heat pipe's wick and wall (unchanged)
TMPC = 700 # temperature of everything else (structure, control rod, shielding,
           # Argon, etc.) (unchanged)
TMPD = 850 # temperature of molten salt pool (unchanged)

#########################################################################
## VARIABLES CHANGES FOR TEMPERATURE REACTIVITY COEFFICIENT SIMULATION ##
#########################################################################
# Temperature (in K)      #                                            ##
TMPX = 950                # temperature of fuel                        ##
TMPY = 900                # temperature of heat pipe's vapor           ##
TMPZ = 700                # temperature of reflectors                  ##
#########################################################################
# Density (in g/cm3)      #                                            ##
DENX = 13.32535           # density of fuel                            ##
DENY = 1.72676E-05        # density of heat pipe's vapor               ##
DENZ = 1.82046716716717   # density of reflectors                      ##
#########################################################################

m100 = openmc.Material(name='Na, Vapor')
m100.set_density('g/cm3', DENY)
m100.add_element('Na', 1.0, 'wo')
m100.temperature = TMPY

m200 = openmc.Material(name='Na, Liquid')
m200.set_density('g/cm3', 0.802312)
m200.add_element('Na', 1.0, 'wo')
m200.temperature = TMPB

m300 = openmc.Material(name='SS 304 Wall')
m300.set_density('g/cm3', 7.47403437848045)
m300.add_element('C', 0.0003, 'wo')
m300.add_element('Si', 0.005, 'wo')
m300.add_element('P', 0.000225, 'wo')
m300.add_element('S', 0.00015, 'wo')
m300.add_element('Cr', 0.1700003, 'wo')
m300.add_element('Mn', 0.01, 'wo')
m300.add_element('Fe', 0.6693237, 'wo')
m300.add_element('Ni', 0.120001, 'wo')
m300.add_element('Mo', 0.025, 'wo')
m300.temperature = TMPB

m400 = openmc.Material(name='Uranium Carbide, UC')
m400.set_density('g/cm3', DENX)
m400.add_nuclide('U234', E24, percent_type='wo')
m400.add_nuclide('U235', E25, percent_type='wo')
m400.add_nuclide('U236', E26, percent_type='wo')
m400.add_nuclide('U238', E28, percent_type='wo')
m400.add_element('C', EC, percent_type='wo')
m400.temperature = TMPX

m500 = openmc.Material(name='SS 304 Cladding')
m500.set_density('g/cm3', 7.44743856292913)
m500.add_element('C', 0.0003, 'wo')
m500.add_element('Si', 0.005, 'wo')
m500.add_element('P', 0.000225, 'wo')
m500.add_element('S', 0.00015, 'wo')
m500.add_element('Cr', 0.1700003, 'wo')
m500.add_element('Mn', 0.01, 'wo')
m500.add_element('Fe', 0.6693237, 'wo')
m500.add_element('Ni', 0.120001, 'wo')
m500.add_element('Mo', 0.025, 'wo')
m500.temperature = TMPA

m600 = openmc.Material(name='Graphite, Reactor Grade, Moderator')
m600.set_density('g/cm3', 1.81221858141858)
m600.add_element('C', 0.999999, 'wo')
m600.add_element('B', 0.000001, 'wo')
m600.add_s_alpha_beta('c_Graphite')
m600.temperature = TMPA

m666 = openmc.Material(name='Graphite, Reactor Grade, Reflector')
m666.set_density('g/cm3', DENZ)
m666.add_element('C', 0.999999, 'wo')
m666.add_element('B', 0.000001, 'wo')
m666.add_s_alpha_beta('c_Graphite')
m666.temperature = TMPZ

m700 = openmc.Material(name='SS 304 Other')
m700.set_density('g/cm3', 7.56421837080509)
m700.add_element('C', 0.0003, 'wo')
m700.add_element('Si', 0.005, 'wo')
m700.add_element('P', 0.000225, 'wo')
m700.add_element('S', 0.00015, 'wo')
m700.add_element('Cr', 0.1700003, 'wo')
m700.add_element('Mn', 0.01, 'wo')
m700.add_element('Fe', 0.6693237, 'wo')
m700.add_element('Ni', 0.120001, 'wo')
m700.add_element('Mo', 0.025, 'wo')
m700.temperature = TMPC

m800 = openmc.Material(name='Argon')
m800.set_density('g/cm3', 0.0006862)
m800.add_element('Ar', 1.0, 'wo')
m800.temperature = TMPC

m900 = openmc.Material(name='B4C Absorber')
m900.set_density('g/cm3', 2.52)
m900.add_element('C', 0.2, 'ao')
m900.add_nuclide('B10', 0.16, 'ao')
m900.add_nuclide('B11', 0.64, 'ao')
m900.temperature = TMPC

m990 = openmc.Material(name='KNO3-NaNO3 50/50')
m990.set_density('g/cm3', 2.176)
m990.add_element('K', 0.1, 'ao')
m990.add_element('N', 0.2, 'ao')
m990.add_element('Na', 0.1, 'ao')
m990.add_element('O', 0.6, 'ao')
m990.temperature = TMPD

m995 = openmc.Material(name='Nitrogen')
m995.set_density('g/cm3', 0.0009855)
m995.add_element('N', 1.0, 'wo')
m995.temperature = TMP0

m999 = openmc.Material(name='Concrete, Ordinary')
m999.set_density('g/cm3', 2.35)
m999.add_element('H', 0.005558, 'wo')
m999.add_element('O', 0.498076, 'wo')
m999.add_element('Na', 0.017101, 'wo')
m999.add_element('Mg', 0.002565, 'wo')
m999.add_element('Al', 0.045746, 'wo')
m999.add_element('Si', 0.315092, 'wo')
m999.add_element('S', 0.001283, 'wo')
m999.add_element('K', 0.019239, 'wo')
m999.add_element('Ca', 0.082941, 'wo')
m999.add_element('Fe', 0.012398, 'wo')
m999.temperature = TMP0

m = openmc.Materials([m100,m200,m300,m400,m500,m600,m666,\
                      m700,m800,m900,m990,m995,m999])
m.cross_sections = '/home/fadilnaufal/OpenMC_CS/endfb80/cross_sections.xml'
m.export_to_xml()

################################################################################
#                                 Geometry                                     #
################################################################################

###########################################################
#                  Control Rod Mechanism                  #
###########################################################
               ##   POS1 is a position of Control Rod    ##
POS1 = -50.0   ##                                        ##
               ##   -50.0 --> bottom of CR are straight  ##
#################   with the bottom of CR upper channel  ##
Z2 = 50 - POS1 ##   ...                                  ##
T0 = 100.0     ##   100.0 --> bottom of CR are straight  ##
Z1 = Z2 + T0   ##   with the bottom of active core       ##
###########################################################

###########################################################
#              Movable Reflector Mechanism                #
###########################################################
               ##   POS2 is a position of Mov. Reflector ##
POS2 = 100.0   ##                                        ##
               ##   0.0 --> top of MR are straight       ##
#################   with the bottom of active core       ##
Z3 = -50+POS2  ##   ...                                  ##
T00 = 100.0    ##   100.0 --> top of MR are straight     ##
Z4 = Z3 - T00  ##   with the top of active core          ##
###########################################################

#################################   VARIABLES  #################################

#fuel element
r1 = 0.55           #radius of vapor
r2 = 0.65           #radius of wick
r3 = 0.75           #radius of wall
r4 = 1.28           #radius of fuel
r5 = 1.38           #radius of cladding
ptch = 2.8          #pitch of fuel element

#movable reflector + space outside movable reflector
rinn = 25.0         #radius of movable reflector
rout = 26.0         #radius of space outside movable reflector
pos1 = 83.0         #coordinate(s) of position of movable reflector
cos45 = 0.70710678118655
pos2 = cos45 * pos1 #coordinate(s) of position of movable reflector

#control rod + space outside control rod
pos3 = 25           #coordinate(s) of position of control rod
pos4 = 60           #coordinate(s) of position of control rod
rcr1 = 3.95         #radius of control rod
rcr2 = 4            #radius of space outside control rod

#upper plenum
rar1 = 20           #radius of space in upper plenum (|| to mov. reflector)
rar2 = 22           #radius of space in the upper plenum
rcr3 = 5            #radius of space in upper plenum (|| to control rod)

#outer part
rref = 111.0        #radius of outer reflector
tves = 5.0          #radius of outer vessel
rtot = rref + tves  #total radius of core
hi = 100.0          #height of active core
rrefb = 55.0        #outer radius of active core
hrefb = 100.0       #height of bottom reflector
hrefa = 50.0        #height of upper reflector
tkss = 5.0          #thickness of SS on the upper plenum and the very bottom
hshld = 20.0        #height of shield above the upper reflector
hic1 = 100.0        #height of heat pipe in the intermediate coolant
hic2 = 50.0         #height of the rest of the intermediate coolant
hic0 = hic1 + hic2

#neutron auxilary
t001 = 2.0          #thickness of base of neutron aux
h001 = 30.0         #height of base of neutron (Ar) aux
h002 = 50.0         #height of base of neutron (B4C) aux
r001 = 18.0         #radius of base of neutron (B4C) aux
r002 = 20.0         #radius of base of neutron (SS) aux
r003 = 22.0         #radius of base of neutron (IC) aux

#very outer part
tN = 2.0            #thickness of Nitrogen vessel
rN = rtot+tN        #radius of Nitrogen vessel

#the most outer part
rC = 200            #outer radius of concrete
z_minC = 250        #-z position of concrete

#################################   SURFACE   ##################################

#fuel element
s01 = openmc.ZCylinder(x0=0.0, y0=0.0, r=r1, boundary_type='transmission')
s02 = openmc.ZCylinder(x0=0.0, y0=0.0, r=r2, boundary_type='transmission')
s03 = openmc.ZCylinder(x0=0.0, y0=0.0, r=r3, boundary_type='transmission')
s04 = openmc.ZCylinder(x0=0.0, y0=0.0, r=r4, boundary_type='transmission')
s05 = openmc.ZCylinder(x0=0.0, y0=0.0, r=r5, boundary_type='transmission')

#movable reflector
s11 = openmc.ZCylinder(x0=pos1, y0=0.0, r=rinn, boundary_type='transmission')
s12 = openmc.ZCylinder(x0=0.0, y0=-pos1, r=rinn, boundary_type='transmission')
s13 = openmc.ZCylinder(x0=-pos1, y0=0.0, r=rinn, boundary_type='transmission')
s14 = openmc.ZCylinder(x0=0.0, y0=pos1, r=rinn, boundary_type='transmission')
s15 = openmc.ZCylinder(x0=pos2, y0=pos2, r=rinn, boundary_type='transmission')
s16 = openmc.ZCylinder(x0=pos2, y0=-pos2, r=rinn, boundary_type='transmission')
s17 = openmc.ZCylinder(x0=-pos2, y0=-pos2, r=rinn, boundary_type='transmission')
s18 = openmc.ZCylinder(x0=-pos2, y0=pos2, r=rinn, boundary_type='transmission')

s21 = openmc.ZCylinder(x0=pos1, y0=0.0, r=rout, boundary_type='transmission')
s22 = openmc.ZCylinder(x0=0.0, y0=-pos1, r=rout, boundary_type='transmission')
s23 = openmc.ZCylinder(x0=-pos1, y0=0.0, r=rout, boundary_type='transmission')
s24 = openmc.ZCylinder(x0=0.0, y0=pos1, r=rout, boundary_type='transmission')
s25 = openmc.ZCylinder(x0=pos2, y0=pos2, r=rout, boundary_type='transmission')
s26 = openmc.ZCylinder(x0=pos2, y0=-pos2, r=rout, boundary_type='transmission')
s27 = openmc.ZCylinder(x0=-pos2, y0=-pos2, r=rout, boundary_type='transmission')
s28 = openmc.ZCylinder(x0=-pos2, y0=pos2, r=rout, boundary_type='transmission')

#upper plenum (parallel to movable reflector)
s31 = openmc.ZCylinder(x0=pos1, y0=0.0, r=rar1, boundary_type='transmission')
s32 = openmc.ZCylinder(x0=0.0, y0=-pos1, r=rar1, boundary_type='transmission')
s33 = openmc.ZCylinder(x0=-pos1, y0=0.0, r=rar1, boundary_type='transmission')
s34 = openmc.ZCylinder(x0=0.0, y0=pos1, r=rar1, boundary_type='transmission')
s35 = openmc.ZCylinder(x0=pos2, y0=pos2, r=rar1, boundary_type='transmission')
s36 = openmc.ZCylinder(x0=pos2, y0=-pos2, r=rar1, boundary_type='transmission')
s37 = openmc.ZCylinder(x0=-pos2, y0=-pos2, r=rar1, boundary_type='transmission')
s38 = openmc.ZCylinder(x0=-pos2, y0=pos2, r=rar1, boundary_type='transmission')

s41 = openmc.ZCylinder(x0=pos1, y0=0.0, r=rar2, boundary_type='transmission')
s42 = openmc.ZCylinder(x0=0.0, y0=-pos1, r=rar2, boundary_type='transmission')
s43 = openmc.ZCylinder(x0=-pos1, y0=0.0, r=rar2, boundary_type='transmission')
s44 = openmc.ZCylinder(x0=0.0, y0=pos1, r=rar2, boundary_type='transmission')
s45 = openmc.ZCylinder(x0=pos2, y0=pos2, r=rar2, boundary_type='transmission')
s46 = openmc.ZCylinder(x0=pos2, y0=-pos2, r=rar2, boundary_type='transmission')
s47 = openmc.ZCylinder(x0=-pos2, y0=-pos2, r=rar2, boundary_type='transmission')
s48 = openmc.ZCylinder(x0=-pos2, y0=pos2, r=rar2, boundary_type='transmission')

#control rod
s51 = openmc.ZCylinder(x0=pos3, y0=pos4, r=rcr1, boundary_type='transmission')
s52 = openmc.ZCylinder(x0=pos4, y0=pos3, r=rcr1, boundary_type='transmission')
s53 = openmc.ZCylinder(x0=pos4, y0=-pos3, r=rcr1, boundary_type='transmission')
s54 = openmc.ZCylinder(x0=pos3, y0=-pos4, r=rcr1, boundary_type='transmission')
s55 = openmc.ZCylinder(x0=-pos3, y0=-pos4, r=rcr1, boundary_type='transmission')
s56 = openmc.ZCylinder(x0=-pos4, y0=-pos3, r=rcr1, boundary_type='transmission')
s57 = openmc.ZCylinder(x0=-pos4, y0=pos3, r=rcr1, boundary_type='transmission')
s58 = openmc.ZCylinder(x0=-pos3, y0=pos4, r=rcr1, boundary_type='transmission')

s61 = openmc.ZCylinder(x0=pos3, y0=pos4, r=rcr2, boundary_type='transmission')
s62 = openmc.ZCylinder(x0=pos4, y0=pos3, r=rcr2, boundary_type='transmission')
s63 = openmc.ZCylinder(x0=pos4, y0=-pos3, r=rcr2, boundary_type='transmission')
s64 = openmc.ZCylinder(x0=pos3, y0=-pos4, r=rcr2, boundary_type='transmission')
s65 = openmc.ZCylinder(x0=-pos3, y0=-pos4, r=rcr2, boundary_type='transmission')
s66 = openmc.ZCylinder(x0=-pos4, y0=-pos3, r=rcr2, boundary_type='transmission')
s67 = openmc.ZCylinder(x0=-pos4, y0=pos3, r=rcr2, boundary_type='transmission')
s68 = openmc.ZCylinder(x0=-pos3, y0=pos4, r=rcr2, boundary_type='transmission')

#upper plenum (parallel to control rod)
s71 = openmc.ZCylinder(x0=pos3, y0=pos4, r=rcr3, boundary_type='transmission')
s72 = openmc.ZCylinder(x0=pos4, y0=pos3, r=rcr3, boundary_type='transmission')
s73 = openmc.ZCylinder(x0=pos4, y0=-pos3, r=rcr3, boundary_type='transmission')
s74 = openmc.ZCylinder(x0=pos3, y0=-pos4, r=rcr3, boundary_type='transmission')
s75 = openmc.ZCylinder(x0=-pos3, y0=-pos4, r=rcr3, boundary_type='transmission')
s76 = openmc.ZCylinder(x0=-pos4, y0=-pos3, r=rcr3, boundary_type='transmission')
s77 = openmc.ZCylinder(x0=-pos4, y0=pos3, r=rcr3, boundary_type='transmission')
s78 = openmc.ZCylinder(x0=-pos3, y0=pos4, r=rcr3, boundary_type='transmission')

#outer part
s91 = openmc.ZCylinder(x0=0.0, y0=0.0, r=rref, boundary_type='transmission')
s92 = openmc.ZCylinder(x0=0.0, y0=0.0, r=rtot, boundary_type='transmission')
s93 = openmc.ZPlane(z0=hi/2, boundary_type='transmission')
s94 = openmc.ZPlane(z0=-hi/2, boundary_type='transmission')

CR93 = openmc.ZPlane(z0=Z1, boundary_type='transmission')
CR94 = openmc.ZPlane(z0=Z2, boundary_type='transmission')

MR93 = openmc.ZPlane(z0=Z3, boundary_type='transmission')
MR94 = openmc.ZPlane(z0=Z4, boundary_type='transmission')

s30 = openmc.ZCylinder(x0=0.0, y0=0.0, r=rrefb, boundary_type='transmission')
s40 = openmc.ZCylinder(x0=0.0, y0=0.0, r=rrefb+2.0, boundary_type='transmission')
s95 = openmc.ZPlane(z0=((-hi/2)-hrefb), boundary_type='transmission')
s96 = openmc.ZPlane(z0=((hi/2)+hrefa), boundary_type='transmission')

s97 = openmc.ZPlane(z0=((hi/2)+hrefa+hshld), boundary_type='transmission')
s98 = openmc.ZPlane(z0=((hi/2)+hrefa+tkss), boundary_type='transmission')

s81 = openmc.ZCylinder(x0=0.0, y0=0.0, r=65.0, boundary_type='transmission')
s82 = openmc.ZCylinder(x0=0.0, y0=0.0, r=67.0, boundary_type='transmission')

s90 = openmc.ZPlane(z0=((hi/2)+hrefa+hshld+hic1), boundary_type='transmission')

s990 = openmc.ZPlane(z0=((-hi/2)-hrefb-tkss), boundary_type='transmission')

#neutron auxilary
s001 = openmc.ZCylinder(x0=pos3, y0=pos4, r=r001, boundary_type='transmission')
s002 = openmc.ZCylinder(x0=-pos3, y0=-pos4, r=r001, boundary_type='transmission')
s003 = openmc.ZCylinder(x0=pos3, y0=pos4, r=r002, boundary_type='transmission')
s004 = openmc.ZCylinder(x0=-pos3, y0=-pos4, r=r002, boundary_type='transmission')
s005 = openmc.ZCylinder(x0=pos3, y0=pos4, r=r003, boundary_type='transmission')
s006 = openmc.ZCylinder(x0=-pos3, y0=-pos4, r=r003, boundary_type='transmission')
#
s101 = openmc.ZPlane(z0=(-hi/2)-hrefb-tkss-t001, boundary_type='transmission')
s102 = openmc.ZPlane(z0=(-hi/2)-hrefb-h001, boundary_type='transmission')
s103 = openmc.ZPlane(z0=(-hi/2)-hrefb-h002, boundary_type='transmission')
s104 = openmc.ZPlane(z0=(-hi/2)-hrefb-h002-t001, boundary_type='transmission')
s105 = openmc.ZPlane(z0=(-hi/2)-hrefb-h002-(2*t001), boundary_type='transmission')

#Nitrogen vessel
s201 = openmc.ZCylinder(x0=0.0, y0=0.0, r=rN, boundary_type='transmission')

#the most outer part
s99 = openmc.ZPlane(z0=((hi/2)+hrefa+hshld+hic0), boundary_type='vacuum')
s999 = openmc.ZPlane(z0=(-z_minC), boundary_type='vacuum')
s9999 = openmc.ZCylinder(x0=0.0, y0=0.0, r=rC, boundary_type='vacuum')

###################################   CELL   ###################################

c01 = openmc.Cell(fill=m100, region=-s01 & -s96 & +s94)
c02 = openmc.Cell(fill=m200, region=+s01 & -s02 & -s96 & +s94)
c03 = openmc.Cell(fill=m300, region=+s02 & -s03 & -s96 & +s94)
c04 = openmc.Cell(fill=m400, region=+s03 & -s04 & -s93 & +s94)
c05 = openmc.Cell(fill=m500, region=+s04 & -s05 & -s93 & +s94)
c06 = openmc.Cell(fill=m600, region=+s05 & -s93 & +s94)
c07 = openmc.Cell(fill=m666, region=+s03 & -s96 & +s93) #Upper Reflector
pin1 = openmc.Universe(cells=(c01, c02, c03, c04, c05, c06, c07))

c08 = openmc.Cell(fill=m100, region=-s01 & -s97 & +s96)
c09 = openmc.Cell(fill=m200, region=+s01 & -s02 & -s97 & +s96)
c10 = openmc.Cell(fill=m300, region=+s02 & -s03 & -s97 & +s96)
c11 = openmc.Cell(fill=m700, region=+s03 & -s97 & +s96)
pin2 = openmc.Universe(cells=(c08, c09, c10, c11))

c12 = openmc.Cell(fill=m100, region=-s01 & -s90 & +s97)
c13 = openmc.Cell(fill=m200, region=+s01 & -s02 & -s90 & +s97)
c14 = openmc.Cell(fill=m300, region=+s02 & -s03 & -s90 & +s97)
c15 = openmc.Cell(fill=m990, region=+s03 & -s90 & +s97)
pin3 = openmc.Universe(cells=(c12, c13, c14, c15))

#lower plenum lattice
cellC = openmc.Cell(fill=m666) #Side Reflector
univC = openmc.Universe(cells=(cellC,))

ring00 = ([univC]*6 + [pin1]*10 + [univC]*5)*6
ring01 = ([univC]*4 + [pin1]*13 + [univC]*3)*6
ring02 = ([univC] + [pin1]*18)*6
ring03 = [pin1]*18*6
ring04 = [pin1]*17*6
ring05 = [pin1]*16*6
ring06 = [pin1]*15*6
ring07 = [pin1]*14*6
ring08 = [pin1]*13*6
ring09 = [pin1]*12*6
ring10 = [pin1]*11*6
ring11 = [pin1]*10*6
ring12 = [pin1]*9*6
ring13 = [pin1]*8*6
ring14 = [pin1]*7*6
ring15 = [pin1]*6*6
ring16 = [pin1]*5*6
ring17 = [pin1]*4*6
ring18 = [pin1]*3*6
ring19 = [pin1]*2*6
ring20 = [pin1]*1*6
ring21 = [pin1]

lat1 = openmc.HexLattice()
lat1.orientation = 'x'
lat1.center = (0., 0.)
lat1.pitch = (ptch,)
lat1.outer = univC
lat1.universes = [ ring00,ring01,ring02,ring03,ring04,ring05,ring06,ring07,
                   ring08,ring09,ring10,ring11,ring12,ring13,ring14,ring15,
                   ring16,ring17,ring18,ring19,ring20,ring21 ]

#middle plenum lattice
cellSS = openmc.Cell(fill=m700)
univSS = openmc.Universe(cells=(cellSS,))

ring30 = ([univSS]*6 + [pin2]*10 + [univSS]*5)*6
ring31 = ([univSS]*4 + [pin2]*13 + [univSS]*3)*6
ring32 = ([univSS] + [pin2]*18)*6
ring33 = [pin2]*18*6
ring34 = [pin2]*17*6
ring35 = [pin2]*16*6
ring36 = [pin2]*15*6
ring37 = [pin2]*14*6
ring38 = [pin2]*13*6
ring39 = [pin2]*12*6
ring40 = [pin2]*11*6
ring41 = [pin2]*10*6
ring42 = [pin2]*9*6
ring43 = [pin2]*8*6
ring44 = [pin2]*7*6
ring45 = [pin2]*6*6
ring46 = [pin2]*5*6
ring47 = [pin2]*4*6
ring48 = [pin2]*3*6
ring49 = [pin2]*2*6
ring50 = [pin2]*1*6
ring51 = [pin2]

lat2 = openmc.HexLattice()
lat2.orientation = 'x'
lat2.center = (0., 0.)
lat2.pitch = (ptch,)
lat2.outer = univSS
lat2.universes = [ ring30,ring31,ring32,ring33,ring34,ring35,ring36,ring37,
                   ring38,ring39,ring40,ring41,ring42,ring43,ring44,ring45,
                   ring46,ring47,ring48,ring49,ring50,ring51 ]

#upper plenum lattice
cellIC = openmc.Cell(fill=m990)
univIC = openmc.Universe(cells=(cellIC,))

ring60 = ([univIC]*6 + [pin3]*10 + [univIC]*5)*6
ring61 = ([univIC]*4 + [pin3]*13 + [univIC]*3)*6
ring62 = ([univIC] + [pin3]*18)*6
ring63 = [pin3]*18*6
ring64 = [pin3]*17*6
ring65 = [pin3]*16*6
ring66 = [pin3]*15*6
ring67 = [pin3]*14*6
ring68 = [pin3]*13*6
ring69 = [pin3]*12*6
ring70 = [pin3]*11*6
ring71 = [pin3]*10*6
ring72 = [pin3]*9*6
ring73 = [pin3]*8*6
ring74 = [pin3]*7*6
ring75 = [pin3]*6*6
ring76 = [pin3]*5*6
ring77 = [pin3]*4*6
ring78 = [pin3]*3*6
ring79 = [pin3]*2*6
ring80 = [pin3]*1*6
ring81 = [pin3]

lat3 = openmc.HexLattice()
lat3.orientation = 'x'
lat3.center = (0., 0.)
lat3.pitch = (ptch,)
lat3.outer = univIC
lat3.universes = [ ring60,ring61,ring62,ring63,ring64,ring65,ring66,ring67,
                   ring68,ring69,ring70,ring71,ring72,ring73,ring74,ring75,
                   ring76,ring77,ring78,ring79,ring80,ring81 ]

#initialize root universe
U = openmc.Universe()

#movable reflector
inn = (-s11|-s12|-s13|-s14|-s15|-s16|-s17|-s18) & -MR93 & +MR94
c21 = openmc.Cell(fill=m666, region=inn)
out = (-s21|-s22|-s23|-s24|-s25|-s26|-s27|-s28) & -s96 & +s95
c22 = openmc.Cell(fill=m800, region=~inn & out)

#control rod
inn1 = (-s51|-s52|-s53|-s54|-s55|-s56|-s57|-s58) & -CR93 & +CR94
c23 = openmc.Cell(fill=m900, region=inn1)
out1 = (-s61|-s62|-s63|-s64|-s65|-s66|-s67|-s68) & -s96 & +s95
c24 = openmc.Cell(fill=m800, region=~inn1 & out1)

#bottom reflector
s_reflb = -s30 & -s94 & +s95
reflb = openmc.Cell(fill=m666, region=s_reflb)

#ss that separates upper and lower plenum
upmr = -s31|-s32|-s33|-s34|-s35|-s36|-s37|-s38
upcr = -s61|-s62|-s63|-s64|-s65|-s66|-s67|-s68
uppl = -s31|-s63|-s36|-s64|-s32|-s65|-s37|-s66| \
       -s33|-s67|-s38|-s68|-s34|-s61|-s35|-s62
c25 = openmc.Cell(fill=m700, region= +s30 & -s91 & +s96 & -s98 & (~uppl))

#ss around shielding
c26 = openmc.Cell(fill=m700, region= +s30 & -s40 & +s98 & -s97)

#intermediate coolant on MIDDLE PLENUM
out2 = -s41|-s42|-s43|-s44|-s45|-s46|-s47|-s48
out3 = -s71|-s72|-s73|-s74|-s75|-s76|-s77|-s78
int1 = -s91 & +s82 & +s98 & -s97 & (~out2) & (~out3)
int2 = +s40 & -s81 & +s98 & -s97 & (~out2) & (~out3)
c27 = openmc.Cell(fill=m990, region= (int1) | (int2) )

#ss that connects cr and mr on MIDDLE PLENUM
x= -s41|-s73|-s46|-s74|-s42|-s75|-s47|-s76| \
   -s43|-s77|-s48|-s78|-s44|-s71|-s45|-s72
#backslash for row over 80 columns
c28 = openmc.Cell(fill=m700, region=(+s81&-s82&+s98&-s97) & (~x))

#ss on MIDDLE PLENUM
ssmr = (~upmr)&(out2)
sscr = (~upcr)&(out3)
c29=openmc.Cell( fill=m700,region=+s98&-s97& ((ssmr)|(sscr)) )

#argon on MIDDLE PLENUM
c30 = openmc.Cell(fill=m800, region= +s96 & -s97 & (upmr | upcr))

#intermediate coolant on UPPER PLENUM
int3 = -s91 & +s82 & +s97 & -s99 & (~out2) & (~out3)
int4 = +s30 & -s81 & +s97 & -s99 & (~out2) & (~out3)
c027 = openmc.Cell(fill=m990, region= (int3) | (int4) )

#ss that connects cr and mr on UPPER PLENUM
c028 = openmc.Cell(fill=m700, region=(+s81&-s82&+s97&-s99) & (~x))

#ss on UPPER PLENUM
c029=openmc.Cell( fill=m700,region=+s97&-s99& ((ssmr)|(sscr)) )

#argon on UPPER PLENUM
c030 = openmc.Cell(fill=m800, region= +s97 & -s99 & (upmr | upcr))

#main cell
main1 = openmc.Cell(fill=lat1)
main1.region = -s91 & ~out & ~out1 & ~s_reflb & -s96 & +s95

up = -s97 & +s96 & -s91 & +s30
main2 = openmc.Cell(fill=lat2)
main2.region = -s91 & ~up & -s97 & +s96

#upper intermediate coolant
ic = -s30 & -s99 & +s90
c40 = openmc.Cell(fill=m990, region=ic)

up2 = -s99 & +s97 & -s91 & +s30
main3 = openmc.Cell(fill=lat3)
main3.region = -s91 & ~ic & ~up2 & -s99 & +s97

#neutron auxilary
rg01 = -s61 & -s95 & +s102
cn01 = openmc.Cell(fill=m800, region=rg01)
rg02 = -s001 & -s95 & +s103
cn02 = openmc.Cell(fill=m900, region=~rg01 & rg02)
rg03 = -s003 & -s95 & +s104
cn03 = openmc.Cell(fill=m700, region=~rg02 & rg03)
rg04 = -s003 & -s101 & +s104
rg05 = -s005 & -s101 & +s105
cn04 = openmc.Cell(fill=m990, region=~rg04 & rg05)
#
rg11 = -s65 & -s95 & +s102
cn11 = openmc.Cell(fill=m800, region=rg11)
rg12 = -s002 & -s95 & +s103
cn12 = openmc.Cell(fill=m900, region=~rg11 & rg12)
rg13 = -s004 & -s95 & +s104
cn13 = openmc.Cell(fill=m700, region=~rg12 & rg13)
rg14 = -s004 & -s101 & +s104
rg15 = -s006 & -s101 & +s105
cn14 = openmc.Cell(fill=m990, region=~rg14 & rg15)

#SS vessel
outpt1 = -s91 & -s99 & +s95
outpt2 = -s92 & -s99 & +s990
o1 = -s003 & -s95 & +s990
o2 = -s004 & -s95 & +s990
o0 = o1 | o2
vessel = openmc.Cell(fill=m700)
vessel.region = (outpt2) & ~outpt1 & ~o0

#N vessel
oN1 = -s92 & -s99 & +s990
oN2 = -s201 & -s99 & +s101
oN3 = -s003 & -s990 & +s101
oN4 = -s004 & -s990 & +s101
oN0 = oN3 | oN3
N_vessel = openmc.Cell(fill=m995)
N_vessel.region = ~oN1 & oN2 & ~oN0

#concrete
oC1 = -s201 & -s99 & +s101
oC2 = -s9999 & -s99 & +s999
oC3 = -s005 & -s101 & +s105
oC4 = -s006 & -s101 & +s105
oC0 = oC3 | oC4
C_vessel = openmc.Cell(fill=m999)
C_vessel.region = ~oC1 & oC2 & ~oC0

U.add_cells((main1,main2,main3,\
             c21,c22,c23,c24,reflb,\
             c25,c26,c27,c28,c29,c30,\
             c40,c027,c028,c029,c030,\
             cn01,cn02,cn03,cn04,\
             cn11,cn12,cn13,cn14,\
             vessel,N_vessel,C_vessel))

################################################################################
#                                   Plots                                      #
################################################################################

pixlxy = (rC+5)*2
pixlxz = 550

colors = {}
colors[m100] = 'green'
colors[m200] = 'orange'
colors[m300] = 'blue'
colors[m400] = 'yellow'
colors[m500] = 'blue'
colors[m600] = 'gray'
colors[m666] = 'gray'
colors[m700] = 'blue'
colors[m800] = 'pink'
colors[m900] = 'fuchsia'
colors[m990] = 'lightblue'
colors[m995] = 'palegreen'
colors[m999] = 'dimgray'

U.plot(
    origin=(0.0, 0.0, 0.0),
    width=(pixlxz, pixlxz),
    pixels=(500, 500),
    basis='xz',
    color_by='material',
    colors=colors)
plt.show()

geom = openmc.Geometry(U)
geom.export_to_xml()

################################################################################
#                                 Settings                                     #
################################################################################

batches = 110
inactive = 10
particles = 10000
settings_file = openmc.Settings()
settings_file.photon_transport = True #initializing neutron-photon coupled
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.temperature = {'method': 'interpolation'}

# Initial Uniform Spatial Source Distribution Over Fissionable Zones
uniform_dist = openmc.stats.Point(xyz=(0.0, 1.0, 0.0))
settings_file.source = openmc.Source(space=uniform_dist)
settings_file.export_to_xml()

################################################################################
#                                  Tally                                       #
################################################################################

tallies = openmc.Tallies()

filter_cell0 = openmc.CellFilter(main1)
tally0 = openmc.Tally()
tally0.filters = [filter_cell0]
tally0.scores = ['heating']
tallies.append(tally0)

filter_cell1 = openmc.DistribcellFilter([c04])
tally1 = openmc.Tally()
tally1.filters = [filter_cell1]
tally1.scores = ['heating']
tallies.append(tally1)

mesh2 = openmc.RegularMesh(mesh_id=1)
mesh2.dimension = [1, 1, 10]
mesh2.lower_left = [-1.28, -1.28, -50]
mesh2.upper_right = [1.28, 1.28, 50]
MESH02 = openmc.MeshFilter(mesh2)

tally2 = openmc.Tally()
tally2.filters = [MESH02]
tally2.scores = ['heating']
tallies.append(tally2)

tallies.export_to_xml()
