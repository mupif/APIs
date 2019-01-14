_____tm.out
comment
TransientTransport nsteps 20 deltat 3600. rtolv 1.e-4 alpha 0.6 exportfields 1 5 nmodules 1
vtkxml tstep_all domain_all primvars 1 6 vars 2 37 56 stype 1
domain HeatTransfer
OutputManager
ndofman 24 nelem 5 ncrosssect 1 nmat 1 nbc 4 nic 0 nltf 1 nset 3
node     1  coords 3     0.4000     0.4000     0.4000
node     2  coords 3     0.4400     0.4000     0.4000
node     3  coords 3     0.4800     0.4000     0.4000
node     4  coords 3     0.5200     0.4000     0.4000
node     5  coords 3     0.5600     0.4000     0.4000
node     6  coords 3     0.6000     0.4000     0.4000
node     7  coords 3     0.4000     0.6000     0.4000
node     8  coords 3     0.4400     0.6000     0.4000
node     9  coords 3     0.4800     0.6000     0.4000
node    10  coords 3     0.5200     0.6000     0.4000
node    11  coords 3     0.5600     0.6000     0.4000
node    12  coords 3     0.6000     0.6000     0.4000
node    13  coords 3     0.4000     0.4000     0.6000
node    14  coords 3     0.4400     0.4000     0.6000
node    15  coords 3     0.4800     0.4000     0.6000
node    16  coords 3     0.5200     0.4000     0.6000
node    17  coords 3     0.5600     0.4000     0.6000
node    18  coords 3     0.6000     0.4000     0.6000
node    19  coords 3     0.4000     0.6000     0.6000
node    20  coords 3     0.4400     0.6000     0.6000
node    21  coords 3     0.4800     0.6000     0.6000
node    22  coords 3     0.5200     0.6000     0.6000
node    23  coords 3     0.5600     0.6000     0.6000
node    24  coords 3     0.6000     0.6000     0.6000

brick1ht     1  nodes 8     1    13    14     2     7    19    20     8
brick1ht     2  nodes 8     2    14    15     3     8    20    21     9
brick1ht     3  nodes 8     3    15    16     4     9    21    22    10
brick1ht     4  nodes 8     4    16    17     5    10    22    23    11
brick1ht     5  nodes 8     5    17    18     6    11    23    24    12

SimpleTransportCS 1 mat 1 set 1
IsoHeat 1 d 2400. k 1.0 c 1000.0

constantsurfaceload 1 loadTimeFunction 1 components 1 -10.0 properties 1 a 20.0 loadtype 3 set 2
constantsurfaceload 2 loadTimeFunction 1 components 1 -10.0 properties 1 a 20.0 loadtype 3 set 3

constantsurfaceload 3 loadTimeFunction 1 components 1 -10.0 properties 1 e 0.8 loadtype 7 set 2
constantsurfaceload 4 loadTimeFunction 1 components 1 -10.0 properties 1 e 0.8 loadtype 7 set 3

ConstantFunction 1 f(t) 1.0
Set 1 elementranges {(1 5)}
Set 2 elementboundaries 2 1 3
Set 3 elementboundaries 2 5 5
