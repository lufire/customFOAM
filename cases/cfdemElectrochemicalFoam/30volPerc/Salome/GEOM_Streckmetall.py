# --- salome initialization
import salome
salome.salome_init()

# --- geom Python interface
import GEOM
from salome.geom import geomBuilder
geompy = geomBuilder.New(salome.myStudy)

# --- Create a box and publish it into the salome study
Vertex_1 = geompy.MakeVertex(-0.05, 0.10, 0.005)
Vertex_2 = geompy.MakeVertex(-0.04877641291, 0.105, 0.005)
Vertex_3 = geompy.MakeVertex(-0.04522542486, 0.11, 0.005)
Vertex_4 = geompy.MakeVertex(-0.03969463131, 0.115, 0.005)
Vertex_5 = geompy.MakeVertex(-0.03272542486, 0.12, 0.005)
Vertex_6 = geompy.MakeVertex(-0.025, 0.125, 0.005)
Vertex_7 = geompy.MakeVertex(-0.01727457514, 0.13, 0.005)
Vertex_8 = geompy.MakeVertex(-0.01030536869, 0.135, 0.005)
Vertex_9 = geompy.MakeVertex(-0.00477457514, 0.14, 0.005)
Vertex_10 = geompy.MakeVertex(-0.00122358709, 0.145, 0.005)
Vertex_11 = geompy.MakeVertex(0, 0.15, 0.005)
Vertex_12 = geompy.MakeVertex(-0.00122358709, 0.155, 0.005)
Vertex_13 = geompy.MakeVertex(-0.00477457514, 0.16, 0.005)
Vertex_14 = geompy.MakeVertex(-0.01030536869, 0.165, 0.005)
Vertex_15 = geompy.MakeVertex(-0.01727457514, 0.17, 0.005)
Vertex_16 = geompy.MakeVertex(-0.025, 0.175, 0.005)
Vertex_17 = geompy.MakeVertex(-0.03272542486,	0.18, 0.005)
Vertex_18 = geompy.MakeVertex(-0.03969463131,	0.185, 0.005)
Vertex_19 = geompy.MakeVertex(-0.04522542486,	0.19, 0.005)
Vertex_20 = geompy.MakeVertex(-0.04877641291,	0.195, 0.005)
Vertex_21 = geompy.MakeVertex(-0.05, 0.2, 0.005)

geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Vertex_10, 'Vertex_10' )
geompy.addToStudy( Vertex_11, 'Vertex_11' )
geompy.addToStudy( Vertex_12, 'Vertex_12' )
geompy.addToStudy( Vertex_13, 'Vertex_13' )
geompy.addToStudy( Vertex_14, 'Vertex_14' )
geompy.addToStudy( Vertex_15, 'Vertex_15' )
geompy.addToStudy( Vertex_16, 'Vertex_16' )
geompy.addToStudy( Vertex_17, 'Vertex_17' )
geompy.addToStudy( Vertex_18, 'Vertex_18' )
geompy.addToStudy( Vertex_19, 'Vertex_19' )
geompy.addToStudy( Vertex_20, 'Vertex_20' )
geompy.addToStudy( Vertex_21, 'Vertex_21' )

# --- update the study object browser
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
