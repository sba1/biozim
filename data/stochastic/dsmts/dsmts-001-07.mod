@model:2.3.1=BirthDeath07 "Birth-death model (001), variant 07"
@units
 substance=item
@compartments
 Cell
@species
 Cell:X=100 s
 Cell:Sink=0 s
@parameters
 Lambda=0.1
 Mu=0.11
@reactions
@r=Birth
 X ->  2X
 Lambda*X
@r=Death
 X -> Sink
 Mu*X
