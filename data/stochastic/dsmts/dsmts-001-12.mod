@model:2.3.1=BirthDeath12 "Birth-death model (001), variant 12"
@units
 substance=item
@compartments
 Cell
@species
 Cell:X=100 s
@parameters
 Lambda=0.1
 Mu=0.11
@reactions
@r=Birth
 X ->  2X
 Lambda*X*0.5*2
@r=Death
 X -> 
 Mu*X
