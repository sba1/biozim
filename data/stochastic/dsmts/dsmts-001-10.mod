@model:2.3.1=BirthDeath10 "Birth-death model (001), variant 10"
@units
 substance=item
@compartments
 Cell=1
@species
 Cell:X=100
@parameters
 Lambda=0.1
 Mu=0.11
@reactions
@r=Birth
 X ->  2X
 Lambda*X
@r=Death
 X -> 
 Mu*X
