@model:2.3.1=BirthDeath18 "Birth-death model (001), variant 18"
@units
 substance=item
@compartments
 Cell=0.5
@species
 Cell:X=100 s
@parameters
 Lambda=0.1
 Mu=0.11
@reactions
@r=Birth
 X ->  2X
 Cell*Lambda*X
@r=Death
 X -> 
 Cell*Mu*X
