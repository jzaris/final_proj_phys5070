# final_proj_phys5070

This is an N-body gravity simulation project. 
The repository consists of the Jupyter Notebook file 'Gravity main' and the python file 'gravity_functions'.  

'Gravity main' includes tests of the gravity code against known results and some interesting applications of the code to more complicated problems.
'gravity_functions' contains functions used to evolve gravitationally interacting systems, and is called many times in 'Gravity main'.

It would be easiest to evaluate the cells of 'Gravity main' in order.  I tried to interleave text and code so that the user
can understand the function of each cell.  

'gravity_functions' would only need to be analyzed to understand the time-stepping algorithm in which the momentum kick on each mass is calculated and the positions 
of the masses are evolved in time.
 
