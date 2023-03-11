# backmapping
wello and helcome!
i calculated about how long this would take me, but boy am i bad at math =>
now go home

# stuff to get done
  - analysis package
    - need a series of functions to analyze the output of the backmapping workflow at each step
      - g(r)
      - bond length
      - dihedral angles
      - bond angles
  - write_backmapping protocol
    - need to pattern out how the current bashmap.sh and in. files on talapas work
    and organize around that

# hopes and dreams
  - taking the analysis package and creating a sort of monte carlo style simulation where we can map how 
  certain thermodynamic factors or other general parameters within lammps affect the structure
    - in less eloquent terms, i have a vague understanding of what a force field is, instead of learning how to 
    effectively make one i would rather brute force my way to victory!
  - one word, proteins
    - best way would be to just have a general program that can orient the specific amino acids in their given order
      - start with a large bead 
        - easy
      - move on using the splitting method we normally do
        - easy
      - onece the full chain is there, split the backbone from the r as beads
        - maybe easy
      - finally atomize the backbone and sidechain
        - easy
