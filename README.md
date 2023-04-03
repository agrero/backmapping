# backmapping
- Below are the files and a general overview of what they do
  - backpack (module)
    - mapanly.py - takes in the dataframes of coordinates/bonds/dihedrals etc and 
    outputs the bond lengths/angles of each chain
    - mapback.py - all functions related to physically backmapping those course grained models
    - writepack.py - all functions regarding reading and writing lammps formatted 
    coordinate files
  - backmap.py - run this to backmap a specified lammps configured coordinate file
  - default_configs.py - default lammps configs as a dictionary 
  - speed_test.py - used to test the speed of the backmapping algorithm with changing monomer no
  - test.py - file for testing individual components
  
# stuff to get done
  - analysis package
    - need a series of functions to analyze the output of the backmapping workflow at each step
      - g(r)
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
