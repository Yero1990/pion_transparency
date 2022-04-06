# CLAS12ROOT on ifarm
To run clas12root on ifarm, do:
>    `source /group/clas12/packages/setup.csh
(or setup.sh for bash)` <br>
`module load clas12/pro`

Then, from a ROOT C++ script, the clas12root, and related tools can be called, and its methods can be invoked as well. For example, one can create class instances: <br>
 
>`clas12root::HipoChain chain` <br>
`chain.Add("path/to/hipo_file")` to combine multiple hipo files of same run <br>
`chain.next()` for event loop <br>

To create a ROOTfile with leafs and branches, one can do:

>`hiporoot::ParticleTree treemaker("file.hipo", "test.root");` <br>
>`treemaker.Branch("my_branch/F");` <br>
>`treemaker.Fill();` 


These lines can be called within your own ROOT macro, for example.