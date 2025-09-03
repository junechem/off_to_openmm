Problem:
    Currently, there exists no framework for copying force field parameters from CRYOFF (Create Your Own Force Field) output (.off files)
        to the proper sections of an OpenMM .xml force field file. Thus, if someone desired to run an AFM (Adaptive Force Matching) model on OpenMM,
        they would have to go through the painstaking process of copying all the parameters to the .xml file by hand, which is
        very challenging.
    There are other difficulties associated with running simulations on OpenMM. One of those is getting the proper .pdb format.
        Usually, .pdb files that will be run on OpenMM will not include the important section related to bonding information
        at the bottom of the .pdb file, which is necessary to run simulations in OpenMM with custom forces like we will be doing.
        Thus, it will be necessary to produce a standalone tool that can add this bonding information in to the bottom of the
        .pdb file.


What are we building:
    1) We are building a set of python tools which can be used to produce .xml force field files for OpenMM simulations
    2) We are also creating a set of tools which make running simulations using OpenMM easier
Who is it for:
    This set of scripts is for the person who is interested in running Adaptive Force Matching based force fields
        in OpenMM, either for enhanced speed capabilities due to using a GPU or for PIMD simulations which can be used
        to compute nuclear quantum effects. It is also to increase the possible audience of people who would like to use
        the AFM method to produce force fields.
Why does it matter now:
    It matters now because I am trying to run OpenMM simulations using AFM based models, but cannot because no current scripts exist to make this happen.
