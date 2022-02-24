.. _reference: https://doi.org/10.1016/S1090-7807(02)00014-9
.. _XPLOR-NIH: https://nmr.cit.nih.gov/xplor-nih

This repository includes the Python code and test dataset for the
PCS-CPMG methodology as described in:

John B. Stiller, Renee Otten, Daniel Häussinger, Pascal S. Rieder, Douglas L. Theobald & Dorothee Kern.
*Structure determination of high-energy states in a dynamic protein ensemble* **(2022)**,
Nature, doi: `10.1038/s41586-022-04468-9 <https://www.nature.com/articles/s41586-022-04468-9>`_

The main file is ``refine_EM_PCS.py``, which is a Python-based script
written for the XPLOR-NIH software [`reference`_]. It is a modified version
of the ``refine.py`` script that is part of the `XPLOR-NIH`_ package.

All the work was done using XPLOR-NIH version 2.46 (i.e., Python 2.7), but
it should work with all versions in the 2.x branch. Please note, that we
are currently working on making sure that it runs with XPLOR-NIH v3.x and
to make it more user-friendly.

The script takes degenerate pseudocontact shift (PCS) data extracted
extracted from Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion
experiments and uses it to solve the high-energy structure. The
methodology relies on an Expectation Maximization algorithm for handling
of the degenerate data.


Test dataset and output of a run for adenylate kinase
-----------------------------------------------------
The following files are used in the example run for adenylate kinase:

- 4AKE_Modeled.pdb: chain A of PDB entry `4AKE <https://www.rcsb.org/structure/4AKE>`_,
  modeled with the *Geobacillus stearothermophilus* sequence. Open Adk structure
  that is a simulated high-energy state in this case.

- 4AKE_Modeled_MinorPCS_True_Co.txt: simulated PCS generated for the 4AKE structure.

- 4AKE_Modeled_FourPCS_Co.txt: four-fold degenerate simulated PCSs generated for
  the ``4AKE`` structure.


Flow-chart for a typical, iterative PCS-CPMG run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. calculate 50-100 structures with XPLOR-NIH (``Run1``) using the protocol in ``refine_EM_PCS.py``
2. extract information from each run:
   - the take each structure and consider its likelihood (using the script ``Make_MC_Files_Indiv.py``)
   - calculate features for each Adk structure (using the script ``Extract_Information_Indiv.py``)
3. start further iterations (``Run2`` and ``Run3``), but use the output from the previous step as the
   starting point. The parameters are very similar, if not the same, but we *automatically*
   use the refined PDB structure and PCS file with optimized probabilities from the previous run.


Peaklists for lanthanide-tagged ubiquitin
-----------------------------------------
Finally, peak assignments for the ubiquitin mutants with diamagnetic/paramagnetic
lanthanide-tags are also provided and can be found in the ``ubiquitin`` directory:

- gNhsqc_ubiquitin-K6C-Lu-M7-PyThz_25C.list
- gNhsqc_ubiquitin-K6C-Tm-M7-PyThz_25C.list
- gNhsqc_ubiquitin-K6C-Lu-M8-4R4S_25C.list
- gNhsqc_ubiquitin-K6C-Tm-M8-4R4S_25C.list
- gNhsqc_ubiquitin-S20C-Lu-M7-PyThz_25C.list
- gNhsqc_ubiquitin-S20C-Tm-M7-PyThz_25C.list

where ``M7`` and ``M8`` refer to the DOTA-M7PyThiazole and DOTA-M8-(4R4S)-SSPy
lanthanide-binding tags, respectively; ``Lu`` is the diamagnetic lutetium and
``Tm`` is the paramagnetic thulium. Side-chains NH2 groups are randomly
assigned as ``SC#?-?``, where ``#`` is a number and ``?`` is ``a`` and ``b``
for the downfield and upfield cross peak, respectively. Cross-peaks originating
from untagged protein are suffixed by a ``u`` (e.g., ``Q2Nu-Hu``).
