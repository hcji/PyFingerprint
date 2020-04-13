## PyFingerprint
***

There are many types of chemical fingerprint for describing the molecule provided by different tools, such as RDKit, CDK and OpenBabel. This package aims to summarize them all.

### Dependencies

 1. Anaconda for python 3.6
 2. Java Runtime Environment 8.0
 3. jpype
 
        pip install jpype1

 4. RDKit

        conda install -c rdkit rdkit
 5. [Pybel](https://open-babel.readthedocs.io/en/latest/UseTheLibrary/PythonInstall.html) (Optional, only for OpenBabel fingerprints)
 
### Install

	pip install git+git://github.com/hcji/PyFingerprint@master

### Usage

	from PyFingerprint.All_Fingerprint import get_fingerprint
	fps = get_fingerprint('CCCCN', fp_type='daylight')
	
### Cite

	@article {Ji2020.03.30.017137,
	author = {Ji, Hongchao and Lu, Hongmei and Zhang, Zhimin},
	title = {Predicting Molecular Fingerprint from Electron-Ionization Mass Spectrum with Deep Neural Networks},
	elocation-id = {2020.03.30.017137},
	year = {2020},
	doi = {10.1101/2020.03.30.017137},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Electron-ionization mass spectrometry (EI-MS) hyphenated gas chromatography (GC) is the workhorse to analyze volatile compounds in complex samples. The spectral matching method can only identify compounds within spectral database. In response, we present a deep-learning-based approach (DeepEI) for structure elucidation of unknown compound with its EI-MS spectrum. DeepEI employs deep neural networks to predict molecular fingerprint from EI-MS spectrum, and searches molecular structure database with the predicted fingerprints. In addition, a convolutional neural network was also trained to filter the structures in database and improve the identification performance. Our method shows improvement on the competing method NEIMS in identification accuracy on both NIST test dataset and MassBank dataset. Furthermore, DeepEI (spectrum to fingerprint) and NEIMS (fingerprint to spectrum) can be combined to improve identification accuracy.},
	URL = {https://www.biorxiv.org/content/early/2020/04/01/2020.03.30.017137},
	eprint = {https://www.biorxiv.org/content/early/2020/04/01/2020.03.30.017137.full.pdf},
	journal = {bioRxiv}
	}

### Support fingerprint types:

	**daylight**: Considers paths of a given length. These are hashed fingerprints, with a default length of 1024.
	**extended**: Similar to the standard type, but takes rings and atomic properties into account into account.
	**graph**: Similar to the standard type by simply considers connectivity.
	**hybridization**: Similar to the standard type, but only consider hybridization state.
	**estate**: 79 bit fingerprints corresponding to the E-State atom types described by Hall and Kier.
	**pubchem**: 881 bit fingerprints defined by PubChem.
	**klekota-roth**: 4860 bit fingerprint defined by Klekota and Roth.
	**shortestpath**: A fingerprint based on the shortest paths between pairs of atoms and takes into account ring systems, charges etc.
	**signature**: A feature,count type of fingerprint, similar in nature to circular fingerprints, but based on the signature descriptor.
	**circular**: An implementation of the ECFP6 fingerprint.
	**lingo**: An implementation of the LINGO fingerprint.
	**rdkit**: Another implementation of a Daylight-like fingerprint by RDKit.
	**maccs**: The popular 166 bit MACCS keys described by MDL.
	**Avalon**: Substructure or similarity Avalon fingerprint.
	**AtomPair**: RDKit Atom-Pair fingerprint.
	**TopologicalTorsion**: RDKit Topological-Torsion Fingerprint.
	**Morgan**: RDKit Morgan fingerprint.
	**FP2**: OpenBabel FP2 fingerprint, which indexes small molecule fragments based on linear segments of up to 7 atoms in length.
	**FP3**: OpenBabel FP3 fingerprint, which is a fingerprint method created from a set of SMARTS patterns defining functional groups.
	**FP4**: OpenBabel FP4 fingerprint, which is a fingerprint method created from a set of SMARTS patterns defining functional groups.
	

