from setuptools import setup, find_packages

setup(name='PyFingerprint',
      version='0.1.1',
      description='Python tool for generate fingerprints of a molecule',
      license='AGPLv3',
      author='Ji Hongchao',
      author_email='ji.hongchao@foxmail.com',
      url='https://github.com/hcji/PyFingerprint',
	  include_package_data = True,
      packages=find_packages(),
      install_requires=[
        "numpy          == 1.16.5",
        "h5py           == 2.9.0",
        "tensorflow-gpu == 2.0.0",
        "tqdm           == 4.35.0",
        "scikit-learn   == 0.21.3",
        "scipy          == 1.3.1",
        "ipykernel      == 5.1.1",
        "ipython",
        "jpype1         == 1.3.0",
        "matplotlib     == 3.1.1",
        "pandas         == 0.25.1",
	"molsets        == 0.2"
    ],
     )