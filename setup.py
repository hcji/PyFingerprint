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
  "jpype1",
  "gensim",
  "scikit-learn",
  "tqdm"
    ],
     )