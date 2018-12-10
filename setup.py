from setuptools import setup, find_packages

setup(name='PyFingerprint',
      version='0.1.0',
      description='Python tool for generate fingerprints of a molecule',
      license='AGPLv3',
      author='Ji Hongchao',
      author_email='ji.hongchao@foxmail.com',
      url='https://github.com/hcji/PyFingerprint',
      data_files = [('', ['PyFingerprint/CDK/cdk-2.2.jar'])],
      packages=find_packages()
     )