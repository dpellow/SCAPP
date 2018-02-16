from setuptools import setup

setup(name='recycler2',
      version='0.1',
      description='Recycler2: an updated algorithm for detecting plasmids from de novo assembly graphs',
      url='https://github.com/dpellow/Recycler2',
      author='David Pellow',
      author_email = '',
      license='BSD-3-Clause',
      scripts = ['bin/recycle.py', 'bin/make_fasta_from_fastg.py', 'bin/get_simple_cycs.py'],
      packages = ['recyclelib'],
      requires=['python (<3.0)'],
      install_requires=[
        'networkx',
        'pysam',
        'nose',
        'numpy']
      )
