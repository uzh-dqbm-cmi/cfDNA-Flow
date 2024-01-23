from setuptools import setup, find_packages

setup(name='cfDNA_pipeline',
      version='2021.18.11',
      description='cfDNA preprocessing pipeline',
      url='https://github.com/uzh-dqbm-cmi/cfDNA_pipeline',
      license='MIT',
      classifiers=[  # https://pypi.org/classifiers/
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.9',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          ],
      keywords='cfDNA ctDNA bioinformatics pipeline',
      python_requires='>3.9.0',
      packages=find_packages(exclude=['test']),
      include_package_data=True,
      # setup_requires=['sphinx', 'sphinx-argparse', 'sphinx-argparse-cli'],
      install_requires=[
          'numpy',
          'pandas',
          'snakemake'],
      entry_points={
          'console_scripts': [
              'cfDNA_pipeline = cfDNA:main'
          ],
      },
      zip_safe=False)
