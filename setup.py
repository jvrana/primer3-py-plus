from distutils.core import setup

setup(name='primer3plus',
      packages=['primer3plus'],
      package_data={
          'primer3plus': ['primer3_params_raw.txt']
      },
      version='0.0.1'
     )
